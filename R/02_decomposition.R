# =============================================================================
# 02: Who carries the male unemployment burden?
# Shift-share / Oaxaca-style decomposition of the male-female unemployment
# rate gap into:
#   COMPOSITION: men are over-weighted in high-unemployment age-industry cells
#   WITHIN:      men have higher unemployment inside the same cell
# using symmetric (two-fold, mean-weighted) weights so the split does not
# depend on which gender is treated as the baseline. The identity is exact:
#   UR_men - UR_women = sum_c comp_c + sum_c within_c
#
# Data: StatCan 14-10-0022-01, monthly UNADJUSTED, Canada, Men+/Women+,
#       ages 15-24 / 25-54 / 55+, all industries, Jan 2017 - Jun 2026.
# Seasonality is handled by aggregating each period to annual-equivalent
# sums (all months of the period), never comparing raw months.
# Run from project root: Rscript R/02_decomposition.R
# =============================================================================

library(lmtest)
library(sandwich)

working_directory <- dirname(rstudioapi::documentPath())
setwd(working_directory)

raw_path <- "../data/1410002201-eng.csv"
dir.create("output", showWarnings = FALSE)

# ---- 1. Parse the 4-level pivot header ---------------------------------------
con <- file(raw_path, encoding = "UTF-8-BOM")
lines <- readLines(con, warn = FALSE); close(con)
# Force the column count instead of letting read.csv sniff it from the narrow
# metadata lines (count.fields is also unreliable on some StatCan rows, so we
# take the max comma count as the width and let short rows pad with blanks).
width <- max(vapply(lines, function(l)
  length(scan(text = l, what = "", sep = ",", quiet = TRUE)), 1L)) 
raw <- read.csv(text = paste(lines, collapse = "\n"), header = FALSE,
                fill = TRUE, col.names = paste0("V", seq_len(width)),
                colClasses = "character")
cells <- trimws(gsub("\u00a0", " ", as.matrix(raw)))

strip_fn <- function(x) trimws(gsub(" [0-9 ]+$", "", x))  # drop footnote refs

stat_row <- which(cells[, 1] == "Labour force characteristics")[1]
gen_row  <- which(startsWith(cells[, 1], "Gender"))[1]
age_row  <- which(cells[, 1] == "Age group")[1]
mon_row  <- which(startsWith(cells[, 1], "North American Industry"))[1]
stopifnot(!anyNA(c(stat_row, gen_row, age_row, mon_row)))

ffill <- function(x) { x[x == ""] <- NA
  for (i in seq_along(x)[-1]) if (is.na(x[i])) x[i] <- x[i - 1]
  x }

stat_f <- strip_fn(ffill(cells[stat_row, ]))
gen_f  <- ffill(cells[gen_row, ])
age_f  <- strip_fn(gsub("\\s*\\(x 1,000\\)", "", ffill(cells[age_row, ])))
mons   <- cells[mon_row, ]

suppressWarnings({ ok <- Sys.setlocale("LC_TIME", "C")
                   if (identical(ok, "")) Sys.setlocale("LC_TIME", "English") })
dates_all <- as.Date(paste0("01 ", mons), format = "%d %B %Y")
data_cols <- which(!is.na(dates_all))

# ---- 2. Tidy long format for the 17 mutually exclusive top-level industries --
industries <- c(
  "Agriculture", "Forestry, fishing, mining, quarrying, oil and gas",
  "Utilities", "Construction", "Manufacturing", "Wholesale and retail trade",
  "Transportation and warehousing",
  "Finance, insurance, real estate, rental and leasing",
  "Professional, scientific and technical services",
  "Business, building and other support services", "Educational services",
  "Health care and social assistance", "Information, culture and recreation",
  "Accommodation and food services",
  "Other services (except public administration)", "Public administration",
  "Unclassified industries")

row_lab  <- strip_fn(cells[, 1])
ind_rows <- match(c(industries, "Total, all industries"), row_lab)
stopifnot(!anyNA(ind_rows))

num <- function(x) suppressWarnings(as.numeric(gsub(",", "", x)))

long <- do.call(rbind, lapply(ind_rows, function(i) data.frame(
  industry = row_lab[i], stat = stat_f[data_cols], gender = gen_f[data_cols],
  age = age_f[data_cols], date = dates_all[data_cols],
  value = num(cells[i, data_cols]))))

# StatCan suppresses cells under ~1,500 persons ('x'); treat as 0 -- the bias
# is bounded by 1.5 (thousand) per suppressed cell and is negligible here.
n_supp <- sum(cells[ind_rows, data_cols] == "x")
long$value[is.na(long$value)] <- 0
cat(sprintf("Parsed %d values (%d suppressed cells set to 0)\n",
            nrow(long), n_supp))

wide <- merge(
  setNames(subset(long, stat == "Employment",
                  c("industry", "gender", "age", "date", "value")),
           c("industry", "gender", "age", "date", "emp")),
  setNames(subset(long, stat == "Unemployment",
                  c("industry", "gender", "age", "date", "value")),
           c("industry", "gender", "age", "date", "unemp")))
wide$lf <- wide$emp + wide$unemp

# ---- 3. Validation against StatCan's published unemployment rates ------------
pub <- subset(long, stat == "Unemployment rate" &
                    industry == "Total, all industries")
tot <- subset(wide, industry == "Total, all industries")
chk <- merge(tot, pub[c("gender", "age", "date", "value")])
chk <- subset(chk, lf > 0 & value > 0)
max_dev <- max(abs(100 * chk$unemp / chk$lf - chk$value))
cat(sprintf("Validation: computed vs published UR, max deviation %.3f pp\n",
            max_dev))
stopifnot(max_dev < 0.15)   # published rates are rounded to 0.1

cell <- subset(wide, industry != "Total, all industries")

# ---- 4. Decomposition --------------------------------------------------------
# For a set of months, sum flows per gender x (age x industry) cell, then:
#   u_gc = cell unemployment rate; s_gc = cell share of gender g labour force
#   comp_c   = (s_mc - s_wc) * mean(u_mc, u_wc)
#   within_c = mean(s_mc, s_wc) * (u_mc - u_wc)
decompose <- function(df, from, to, label) {
  d <- subset(df, date >= from & date <= to)
  a <- aggregate(cbind(emp, unemp, lf) ~ gender + age + industry, d, sum)
  m <- subset(a, gender == "Men+");  w <- subset(a, gender == "Women+")
  key <- c("age", "industry")
  mw  <- merge(m[c(key, "unemp", "lf")], w[c(key, "unemp", "lf")],
               by = key, suffixes = c("_m", "_w"))
  mw$u_m <- 100 * mw$unemp_m / mw$lf_m
  mw$u_w <- 100 * mw$unemp_w / mw$lf_w
  mw$u_m[!is.finite(mw$u_m)] <- 0;  mw$u_w[!is.finite(mw$u_w)] <- 0
  mw$s_m <- mw$lf_m / sum(mw$lf_m); mw$s_w <- mw$lf_w / sum(mw$lf_w)
  mw$composition <- (mw$s_m - mw$s_w) * (mw$u_m + mw$u_w) / 2
  mw$within      <- (mw$s_m + mw$s_w) / 2 * (mw$u_m - mw$u_w)
  mw$total_contrib <- mw$composition + mw$within
  ur_m <- 100 * sum(mw$unemp_m) / sum(mw$lf_m)
  ur_w <- 100 * sum(mw$unemp_w) / sum(mw$lf_w)
  cat(sprintf("\n===== %s =====\n", label))
  cat(sprintf("UR men %.2f | UR women %.2f | gap %+.2f pp\n", ur_m, ur_w,
              ur_m - ur_w))
  cat(sprintf("  composition (industry/age mix): %+.2f pp (%.0f%% of gap)\n",
              sum(mw$composition), 100 * sum(mw$composition) / (ur_m - ur_w)))
  cat(sprintf("  within-cell (rate differences): %+.2f pp (%.0f%% of gap)\n",
              sum(mw$within), 100 * sum(mw$within) / (ur_m - ur_w)))
  top <- mw[order(-mw$total_contrib),
            c("age", "industry", "u_m", "u_w", "s_m", "s_w",
              "composition", "within", "total_contrib")]
  cat("\nTop 10 cells driving the male excess:\n")
  print(head(data.frame(lapply(top, function(x)
    if (is.numeric(x)) round(x, 3) else x)), 10), row.names = FALSE)
  mw$period <- label
  invisible(mw)
}

pre  <- decompose(cell, as.Date("2017-01-01"), as.Date("2019-12-01"),
                  "Pre-COVID (2017-2019)")
post <- decompose(cell, as.Date("2025-07-01"), as.Date("2026-06-01"),
                  "Latest 12 months (Jul 2025 - Jun 2026)")

both <- rbind(pre, post)
write.csv(both[order(both$period, -both$total_contrib),
               c("period", "age", "industry", "u_m", "u_w", "s_m", "s_w",
                 "composition", "within", "total_contrib")],
          "output/decomposition_cells.csv", row.names = FALSE)

# ---- 5. Plot: top contributing cells, both periods ---------------------------
lab <- function(d) paste0(sub(" years.*", "", d$age), " | ",
                          substr(d$industry, 1, 32))
top_cells <- unique(c(head(pre$industry[order(-pre$total_contrib)], 0),
                      lab(head(pre[order(-pre$total_contrib), ], 8)),
                      lab(head(post[order(-post$total_contrib), ], 8))))
pre$cell  <- lab(pre);  post$cell <- lab(post)
p1 <- pre[match(top_cells, pre$cell), "total_contrib"]
p2 <- post[match(top_cells, post$cell), "total_contrib"]

png("output/decomposition_top_cells.png", width = 1500, height = 700, res = 110)
par(mar = c(4, 16, 3, 1))
o <- order(p2)
barplot(t(cbind(p1[o], p2[o])), beside = TRUE, horiz = TRUE,
        names.arg = top_cells[o], las = 1, cex.names = 0.75,
        col = c("grey70", "#2166ac"),
        xlab = "Contribution to male-female UR gap (pp)",
        main = "Which age x industry cells drive men's excess unemployment")
legend("bottomright", c("2017-2019", "Jul 2025 - Jun 2026"),
       fill = c("grey70", "#2166ac"), bty = "n")
abline(v = 0)
dev.off()
cat("\nWrote output/decomposition_cells.csv and output/decomposition_top_cells.png\n")

# ---- Top 10 cells on the women's side ---------------------------------------
for (p in unique(both$period)) {
  d <- both[both$period == p, ]
  d <- d[order(d$total_contrib), ]          # most negative first
  cat(sprintf("\n===== %s -- top 10 cells offsetting the male excess =====\n", p))
  print(head(data.frame(lapply(d[c("age", "industry", "u_m", "u_w",
                                   "s_m", "s_w", "composition", "within", "total_contrib")],
                               function(x) if (is.numeric(x)) round(x, 3) else x)), 10),
        row.names = FALSE)
}

# ---- Diverging plot: both sides of the gap ----------------------------------
sel <- c(head(post$cell[order(-post$total_contrib)], 8),
         rev(head(post$cell[order(post$total_contrib)], 8)))
p1 <- pre$total_contrib[match(sel, pre$cell)]
p2 <- post$total_contrib[match(sel, post$cell)]
o  <- rev(seq_along(sel))
colmat <- rbind(rep("grey75", length(sel)),
                ifelse(p2[o] >= 0, "#2166ac", "#b2182b"))
png("output/decomposition_diverging.png", width = 1500, height = 850, res = 110)
par(mar = c(4.5, 15.5, 4, 2))
barplot(t(cbind(p1, p2))[, o], beside = TRUE, horiz = TRUE,
        names.arg = sel[o], las = 1, cex.names = 0.72,
        xlim = range(c(p1, p2, 0)) * 1.15, col = colmat, border = "grey30",
        xlab = "Contribution to male-female unemployment rate gap (pp)",
        main = "Both sides of the gender unemployment gap")
mtext("Right of zero: pushes men's rate above women's. Left: offsets it. Grey = 2017-2019, colour = Jul 2025 - Jun 2026",
      side = 3, line = 0.3, cex = 0.78)
abline(v = 0, lwd = 1.5)
legend("bottomright", c("2017-2019", "Latest: male excess", "Latest: female offset"),
       fill = c("grey75", "#2166ac", "#b2182b"), bty = "n", cex = 0.8)
dev.off()
