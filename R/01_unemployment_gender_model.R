# =============================================================================
# Gender differentials in Canadian unemployment: velocity & acceleration model
#
#   UR_it = b0 + b1*Time + b2*Time^2 + b3*D + b4*(D x Time) + b5*(D x Time^2) + e
#
#   D = 1 for Men+, 0 for Women+ (women are the reference group, so POSITIVE
#   b3/b4/b5 directly support the hypothesis that men's unemployment level /
#   velocity / acceleration exceeds women's).
#
# Data: StatCan Table 14-10-0287-01, monthly SA, Canada, 15+, Jan 2017-Jun 2026
# Deps: lmtest, sandwich (install.packages(c("lmtest","sandwich")) if needed)
# Run from the project root: Rscript R/01_unemployment_gender_model.R
# =============================================================================

library(lmtest)
library(sandwich)

working_directory <- dirname(rstudioapi::documentPath())
setwd(working_directory)

raw_path <- "../data/1410028701-eng.csv"
dir.create("output", showWarnings = FALSE)

# ---- Helper Function ---------------------------------------------------------
save_results <- function(m, file) {
  ct  <- hac(m)
  out <- data.frame(term = rownames(ct), unclass(ct), check.names = FALSE)
  names(out) <- c("term", "estimate", "hac_se", "t_value", "p_value")
  out$signif <- cut(out$p_value, c(0, .001, .01, .05, .1, 1),
                    labels = c("***", "**", "*", ".", ""), include.lowest = TRUE)
  out$n_obs  <- c(nobs(m), rep(NA, nrow(out) - 1))
  out$adj_r2 <- c(round(summary(m)$adj.r.squared, 4), rep(NA, nrow(out) - 1))
  write.csv(out, file.path("output", file), row.names = FALSE)
}

# ---- 1. Clean the StatCan pivot export ---------------------------------------
# Layout: row of month labels; each stat is one row; cols 3..116 = Men+,
# cols 117..230 = Women+ (114 months each). We locate rows by label to be safe.

# StatCan exports have a BOM and a footnote block with unbalanced quotes that
# breaks read.csv, so: read lines, keep only the wide data block, then parse.
con <- file(raw_path, encoding = "UTF-8-BOM")
lines <- readLines(con, warn = FALSE)
close(con)
# Keep only the full-width rows: the narrow metadata/footnote lines would make
# read.csv infer a 1-column layout and wrap the data (it sniffs the first 5
# lines to set the column count). Everything we need is in the wide rows.
nf   <- suppressWarnings(count.fields(textConnection(lines), sep = ","))
keep <- which(!is.na(nf) & nf == max(nf, na.rm = TRUE))
raw  <- read.csv(text = paste(lines[keep], collapse = "\n"),
                 header = FALSE, stringsAsFactors = FALSE)

date_row <- which(raw[[1]] == "Labour force characteristics")[1]
ur_row   <- grep("^Unemployment rate", raw[[1]])[1]
gen_row  <- grep("^Gender", raw[[2]])[1]

men_start   <- which(raw[gen_row, ] == "Men+")[1]
women_start <- which(raw[gen_row, ] == "Women+")[1]
n_months    <- women_start - men_start

dates_chr <- as.character(raw[date_row, men_start:(men_start + n_months - 1)])
dates     <- as.Date(paste0("01 ", dates_chr), format = "%d %B %Y")
num <- function(x) as.numeric(gsub(",", "", as.character(x)))

lfs <- data.frame(
  date   = rep(dates, 2),
  gender = rep(c("Women", "Men"), each = n_months),
  ur     = c(num(raw[ur_row, women_start:(women_start + n_months - 1)]),
             num(raw[ur_row, men_start:(men_start + n_months - 1)]))
)

# Model variables: Time = months since Jan 2017 (0-indexed), D = 1 for men
lfs$time  <- rep(seq_len(n_months) - 1, 2)
lfs$time2 <- lfs$time^2
lfs$D     <- as.integer(lfs$gender == "Men")
# COVID shock indicator (Mar 2020 - Dec 2021) for the robustness spec
lfs$covid <- as.integer(lfs$date >= as.Date("2020-03-01") &
                        lfs$date <= as.Date("2021-12-31"))

write.csv(lfs, "output/lfs_tidy.csv", row.names = FALSE)
stopifnot(nrow(lfs) == 2 * n_months, !anyNA(lfs$ur))
cat(sprintf("Cleaned: %d months x 2 genders (%s to %s)\n\n",
            n_months, format(min(dates), "%b %Y"), format(max(dates), "%b %Y")))

# ---- Dataset summary: men vs women levels and shares ------------------------
# script: it only needs the parsed StatCan sheet (`cells` or `raw`) in memory.

if (!exists("cells")) cells <- trimws(gsub("\u00a0", " ", as.matrix(raw)))

.loc <- function(labels) {
  hit <- which(cells %in% labels)[1]
  c(row = (hit - 1) %% nrow(cells) + 1, col = (hit - 1) %/% nrow(cells) + 1)
}
.mcol <- .loc(c("Men+", "Men"))["col"]
.wcol <- .loc(c("Women+", "Women"))["col"]
.n    <- abs(.wcol - .mcol)

.pull <- function(label) {
  # first row whose label matches AND whose data type is seasonally adjusted
  i <- which(startsWith(cells[, 1], label) &
               grepl("adjusted", cells[, 2], ignore.case = TRUE))[1]
  num <- function(x) as.numeric(gsub(",", "", x))
  list(men   = num(cells[i, .mcol:(.mcol + .n - 1)]),
       women = num(cells[i, .wcol:(.wcol + .n - 1)]))
}

.stats <- c("Population", "Labour force", "Employment", "Full-time employment",
            "Part-time employment", "Unemployment", "Unemployment rate",
            "Participation rate", "Employment rate")

summary_tbl <- do.call(rbind, lapply(.stats, function(s) {
  x <- .pull(s)
  data.frame(statistic  = s,
             men_avg    = round(mean(x$men),   1),
             women_avg  = round(mean(x$women), 1),
             men_last   = tail(x$men,   1),
             women_last = tail(x$women, 1))
}))
cat("\n---- Summary (avg over sample; 'last' = most recent month) ----\n")
print(summary_tbl, row.names = FALSE)

cat("\n---- Men's share of each group (%) ----\n")
for (s in c("Population", "Labour force", "Employment", "Unemployment")) {
  x <- .pull(s)
  cat(sprintf("%-14s avg: %5.2f   latest: %5.2f\n", s,
              100 * mean(x$men) / (mean(x$men) + mean(x$women)),
              100 * tail(x$men, 1) / (tail(x$men, 1) + tail(x$women, 1))))
}
cat("\nNote: values are weighted population estimates (x 1,000), not survey\n")
cat("respondent counts; the LFS file does not include raw sample sizes.\n\n")
# -----------------------------------------------------------------------------


# ---- 2. Fit the models -------------------------------------------------------
# Newey-West (HAC) standard errors, 12 lags: monthly macro series are strongly
# autocorrelated and plain OLS SEs would wildly overstate significance.
hac <- function(m) coeftest(m, vcov = NeweyWest(m, lag = 12, prewhite = FALSE))

report <- function(m, label) {
  cat("=====", label, "=====\n")
  print(hac(m))
  cat(sprintf("Adj. R-squared: %.3f\n\n", summary(m)$adj.r.squared))
}

# (a) Your model, exactly as specified
m_main <- lm(ur ~ time + time2 + D + D:time + D:time2, data = lfs)
report(m_main, "Main model (full sample 2017-2026)")
save_results(m_main,  "model_main.csv")

# (b) Robustness: soak up the pandemic spike with a COVID indicator
m_covid <- lm(ur ~ time + time2 + D + D:time + D:time2 + covid + covid:D,
              data = lfs)
report(m_covid, "Robustness: + COVID dummy (Mar 2020 - Dec 2021)")
save_results(m_covid, "model_covid_dummy.csv")

# (c) Robustness: post-pandemic subsample only (Jan 2022+), time re-indexed
post <- subset(lfs, date >= as.Date("2022-01-01"))
post$time  <- post$time - min(post$time)
post$time2 <- post$time^2
m_post <- lm(ur ~ time + time2 + D + D:time + D:time2, data = post)
report(m_post, "Robustness: post-COVID subsample (2022+)")
save_results(m_post,  "model_post2022.csv")

# (d) The gap directly: (men - women) each month, regressed on time + time^2.
#     Equivalent to b3/b4/b5 but often the cleanest single readout.
gap <- data.frame(
  date  = dates,
  gap   = lfs$ur[lfs$gender == "Men"] - lfs$ur[lfs$gender == "Women"],
  time  = seq_len(n_months) - 1
)
gap$time2 <- gap$time^2
m_gap <- lm(gap ~ time + time2, data = gap)
report(m_gap, "Gap model: (Men - Women) ~ time + time^2")
save_results(m_gap,   "model_gap.csv")

# ---- 3. Plot -----------------------------------------------------------------
png("output/unemployment_by_gender.png", width = 1100, height = 900, res = 120)
par(mfrow = c(2, 1), mar = c(3, 4, 2.5, 1))

men_ur   <- lfs$ur[lfs$gender == "Men"]
women_ur <- lfs$ur[lfs$gender == "Women"]
plot(dates, men_ur, type = "l", col = "#2166ac", lwd = 2,
     ylim = range(lfs$ur), ylab = "Unemployment rate (%)", xlab = "",
     main = "Unemployment rate by gender, Canada 15+, SA")
lines(dates, women_ur, col = "#b2182b", lwd = 2)
lines(dates, predict(m_main)[lfs$gender == "Men"],   col = "#2166ac", lty = 2)
lines(dates, predict(m_main)[lfs$gender == "Women"], col = "#b2182b", lty = 2)
legend("topright", c("Men", "Women", "Quadratic fit"),
       col = c("#2166ac", "#b2182b", "grey40"),
       lty = c(1, 1, 2), lwd = c(2, 2, 1), bty = "n")

plot(gap$date, gap$gap, type = "h", col = ifelse(gap$gap >= 0, "#2166ac", "#b2182b"),
     ylab = "Men - Women (pp)", xlab = "", main = "Gender unemployment gap")
lines(gap$date, predict(m_gap), lwd = 2)
abline(h = 0, col = "grey50")
dev.off()

cat("Wrote output/lfs_tidy.csv and output/unemployment_by_gender.png\n")
