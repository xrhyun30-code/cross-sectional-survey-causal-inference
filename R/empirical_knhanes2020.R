## Empirical illustration: weekday-weekend sleep irregularity and obesity, KNHANES 2020 (adults >=19 y)
## Exposure: |weekend - weekday sleep hours| > 2 h (BP16_1, BP16_2; whole hours in the public file)
## Outcome: measured BMI >= 25 kg/m2 (HE_BMI)
## Codes below follow the labels/value patterns in hn20_all; verify against the official codebook.
suppressPackageStartupMessages({library(haven); library(survey)})
options(survey.lonely.psu = "adjust")

DATA <- Sys.getenv("KNHANES2020_SAS", "data/hn20_all.sas7bdat")   # public-use file from KDCA
OUT  <- "output"
dir.create(OUT, showWarnings = FALSE)

d0 <- as.data.frame(read_sas(DATA))
d0[] <- lapply(d0, function(x) if (inherits(x, "haven_labelled")) as.vector(x) else x)

## ---- sample flow ----
flow <- data.frame(step = character(), n = integer(), stringsAsFactors = FALSE)
add <- function(s, n) flow[nrow(flow) + 1, ] <<- list(s, n)
add("KNHANES 2020 participants", nrow(d0))
d <- d0[d0$age >= 19, ];                                   add("Aged >=19 years", nrow(d))
d <- d[!is.na(d$wt_itvex) & d$wt_itvex > 0, ];            add("With health interview/examination weight (wt_itvex)", nrow(d))
slp_ok <- function(x) !is.na(x) & x < 88
d <- d[slp_ok(d$BP16_1) & slp_ok(d$BP16_2), ];            add("Valid weekday and weekend sleep duration", nrow(d))
d <- d[!is.na(d$HE_BMI), ];                               add("Measured BMI available", nrow(d))

## ---- variables ----
d$A <- as.integer(abs(d$BP16_2 - d$BP16_1) > 2)
d$Y <- as.integer(d$HE_BMI >= 25)
d$sex <- factor(d$sex, 1:2, c("Men", "Women"))
d$agegrp <- cut(d$age, c(19, 30, 40, 50, 60, 70, Inf), right = FALSE,
                labels = c("19-29", "30-39", "40-49", "50-59", "60-69", "70+"))
d$income <- factor(d$incm, 1:4, c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"))
d$educ4 <- factor(d$edu, 1:4, c("<=Elementary", "Middle school", "High school", ">=College"))
d$smoke <- factor(ifelse(d$BS3_1 %in% 1:2, "Current", ifelse(d$BS3_1 == 3, "Former",
                  ifelse(d$BS3_1 == 8, "Never", NA))), c("Never", "Former", "Current"))
d$alcohol <- factor(ifelse(d$BD1_11 %in% c(1, 2, 8), "<1/month or never",
                    ifelse(d$BD1_11 %in% 3:4, "1-4/month", ifelse(d$BD1_11 %in% 5:6, ">=2/week", NA))),
                    c("<1/month or never", "1-4/month", ">=2/week"))
d$walk <- factor(ifelse(d$BE3_31 == 1, "None", ifelse(d$BE3_31 %in% 2:5, "1-4 days/wk",
                 ifelse(d$BE3_31 %in% 6:8, ">=5 days/wk", NA))), c("None", "1-4 days/wk", ">=5 days/wk"))
## BE3_31: 1 = no walking (BE3_32 is "not applicable" (88) for exactly these respondents), 2-8 = 1-7 days/week, 99 = unknown.
## Verified from within-file consistency (BE3_31 x BE3_32); confirm against the official codebook.

covs <- c("agegrp", "sex", "income", "educ4", "smoke", "alcohol", "walk")
miss <- sapply(covs, function(v) sum(is.na(d[[v]])))
write.csv(data.frame(variable = covs, n_missing = miss), file.path(OUT, "missing_covariates.csv"), row.names = FALSE)
cc <- complete.cases(d[, c(covs, "A", "Y")])
add("Complete covariate data (analytic sample)", sum(cc))
write.csv(flow, file.path(OUT, "flow.csv"), row.names = FALSE)
d$cc <- cc
d$agec <- (d$age - 50) / 10   # continuous age (for modelling)

## ---- design ----
des_full <- svydesign(ids = ~psu, strata = ~kstrata, weights = ~wt_itvex, nest = TRUE,
                      data = rbind(d[, names(d)], NULL))
rep_full <- as.svrepdesign(des_full, type = "JKn")
des <- subset(rep_full, cc)           # domain analysis keeps design-based variance
dat <- d[d$cc, ]

## ---- Table: weighted characteristics ----
cat_vars <- c("agegrp", "sex", "income", "educ4", "smoke", "alcohol", "walk")
rows <- list()
tot_row <- function(sub, lab) {
  data.frame(variable = "N", level = lab, n_A0 = NA, pct_A0 = NA, n_A1 = NA, pct_A1 = NA,
             n_all = NA, pct_all = NA, p = NA)
}
rows[[1]] <- data.frame(variable = "Participants (n)", level = "", n_A0 = sum(dat$A == 0), pct_A0 = NA,
                        n_A1 = sum(dat$A == 1), pct_A1 = NA, n_all = nrow(dat), pct_all = NA, p = NA)
for (v in cat_vars) {
  f <- as.formula(paste0("~", v))
  tab <- svyby(f, ~A, des, svymean, vartype = NULL)
  all <- svymean(f, des)
  p <- tryCatch(svychisq(as.formula(paste0("~", v, "+A")), des)$p.value, error = function(e) NA)
  lv <- levels(dat[[v]])
  for (l in lv) {
    nm <- paste0(v, l)
    rows[[length(rows) + 1]] <- data.frame(
      variable = v, level = l,
      n_A0 = sum(dat[[v]] == l & dat$A == 0), pct_A0 = 100 * tab[tab$A == 0, nm],
      n_A1 = sum(dat[[v]] == l & dat$A == 1), pct_A1 = 100 * tab[tab$A == 1, nm],
      n_all = sum(dat[[v]] == l), pct_all = 100 * unname(all[nm]),
      p = if (l == lv[1]) p else NA)
  }
}
## age as mean (SE)
ag <- svyby(~age, ~A, des, svymean); agall <- svymean(~age, des)
rows[[length(rows) + 1]] <- data.frame(variable = "age_mean", level = "mean", n_A0 = NA, pct_A0 = ag$age[1],
  n_A1 = NA, pct_A1 = ag$age[2], n_all = NA, pct_all = unname(agall["age"]),
  p = tryCatch(svyttest(age ~ A, des)$p.value, error = function(e) NA))
chars <- do.call(rbind, rows)
write.csv(chars, file.path(OUT, "table_characteristics.csv"), row.names = FALSE)

## ---- prevalence ----
prev <- svyby(~Y, ~A, des, svymean); prevall <- svymean(~Y, des); pA <- svymean(~A, des)
write.csv(data.frame(what = c("Y|A=0", "Y|A=1", "Y overall", "A overall"),
  est = c(prev$Y, coef(prevall), coef(pA)),
  se = c(prev$se, SE(prevall), SE(pA)), n = c(sum(dat$A == 0), sum(dat$A == 1), nrow(dat), nrow(dat))),
  file.path(OUT, "prevalence.csv"), row.names = FALSE)

## ---- estimators ----
fml_x <- paste(c("agec", "I(agec^2)", "sex", "income", "educ4", "smoke", "alcohol", "walk"), collapse = " + ")
f_ps <- as.formula(paste("A ~", fml_x)); f_y <- as.formula(paste("Y ~ A +", fml_x))
wmean <- function(x, w) sum(w * x) / sum(w)

fit_all <- function(w, data, trunc = NULL) {
  ok <- w > 0; dd <- data[ok, ]; dd$ww_ <- w[ok] / mean(w[ok]); ww <- dd$ww_   # rescale: estimators are scale-invariant, glm is not
  psf <- suppressWarnings(glm(f_ps, quasibinomial, dd, weights = ww_))
  ps <- predict(psf, type = "response")
  yf <- suppressWarnings(glm(f_y, quasibinomial, dd, weights = ww_))
  m1 <- predict(yf, transform(dd, A = 1), type = "response")
  m0 <- predict(yf, transform(dd, A = 0), type = "response")
  if (!is.null(trunc)) {   # truncate the stabilised-free IPW at the given percentiles
    iw <- ifelse(dd$A == 1, 1 / ps, 1 / (1 - ps)); q <- quantile(iw, trunc)
    iw <- pmin(pmax(iw, q[1]), q[2])
    ps <- ifelse(dd$A == 1, 1 / iw, 1 - 1 / iw)       # implied P(A=a|X) for the observed arm
  }
  ipw_w <- ifelse(dd$A == 1, 1 / ps, 1 / (1 - ps))
  p1c <- wmean(dd$Y[dd$A == 1], ww[dd$A == 1]); p0c <- wmean(dd$Y[dd$A == 0], ww[dd$A == 0])
  p1g <- wmean(m1, ww); p0g <- wmean(m0, ww)
  p1i <- wmean(dd$Y[dd$A == 1], (ww * ipw_w)[dd$A == 1]); p0i <- wmean(dd$Y[dd$A == 0], (ww * ipw_w)[dd$A == 0])
  mu1 <- dd$A * (dd$Y - m1) / ps + m1; mu0 <- (1 - dd$A) * (dd$Y - m0) / (1 - ps) + m0
  p1a <- wmean(mu1, ww); p0a <- wmean(mu0, ww)
  c(crude_rd = p1c - p0c, crude_rr = p1c / p0c, logit_or = unname(exp(coef(yf)["A"])),
    gcomp_rd = p1g - p0g, gcomp_rr = p1g / p0g,
    ipw_rd = p1i - p0i, ipw_rr = p1i / p0i,
    aipw_rd = p1a - p0a, aipw_rr = p1a / p0a)
}

main <- withReplicates(des, function(w, data) fit_all(w, data), return.replicates = TRUE)
est <- coef(main); se <- SE(main)
## RR / OR CIs on log scale
is_ratio <- grepl("_rr$|_or$", names(est))
lo <- ifelse(is_ratio, exp(log(est) - 1.96 * se / est), est - 1.96 * se)
hi <- ifelse(is_ratio, exp(log(est) + 1.96 * se / est), est + 1.96 * se)
res <- data.frame(analysis = "Primary (untruncated)", measure = names(est), est = est, se = se, lo = lo, hi = hi)

trn <- withReplicates(des, function(w, data) fit_all(w, data, trunc = c(0.01, 0.99)), return.replicates = TRUE)
e2 <- coef(trn)[c("ipw_rd", "ipw_rr", "aipw_rd", "aipw_rr")]; s2 <- SE(trn)[c("ipw_rd", "ipw_rr", "aipw_rd", "aipw_rr")]
r2 <- grepl("_rr$", names(e2))
res <- rbind(res, data.frame(analysis = "IPW weights truncated at 1st/99th pct", measure = names(e2), est = e2, se = s2,
  lo = ifelse(r2, exp(log(e2) - 1.96 * s2 / e2), e2 - 1.96 * s2),
  hi = ifelse(r2, exp(log(e2) + 1.96 * s2 / e2), e2 + 1.96 * s2)))
write.csv(res, file.path(OUT, "estimates.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---- balance / overlap / weights ----
w0 <- dat$wt_itvex
psf <- glm(f_ps, quasibinomial, transform(dat, w0_ = w0 / mean(w0)), weights = w0_); ps <- predict(psf, type = "response")
ipw <- ifelse(dat$A == 1, 1 / ps, 1 / (1 - ps)); wc <- w0 * ipw
X <- model.matrix(as.formula(paste("~", fml_x)), dat)[, -1, drop = FALSE]
wsd <- function(x, w) sqrt(sum(w * (x - wmean(x, w))^2) / sum(w))
smd <- function(x, w) {
  m1 <- wmean(x[dat$A == 1], w[dat$A == 1]); m0 <- wmean(x[dat$A == 0], w[dat$A == 0])
  s <- sqrt((wsd(x[dat$A == 1], w[dat$A == 1])^2 + wsd(x[dat$A == 0], w[dat$A == 0])^2) / 2)
  (m1 - m0) / s
}
bal <- data.frame(covariate = colnames(X),
  smd_survey_weighted = apply(X, 2, smd, w = w0),
  smd_after_ipw = apply(X, 2, smd, w = wc))
write.csv(bal, file.path(OUT, "balance.csv"), row.names = FALSE)
ess <- function(w) sum(w)^2 / sum(w^2)
wd <- data.frame(
  what = c("Survey weight", "IPW (unstabilized)", "Combined weight (survey x IPW)"),
  min = c(min(w0), min(ipw), min(wc)), p1 = c(quantile(w0, .01), quantile(ipw, .01), quantile(wc, .01)),
  median = c(median(w0), median(ipw), median(wc)), p99 = c(quantile(w0, .99), quantile(ipw, .99), quantile(wc, .99)),
  max = c(max(w0), max(ipw), max(wc)), ess = c(ess(w0), ess(ipw), ess(wc)),
  ess_by_A1 = c(ess(w0[dat$A == 1]), ess((ipw)[dat$A == 1]), ess(wc[dat$A == 1])),
  ess_by_A0 = c(ess(w0[dat$A == 0]), ess((ipw)[dat$A == 0]), ess(wc[dat$A == 0])))
write.csv(wd, file.path(OUT, "weight_diagnostics.csv"), row.names = FALSE)
ps_sum <- data.frame(group = c("A=0", "A=1"),
  min = tapply(ps, dat$A, min), p5 = tapply(ps, dat$A, quantile, .05), median = tapply(ps, dat$A, median),
  p95 = tapply(ps, dat$A, quantile, .95), max = tapply(ps, dat$A, max))
write.csv(ps_sum, file.path(OUT, "ps_overlap.csv"), row.names = FALSE)
cat("max |SMD| before:", round(max(abs(bal$smd_survey_weighted)), 3), " after IPW:", round(max(abs(bal$smd_after_ipw)), 3), "\n")
print(res, digits = 4); print(flow); print(prev)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
