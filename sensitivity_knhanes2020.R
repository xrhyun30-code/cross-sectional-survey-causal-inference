## Exploratory (post hoc) sensitivity analyses for the KNHANES 2020 illustration.
## Re-uses data preparation and estimators from empirical_knhanes2020.R (evaluated, not copied).
src <- readLines("R/empirical_knhanes2020.R")
i_des <- grep("^## ---- design ----", src)[1]
i_est <- grep("^## ---- estimators ----", src)[1]
i_main <- grep("^main <- withReplicates", src)[1]
eval(parse(text = src[1:(i_des - 1)]))                 # data prep (d, covs, cc, A, Y)
eval(parse(text = src[(i_est + 1):(i_main - 1)]))      # f_ps, f_y, fit_all
OUT <- "output"

dat <- d[d$cc, ]
dat$diff <- abs(dat$BP16_2 - dat$BP16_1)
dat$BMI <- dat$HE_BMI
des_full <- svydesign(ids = ~psu, strata = ~kstrata, weights = ~wt_itvex, nest = TRUE, data = dat)
rep_full <- as.svrepdesign(des_full, type = "JKn")

specs <- list(
  list(id = "Primary", label = "Primary: difference > 2 h; BMI >= 25", A = quote(diff > 2), Y = quote(BMI >= 25), sub = NULL),
  list(id = "E1", label = "Exposure: difference > 1 h", A = quote(diff > 1), Y = quote(BMI >= 25), sub = NULL),
  list(id = "E2", label = "Exposure: difference > 3 h", A = quote(diff > 3), Y = quote(BMI >= 25), sub = NULL),
  list(id = "O1", label = "Outcome: BMI >= 30", A = quote(diff > 2), Y = quote(BMI >= 30), sub = NULL),
  list(id = "P1", label = "Adults aged 19-39 years", A = quote(diff > 2), Y = quote(BMI >= 25), sub = quote(age < 40)),
  list(id = "P2", label = "Adults aged 40 years or older", A = quote(diff > 2), Y = quote(BMI >= 25), sub = quote(age >= 40))
)
res <- list()
for (sp in specs) {
  dd <- dat; dd$A <- as.integer(eval(sp$A, dd)); dd$Y <- as.integer(eval(sp$Y, dd))
  des_sp <- as.svrepdesign(svydesign(ids = ~psu, strata = ~kstrata, weights = ~wt_itvex, nest = TRUE, data = dd), type = "JKn")
  if (!is.null(sp$sub)) des_sp <- subset(des_sp, eval(sp$sub, dd))
  nA <- sum(dd$A[if (is.null(sp$sub)) TRUE else eval(sp$sub, dd)])
  n <- if (is.null(sp$sub)) nrow(dd) else sum(eval(sp$sub, dd))
  r <- tryCatch(withReplicates(des_sp, function(w, data) fit_all(w, data), return.replicates = FALSE), error = function(e) NULL)
  if (is.null(r)) { cat("failed", sp$id, "\n"); next }
  est <- coef(r); se <- SE(r)
  pA <- coef(svymean(~A, des_sp)); pY <- coef(svymean(~Y, des_sp))
  ratio <- grepl("_rr$|_or$", names(est))
  res[[sp$id]] <- data.frame(id = sp$id, label = sp$label, n = n, n_exposed = nA, wprev_exposed = pA, wprev_outcome = pY,
    measure = names(est), est = est, se = se,
    lo = ifelse(ratio, exp(log(est) - 1.96 * se / est), est - 1.96 * se),
    hi = ifelse(ratio, exp(log(est) + 1.96 * se / est), est + 1.96 * se))
  cat(sp$id, "done; n exposed", nA, "\n")
}
write.csv(do.call(rbind, res), file.path(OUT, "sensitivity_estimates.csv"), row.names = FALSE)
