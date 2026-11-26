#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(mice)
  library(parallel)
})

options(survey.lonely.psu = "adjust")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}

N_SIM <- as.integer(arg_value("--sim-reps", "1000"))
M_IMP <- as.integer(arg_value("--mice-m", "20"))
MAXIT <- as.integer(arg_value("--mice-maxit", "10"))
SCENARIOS <- strsplit(arg_value("--scenarios", "S3a,S3b,S6"), ",", fixed = TRUE)[[1L]]
RUN_KNHANES <- !("--skip-knhanes" %in% args)
RUN_SIMULATIONS <- !("--skip-simulations" %in% args)
N_CORES <- min(4L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
SEED <- 20261126L

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1L]])
script_path <- normalizePath(gsub("~\\+~", " ", script_path))
submission_dir <- dirname(script_path)
project_dir <- dirname(submission_dir)
output_dir <- file.path(submission_dir, "analysis_outputs")
dir.create(output_dir, showWarnings = FALSE)

data_path <- Sys.getenv(
  "KNHANES2020_SAS",
  unset = file.path("data", "hn20_all.sas7bdat")
)
if (RUN_KNHANES && !file.exists(data_path)) {
  stop("KNHANES data file not found: ", data_path)
}

## ---------------------- S3a: outcome-model-only misspecification
N_PSU <- 400L
PSU_SIZE <- 250L
N_STRATA <- 20L
PSU_SAMPLE_PER_STRATUM <- 3L

make_s3a_population <- function(seed) {
  set.seed(seed)
  n <- N_PSU * PSU_SIZE
  psu <- rep(seq_len(N_PSU), each = PSU_SIZE)
  strata <- rep(seq_len(N_STRATA), each = N_PSU / N_STRATA * PSU_SIZE)
  u_psu <- rnorm(N_PSU, 0, 0.3)[psu]
  x1 <- rnorm(n) + u_psu
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)
  L <- 0.5 * x1 + 0.5 * x2 + 0.3 * x3
  h <- 0.8 * (x1^2 - 1) + 0.5 * x2 * x3
  p_a <- plogis(-0.5 + 0.5 * L)
  a <- rbinom(n, 1, p_a)
  p_y1 <- plogis(-1 + 0.5 * L + h + 0.5)
  p_y0 <- plogis(-1 + 0.5 * L + h)
  y <- rbinom(n, 1, ifelse(a == 1, p_y1, p_y0))
  pi_i <- pmin(0.9, 0.2 * exp(0.4 * x1))
  pi_i <- pi_i / mean(pi_i) * 0.2
  list(
    population = data.frame(psu, strata, X1 = x1, X2 = x2, X3 = x3,
                            A = a, Y = y, pi_i),
    truth = mean(p_y1 - p_y0),
    prev_a = mean(a),
    prev_y = mean(y)
  )
}

make_s3b_population <- function(seed) {
  set.seed(seed)
  n <- N_PSU * PSU_SIZE
  psu <- rep(seq_len(N_PSU), each = PSU_SIZE)
  strata <- rep(seq_len(N_STRATA), each = N_PSU / N_STRATA * PSU_SIZE)
  u_psu <- rnorm(N_PSU, 0, 0.3)[psu]
  x1 <- rnorm(n) + u_psu
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)
  L <- 0.5 * x1 + 0.5 * x2 + 0.3 * x3
  h <- 0.8 * (x1^2 - 1) + 0.5 * x2 * x3
  p_a <- plogis(-0.5 + 0.5 * L + h)
  a <- rbinom(n, 1, p_a)
  p_y1 <- plogis(-1 + 0.5 * L + h + 0.5)
  p_y0 <- plogis(-1 + 0.5 * L + h)
  y <- rbinom(n, 1, ifelse(a == 1, p_y1, p_y0))
  pi_i <- pmin(0.9, 0.2 * exp(0.4 * x1))
  pi_i <- pi_i / mean(pi_i) * 0.2
  list(
    population = data.frame(psu, strata, X1 = x1, X2 = x2, X3 = x3,
                            A = a, Y = y, pi_i),
    truth = mean(p_y1 - p_y0),
    prev_a = mean(a),
    prev_y = mean(y)
  )
}

make_s6_population <- function(seed) {
  set.seed(seed)
  n <- N_PSU * PSU_SIZE
  psu <- rep(seq_len(N_PSU), each = PSU_SIZE)
  strata <- rep(seq_len(N_STRATA), each = N_PSU / N_STRATA * PSU_SIZE)
  u_psu <- rnorm(N_PSU, 0, 0.3)[psu]
  x1 <- rnorm(n) + u_psu
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)
  L <- 0.5 * x1 + 0.5 * x2 + 0.3 * x3
  y0 <- 0.6 * x1 + rnorm(n)
  p_a <- plogis(-0.5 + 0.5 * L + 1.5 * y0)
  a <- rbinom(n, 1, p_a)
  p_y1 <- plogis(-1 + 0.5 * L + 0.5 + 1.5 * y0)
  p_y0 <- plogis(-1 + 0.5 * L + 1.5 * y0)
  y <- rbinom(n, 1, ifelse(a == 1, p_y1, p_y0))
  pi_i <- pmin(0.9, 0.2 * exp(0.4 * x1))
  pi_i <- pi_i / mean(pi_i) * 0.2
  list(
    population = data.frame(psu, strata, X1 = x1, X2 = x2, X3 = x3,
                            A = a, Y = y, pi_i),
    truth = mean(p_y1 - p_y0),
    prev_a = mean(a),
    prev_y = mean(y)
  )
}

draw_survey_sample <- function(population) {
  psu_tab <- unique(population[, c("psu", "strata")])
  selected <- unlist(lapply(
    split(psu_tab$psu, psu_tab$strata),
    function(z) sample(z, PSU_SAMPLE_PER_STRATUM)
  ), use.names = FALSE)
  sample <- population[population$psu %in% selected, , drop = FALSE]
  sample <- sample[runif(nrow(sample)) < sample$pi_i, , drop = FALSE]
  pi_psu <- PSU_SAMPLE_PER_STRATUM / (N_PSU / N_STRATA)
  sample$w <- 1 / (pi_psu * sample$pi_i)
  sample
}

wmean <- function(x, w) sum(w * x) / sum(w)

estimate_four <- function(d) {
  d$w <- d$w / mean(d$w)
  w <- d$w
  crude <- wmean(d$Y[d$A == 1], w[d$A == 1]) -
    wmean(d$Y[d$A == 0], w[d$A == 0])
  ps_fit <- suppressWarnings(
    glm(A ~ X1 + X2 + X3, quasibinomial(), data = d, weights = w)
  )
  ps <- predict(ps_fit, type = "response")
  y_fit <- suppressWarnings(
    glm(Y ~ A + X1 + X2 + X3, quasibinomial(), data = d, weights = w)
  )
  m1 <- predict(y_fit, transform(d, A = 1), type = "response")
  m0 <- predict(y_fit, transform(d, A = 0), type = "response")
  standardization <- wmean(m1 - m0, w)
  ipw <- wmean(d$Y[d$A == 1], (w / ps)[d$A == 1]) -
    wmean(d$Y[d$A == 0], (w / (1 - ps))[d$A == 0])
  psi <- (d$A * (d$Y - m1) / ps + m1) -
    ((1 - d$A) * (d$Y - m0) / (1 - ps) + m0)
  aipw <- wmean(psi, w)
  c(crude = crude, standardization = standardization,
    ipw = ipw, aipw = aipw)
}

jkn_estimate <- function(sample) {
  full <- estimate_four(sample)
  psus <- unique(sample[, c("psu", "strata")])
  variance <- setNames(numeric(length(full)), names(full))
  for (h in unique(psus$strata)) {
    ph <- psus$psu[psus$strata == h]
    nh <- length(ph)
    if (nh < 2L) stop("Fewer than two sampled PSUs remain in a stratum.")
    for (j in ph) {
      replicate <- sample[sample$psu != j, , drop = FALSE]
      idx <- replicate$strata == h
      replicate$w[idx] <- replicate$w[idx] * nh / (nh - 1)
      estimate <- estimate_four(replicate)
      variance <- variance + (nh - 1) / nh * (estimate - full)^2
    }
  }
  cbind(estimate = full, se = sqrt(variance))
}

run_simulation_scenario <- function(scenario, population, seed_offset) {
  cat(scenario, "true RD:", format(population$truth, digits = 8), "\n")
  run_one <- function(i) {
    set.seed(SEED + seed_offset + i)
    sample <- draw_survey_sample(population$population)
    result <- jkn_estimate(sample)
    data.frame(
      replicate = i,
      estimator = rownames(result),
      estimate = result[, "estimate"],
      se = result[, "se"],
      truth = population$truth,
      n = nrow(sample)
    )
  }
  res <- parallel::mclapply(seq_len(N_SIM), run_one, mc.cores = N_CORES)
  if (any(vapply(res, inherits, logical(1), "try-error"))) {
    stop("At least one ", scenario, " replicate failed.")
  }
  draws <- do.call(rbind, res)
  draws$bias <- draws$estimate - draws$truth
  draws$covered <- abs(draws$estimate - draws$truth) <= 1.96 * draws$se
  summary <- do.call(rbind, lapply(split(draws, draws$estimator), function(x) {
    data.frame(
      scenario = scenario,
      estimator = x$estimator[[1L]],
      truth = unique(x$truth),
      mean_estimate = mean(x$estimate),
      bias = mean(x$bias),
      empirical_se = sd(x$estimate),
      mean_jkn_se = mean(x$se),
      rmse = sqrt(mean(x$bias^2)),
      coverage = mean(x$covered),
      mean_n = mean(x$n),
      mcse_coverage = sqrt(mean(x$covered) * (1 - mean(x$covered)) / nrow(x))
    )
  }))
  write.csv(draws, file.path(output_dir, paste0(scenario, "_replicates.csv")),
            row.names = FALSE)
  write.csv(summary, file.path(output_dir, paste0(scenario, "_summary.csv")),
            row.names = FALSE)
  write.csv(data.frame(
    scenario = scenario, truth = population$truth,
    exposure_prevalence = population$prev_a,
    outcome_prevalence = population$prev_y,
    replicates = N_SIM
  ), file.path(output_dir, paste0(scenario, "_population.csv")), row.names = FALSE)
  summary
}

## ---------------------- KNHANES analytic data and replicate design
make_knhanes_data <- function() {
  d0 <- as.data.frame(haven::read_sas(data_path))
  d0[] <- lapply(d0, function(x) {
    if (inherits(x, "haven_labelled")) as.vector(x) else x
  })
  valid_sleep <- !is.na(d0$BP16_1) & d0$BP16_1 < 88 &
    !is.na(d0$BP16_2) & d0$BP16_2 < 88
  keep <- !is.na(d0$age) & d0$age >= 19 &
    !is.na(d0$wt_itvex) & d0$wt_itvex > 0 &
    valid_sleep & !is.na(d0$HE_BMI)
  d <- d0[keep, , drop = FALSE]
  d$A <- as.integer(abs(d$BP16_2 - d$BP16_1) > 2)
  d$Y <- as.integer(d$HE_BMI >= 25)
  d$sex <- factor(d$sex, 1:2, c("Men", "Women"))
  d$agec <- (d$age - 50) / 10
  d$agec2 <- d$agec^2
  d$income <- factor(d$incm, 1:4, c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"))
  d$educ4 <- factor(d$edu, 1:4,
                    c("<=Elementary", "Middle school", "High school", ">=College"))
  d$smoke <- factor(
    ifelse(d$BS3_1 %in% 1:2, "Current",
           ifelse(d$BS3_1 == 3, "Former",
                  ifelse(d$BS3_1 == 8, "Never", NA))),
    c("Never", "Former", "Current")
  )
  d$alcohol <- factor(
    ifelse(d$BD1_11 %in% c(1, 2, 8), "<1/month or never",
           ifelse(d$BD1_11 %in% 3:4, "1-4/month",
                  ifelse(d$BD1_11 %in% 5:6, ">=2/week", NA))),
    c("<1/month or never", "1-4/month", ">=2/week")
  )
  d$walk <- factor(
    ifelse(d$BE3_31 == 1, "None",
           ifelse(d$BE3_31 %in% 2:5, "1-4 days/wk",
                  ifelse(d$BE3_31 %in% 6:8, ">=5 days/wk", NA))),
    c("None", "1-4 days/wk", ">=5 days/wk")
  )
  d$kstrata <- factor(d$kstrata)
  d
}

covariates <- c("sex", "income", "educ4", "smoke", "alcohol", "walk")
model_terms <- "agec + agec2 + sex + income + educ4 + smoke + alcohol + walk"
ps_formula <- as.formula(paste("A ~", model_terms))
outcome_formula <- as.formula(paste("Y ~ A +", model_terms))

fit_survey_estimators <- function(rep_weights, data, ps_weighted = TRUE) {
  keep <- is.finite(rep_weights) & rep_weights > 0
  d <- data[keep, , drop = FALSE]
  w <- rep_weights[keep]
  w <- w / mean(w)
  d$analysis_weight <- w
  ps_fit <- if (ps_weighted) {
    suppressWarnings(
      glm(ps_formula, quasibinomial(), data = d, weights = analysis_weight)
    )
  } else {
    suppressWarnings(glm(ps_formula, quasibinomial(), data = d))
  }
  ps <- predict(ps_fit, type = "response")
  y_fit <- suppressWarnings(
    glm(outcome_formula, quasibinomial(), data = d, weights = analysis_weight)
  )
  m1 <- predict(y_fit, transform(d, A = 1), type = "response")
  m0 <- predict(y_fit, transform(d, A = 0), type = "response")
  crude1 <- wmean(d$Y[d$A == 1], w[d$A == 1])
  crude0 <- wmean(d$Y[d$A == 0], w[d$A == 0])
  std1 <- wmean(m1, w)
  std0 <- wmean(m0, w)
  ipw1 <- wmean(d$Y[d$A == 1], (w / ps)[d$A == 1])
  ipw0 <- wmean(d$Y[d$A == 0], (w / (1 - ps))[d$A == 0])
  aipw1 <- wmean(d$A * (d$Y - m1) / ps + m1, w)
  aipw0 <- wmean((1 - d$A) * (d$Y - m0) / (1 - ps) + m0, w)
  c(
    crude_rd = crude1 - crude0, crude_rr = crude1 / crude0,
    standardization_rd = std1 - std0, standardization_rr = std1 / std0,
    ipw_rd = ipw1 - ipw0, ipw_rr = ipw1 / ipw0,
    aipw_rd = aipw1 - aipw0, aipw_rr = aipw1 / aipw0
  )
}

summarize_rep_estimates <- function(design, ps_weighted) {
  fit <- survey::withReplicates(
    design,
    function(w, data) fit_survey_estimators(w, data, ps_weighted),
    return.replicates = TRUE
  )
  estimate <- coef(fit)
  se <- SE(fit)
  data.frame(
    estimator = names(estimate),
    estimate = unname(estimate),
    se = unname(se),
    lower = unname(estimate - 1.96 * se),
    upper = unname(estimate + 1.96 * se)
  )
}

run_knhanes_ps <- function(data) {
  cc <- complete.cases(data[, covariates, drop = FALSE])
  cc_data <- data[cc, , drop = FALSE]
  full_design <- svydesign(
    ids = ~psu, strata = ~kstrata, weights = ~wt_itvex,
    nest = TRUE, data = cc_data
  )
  rep_design <- as.svrepdesign(full_design, type = "JKn")
  weighted <- summarize_rep_estimates(rep_design, TRUE)
  unweighted <- summarize_rep_estimates(rep_design, FALSE)
  weighted$ps_model <- "survey-weighted"
  unweighted$ps_model <- "unweighted"
  estimates <- rbind(weighted, unweighted)

  ps_unweighted <- suppressWarnings(
    glm(ps_formula, quasibinomial(), data = cc_data)
  )
  ps <- predict(ps_unweighted, type = "response")
  ipw <- ifelse(cc_data$A == 1, 1 / ps, 1 / (1 - ps))
  combined <- cc_data$wt_itvex * ipw
  ess <- function(w) sum(w)^2 / sum(w^2)
  diagnostics <- data.frame(
    ps_model = "unweighted",
    n = nrow(cc_data),
    max_ps = max(ps),
    min_ps = min(ps),
    p99_ipw = unname(quantile(ipw, 0.99)),
    max_ipw = max(ipw),
    p99_combined_weight = unname(quantile(combined, 0.99)),
    max_combined_weight = max(combined),
    ess_survey_weight = ess(cc_data$wt_itvex),
    ess_combined_weight = ess(combined)
  )
  write.csv(estimates, file.path(output_dir, "KNHANES_PS_models.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(output_dir, "KNHANES_unweighted_PS_diagnostics.csv"),
            row.names = FALSE)
  estimates
}

## ---------------------- S1 finite-population truth diagnostic
run_s1_diagnostic <- function() {
  beta_a <- 0
  x <- seq(-3, 3, length.out = 1001)
  p1 <- plogis(-1 + 1 * (0.5 * x) + beta_a)
  p0 <- plogis(-1 + 1 * (0.5 * x))
  structural_rd <- p1 - p0
  if (any(structural_rd != 0)) stop("S1 beta_A=0 did not imply p1=p0.")

  existing_file <- file.path(project_dir, "sim", "sim_results.csv")
  summary <- NULL
  if (file.exists(existing_file)) {
    old <- read.csv(existing_file)
    old <- old[old$scenario == "S1", , drop = FALSE]
    if (nrow(old)) {
      summary <- aggregate(
        cbind(est, se) ~ estimator,
        data = old,
        FUN = mean
      )
      names(summary)[names(summary) == "est"] <- "mean_estimate"
      names(summary)[names(summary) == "se"] <- "mean_jkn_se"
      summary$bias_vs_exact_truth <- summary$mean_estimate
      summary$truth <- 0
    }
  }
  write.csv(data.frame(
    beta_A = beta_a,
    finite_population_truth = 0,
    superpopulation_truth = 0,
    conclusion = "The S1 offset is estimator/sampling bias, not finite-population truth draw error."
  ), file.path(output_dir, "S1_truth_diagnostic.csv"), row.names = FALSE)
  if (!is.null(summary)) {
    write.csv(summary, file.path(output_dir, "S1_existing_simulation_summary.csv"),
              row.names = FALSE)
  }
  summary
}

## ---------------------- multiple imputation and Rubin pooling
run_mice <- function(data) {
  imp_data <- data[, c(
    "agec", "agec2", "A", "Y", "sex", "income", "educ4",
    "smoke", "alcohol", "walk", "wt_itvex", "kstrata"
  )]
  imp_data$log_wt <- log(imp_data$wt_itvex)
  init <- mice::mice(imp_data, maxit = 0, printFlag = FALSE)
  method <- init$method
  predictor <- init$predictorMatrix
  method[] <- ""
  method[covariates] <- "polyreg"
  predictor[, "wt_itvex"] <- 0
  predictor[, "kstrata"] <- 1
  predictor[, "log_wt"] <- 1
  predictor["kstrata", ] <- 0
  predictor["log_wt", ] <- 0
  predictor["wt_itvex", ] <- 0
  imp <- mice::mice(
    imp_data,
    m = M_IMP,
    maxit = MAXIT,
    method = method,
    predictorMatrix = predictor,
    seed = SEED + 5000L,
    printFlag = FALSE
  )
  if (!is.null(imp$loggedEvents) && nrow(imp$loggedEvents)) {
    write.csv(imp$loggedEvents, file.path(output_dir, "MICE_logged_events.csv"),
              row.names = FALSE)
  }

  m <- imp$m
  q <- list()
  u <- list()
  design_df <- NULL
  for (k in seq_len(m)) {
    completed <- mice::complete(imp, k)
    completed$psu <- data$psu
    completed$sex <- factor(completed$sex, levels = levels(data$sex))
    completed$income <- factor(completed$income, levels = levels(data$income))
    completed$educ4 <- factor(completed$educ4, levels = levels(data$educ4))
    completed$smoke <- factor(completed$smoke, levels = levels(data$smoke))
    completed$alcohol <- factor(completed$alcohol, levels = levels(data$alcohol))
    completed$walk <- factor(completed$walk, levels = levels(data$walk))
    completed$kstrata <- factor(completed$kstrata, levels = levels(data$kstrata))
    des <- svydesign(
      ids = ~psu, strata = ~kstrata, weights = ~wt_itvex,
      nest = TRUE, data = completed
    )
    rep_des <- as.svrepdesign(des, type = "JKn")
    stat <- survey::withReplicates(
      rep_des,
      function(w, d) fit_survey_estimators(w, d, TRUE),
      return.replicates = TRUE
    )
    if (is.null(design_df)) design_df <- survey::degf(rep_des)
    q[[k]] <- coef(stat)
    u[[k]] <- diag(vcov(stat))
    names(u[[k]]) <- names(q[[k]])
  }

  measures <- names(q[[1L]])
  pooled <- do.call(rbind, lapply(measures, function(measure) {
    q_i <- vapply(q, `[[`, numeric(1), measure)
    u_i <- vapply(u, `[[`, numeric(1), measure)
    pooled <- mice::pool.scalar(
      Q = q_i, U = u_i, n = design_df + 1, k = 1
    )
    crit <- qt(0.975, df = pooled$df)
    data.frame(
      measure = measure,
      m = m,
      estimate = pooled$qbar,
      within_variance = pooled$ubar,
      between_variance = pooled$b,
      total_se = sqrt(pooled$t),
      df = pooled$df,
      lower = pooled$qbar - crit * sqrt(pooled$t),
      upper = pooled$qbar + crit * sqrt(pooled$t)
    )
  }))
  write.csv(pooled, file.path(output_dir, "KNHANES_MICE_pooled.csv"), row.names = FALSE)
  write.csv(data.frame(m = m, maxit = MAXIT, design_df = design_df,
                       eligible_n = nrow(data)),
            file.path(output_dir, "KNHANES_MICE_run_info.csv"), row.names = FALSE)
  pooled
}

if (RUN_KNHANES) {
  cat("Preparing KNHANES data...\n")
  knhanes <- make_knhanes_data()
  cc <- complete.cases(knhanes[, covariates, drop = FALSE])
  cat("Eligible n =", nrow(knhanes), "; complete-case n =", sum(cc), "\n")
  if (nrow(knhanes) != 5819L || sum(cc) != 5314L) {
    stop("Sample counts do not match the documented KNHANES flow.")
  }

  cat("Running S1 truth diagnostic...\n")
  print(run_s1_diagnostic())

  cat("Running unweighted-PS comparison...\n")
  print(run_knhanes_ps(knhanes))

  cat("Running multiple imputation...\n")
  print(run_mice(knhanes))
}

if (RUN_SIMULATIONS) {
  scenario_specs <- list(
    S3a = list(make = make_s3a_population, seed = SEED + 3000L, offset = 100000L),
    S3b = list(make = make_s3b_population, seed = SEED + 4000L, offset = 200000L),
    S6 = list(make = make_s6_population, seed = SEED + 6000L, offset = 300000L)
  )
  for (scenario in SCENARIOS) {
    if (!scenario %in% names(scenario_specs)) {
      stop("Unknown simulation scenario: ", scenario)
    }
    spec <- scenario_specs[[scenario]]
    population <- spec$make(spec$seed)
    cat("Running ", scenario, " simulation with ", N_SIM,
        " replicates using ", N_CORES, " cores...\n", sep = "")
    print(run_simulation_scenario(scenario, population, spec$offset))
  }
}

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))
cat("All requested analyses completed. Results: ", output_dir, "\n", sep = "")
