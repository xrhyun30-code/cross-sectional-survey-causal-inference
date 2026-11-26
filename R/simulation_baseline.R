## Monte Carlo simulation: survey-weighted crude, G-computation, IPW, AIPW
## Estimand: population marginal risk difference (RD) for binary exposure A, binary outcome Y.
## Truth = finite-population mean(p1 - p0) from the structural outcome model.
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
R_REP <- if (length(args) >= 1) as.integer(args[1]) else 1000L
OUT   <- if (length(args) >= 2) args[2] else "sim_results.csv"
SEED  <- 20261112L

expit <- plogis

## scenario parameters
scen <- data.frame(
  id   = paste0("S", 1:6),
  label = c("Null effect, strong measured confounding",
            "Non-null effect, moderate confounding",
            "Non-linear confounding (misspecified working models)",
            "Limited positivity (strong exposure model)",
            "Unmeasured confounding",
            "Reverse causation (prior outcome state drives exposure)"),
  beta_A  = c(0, 0.5, 0.5, 0.5, 0.5, 0.5),   # log-odds effect of A on Y
  gamma   = c(1.0, 0.5, 0.5, 2.5, 0.5, 0.5), # strength of X -> A
  delta   = c(1.0, 0.5, 0.5, 1.0, 0.5, 0.5), # strength of X -> Y
  nonlin  = c(0, 0, 1, 0, 0, 0),
  unmeas  = c(0, 0, 0, 0, 1, 0),
  reverse = c(0, 0, 0, 0, 0, 1),
  stringsAsFactors = FALSE
)
scen <- scen[scen$id %in% c("S1", "S2", "S4", "S5"), , drop = FALSE]
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

N_PSU <- 400; PSU_SIZE <- 250   # population 100,000
N_STRATA <- 20; PSU_SAMPLE_PER_STRATUM <- 3   # 60 sampled PSUs

make_pop <- function(p, seed) {
  set.seed(seed)
  N <- N_PSU * PSU_SIZE
  psu <- rep(seq_len(N_PSU), each = PSU_SIZE)
  strata <- rep(seq_len(N_STRATA), each = N_PSU / N_STRATA * PSU_SIZE)
  u_psu <- rnorm(N_PSU, 0, 0.3)[psu]
  X1 <- rnorm(N) + u_psu; X2 <- rnorm(N); X3 <- rbinom(N, 1, 0.5)
  Uc <- rnorm(N)                       # unmeasured confounder
  Yprev <- rbinom(N, 1, expit(-0.5 + 0.5 * X1))  # prior outcome state
  f <- if (p$nonlin == 1) (X1^2 - 1) + 0.5 * X2 * X3 + 0.5 * X1 else 0.5 * X1 + 0.5 * X2 + 0.3 * X3
  lin_A <- -0.5 + p$gamma * (0.5 * X1 + 0.5 * X2 + 0.3 * X3) +
    p$gamma * (if (p$nonlin == 1) 0.5 * (X1^2 - 1) else 0) +
    p$unmeas * 1.0 * Uc + p$reverse * 1.5 * Yprev
  A <- rbinom(N, 1, expit(lin_A))
  base <- -1 + p$delta * f + p$unmeas * 1.0 * Uc + p$reverse * 1.5 * Yprev
  p1 <- expit(base + p$beta_A); p0 <- expit(base)
  Y <- rbinom(N, 1, ifelse(A == 1, p1, p0))
  # selection: unequal probabilities depending on X1 (PSU sampled at stratum level)
  pi_i <- pmin(0.9, 0.2 * exp(0.4 * X1)); pi_i <- pi_i / mean(pi_i) * 0.2
  list(pop = data.frame(psu, strata, X1, X2, X3, A, Y, pi_i),
       truth = mean(p1 - p0), n_ps = mean(expit(lin_A)),
       min_ps = min(expit(lin_A)), max_ps = max(expit(lin_A)),
       prev_A = mean(A), prev_Y = mean(Y))
}

draw_sample <- function(pop) {
  # stratified one-stage cluster-ish: sample PSUs within strata, then Poisson within PSU
  psu_tab <- unique(pop[, c("psu", "strata")])
  sel <- unlist(lapply(split(psu_tab$psu, psu_tab$strata),
                       function(z) sample(z, PSU_SAMPLE_PER_STRATUM)))
  s <- pop[pop$psu %in% sel, ]
  s <- s[runif(nrow(s)) < s$pi_i, ]
  pi_psu <- PSU_SAMPLE_PER_STRATUM / (N_PSU / N_STRATA)
  s$w <- 1 / (pi_psu * s$pi_i)
  s
}

wmean <- function(x, w) sum(w * x) / sum(w)

## point estimators (all weighted by survey weight)
est_all <- function(d) {
  d$w <- d$w / mean(d$w)   # estimators are scale-invariant; glm convergence is not
  w <- d$w
  # crude
  crude <- wmean(d$Y[d$A == 1], w[d$A == 1]) - wmean(d$Y[d$A == 0], w[d$A == 0])
  # exposure model (survey weighted)
  ps_fit <- suppressWarnings(glm(A ~ X1 + X2 + X3, quasibinomial, d, weights = w))
  ps <- predict(ps_fit, type = "response")
  # outcome model
  y_fit <- suppressWarnings(glm(Y ~ A + X1 + X2 + X3, quasibinomial, d, weights = w))
  d1 <- transform(d, A = 1); d0 <- transform(d, A = 0)
  m1 <- predict(y_fit, d1, type = "response"); m0 <- predict(y_fit, d0, type = "response")
  gcomp <- wmean(m1 - m0, w)
  # IPW (Hajek), untruncated; combined weight = survey weight x IPW
  ipw <- wmean(d$Y[d$A == 1], (w / ps)[d$A == 1]) - wmean(d$Y[d$A == 0], (w / (1 - ps))[d$A == 0])
  # AIPW
  psi <- (d$A * (d$Y - m1) / ps + m1) - ((1 - d$A) * (d$Y - m0) / (1 - ps) + m0)
  aipw <- wmean(psi, w)
  c(crude = crude, gcomp = gcomp, ipw = ipw, aipw = aipw)
}

## delete-one-PSU jackknife within strata (JKn), ignoring fpc (conservative)
jk_se <- function(s) {
  full <- est_all(s)
  psus <- unique(s[, c("psu", "strata")])
  v <- rep(0, length(full))
  for (h in unique(psus$strata)) {
    ph <- psus$psu[psus$strata == h]; nh <- length(ph)
    for (j in ph) {
      sj <- s[!(s$psu == j), ]
      # reweight remaining PSUs in stratum h
      idx <- sj$strata == h
      sj$w[idx] <- sj$w[idx] * nh / (nh - 1)
      e <- est_all(sj)
      v <- v + (nh - 1) / nh * (e - full)^2
    }
  }
  list(est = full, se = sqrt(v))
}

one_rep <- function(i, pop, truth) {
  set.seed(SEED + i)
  s <- draw_sample(pop)
  out <- tryCatch(jk_se(s), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  data.frame(rep = i, estimator = names(out$est), est = out$est, se = out$se,
             truth = truth, n = nrow(s),
             ess = sum(s$w)^2 / sum(s$w^2))
}

res <- list(); info <- list()
for (k in seq_len(nrow(scen))) {
  p <- scen[k, ]
  scenario_index <- match(p$id, paste0("S", 1:6))
  P <- make_pop(p, SEED + 1000 * scenario_index)
  info[[k]] <- data.frame(id = p$id, truth = P$truth, prev_A = P$prev_A, prev_Y = P$prev_Y,
                          min_ps = P$min_ps, max_ps = P$max_ps)
  cat("Scenario", p$id, "truth RD =", round(P$truth, 4), "\n"); flush.console()
  r <- mclapply(seq_len(R_REP), one_rep, pop = P$pop, truth = P$truth,
                mc.cores = 7, mc.set.seed = FALSE)
  r <- do.call(rbind, r); r$scenario <- p$id
  res[[k]] <- r
}
res <- do.call(rbind, res)
write.csv(res, OUT, row.names = FALSE)
write.csv(do.call(rbind, info), sub("\\.csv$", "_popinfo.csv", OUT), row.names = FALSE)
write.csv(scen, sub("\\.csv$", "_scenarios.csv", OUT), row.names = FALSE)
cat("done\n")
