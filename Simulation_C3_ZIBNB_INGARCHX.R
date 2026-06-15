# ------------------------------
# Load Packages and Set Seed
# ------------------------------
library(mvtnorm)
library(MASS)
library(coda)
library(foreach)
library(doParallel)
set.seed(7)

# ==============================================================================
# SECTION 1: BNB Distribution Utilities
# ==============================================================================

log_dbnb <- function(y, r, gamma_t, phi) {
  lgamma(y + r) - lgamma(y + 1) - lgamma(r) +
    lbeta(phi + r, gamma_t + y) - lbeta(phi, gamma_t)
}

log_fbnb_zero <- function(r, gamma_t, phi) {
  lbeta(phi + r, gamma_t) - lbeta(phi, gamma_t)
}

rbnb <- function(n, r, gamma_t, phi) {
  p <- rbeta(n, phi, gamma_t)
  y <- rnbinom(n, size = r, prob = p)
  return(y)
}

# ==============================================================================
# SECTION 2: Beta Constraint
# ==============================================================================
check_beta_constraints <- function(beta_prop) {
  return(beta_prop[1] > 0 &&
           beta_prop[2] > 0 &&
           beta_prop[3] >= 0 &&
           (beta_prop[2] + beta_prop[3]) < 1)
}

# ==============================================================================
# SECTION 3: Data Generation — k = 2 exogenous covariates
# ==============================================================================
ZIBNB_data <- function(n, rho, beta, omega, r, phi, X_list, b, b0,
                       lambda_init, k) {
  # omega : numeric vector of length k
  # b     : integer vector of length k  (lag indices)
  y        <- rep(0, n)
  lambda_t <- rep(lambda_init, n)
  
  gamma_1 <- (phi - 1) / r * lambda_init
  if (runif(1) < rho) {
    y[1] <- 0L
  } else {
    y[1] <- rbnb(1, r, gamma_1, phi)
  }
  
  for (t in (b0 + 1):n) {
    lambda_t[t] <- beta[1] + beta[2] * y[t - 1] + beta[3] * lambda_t[t - 1]
    if (k >= 1) {
      for (i in 1:k) {
        lambda_t[t] <- lambda_t[t] + omega[i] * X_list[[i]][1, t - b[i]]
      }
    }
    lambda_t[t] <- max(lambda_t[t], 1e-6)
    gamma_t <- (phi - 1) / r * lambda_t[t]
    
    if (runif(1) < rho) {
      y[t] <- 0L
    } else {
      y[t] <- rbnb(1, r, gamma_t, phi)
    }
  }
  
  return(list(y = y, lambda = lambda_t))
}

# ==============================================================================
# SECTION 4: Log-Likelihood — k = 2
# ==============================================================================
log_likelihood <- function(Y, rho, beta, omega, r, phi, X_list, b, b0, k) {
  # omega : numeric vector of length k
  # b     : integer vector of length k
  n        <- length(Y)
  log_like <- 0
  lambda_t <- rep(0.1, n)
  
  for (t in (b0 + 1):n) {
    lambda_t[t] <- beta[1] + beta[2] * Y[t - 1] + beta[3] * lambda_t[t - 1]
    if (k >= 1) {
      for (i in 1:k) {
        lambda_t[t] <- lambda_t[t] + omega[i] * X_list[[i]][1, t - b[i]]
      }
    }
    lambda_t[t] <- max(lambda_t[t], 1e-6)
    gamma_t <- (phi - 1) / r * lambda_t[t]
    
    log_fbnb_y <- lgamma(Y[t] + r) - lgamma(Y[t] + 1) - lgamma(r) +
      lbeta(phi + r, gamma_t + Y[t]) - lbeta(phi, gamma_t)
    log_fbnb_0 <- lbeta(phi + r, gamma_t) - lbeta(phi, gamma_t)
    fbnb_0     <- exp(log_fbnb_0)
    
    if (Y[t] == 0) {
      log_like <- log_like + log(rho + (1 - rho) * fbnb_0)
    } else {
      log_like <- log_like + log(1 - rho) + log_fbnb_y
    }
  }
  
  return(log_like)
}

# ==============================================================================
# SECTION 5: MCMC — k = 2
#
# Two-phase Metropolis-Hastings:
#   Phase 1  (iterations 1 … burn_in)          : Random-Walk MH
#   Phase 2  (iterations burn_in+1 … N_total)  : Independent-Kernel MH
#     using the empirical moments of Phase 1 as the proposal distribution.
#
# MH acceptance rule (both phases):
#   log_accept_ratio = r_prop_lposterior - r_curr_lposterior
#
#   Phase 1 (symmetric RW proposal — proposal terms cancel):
#     r_prop_lposterior = lp_*(theta')
#     r_curr_lposterior = lp_*(theta)
#
#   Phase 2 (independent proposal q):
#     r_prop_lposterior = lp_*(theta') + log q(theta)   [lq_rev]
#     r_curr_lposterior = lp_*(theta)  + log q(theta')  [lq_fwd]
#
# Parameter vector (10 columns):
#   [rho, beta1, beta2, beta3, omega1, omega2, r, phi, b1, b2]
# ==============================================================================
run_mcmc <- function(Y, X_list, k = 2,
                     N_total   = 40000,
                     burn_in   = 16000,
                     thin      = 4,
                     b0_mcmc   = 3,
                     prior_hyp = list(e1 = 1, e2 = 1,
                                      c1 = 1, c2 = 1,
                                      a1 = 3, a2 = 1,
                                      d1 = 3, d2 = 1))
{
  e1 <- prior_hyp$e1; e2 <- prior_hyp$e2   # rho    ~ Beta(e1, e2)
  c1 <- prior_hyp$c1; c2 <- prior_hyp$c2   # omega  ~ Gamma(c1, rate=c2)
  a1 <- prior_hyp$a1; a2 <- prior_hyp$a2   # r      ~ Gamma(a1, rate=a2)
  d1 <- prior_hyp$d1; d2 <- prior_hyp$d2   # phi    ~ Gamma-like (d1,d2), phi>2
  
  N_phase2 <- N_total - burn_in
  
  # ---- Random-walk stepfburn sizes (Phase 1) ----
  step_rho   <- 0.05
  step_beta  <- c(0.12, 0.12, 0.12)
  step_omega <- c(0.25, 0.25)   # one per covariate
  step_r     <- 3.0
  step_phi   <- 2.0
  
  # ---- Initial values ----
  rho_curr   <- 0.3
  beta_curr  <- c(0.1, 0.1, 0.1)
  omega_curr <- c(0.1, 0.1)    # length-2
  r_curr     <- 2
  phi_curr   <- 3
  b          <- c(1, 1)        # length-2 lag vector
  
  # Column layout (10 parameters):
  #  1      : rho
  #  2-4    : beta1, beta2, beta3
  #  5-6    : omega1, omega2
  #  7      : r
  #  8      : phi
  #  9-10   : b1, b2
  n_param   <- 10
  idx_rho   <- 1
  idx_beta  <- 2:4
  idx_omega <- 5:6
  idx_r     <- 7
  idx_phi   <- 8
  idx_b     <- 9:10
  
  all_samples     <- matrix(0, nrow = N_total, ncol = n_param)
  burn_in_samples <- matrix(0, nrow = burn_in,  ncol = n_param)
  
  pack_params <- function()
    c(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
  
  # ============================================================
  # Log-posterior functions
  # lp_*(focal_param | rest, Y) = log-likelihood + log-prior
  # ============================================================
  
  lp_rho <- function(rho, beta, omega, r, phi, b) {
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k) +
      dbeta(rho, e1, e2, log = TRUE)
  }
  
  lp_beta <- function(rho, beta, omega, r, phi, b) {
    # Flat prior on beta (support enforced by proposal constraint)
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k)
  }
  
  lp_omega1 <- function(rho, beta, omega, r, phi, b) {
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k) +
      dgamma(omega[1], c1, rate = c2, log = TRUE)
  }
  
  lp_omega2 <- function(rho, beta, omega, r, phi, b) {
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k) +
      dgamma(omega[2], c1, rate = c2, log = TRUE)
  }
  
  lp_r <- function(rho, beta, omega, r, phi, b) {
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k) +
      dgamma(r, a1, rate = a2, log = TRUE)
  }
  
  lp_phi <- function(rho, beta, omega, r, phi, b) {
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi) - d2 * phi
  }
  
  # ============================================================
  # MH accept / reject helper
  # ============================================================
  mh_accept <- function(log_accept_ratio) {
    if (!is.finite(log_accept_ratio)) return(FALSE)
    log(runif(1)) < min(0, log_accept_ratio)
  }
  
  # ---- b1 and b2: exact discrete samplers ----
  # Each b_i is sampled from its full conditional with the other fixed.
  sample_b1 <- function() {
    log_lik_b1 <- sapply(1:b0_mcmc, function(j) {
      log_likelihood(Y, rho_curr, beta_curr, omega_curr,
                     r_curr, phi_curr, X_list, c(j, b[2]), b0_mcmc, k)
    })
    w <- exp(log_lik_b1 - max(log_lik_b1))
    sample(1:b0_mcmc, 1, prob = w / sum(w))
  }
  
  sample_b2 <- function() {
    log_lik_b2 <- sapply(1:b0_mcmc, function(j) {
      log_likelihood(Y, rho_curr, beta_curr, omega_curr,
                     r_curr, phi_curr, X_list, c(b[1], j), b0_mcmc, k)
    })
    w <- exp(log_lik_b2 - max(log_lik_b2))
    sample(1:b0_mcmc, 1, prob = w / sum(w))
  }
  
  acc_rho_p1   <- 0
  acc_beta_p1  <- 0
  acc_omega_p1 <- c(0, 0)
  acc_r_p1     <- 0
  acc_phi_p1   <- 0
  
  # ============================================================
  # PHASE 1: Random-Walk MH  (iterations 1 … burn_in)
  # ============================================================
  for (iter in 1:burn_in) {
    if (iter %% 1000 == 0) cat("  Burn-in Iteration:", iter, "\n")
    
    # ---------- rho ----------
    repeat {
      rho_prop <- rho_curr + rnorm(1, 0, step_rho)
      if (rho_prop > 0 && rho_prop < 1) break
    }
    log_ar <- lp_rho(rho_prop, beta_curr, omega_curr, r_curr, phi_curr, b) -
      lp_rho(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    if (mh_accept(log_ar)) { rho_curr <- rho_prop; acc_rho_p1 <- acc_rho_p1 + 1 }
    
    # ---------- beta ----------
    repeat {
      beta_prop <- beta_curr + rnorm(3, 0, step_beta)
      if (check_beta_constraints(beta_prop)) break
    }
    log_ar <- lp_beta(rho_curr, beta_prop, omega_curr, r_curr, phi_curr, b) -
      lp_beta(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    if (mh_accept(log_ar)) { beta_curr <- beta_prop; acc_beta_p1 <- acc_beta_p1 + 1 }
    
    # ---------- omega1 ----------
    repeat {
      omega_prop1 <- omega_curr[1] + rnorm(1, 0, step_omega[1])
      if (omega_prop1 > 0) break
    }
    omega_test1 <- c(omega_prop1, omega_curr[2])
    log_ar <- lp_omega1(rho_curr, beta_curr, omega_test1, r_curr, phi_curr, b) -
      lp_omega1(rho_curr, beta_curr, omega_curr,  r_curr, phi_curr, b)
    if (mh_accept(log_ar)) {
      omega_curr[1] <- omega_prop1; acc_omega_p1[1] <- acc_omega_p1[1] + 1
    }
    
    # ---------- omega2 ----------
    repeat {
      omega_prop2 <- omega_curr[2] + rnorm(1, 0, step_omega[2])
      if (omega_prop2 > 0) break
    }
    omega_test2 <- c(omega_curr[1], omega_prop2)
    log_ar <- lp_omega2(rho_curr, beta_curr, omega_test2, r_curr, phi_curr, b) -
      lp_omega2(rho_curr, beta_curr, omega_curr,  r_curr, phi_curr, b)
    if (mh_accept(log_ar)) {
      omega_curr[2] <- omega_prop2; acc_omega_p1[2] <- acc_omega_p1[2] + 1
    }
    
    # ---------- r ----------
    repeat {
      r_prop <- r_curr + rnorm(1, 0, step_r)
      if (r_prop > 0) break
    }
    log_ar <- lp_r(rho_curr, beta_curr, omega_curr, r_prop, phi_curr, b) -
      lp_r(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    if (mh_accept(log_ar)) { r_curr <- r_prop; acc_r_p1 <- acc_r_p1 + 1 }
    
    # ---------- phi (phi > 2) ----------
    repeat {
      phi_prop <- phi_curr + rnorm(1, 0, step_phi)
      if (phi_prop > 2) break
    }
    log_ar <- lp_phi(rho_curr, beta_curr, omega_curr, r_curr, phi_prop, b) -
      lp_phi(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    if (mh_accept(log_ar)) { phi_curr <- phi_prop; acc_phi_p1 <- acc_phi_p1 + 1 }
    
    # ---------- b1, b2 (exact discrete sampling) ----------
    b[1] <- sample_b1()
    b[2] <- sample_b2()
    
    burn_in_samples[iter, ] <- pack_params()
    all_samples[iter, ]     <- pack_params()
  }
  
  cat("\n  Phase 1 Acceptance Rates:\n")
  cat("    rho:   ", round(acc_rho_p1      / burn_in, 4), "\n")
  cat("    beta:  ", round(acc_beta_p1     / burn_in, 4), "\n")
  cat("    omega1:", round(acc_omega_p1[1] / burn_in, 4), "\n")
  cat("    omega2:", round(acc_omega_p1[2] / burn_in, 4), "\n")
  cat("    r:     ", round(acc_r_p1        / burn_in, 4), "\n")
  cat("    phi:   ", round(acc_phi_p1      / burn_in, 4), "\n")
  
  # ============================================================
  # Empirical proposal moments from Phase 1 (drop first 1,000)
  # ============================================================
  warmup <- 1:1000
  
  mu_rho    <- mean(burn_in_samples[-warmup, idx_rho])
  sd_rho    <- sd(  burn_in_samples[-warmup, idx_rho])
  
  mu_beta   <- colMeans(burn_in_samples[-warmup, idx_beta])
  cov_beta  <- cov(    burn_in_samples[-warmup, idx_beta])
  
  mu_omega  <- colMeans(burn_in_samples[-warmup, idx_omega, drop = FALSE])
  sd_omega1 <- sd(burn_in_samples[-warmup, idx_omega[1]])
  sd_omega2 <- sd(burn_in_samples[-warmup, idx_omega[2]])
  
  mu_r      <- mean(burn_in_samples[-warmup, idx_r])
  sd_r      <- sd(  burn_in_samples[-warmup, idx_r])
  
  mu_phi    <- mean(burn_in_samples[-warmup, idx_phi])
  sd_phi    <- sd(  burn_in_samples[-warmup, idx_phi])
  
  acc_rho_p2   <- 0
  acc_beta_p2  <- 0
  acc_omega_p2 <- c(0, 0)
  acc_r_p2     <- 0
  acc_phi_p2   <- 0
  
  # ============================================================
  # PHASE 2: Independent-Kernel MH  (iterations burn_in+1 … N_total)
  # ============================================================
  for (iter in 1:N_phase2) {
    if (iter %% 1000 == 0) cat("  Independent-Kernel Iteration:", iter, "\n")
    
    # ---------- rho ----------
    repeat {
      rho_prop <- rnorm(1, mu_rho, sd_rho)
      if (rho_prop > 0 && rho_prop < 1) break
    }
    r_prop_lpost <- lp_rho(rho_prop, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(rho_curr, mu_rho, sd_rho, log = TRUE)    # lq_rev
    r_curr_lpost <- lp_rho(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(rho_prop, mu_rho, sd_rho, log = TRUE)    # lq_fwd
    if (mh_accept(r_prop_lpost - r_curr_lpost)) {
      rho_curr <- rho_prop; acc_rho_p2 <- acc_rho_p2 + 1
    }
    
    # ---------- beta ----------
    repeat {
      beta_prop <- as.numeric(mvrnorm(1, mu_beta, cov_beta))
      if (check_beta_constraints(beta_prop)) break
    }
    r_prop_lpost <- lp_beta(rho_curr, beta_prop, omega_curr, r_curr, phi_curr, b) +
      mvtnorm::dmvnorm(beta_curr, mu_beta, cov_beta, log = TRUE)  # lq_rev
    r_curr_lpost <- lp_beta(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      mvtnorm::dmvnorm(beta_prop, mu_beta, cov_beta, log = TRUE)  # lq_fwd
    if (mh_accept(r_prop_lpost - r_curr_lpost)) {
      beta_curr <- beta_prop; acc_beta_p2 <- acc_beta_p2 + 1
    }
    
    # ---------- omega1 ----------
    repeat {
      omega_prop1 <- rnorm(1, mu_omega[1], sd_omega1)
      if (omega_prop1 > 0) break
    }
    omega_test1  <- c(omega_prop1, omega_curr[2])
    r_prop_lpost <- lp_omega1(rho_curr, beta_curr, omega_test1, r_curr, phi_curr, b) +
      dnorm(omega_curr[1], mu_omega[1], sd_omega1, log = TRUE)    # lq_rev
    r_curr_lpost <- lp_omega1(rho_curr, beta_curr, omega_curr,  r_curr, phi_curr, b) +
      dnorm(omega_prop1,   mu_omega[1], sd_omega1, log = TRUE)    # lq_fwd
    if (mh_accept(r_prop_lpost - r_curr_lpost)) {
      omega_curr[1] <- omega_prop1; acc_omega_p2[1] <- acc_omega_p2[1] + 1
    }
    
    # ---------- omega2 ----------
    repeat {
      omega_prop2 <- rnorm(1, mu_omega[2], sd_omega2)
      if (omega_prop2 > 0) break
    }
    omega_test2  <- c(omega_curr[1], omega_prop2)
    r_prop_lpost <- lp_omega2(rho_curr, beta_curr, omega_test2, r_curr, phi_curr, b) +
      dnorm(omega_curr[2], mu_omega[2], sd_omega2, log = TRUE)    # lq_rev
    r_curr_lpost <- lp_omega2(rho_curr, beta_curr, omega_curr,  r_curr, phi_curr, b) +
      dnorm(omega_prop2,   mu_omega[2], sd_omega2, log = TRUE)    # lq_fwd
    if (mh_accept(r_prop_lpost - r_curr_lpost)) {
      omega_curr[2] <- omega_prop2; acc_omega_p2[2] <- acc_omega_p2[2] + 1
    }
    
    # ---------- r ----------
    repeat {
      r_prop <- rnorm(1, mu_r, sd_r)
      if (r_prop > 0) break
    }
    r_prop_lpost <- lp_r(rho_curr, beta_curr, omega_curr, r_prop, phi_curr, b) +
      dnorm(r_curr, mu_r, sd_r, log = TRUE)    # lq_rev
    r_curr_lpost <- lp_r(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(r_prop, mu_r, sd_r, log = TRUE)    # lq_fwd
    if (mh_accept(r_prop_lpost - r_curr_lpost)) {
      r_curr <- r_prop; acc_r_p2 <- acc_r_p2 + 1
    }
    
    # ---------- phi (phi > 2) ----------
    repeat {
      phi_prop <- rnorm(1, mu_phi, sd_phi)
      if (phi_prop > 2) break
    }
    r_prop_lpost <- lp_phi(rho_curr, beta_curr, omega_curr, r_curr, phi_prop, b) +
      dnorm(phi_curr, mu_phi, sd_phi, log = TRUE)    # lq_rev
    r_curr_lpost <- lp_phi(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(phi_prop, mu_phi, sd_phi, log = TRUE)    # lq_fwd
    if (mh_accept(r_prop_lpost - r_curr_lpost)) {
      phi_curr <- phi_prop; acc_phi_p2 <- acc_phi_p2 + 1
    }
    
    # ---------- b1, b2 (exact discrete sampling) ----------
    b[1] <- sample_b1()
    b[2] <- sample_b2()
    
    all_samples[burn_in + iter, ] <- pack_params()
  }
  
  cat("\n  Phase 2 Acceptance Rates:\n")
  cat("    rho:   ", round(acc_rho_p2      / N_phase2, 4), "\n")
  cat("    beta:  ", round(acc_beta_p2     / N_phase2, 4), "\n")
  cat("    omega1:", round(acc_omega_p2[1] / N_phase2, 4), "\n")
  cat("    omega2:", round(acc_omega_p2[2] / N_phase2, 4), "\n")
  cat("    r:     ", round(acc_r_p2        / N_phase2, 4), "\n")
  cat("    phi:   ", round(acc_phi_p2      / N_phase2, 4), "\n")
  
  # ============================================================
  # THINNING
  #   1. Extract every `thin`-th row from the full N_total chain.
  #   2. Discard the first (burn_in / thin) thinned rows.
  #   3. Retain the remainder as posterior samples.
  # ============================================================
  thin_idx          <- seq(thin, N_total, by = thin)
  thinned_all       <- all_samples[thin_idx, , drop = FALSE]
  n_burnin_thinned  <- burn_in / thin
  posterior_samples <- thinned_all[(n_burnin_thinned + 1):nrow(thinned_all), ]
  
  cat("\n  Thinning summary:\n")
  cat("    Total raw iterates     :", N_total, "\n")
  cat("    Thinning interval      :", thin, "\n")
  cat("    Thinned rows (total)   :", length(thin_idx), "\n")
  cat("    Burn-in thinned rows   :", n_burnin_thinned, "\n")
  cat("    Posterior samples kept :", nrow(posterior_samples), "\n")
  
  return(list(
    samples       = posterior_samples,
    n_param       = n_param,
    k             = k,
    idx_rho       = idx_rho,
    idx_beta      = idx_beta,
    idx_omega     = idx_omega,
    idx_r         = idx_r,
    idx_phi       = idx_phi,
    idx_b         = idx_b,
    accept_rho    = acc_rho_p2      / N_phase2,
    accept_beta   = acc_beta_p2     / N_phase2,
    accept_omega1 = acc_omega_p2[1] / N_phase2,
    accept_omega2 = acc_omega_p2[2] / N_phase2,
    accept_r      = acc_r_p2        / N_phase2,
    accept_phi    = acc_phi_p2      / N_phase2
  ))
}

# ==============================================================================
# SECTION 6: Summary Table Helper
# ==============================================================================
make_summary_table <- function(post_means, param_names, true_vals) {
  n_p      <- length(param_names)
  means_mx <- post_means[, 1:n_p, drop = FALSE]
  
  col_mean   <- round(colMeans(means_mx), 4)
  col_median <- round(apply(means_mx, 2, median), 4)
  col_sd     <- round(apply(means_mx, 2, sd), 4)
  col_p025   <- round(apply(means_mx, 2, quantile, probs = 0.025), 4)
  col_p975   <- round(apply(means_mx, 2, quantile, probs = 0.975), 4)
  col_bias   <- round(col_mean - true_vals, 4)
  col_rmse   <- round(sqrt(colMeans(
    (means_mx - matrix(true_vals, nrow = nrow(means_mx),
                       ncol = n_p, byrow = TRUE))^2)), 4)
  
  data.frame(
    Parameter = param_names,
    True      = true_vals,
    Mean      = col_mean,
    Median    = col_median,
    Std       = col_sd,
    P0.025    = col_p025,
    P0.975    = col_p975,
    Bias      = col_bias,
    RMSE      = col_rmse,
    row.names = NULL
  )
}

# ==============================================================================
# SECTION 7: Common Settings
# ==============================================================================
n             <- 1000
rho_true      <- 0.30
beta_true     <- c(0.50, 0.30, 0.15)
r_true        <- 3
phi_true      <- 4
omega_true_c3 <- c(0.12, 0.10)   # two omega true values
b_true_c3     <- c(1, 2)         # two lag true values
lambda_init   <- 0.1

# ==============================================================================
# SECTION 8: REPLICATION BLOCK — PARALLEL
# ==============================================================================
prime_seeds <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
                 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
                 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
                 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
                 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
                 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
                 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
                 467, 479, 487, 491, 499, 503, 509, 521, 523, 541)

M           <- 100
n_cont_c3   <- 8   # rho, beta1, beta2, beta3, omega1, omega2, r, phi
# (b1, b2 are discrete — reported separately via mode)

param_names_c3 <- c("rho", "beta1", "beta2", "beta3",
                    "omega1", "omega2", "r", "phi")

col_labels_c3 <- c(
  paste0("mean_",   param_names_c3),
  paste0("median_", param_names_c3),
  paste0("sd_",     param_names_c3),
  paste0("p025_",   param_names_c3),
  paste0("p975_",   param_names_c3),
  "mode_b1", "mode_b2",
  "accept_rho", "accept_beta",
  "accept_omega1", "accept_omega2",
  "accept_r", "accept_phi"
)

mode_fn <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

checkpoint_file_c3 <- "checkpoint_case3_k2_thinned.csv"

if (file.exists(checkpoint_file_c3)) {
  post_means_c3 <- as.matrix(read.csv(checkpoint_file_c3, check.names = FALSE))
  start_m_c3    <- nrow(post_means_c3) + 1
  cat("Case 3 (k=2): Resuming from replication", start_m_c3, "\n")
} else {
  post_means_c3 <- matrix(NA, nrow = 0, ncol = length(col_labels_c3))
  start_m_c3    <- 1
}

# ---- Parallel replication loop ----
if (start_m_c3 <= M) {
  
  ncores <- max(1L, parallel::detectCores() - 1L)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  cat("Running Case 3 (k=2) replications", start_m_c3, "to", M, "on", ncores, "cores...\n")
  t0 <- Sys.time()
  
  sim_results <- foreach(m = start_m_c3:M,
                         .combine = rbind,
                         .packages = c("mvtnorm", "MASS", "coda"),
                         .export = c("ZIBNB_data", "run_mcmc", "log_likelihood",
                                     "rbnb", "log_dbnb", "check_beta_constraints", "log_fbnb_zero",
                                     "prime_seeds", "n", "rho_true", "beta_true", "omega_true_c3",
                                     "r_true", "phi_true", "b_true_c3", "lambda_init", "n_cont_c3", 
                                     "param_names_c3", "mode_fn")) %dopar% {
                                       
                                       set.seed(prime_seeds[m])
                                       
                                       # Generate two covariates
                                       X1_rep <- matrix(rnorm(n, 0, 1), 1, ncol = n)
                                       X2_rep <- matrix(rnorm(n, 0, 1), 1, ncol = n)
                                       
                                       sim_rep <- ZIBNB_data(n           = n,
                                                             rho         = rho_true,
                                                             beta        = beta_true,
                                                             omega       = omega_true_c3,
                                                             r           = r_true,
                                                             phi         = phi_true,
                                                             X_list      = list(X1_rep, X2_rep),
                                                             b           = b_true_c3,
                                                             b0          = max(b_true_c3),
                                                             lambda_init = lambda_init,
                                                             k           = 2)
                                       Y_rep <- sim_rep$y
                                       
                                       res <- run_mcmc(Y       = Y_rep,
                                                       X_list  = list(X1_rep, X2_rep),
                                                       k       = 2,
                                                       N_total = 40000,
                                                       burn_in = 16000,
                                                       thin    = 4,
                                                       b0_mcmc = 3)
                                       
                                       # Traceplots and ACF for replication 1 only
                                       if (m == 1) {
                                         true_vals_plot <- c(rho_true, beta_true, omega_true_c3, r_true, phi_true)
                                         
                                         pdf("traceplot_case3_k2_thinned.pdf")
                                         par(mfrow = c(3, 3))
                                         for (i in 1:n_cont_c3) {
                                           plot(res$samples[, i], type = "l",
                                                main = param_names_c3[i],
                                                ylab = "Value", xlab = "Thinned Iteration")
                                           abline(h = true_vals_plot[i], col = "red", lwd = 2)
                                         }
                                         dev.off()
                                         
                                         pdf("acf_case3_k2_thinned.pdf")
                                         par(mfrow = c(3, 3))
                                         for (i in 1:n_cont_c3) {
                                           acf(res$samples[, i], main = param_names_c3[i])
                                         }
                                         dev.off()
                                       }
                                       
                                       mode_b1 <- mode_fn(res$samples[, res$idx_b[1]])
                                       mode_b2 <- mode_fn(res$samples[, res$idx_b[2]])
                                       
                                       new_row <- c(
                                         apply(res$samples[, 1:n_cont_c3], 2, mean),
                                         apply(res$samples[, 1:n_cont_c3], 2, median),
                                         apply(res$samples[, 1:n_cont_c3], 2, sd),
                                         apply(res$samples[, 1:n_cont_c3], 2, quantile, probs = 0.025),
                                         apply(res$samples[, 1:n_cont_c3], 2, quantile, probs = 0.975),
                                         mode_b1, mode_b2,
                                         res$accept_rho, res$accept_beta,
                                         res$accept_omega1, res$accept_omega2,
                                         res$accept_r, res$accept_phi
                                       )
                                       
                                       return(new_row)
                                     }
  
  stopCluster(cl)
  
  # Combine with any existing checkpoints and save
  post_means_c3           <- rbind(post_means_c3, sim_results)
  colnames(post_means_c3) <- col_labels_c3
  
  write.csv(post_means_c3, checkpoint_file_c3, row.names = FALSE)
  
  cat("\nAll replications complete in:",
      round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes.\n")
}

# ==============================================================================
# SECTION 9: SIMULATION SUMMARY TABLE — Case 3 (k = 2, thinned)
# ==============================================================================
true_cont_c3 <- c(rho_true, beta_true, omega_true_c3, r_true, phi_true)
result_c3    <- make_summary_table(post_means_c3, param_names_c3, true_cont_c3)

cat("\n========== Case 3 (k=2, thinned) Summary across",
    nrow(post_means_c3), "replications ==========\n")
print(result_c3)
write.csv(result_c3, "result_case3_k2_thinned.csv", row.names = FALSE)

cat("\nAverage Phase 2 Acceptance Rates:\n")
cat("  rho:   ", round(mean(post_means_c3[, "accept_rho"]),    4), "\n")
cat("  beta:  ", round(mean(post_means_c3[, "accept_beta"]),   4), "\n")
cat("  omega1:", round(mean(post_means_c3[, "accept_omega1"]), 4), "\n")
cat("  omega2:", round(mean(post_means_c3[, "accept_omega2"]), 4), "\n")
cat("  r:     ", round(mean(post_means_c3[, "accept_r"]),      4), "\n")
cat("  phi:   ", round(mean(post_means_c3[, "accept_phi"]),    4), "\n")

cat("\nMode of b1 and b2 across replications:\n")
cat("  Mode(b1):", mode_fn(post_means_c3[, "mode_b1"]),
    "  (True b1 =", b_true_c3[1], ")\n")
cat("  Mode(b2):", mode_fn(post_means_c3[, "mode_b2"]),
    "  (True b2 =", b_true_c3[2], ")\n")
