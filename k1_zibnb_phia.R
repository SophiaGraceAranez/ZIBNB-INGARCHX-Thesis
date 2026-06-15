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
# SECTION 3: Data Generation — k = 1 exogenous covariate
# ==============================================================================
ZIBNB_data <- function(n, rho, beta, omega, r, phi, X_list, b, b0,
                       lambda_init, k) {
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
# SECTION 4: Log-Likelihood — k = 1
# ==============================================================================
log_likelihood <- function(Y, rho, beta, omega, r, phi, X_list, b, b0, k) {
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
# SECTION 5: MCMC — k = 1
#
# Two-phase Metropolis–Hastings:
#   Phase 1  (iterations 1 … burn_in)     : Random-Walk MH
#   Phase 2  (iterations burn_in+1 … N_total): Independent-Kernel MH
#     using the empirical moments of Phase 1 as the proposal distribution.
#
# MH acceptance rule (both phases):
#
#   log_accept_ratio = r_prop_lposterior − r_curr_lposterior
#
#   where
#     Phase 1 (symmetric RW proposal):
#       r_prop_lposterior = lp(θ')          [proposal terms cancel]
#       r_curr_lposterior = lp(θ)
#
#     Phase 2 (independent proposal q):
#       r_prop_lposterior = lp(θ') + log q(θ  | θ')   [= lp(θ') + log q(θ)]
#       r_curr_lposterior = lp(θ)  + log q(θ' | θ)    [= lp(θ)  + log q(θ')]
#
# Parameter vector (8 columns):
#   [rho, beta1, beta2, beta3, omega1, r, phi, b1]
# ==============================================================================
run_mcmc <- function(Y, X_list, k = 1,
                     N_total   = 40000,
                     burn_in   = 16000,
                     thin      = 4,
                     b0_mcmc   = 3,
                     prior_hyp = list(e1 = 1, e2 = 1,
                                      c1 = 1, c2 = 1,
                                      a1 = 3, a2 = 1,
                                      d1 = 3, d2 = 1))
{
  e1 <- prior_hyp$e1; e2 <- prior_hyp$e2   # rho   ~ Beta(e1, e2)
  c1 <- prior_hyp$c1; c2 <- prior_hyp$c2   # omega ~ Gamma(c1, rate=c2)
  a1 <- prior_hyp$a1; a2 <- prior_hyp$a2   # r     ~ Gamma(a1, rate=a2)
  d1 <- prior_hyp$d1; d2 <- prior_hyp$d2   # phi   ~ Gamma-like (d1,d2), phi>2
  
  N_phase2 <- N_total - burn_in
  
  # ---- Random-walk stepfburn sizes (Phase 1) ----
  step_rho   <- 0.05
  step_beta  <- c(0.12, 0.12, 0.12)
  step_omega <- 0.25
  step_r     <- 3.0
  step_phi   <- 2.0
  
  # ---- Initial values ----
  rho_curr   <- 0.3
  beta_curr  <- c(0.1, 0.1, 0.1)
  omega_curr <- 0.1
  r_curr     <- 2
  phi_curr   <- 3
  b          <- 1
  
  n_param <- 8   # rho, beta(3), omega, r, phi, b
  
  all_samples     <- matrix(0, nrow = N_total, ncol = n_param)
  burn_in_samples <- matrix(0, nrow = burn_in, ncol = n_param)
  
  idx_rho   <- 1
  idx_beta  <- 2:4
  idx_omega <- 5
  idx_r     <- 6
  idx_phi   <- 7
  idx_b     <- 8
  
  pack_params <- function()
    c(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
  
  # ============================================================
  # Log-posterior functions
  # lp_*(focal_param | rest, Y) = log-likelihood + log-prior
  #
  # Each function takes the full parameter set explicitly so it
  # can be called with either the current value or a proposal.
  # Beta has a flat prior (constraints are enforced by the
  # proposal mechanism), so lp_beta returns only the likelihood.
  # ============================================================
  
  lp_rho <- function(rho, beta, omega, r, phi, b) {
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k) +
      dbeta(rho, e1, e2, log = TRUE)
  }
  
  lp_beta <- function(rho, beta, omega, r, phi, b) {
    # Flat prior on beta (support enforced by proposal constraint)
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k)
  }
  
  lp_omega <- function(rho, beta, omega, r, phi, b) {
    log_likelihood(Y, rho, beta, omega, r, phi, X_list, b, b0_mcmc, k) +
      dgamma(omega, c1, rate = c2, log = TRUE)
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
  #
  #   Receives the scalar log_accept_ratio and returns TRUE/FALSE.
  #   A non-finite ratio (NA, NaN, ±Inf) is treated as rejection.
  #
  #   Caller is responsible for forming:
  #     Phase 1:  log_ar = lp_*(prop) − lp_*(curr)
  #     Phase 2:  log_ar = [lp_*(prop) + lq_rev] − [lp_*(curr) + lq_fwd]
  # ============================================================
  mh_accept <- function(log_accept_ratio) {
    if (!is.finite(log_accept_ratio)) return(FALSE)
    log(runif(1)) < min(0, log_accept_ratio)
  }
  
  # ---- b: exact discrete sampler (used in both phases) ----
  sample_b <- function() {
    log_lik_b <- sapply(
      1:b0_mcmc,
      function(j) log_likelihood(Y, rho_curr, beta_curr, omega_curr,
                                 r_curr, phi_curr, X_list, j, b0_mcmc, k)
    )
    w <- exp(log_lik_b - max(log_lik_b))  # numerically stable softmax
    sample(1:b0_mcmc, 1, prob = w / sum(w))
  }
  
  acc_rho_p1 <- acc_beta_p1 <- acc_omega_p1 <- acc_r_p1 <- acc_phi_p1 <- 0
  
  # ============================================================
  # PHASE 1: Random-Walk MH  (iterations 1 … burn_in)
  #
  # Symmetric Gaussian proposal ⟹ proposal terms cancel:
  #   log_accept_ratio = r_prop_lposterior − r_curr_lposterior
  #                    = lp_*(prop) − lp_*(curr)
  # ============================================================
  for (iter in 1:burn_in) {
    if (iter %% 1000 == 0) cat("  Burn-in Iteration:", iter, "\n")
    
    # ---------- rho ----------
    repeat {
      rho_prop <- rho_curr + rnorm(1, 0, step_rho)
      if (rho_prop > 0 && rho_prop < 1) break
    }
    r_prop_lposterior <- lp_rho(rho_prop, beta_curr, omega_curr, r_curr, phi_curr, b)
    r_curr_lposterior <- lp_rho(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      rho_curr    <- rho_prop
      acc_rho_p1  <- acc_rho_p1 + 1
    }
    
    # ---------- beta ----------
    repeat {
      beta_prop <- beta_curr + rnorm(3, 0, step_beta)
      if (check_beta_constraints(beta_prop)) break
    }
    r_prop_lposterior <- lp_beta(rho_curr, beta_prop, omega_curr, r_curr, phi_curr, b)
    r_curr_lposterior <- lp_beta(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      beta_curr   <- beta_prop
      acc_beta_p1 <- acc_beta_p1 + 1
    }
    
    # ---------- omega ----------
    repeat {
      omega_prop <- omega_curr + rnorm(1, 0, step_omega)
      if (omega_prop > 0) break
    }
    r_prop_lposterior <- lp_omega(rho_curr, beta_curr, omega_prop, r_curr, phi_curr, b)
    r_curr_lposterior <- lp_omega(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      omega_curr    <- omega_prop
      acc_omega_p1  <- acc_omega_p1 + 1
    }
    
    # ---------- r ----------
    repeat {
      r_prop <- r_curr + rnorm(1, 0, step_r)
      if (r_prop > 0) break
    }
    r_prop_lposterior <- lp_r(rho_curr, beta_curr, omega_curr, r_prop, phi_curr, b)
    r_curr_lposterior <- lp_r(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      r_curr    <- r_prop
      acc_r_p1  <- acc_r_p1 + 1
    }
    
    # ---------- phi  (phi > 2) ----------
    repeat {
      phi_prop <- phi_curr + rnorm(1, 0, step_phi)
      if (phi_prop > 2) break
    }
    r_prop_lposterior <- lp_phi(rho_curr, beta_curr, omega_curr, r_curr, phi_prop, b)
    r_curr_lposterior <- lp_phi(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b)
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      phi_curr    <- phi_prop
      acc_phi_p1  <- acc_phi_p1 + 1
    }
    
    # ---------- b (discrete — exact posterior sampling) ----------
    b <- sample_b()
    
    burn_in_samples[iter, ] <- pack_params()
    all_samples[iter, ]     <- pack_params()
  }
  
  cat("\n  Phase 1 Acceptance Rates:\n")
  cat("    rho:  ", round(acc_rho_p1   / burn_in, 4), "\n")
  cat("    beta: ", round(acc_beta_p1  / burn_in, 4), "\n")
  cat("    omega:", round(acc_omega_p1 / burn_in, 4), "\n")
  cat("    r:    ", round(acc_r_p1     / burn_in, 4), "\n")
  cat("    phi:  ", round(acc_phi_p1   / burn_in, 4), "\n")
  
  # ============================================================
  # Empirical proposal moments from Phase 1 (drop first 1,000)
  # These parameterise the independent kernel used in Phase 2.
  # ============================================================
  warmup <- 1:1000
  
  mu_rho    <- mean(burn_in_samples[-warmup, idx_rho])
  sd_rho    <- sd(  burn_in_samples[-warmup, idx_rho])
  
  mu_beta   <- colMeans(burn_in_samples[-warmup, idx_beta])
  cov_beta  <- cov(    burn_in_samples[-warmup, idx_beta])
  
  mu_omega  <- mean(burn_in_samples[-warmup, idx_omega])
  sd_omega  <- sd(  burn_in_samples[-warmup, idx_omega])
  
  mu_r      <- mean(burn_in_samples[-warmup, idx_r])
  sd_r      <- sd(  burn_in_samples[-warmup, idx_r])
  
  mu_phi    <- mean(burn_in_samples[-warmup, idx_phi])
  sd_phi    <- sd(  burn_in_samples[-warmup, idx_phi])
  
  acc_rho_p2 <- acc_beta_p2 <- acc_omega_p2 <- acc_r_p2 <- acc_phi_p2 <- 0
  
  # ============================================================
  # PHASE 2: Independent-Kernel MH  (iterations burn_in+1 … N_total)
  #
  # Asymmetric independent proposal q(θ) fitted in Phase 1:
  #   log_accept_ratio =  r_prop_lposterior − r_curr_lposterior
  #
  #   where
  #     r_prop_lposterior = lp_*(θ')   + log q(θ  | θ')
  #                       = lp_*(θ')   + log q(θ)       [independent q]
  #
  #     r_curr_lposterior = lp_*(θ)    + log q(θ' | θ)
  #                       = lp_*(θ)    + log q(θ')      [independent q]
  #
  # The lq_rev / lq_fwd naming is:
  #   lq_fwd = log q(θ' | θ) = log q(θ')   [forward  proposal density]
  #   lq_rev = log q(θ  | θ') = log q(θ)   [reverse proposal density]
  # ============================================================
  for (iter in 1:N_phase2) {
    if (iter %% 1000 == 0) cat("  Independent-Kernel Iteration:", iter, "\n")
    
    # ---------- rho ----------
    repeat {
      rho_prop <- rnorm(1, mu_rho, sd_rho)
      if (rho_prop > 0 && rho_prop < 1) break
    }
    r_prop_lposterior <- lp_rho(rho_prop, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(rho_curr, mu_rho, sd_rho, log = TRUE)   # lq_rev
    r_curr_lposterior <- lp_rho(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(rho_prop, mu_rho, sd_rho, log = TRUE)   # lq_fwd
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      rho_curr    <- rho_prop
      acc_rho_p2  <- acc_rho_p2 + 1
    }
    
    # ---------- beta ----------
    repeat {
      beta_prop <- as.numeric(mvrnorm(1, mu_beta, cov_beta))
      if (check_beta_constraints(beta_prop)) break
    }
    r_prop_lposterior <- lp_beta(rho_curr, beta_prop, omega_curr, r_curr, phi_curr, b) +
      mvtnorm::dmvnorm(beta_curr, mu_beta, cov_beta, log = TRUE)   # lq_rev
    r_curr_lposterior <- lp_beta(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      mvtnorm::dmvnorm(beta_prop, mu_beta, cov_beta, log = TRUE)   # lq_fwd
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      beta_curr   <- beta_prop
      acc_beta_p2 <- acc_beta_p2 + 1
    }
    
    # ---------- omega ----------
    repeat {
      omega_prop <- rnorm(1, mu_omega, sd_omega)
      if (omega_prop > 0) break
    }
    r_prop_lposterior <- lp_omega(rho_curr, beta_curr, omega_prop, r_curr, phi_curr, b) +
      dnorm(omega_curr, mu_omega, sd_omega, log = TRUE)   # lq_rev
    r_curr_lposterior <- lp_omega(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(omega_prop, mu_omega, sd_omega, log = TRUE)   # lq_fwd
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      omega_curr    <- omega_prop
      acc_omega_p2  <- acc_omega_p2 + 1
    }
    
    # ---------- r ----------
    repeat {
      r_prop <- rnorm(1, mu_r, sd_r)
      if (r_prop > 0) break
    }
    r_prop_lposterior <- lp_r(rho_curr, beta_curr, omega_curr, r_prop, phi_curr, b) +
      dnorm(r_curr, mu_r, sd_r, log = TRUE)   # lq_rev
    r_curr_lposterior <- lp_r(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(r_prop, mu_r, sd_r, log = TRUE)   # lq_fwd
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      r_curr    <- r_prop
      acc_r_p2  <- acc_r_p2 + 1
    }
    
    # ---------- phi  (phi > 2) ----------
    repeat {
      phi_prop <- rnorm(1, mu_phi, sd_phi)
      if (phi_prop > 2) break
    }
    r_prop_lposterior <- lp_phi(rho_curr, beta_curr, omega_curr, r_curr, phi_prop, b) +
      dnorm(phi_curr, mu_phi, sd_phi, log = TRUE)   # lq_rev
    r_curr_lposterior <- lp_phi(rho_curr, beta_curr, omega_curr, r_curr, phi_curr, b) +
      dnorm(phi_prop, mu_phi, sd_phi, log = TRUE)   # lq_fwd
    log_accept_ratio  <- r_prop_lposterior - r_curr_lposterior
    if (mh_accept(log_accept_ratio)) {
      phi_curr    <- phi_prop
      acc_phi_p2  <- acc_phi_p2 + 1
    }
    
    # ---------- b (discrete — exact posterior sampling) ----------
    b <- sample_b()
    
    all_samples[burn_in + iter, ] <- pack_params()
  }
  
  cat("\n  Phase 2 Acceptance Rates:\n")
  cat("    rho:  ", round(acc_rho_p2   / N_phase2, 4), "\n")
  cat("    beta: ", round(acc_beta_p2  / N_phase2, 4), "\n")
  cat("    omega:", round(acc_omega_p2 / N_phase2, 4), "\n")
  cat("    r:    ", round(acc_r_p2     / N_phase2, 4), "\n")
  cat("    phi:  ", round(acc_phi_p2   / N_phase2, 4), "\n")
  
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
    samples      = posterior_samples,
    n_param      = n_param,
    k            = k,
    idx_rho      = idx_rho,
    idx_beta     = idx_beta,
    idx_omega    = idx_omega,
    idx_r        = idx_r,
    idx_phi      = idx_phi,
    idx_b        = idx_b,
    accept_rho   = acc_rho_p2   / N_phase2,
    accept_beta  = acc_beta_p2  / N_phase2,
    accept_omega = acc_omega_p2 / N_phase2,
    accept_r     = acc_r_p2     / N_phase2,
    accept_phi   = acc_phi_p2   / N_phase2
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
omega_true_c2 <- 0.12
b_true_c2     <- 1
lambda_init   <- 0.1

# ==============================================================================
# SECTION 8: REPLICATION BLOCK — 100 replications, parallel with checkpoint
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

M         <- 100
n_cont_c2 <- 7      # rho, beta1, beta2, beta3, omega1, r, phi

param_names_c2 <- c("rho", "beta1", "beta2", "beta3", "omega1", "r", "phi")

col_labels_c2 <- c(
  paste0("mean_",   param_names_c2),
  paste0("median_", param_names_c2),
  paste0("sd_",     param_names_c2),
  paste0("p025_",   param_names_c2),
  paste0("p975_",   param_names_c2),
  "mode_b1",
  "accept_rho", "accept_beta", "accept_omega", "accept_r", "accept_phi"
)

mode_fn <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

checkpoint_file_c2 <- "checkpoint_case2_k1_thinned.csv"

if (file.exists(checkpoint_file_c2)) {
  post_means_c2 <- as.matrix(read.csv(checkpoint_file_c2, check.names = FALSE))
  start_m_c2    <- nrow(post_means_c2) + 1
  cat("Case 2 (k=1): Resuming from replication", start_m_c2, "\n")
} else {
  post_means_c2 <- matrix(NA, nrow = 0, ncol = length(col_labels_c2))
  start_m_c2    <- 1
}

if (start_m_c2 <= M) {
  
  ncores <- max(1L, parallel::detectCores() - 1L)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  cat("Running Case 2 (k=1) replications", start_m_c2, "to", M,
      "on", ncores, "cores...\n")
  t0 <- Sys.time()
  
  sim_results_c2 <- foreach(
    m = start_m_c2:M,
    .combine  = rbind,
    .packages = c("mvtnorm", "MASS", "coda"),
    .export   = c("ZIBNB_data", "run_mcmc", "log_likelihood",
                  "rbnb", "log_dbnb", "check_beta_constraints", "log_fbnb_zero",
                  "prime_seeds", "n", "rho_true", "beta_true",
                  "r_true", "phi_true", "lambda_init",
                  "omega_true_c2", "b_true_c2",
                  "n_cont_c2", "param_names_c2", "mode_fn")
  ) %dopar% {
    
    set.seed(prime_seeds[m])
    
    X1_rep <- matrix(rnorm(n, 0, 1), 1, ncol = n)
    
    sim_rep <- ZIBNB_data(n           = n,
                          rho         = rho_true,
                          beta        = beta_true,
                          omega       = omega_true_c2,
                          r           = r_true,
                          phi         = phi_true,
                          X_list      = list(X1_rep),
                          b           = b_true_c2,
                          b0          = 1,
                          lambda_init = lambda_init,
                          k           = 1)
    Y_rep <- sim_rep$y
    
    res <- run_mcmc(Y        = Y_rep,
                    X_list   = list(X1_rep),
                    k        = 1,
                    N_total  = 40000,
                    burn_in  = 16000,
                    thin     = 4,
                    b0_mcmc  = 3)
    
    if (m == 1) {
      true_vals_c2_plot <- c(rho_true, beta_true, omega_true_c2, r_true, phi_true)
      
      pdf("traceplot_case2_k1_thinned.pdf")
      par(mfrow = c(2, 4))
      for (i in 1:n_cont_c2) {
        plot(res$samples[, i], type = "l",
             main = param_names_c2[i],
             ylab = "Value", xlab = "Thinned Iteration")
        abline(h = true_vals_c2_plot[i], col = "red", lwd = 2)
      }
      dev.off()
      
      pdf("acf_case2_k1_thinned.pdf")
      par(mfrow = c(2, 4))
      for (i in 1:n_cont_c2) {
        acf(res$samples[, i], main = param_names_c2[i])
      }
      dev.off()
      
      cat("  Plots saved: traceplot_case2_k1_thinned.pdf and acf_case2_k1_thinned.pdf\n")
    }
    
    mode_b1 <- mode_fn(res$samples[, res$idx_b])
    
    new_row <- c(
      apply(res$samples[, 1:n_cont_c2], 2, mean),
      apply(res$samples[, 1:n_cont_c2], 2, median),
      apply(res$samples[, 1:n_cont_c2], 2, sd),
      apply(res$samples[, 1:n_cont_c2], 2, quantile, probs = 0.025),
      apply(res$samples[, 1:n_cont_c2], 2, quantile, probs = 0.975),
      mode_b1,
      res$accept_rho, res$accept_beta, res$accept_omega,
      res$accept_r,   res$accept_phi
    )
    
    return(new_row)
  }
  
  stopCluster(cl)
  cat("Simulations complete in:",
      round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes.\n")
  
  post_means_c2           <- rbind(post_means_c2, sim_results_c2)
  colnames(post_means_c2) <- col_labels_c2
  
  write.csv(post_means_c2, checkpoint_file_c2, row.names = FALSE)
  cat("  Checkpoint saved:", checkpoint_file_c2, "\n")
}

# ==============================================================================
# SECTION 9: SIMULATION SUMMARY TABLE
# ==============================================================================
true_cont_c2 <- c(rho_true, beta_true, omega_true_c2, r_true, phi_true)
result_c2    <- make_summary_table(post_means_c2, param_names_c2, true_cont_c2)

cat("\n========== Case 2 (k=1, thinned) Summary across",
    nrow(post_means_c2), "replications ==========\n")
print(result_c2)
write.csv(result_c2, "result_case2_k1_thinned.csv", row.names = FALSE)

cat("\nAverage Phase 2 Acceptance Rates:\n")
cat("  rho:  ", round(mean(post_means_c2[, "accept_rho"]),   4), "\n")
cat("  beta: ", round(mean(post_means_c2[, "accept_beta"]),  4), "\n")
cat("  omega:", round(mean(post_means_c2[, "accept_omega"]), 4), "\n")
cat("  r:    ", round(mean(post_means_c2[, "accept_r"]),     4), "\n")
cat("  phi:  ", round(mean(post_means_c2[, "accept_phi"]),   4), "\n")

cat("\nMode of b1 across replications:\n")
cat("  Mode(b1):", mode_fn(post_means_c2[, "mode_b1"]),
    "  (True b1 =", b_true_c2, ")\n")