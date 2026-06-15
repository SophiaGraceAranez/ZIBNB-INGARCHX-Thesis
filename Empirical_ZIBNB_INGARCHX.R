# ==============================================================================
#  EMPIRICAL STUDY: Zero-Inflated Beta-Negative Binomial INGARCHX(1,1) Model
#  (ZIBNB-INGARCHX) via Two-Phase Adaptive MCMC (Chen & So, 2006)
#  Applied to Cenluz Weekly Leptospirosis Cases (2014-2025)
#  Exogenous variables: Weekly Total Rainfall (X1) and Weekly Max Temperature (X2)
#
#  Model structure (ZIBNB):
#    P(Y_t = 0)   = rho + (1 - rho) * f_BNB(0; r, gamma_t, phi)
#    P(Y_t = y)   = (1 - rho) * f_BNB(y; r, gamma_t, phi),  y > 0
#    lambda_t     = beta1 + beta2 * Y_{t-1} + beta3 * lambda_{t-1}
#                   + sum_i omega_i * X_i[t - b_i]          (if k >= 1)
#    gamma_t      = (phi - 1) / r * lambda_t
#  Parameters: rho, beta1, beta2, beta3, [omega1, omega2], r, phi, [b1, b2]
# ==============================================================================

library(mvtnorm)
library(MASS)
library(coda)
library(foreach)
library(doParallel)
library(parallel)
library(readr)
library(dplyr)

# ==============================================================================
# SECTION 1: Load and Prepare Data
# ==============================================================================
raw_data <- read_csv("cenluz2014_2025_real_datset.csv")

# Parse year-week, sort chronologically (data arrives newest-first)
raw_data <- raw_data %>%
  mutate(
    YEAR = as.integer(sub("-.*", "", week)),
    WEEK = as.integer(sub(".*-", "", week))
  ) %>%
  arrange(YEAR, WEEK)

Y  <- as.integer(raw_data$lepto_count)
X1 <- raw_data$total_rain
X2 <- raw_data$max_temp
n  <- length(Y)

cat("=== Descriptive Statistics ===\n")
cat("Total observations:", n, "\n")
cat("Number of zeros:   ", sum(Y == 0), "\n")
cat("Proportion zeros:  ", round(sum(Y == 0) / n, 4), "\n")
cat("Mean:              ", round(mean(Y), 4), "\n")
cat("Variance:          ", round(var(Y), 4), "\n")
cat("Var/Mean ratio:    ", round(var(Y) / mean(Y), 4), "\n")
cat("Min:", min(Y), " Max:", max(Y), "\n\n")

extreme_idx <- which(Y > 100)
if (length(extreme_idx) > 0) {
  cat("Observations above 100:\n")
  cat("Indices:", extreme_idx, "\n")
  cat("Values: ", Y[extreme_idx], "\n")
  cat("Weeks:  ", raw_data$week[extreme_idx], "\n\n")
}

# Scaling: x*_{i,t} = (X_{i,t} - min(X_{i,t})) / S_x^(i)
# where S_x^(i) is the standard deviation of the original covariate X_i
X1_scaled <- (X1 - min(X1)) / sd(X1)
X2_scaled <- (X2 - min(X2)) / sd(X2)

X1_mat <- matrix(X1_scaled, 1, ncol = n)   # total_rain
X2_mat <- matrix(X2_scaled, 1, ncol = n)   # max_temp

# ==============================================================================
# SECTION 2: Parameter Constraints
# ==============================================================================
check_beta_constraints <- function(beta_prop) {
  beta_prop[1] > 0 &&
    beta_prop[2] > 0 &&
    beta_prop[3] >= 0 &&
    (beta_prop[2] + beta_prop[3]) < 1
}

# ==============================================================================
# SECTION 3: ZIBNB Log-Likelihood
# ==============================================================================
#  rho   : structural zero-inflation probability in (0, 1)
#  beta  : c(beta1, beta2, beta3) — INGARCH conditional mean parameters
#  omega : vector length k — exogenous coefficients (> 0)
#  r     : dispersion (> 0)
#  phi   : tail parameter (> 2 for finite variance)
#  b     : integer vector length k — lag indices for exogenous variables
# ==============================================================================
log_likelihood <- function(Y, rho, beta, omega, r, phi, X_list, b, b0, k) {
  n        <- length(Y)
  log_like <- 0
  lambda_t <- rep(mean(Y[Y > 0]), n)
  
  for (t in (b0 + 1):n) {
    
    # --- Conditional mean (INGARCHX) ---
    lambda_t[t] <- beta[1] + beta[2] * Y[t - 1] + beta[3] * lambda_t[t - 1]
    if (k >= 1) {
      for (i in 1:k) {
        lambda_t[t] <- lambda_t[t] + omega[i] * X_list[[i]][1, t - b[i]]
      }
    }
    lambda_t[t] <- max(lambda_t[t], 1e-6)
    
    # --- BNB scale parameter ---
    gamma_t <- (phi - 1) / r * lambda_t[t]
    
    # --- BNB point mass at zero ---
    log_fbnb_0 <- lbeta(phi + r, gamma_t) - lbeta(phi, gamma_t)
    fbnb_0     <- exp(log_fbnb_0)
    
    if (Y[t] == 0) {
      # P(Y_t = 0) = rho + (1 - rho) * f_BNB(0)
      mix_zero <- rho + (1 - rho) * fbnb_0
      mix_zero <- max(mix_zero, 1e-300)
      log_like <- log_like + log(mix_zero)
    } else {
      # P(Y_t = y) = (1 - rho) * f_BNB(y),   y > 0
      log_fbnb_y <- lgamma(Y[t] + r) - lgamma(Y[t] + 1) - lgamma(r) +
        lbeta(phi + r, gamma_t + Y[t]) - lbeta(phi, gamma_t)
      log_like   <- log_like + log(1 - rho) + log_fbnb_y
    }
  }
  return(log_like)
}

# ==============================================================================
# SECTION 4: MCMC — Two-Phase Adaptive (Chen & So, 2006)
# ==============================================================================
run_mcmc_empirical <- function(Y, X_list, k, N = 20000, burn_in = 8000,
                               b0_mcmc = 3,
                               prior_hyp = list(
                                 c1 = 1, c2 = 1,   # Gamma prior on omega
                                 a1 = 3, a2 = 1,   # Gamma prior on r
                                 d1 = 3, d2 = 1,   # Gamma-like prior on phi
                                 e1 = 1, e2 = 1    # Beta prior on rho
                               )) {
  
  c1 <- prior_hyp$c1;  c2 <- prior_hyp$c2
  a1 <- prior_hyp$a1;  a2 <- prior_hyp$a2
  d1 <- prior_hyp$d1;  d2 <- prior_hyp$d2
  e1 <- prior_hyp$e1;  e2 <- prior_hyp$e2
  
  # --- Phase 1 step sizes ---
  # rho bounded (0,1): small step; zero proportion ~2.82% suggests very small rho
  step_size_rho   <- 0.03
  step_size_beta  <- c(0.05, 0.05, 0.05)
  step_size_r     <- 0.80
  step_size_phi   <- 0.60
  step_size_omega <- if (k > 0) rep(0.50, k) else NULL
  
  # --- Initialization ---
  rho_current   <- 0.03          # small structural zero probability
  beta_current  <- c(0.5, 0.1, 0.1)
  omega_current <- if (k > 0) rep(0.1, k) else NULL
  r_current     <- 2
  phi_current   <- 3
  b             <- if (k > 0) rep(1, k) else NULL
  
  # --- Parameter vector layout: [rho, beta(3), omega(k), r, phi, b(k)] ---
  n_param         <- 1 + 3 + k + 1 + 1 + k
  samples         <- matrix(0, nrow = N - burn_in, ncol = n_param)
  burn_in_samples <- matrix(0, nrow = burn_in,     ncol = n_param)
  
  rho_accept_count   <- 0L
  beta_accept_count  <- 0L
  omega_accept_count <- 0L
  r_accept_count     <- 0L
  phi_accept_count   <- 0L
  
  idx_rho   <- 1L
  idx_beta  <- 2:4
  idx_omega <- if (k > 0) (5L):(4L + k)                else integer(0)
  idx_r     <- 4L + k + 1L
  idx_phi   <- 4L + k + 2L
  idx_b     <- if (k > 0) (4L + k + 3L):(4L + 2L*k + 2L) else integer(0)
  
  pack_params <- function() {
    c(rho_current, beta_current,
      if (k > 0) omega_current else numeric(0),
      r_current, phi_current,
      if (k > 0) b else numeric(0))
  }
  
  log_prior_rho <- function(rho) {
    (e1 - 1) * log(rho) + (e2 - 1) * log(1 - rho)
  }
  
  # ============================================================
  # PHASE 1: Random-Walk Metropolis-Hastings (burn-in)
  # ============================================================
  for (iter in 1:burn_in) {
    if (iter %% 2000 == 0) cat("  [k=", k, "] Burn-in Iteration:", iter, "\n")
    
    # ---- rho ----
    repeat {
      rho_proposal <- rho_current + rnorm(1, 0, step_size_rho)
      if (rho_proposal > 0 && rho_proposal < 1) break
    }
    log_ar <- (log_likelihood(Y, rho_proposal, beta_current, omega_current,
                              r_current, phi_current, X_list, b, b0_mcmc, k) +
                 log_prior_rho(rho_proposal)) -
      (log_likelihood(Y, rho_current, beta_current, omega_current,
                      r_current, phi_current, X_list, b, b0_mcmc, k) +
         log_prior_rho(rho_current))
    if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
      rho_current <- rho_proposal
      rho_accept_count <- rho_accept_count + 1L
    }
    
    # ---- beta ----
    repeat {
      beta_proposal <- beta_current + rnorm(3, 0, step_size_beta)
      if (check_beta_constraints(beta_proposal)) break
    }
    log_ar <- log_likelihood(Y, rho_current, beta_proposal, omega_current,
                             r_current, phi_current, X_list, b, b0_mcmc, k) -
      log_likelihood(Y, rho_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k)
    if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
      beta_current <- beta_proposal
      beta_accept_count <- beta_accept_count + 1L
    }
    
    # ---- omega ----
    if (k > 0) {
      omega_proposal <- numeric(k)
      repeat {
        for (i in 1:k) omega_proposal[i] <- omega_current[i] + rnorm(1, 0, step_size_omega[i])
        if (all(omega_proposal > 0)) break
      }
      log_ar <- (log_likelihood(Y, rho_current, beta_current, omega_proposal,
                                r_current, phi_current, X_list, b, b0_mcmc, k) +
                   sum(dgamma(omega_proposal, shape = c1, rate = c2, log = TRUE))) -
        (log_likelihood(Y, rho_current, beta_current, omega_current,
                        r_current, phi_current, X_list, b, b0_mcmc, k) +
           sum(dgamma(omega_current, shape = c1, rate = c2, log = TRUE)))
      if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
        omega_current <- omega_proposal
        omega_accept_count <- omega_accept_count + 1L
      }
    }
    
    # ---- r ----
    repeat {
      r_proposal <- r_current + rnorm(1, 0, step_size_r)
      if (r_proposal > 0) break
    }
    log_ar <- (log_likelihood(Y, rho_current, beta_current, omega_current,
                              r_proposal, phi_current, X_list, b, b0_mcmc, k) +
                 dgamma(r_proposal, shape = a1, rate = a2, log = TRUE)) -
      (log_likelihood(Y, rho_current, beta_current, omega_current,
                      r_current, phi_current, X_list, b, b0_mcmc, k) +
         dgamma(r_current, shape = a1, rate = a2, log = TRUE))
    if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
      r_current <- r_proposal
      r_accept_count <- r_accept_count + 1L
    }
    
    # ---- phi ----
    repeat {
      phi_proposal <- phi_current + rnorm(1, 0, step_size_phi)
      if (phi_proposal > 2) break
    }
    log_ar <- ((d1 - 1) * log(phi_proposal) - d2 * phi_proposal +
                 log_likelihood(Y, rho_current, beta_current, omega_current,
                                r_current, phi_proposal, X_list, b, b0_mcmc, k)) -
      ((d1 - 1) * log(phi_current) - d2 * phi_current +
         log_likelihood(Y, rho_current, beta_current, omega_current,
                        r_current, phi_current, X_list, b, b0_mcmc, k))
    if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
      phi_current <- phi_proposal
      phi_accept_count <- phi_accept_count + 1L
    }
    
    # ---- b (discrete, Gibbs-style) ----
    if (k > 0) {
      for (bi_idx in 1:k) {
        b_temp <- b
        lik_b  <- sapply(1:b0_mcmc, function(j) {
          b_cand           <- b_temp
          b_cand[bi_idx]   <- j
          log_likelihood(Y, rho_current, beta_current, omega_current,
                         r_current, phi_current, X_list, b_cand, b0_mcmc, k)
        })
        max_lik  <- max(lik_b)
        prob     <- exp(lik_b - max_lik)
        prob     <- prob / sum(prob)
        cum_prob <- cumsum(prob)
        U        <- runif(1)
        if (!is.na(cum_prob[1]) && !is.nan(cum_prob[1])) {
          if (U < cum_prob[1]) {
            b[bi_idx] <- 1L
          } else {
            done <- FALSE; I <- 1L
            while (!done) {
              if (U > cum_prob[I] && U < cum_prob[I + 1]) {
                b[bi_idx] <- I + 1L; done <- TRUE
              } else { I <- I + 1L }
            }
          }
        }
      }
    }
    burn_in_samples[iter, ] <- pack_params()
  }
  
  cat("\n  [k=", k, "] Phase 1 Acceptance Rates:\n")
  cat("    rho:  ", round(rho_accept_count   / burn_in, 4), "\n")
  cat("    beta: ", round(beta_accept_count  / burn_in, 4), "\n")
  if (k > 0) cat("    omega:", round(omega_accept_count / burn_in, 4), "\n")
  cat("    r:    ", round(r_accept_count   / burn_in, 4), "\n")
  cat("    phi:  ", round(phi_accept_count / burn_in, 4), "\n")
  
  # --- Empirical moments from Phase 1 (discard first half) ---
  warmup   <- (ceiling(burn_in / 2) + 1):burn_in
  
  mu_rho   <- mean(burn_in_samples[warmup, idx_rho])
  cov_rho  <- var( burn_in_samples[warmup, idx_rho])
  mu_beta  <- colMeans(burn_in_samples[warmup, idx_beta, drop = FALSE])
  cov_beta <- cov(    burn_in_samples[warmup, idx_beta, drop = FALSE])
  
  if (k > 0) {
    if (k == 1L) {
      mu_omega  <- mean(burn_in_samples[warmup, idx_omega])
      cov_omega <- var( burn_in_samples[warmup, idx_omega])
    } else {
      mu_omega  <- colMeans(burn_in_samples[warmup, idx_omega, drop = FALSE])
      cov_omega <- cov(    burn_in_samples[warmup, idx_omega, drop = FALSE])
    }
  }
  
  mu_r    <- mean(burn_in_samples[warmup, idx_r])
  cov_r   <- var( burn_in_samples[warmup, idx_r])
  mu_phi  <- mean(burn_in_samples[warmup, idx_phi])
  cov_phi <- var( burn_in_samples[warmup, idx_phi])
  
  # Reset accept counters for Phase 2
  rho_accept_count <- beta_accept_count <- omega_accept_count <-
    r_accept_count <- phi_accept_count <- 0L
  
  # ============================================================
  # PHASE 2: Independent-Kernel Metropolis-Hastings
  # ============================================================
  for (iter in 1:(N - burn_in)) {
    if (iter %% 2000 == 0) cat("  [k=", k, "] IK Iteration:", iter, "\n")
    
    # ---- rho ----
    repeat {
      rho_proposal <- rnorm(1, mu_rho, sqrt(cov_rho))
      if (rho_proposal > 0 && rho_proposal < 1) break
    }
    log_g_curr <- dnorm(rho_current,  mu_rho, sqrt(cov_rho), log = TRUE)
    log_g_prop <- dnorm(rho_proposal, mu_rho, sqrt(cov_rho), log = TRUE)
    log_ar <- (log_likelihood(Y, rho_proposal, beta_current, omega_current,
                              r_current, phi_current, X_list, b, b0_mcmc, k) +
                 log_prior_rho(rho_proposal) + log_g_curr) -
      (log_likelihood(Y, rho_current, beta_current, omega_current,
                      r_current, phi_current, X_list, b, b0_mcmc, k) +
         log_prior_rho(rho_current) + log_g_prop)
    if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
      rho_current <- rho_proposal
      rho_accept_count <- rho_accept_count + 1L
    }
    
    # ---- beta ----
    repeat {
      beta_proposal <- as.numeric(mvrnorm(1, mu = mu_beta, Sigma = cov_beta))
      if (check_beta_constraints(beta_proposal)) break
    }
    log_g_curr <- mvtnorm::dmvnorm(beta_current,  mean = mu_beta, sigma = cov_beta, log = TRUE)
    log_g_prop <- mvtnorm::dmvnorm(beta_proposal, mean = mu_beta, sigma = cov_beta, log = TRUE)
    log_ar <- (log_likelihood(Y, rho_current, beta_proposal, omega_current,
                              r_current, phi_current, X_list, b, b0_mcmc, k) +
                 log_g_curr) -
      (log_likelihood(Y, rho_current, beta_current, omega_current,
                      r_current, phi_current, X_list, b, b0_mcmc, k) +
         log_g_prop)
    if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
      beta_current <- beta_proposal
      beta_accept_count <- beta_accept_count + 1L
    }
    
    # ---- omega ----
    if (k > 0) {
      if (k == 1L) {
        repeat {
          omega_proposal <- rnorm(1, mu_omega, sqrt(cov_omega))
          if (omega_proposal > 0) break
        }
        log_g_curr <- dnorm(omega_current,  mu_omega, sqrt(cov_omega), log = TRUE)
        log_g_prop <- dnorm(omega_proposal, mu_omega, sqrt(cov_omega), log = TRUE)
      } else {
        repeat {
          omega_proposal <- as.numeric(mvrnorm(1, mu = mu_omega, Sigma = cov_omega))
          if (all(omega_proposal > 0)) break
        }
        log_g_curr <- mvtnorm::dmvnorm(omega_current,  mean = mu_omega, sigma = cov_omega, log = TRUE)
        log_g_prop <- mvtnorm::dmvnorm(omega_proposal, mean = mu_omega, sigma = cov_omega, log = TRUE)
      }
      log_ar <- (log_likelihood(Y, rho_current, beta_current, omega_proposal,
                                r_current, phi_current, X_list, b, b0_mcmc, k) +
                   sum(dgamma(omega_proposal, shape = c1, rate = c2, log = TRUE)) +
                   log_g_curr) -
        (log_likelihood(Y, rho_current, beta_current, omega_current,
                        r_current, phi_current, X_list, b, b0_mcmc, k) +
           sum(dgamma(omega_current, shape = c1, rate = c2, log = TRUE)) +
           log_g_prop)
      if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
        omega_current <- omega_proposal
        omega_accept_count <- omega_accept_count + 1L
      }
    }
    
    # ---- r ----
    repeat {
      r_proposal <- rnorm(1, mu_r, sqrt(cov_r))
      if (r_proposal > 0) break
    }
    log_g_curr <- dnorm(r_current,  mu_r, sqrt(cov_r), log = TRUE)
    log_g_prop <- dnorm(r_proposal, mu_r, sqrt(cov_r), log = TRUE)
    log_ar <- (log_likelihood(Y, rho_current, beta_current, omega_current,
                              r_proposal, phi_current, X_list, b, b0_mcmc, k) +
                 dgamma(r_proposal, shape = a1, rate = a2, log = TRUE) +
                 log_g_curr) -
      (log_likelihood(Y, rho_current, beta_current, omega_current,
                      r_current, phi_current, X_list, b, b0_mcmc, k) +
         dgamma(r_current, shape = a1, rate = a2, log = TRUE) +
         log_g_prop)
    if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
      r_current <- r_proposal
      r_accept_count <- r_accept_count + 1L
    }
    
    # ---- phi ----
    repeat {
      phi_proposal <- rnorm(1, mu_phi, sqrt(cov_phi))
      if (phi_proposal > 2) break
    }
    log_g_curr <- dnorm(phi_current,  mu_phi, sqrt(cov_phi), log = TRUE)
    log_g_prop <- dnorm(phi_proposal, mu_phi, sqrt(cov_phi), log = TRUE)
    log_ar <- ((d1 - 1) * log(phi_proposal) - d2 * phi_proposal +
                 log_likelihood(Y, rho_current, beta_current, omega_current,
                                r_current, phi_proposal, X_list, b, b0_mcmc, k) +
                 log_g_curr) -
      ((d1 - 1) * log(phi_current) - d2 * phi_current +
         log_likelihood(Y, rho_current, beta_current, omega_current,
                        r_current, phi_current, X_list, b, b0_mcmc, k) +
         log_g_prop)
    if (!is.na(log_ar) && !is.nan(log_ar) && runif(1) < min(1, exp(log_ar))) {
      phi_current <- phi_proposal
      phi_accept_count <- phi_accept_count + 1L
    }
    
    # ---- b (discrete) ----
    if (k > 0) {
      for (bi_idx in 1:k) {
        b_temp <- b
        lik_b  <- sapply(1:b0_mcmc, function(j) {
          b_cand         <- b_temp
          b_cand[bi_idx] <- j
          log_likelihood(Y, rho_current, beta_current, omega_current,
                         r_current, phi_current, X_list, b_cand, b0_mcmc, k)
        })
        max_lik  <- max(lik_b)
        prob     <- exp(lik_b - max_lik)
        prob     <- prob / sum(prob)
        cum_prob <- cumsum(prob)
        U        <- runif(1)
        if (!is.na(cum_prob[1]) && !is.nan(cum_prob[1])) {
          if (U < cum_prob[1]) {
            b[bi_idx] <- 1L
          } else {
            done <- FALSE; I <- 1L
            while (!done) {
              if (U > cum_prob[I] && U < cum_prob[I + 1]) {
                b[bi_idx] <- I + 1L; done <- TRUE
              } else { I <- I + 1L }
            }
          }
        }
      }
    }
    samples[iter, ] <- pack_params()
  }
  
  n_phase2   <- N - burn_in
  n_blocks   <- if (k > 0) 5L else 4L
  tot_accept <- rho_accept_count + beta_accept_count +
    r_accept_count + phi_accept_count +
    if (k > 0) omega_accept_count else 0L
  accept_rate <- (tot_accept / n_blocks) / n_phase2
  
  # Per-parameter Phase 2 acceptance rates
  accept_rates_phase2 <- list(
    rho   = round(rho_accept_count   / n_phase2, 4),
    beta  = round(beta_accept_count  / n_phase2, 4),
    omega = if (k > 0) round(omega_accept_count / n_phase2, 4) else NA_real_,
    r     = round(r_accept_count     / n_phase2, 4),
    phi   = round(phi_accept_count   / n_phase2, 4),
    overall = round(accept_rate, 4)
  )
  
  cat("\n  [k=", k, "] Phase 2 Acceptance Rates:\n")
  cat("    rho:  ", round(rho_accept_count   / n_phase2, 4), "\n")
  cat("    beta: ", round(beta_accept_count  / n_phase2, 4), "\n")
  if (k > 0) cat("    omega:", round(omega_accept_count / n_phase2, 4), "\n")
  cat("    r:    ", round(r_accept_count   / n_phase2, 4), "\n")
  cat("    phi:  ", round(phi_accept_count / n_phase2, 4), "\n")
  cat("    Overall:", round(accept_rate, 4), "\n")
  
  return(list(
    samples              = samples,   n_param   = n_param,   k         = k,
    idx_rho              = idx_rho,   idx_beta  = idx_beta,  idx_omega = idx_omega,
    idx_r                = idx_r,     idx_phi   = idx_phi,   idx_b     = idx_b,
    accept_rates_phase2  = accept_rates_phase2
  ))
}

# ==============================================================================
# SECTION 5: Run Four Models in Parallel
# ==============================================================================
model_specs <- list(
  list(k = 0, X_list = NULL,                 seed = 42, label = "Model 1: No Exo"),
  list(k = 1, X_list = list(X1_mat),         seed = 42, label = "Model 2: Rainfall Only"),
  list(k = 1, X_list = list(X2_mat),         seed = 42, label = "Model 3: Max Temp Only"),
  list(k = 2, X_list = list(X1_mat, X2_mat), seed = 42, label = "Model 4: Rain + Max Temp")
)

ncores <- min(4L, max(1L, parallel::detectCores() - 1L))
cl     <- makeCluster(ncores)
registerDoParallel(cl)

cat("Running 4 ZIBNB-INGARCHX models on", ncores, "cores...\n")
t0 <- Sys.time()

all_results <- foreach(
  spec      = model_specs,
  .packages = c("mvtnorm", "MASS"),
  .export   = c("run_mcmc_empirical", "log_likelihood",
                "check_beta_constraints",
                "Y", "X1_mat", "X2_mat", "n")
) %dopar% {
  set.seed(spec$seed)
  cat("Starting:", spec$label, "\n")
  res <- run_mcmc_empirical(Y = Y, X_list = spec$X_list, k = spec$k)
  cat("Completed:", spec$label, "\n")
  res
}

stopCluster(cl)
cat("All models completed in",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes.\n")

res_m1 <- all_results[[1]]
res_m2 <- all_results[[2]]
res_m3 <- all_results[[3]]
res_m4 <- all_results[[4]]

# ==============================================================================
# SECTION 6: Summary Function
# ==============================================================================
mode_fn <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

summarize_empirical <- function(res, param_names, b_names = NULL) {
  samples <- res$samples
  k       <- res$k
  n_cont  <- 1 + 3 + k + 1 + 1    # rho + beta(3) + omega(k) + r + phi
  
  cont_smp   <- samples[, 1:n_cont, drop = FALSE]
  mean_est   <- round(apply(cont_smp, 2, mean),              4)
  median_est <- round(apply(cont_smp, 2, median),            4)
  sd_est     <- round(apply(cont_smp, 2, sd),                4)
  p025_est   <- round(apply(cont_smp, 2, quantile, 0.025),   4)
  p975_est   <- round(apply(cont_smp, 2, quantile, 0.975),   4)
  
  result <- data.frame(
    Parameter = param_names,
    Mean      = mean_est,
    Median    = median_est,
    Std       = sd_est,
    P0.025    = p025_est,
    P0.975    = p975_est,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  if (k > 0 && !is.null(b_names)) {
    for (i in 1:k) {
      mode_b <- mode_fn(samples[, res$idx_b[i]])
      result <- rbind(result, data.frame(
        Parameter = b_names[i], Mean = mode_b,
        Median = NA, Std = NA, P0.025 = NA, P0.975 = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(result)
}

# ==============================================================================
# SECTION 7: Geweke Z-Score and Inefficiency Factor
# ==============================================================================
compute_geweke_IF <- function(res, param_names, model_label) {
  samples  <- res$samples
  k        <- res$k
  n_cont   <- 1 + 3 + k + 1 + 1
  cont_smp <- samples[, 1:n_cont, drop = FALSE]
  
  geweke_z <- numeric(n_cont)
  IF_vals  <- numeric(n_cont)
  
  for (i in 1:n_cont) {
    chain_i     <- mcmc(cont_smp[, i])
    geweke_z[i] <- round(geweke.diag(chain_i, frac1 = 0.1, frac2 = 0.5)$z, 4)
    chain_var   <- var(cont_smp[, i])
    spec0       <- spectrum0.ar(chain_i)$spec
    IF_vals[i]  <- round(spec0 / chain_var, 4)
  }
  
  # Build per-parameter acceptance rates vector matching param_names order
  ar <- res$accept_rates_phase2
  accept_vec <- c(
    ar$rho,
    rep(ar$beta, 3),                              # beta1, beta2, beta3 share one block
    if (k > 0) rep(ar$omega, k) else numeric(0),  # omega(s)
    ar$r,
    ar$phi
  )
  
  diag_table <- data.frame(
    Parameter   = param_names,
    Geweke_Z    = geweke_z,
    IF          = IF_vals,
    AcceptRate  = accept_vec,
    Converged   = ifelse(abs(geweke_z) <= 2, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  cat("\n===", model_label, "— Geweke Z-Score and Inefficiency Factor ===\n")
  print(diag_table)
  cat("Convergence:", sum(abs(geweke_z) <= 2), "of", n_cont,
      "parameters within [-2, 2].\n")
  cat("IF range: [", min(IF_vals), ",", max(IF_vals), "]\n")
  cat("IF mean: ",   round(mean(IF_vals), 4), "\n")
  
  return(diag_table)
}

# ==============================================================================
# SECTION 7b: Comprehensive Diagnostics (Geweke, IF, Bias, Rel.Bias, AR, RMSE)
# ==============================================================================
compute_comprehensive_diag <- function(res, param_names, pred_obj,
                                       Y, model_label, b0_mcmc = 3) {
  samples  <- res$samples
  k        <- res$k
  n_cont   <- 1 + 3 + k + 1 + 1
  cont_smp <- samples[, 1:n_cont, drop = FALSE]
  
  # --- Geweke Z and IF per parameter ---
  geweke_z <- numeric(n_cont)
  IF_vals  <- numeric(n_cont)
  for (i in 1:n_cont) {
    chain_i     <- mcmc(cont_smp[, i])
    geweke_z[i] <- round(geweke.diag(chain_i, frac1 = 0.1, frac2 = 0.5)$z, 4)
    chain_var   <- var(cont_smp[, i])
    spec0       <- spectrum0.ar(chain_i)$spec
    IF_vals[i]  <- round(spec0 / chain_var, 4)
  }
  
  # --- Acceptance rates per parameter (Phase 2) ---
  ar <- res$accept_rates_phase2
  accept_vec <- c(
    ar$rho,
    rep(ar$beta, 3),
    if (k > 0) rep(ar$omega, k) else numeric(0),
    ar$r,
    ar$phi
  )
  
  # --- Bias, Relative Bias, RMSE from predictions ---
  Y_pred    <- pred_obj$Y_pred
  valid_idx <- (b0_mcmc + 1):length(Y)
  resid_vec <- Y[valid_idx] - Y_pred[valid_idx]
  
  bias      <- round(mean(resid_vec), 6)
  rel_bias  <- round(bias / mean(Y[valid_idx]), 6)
  rmse      <- round(sqrt(mean(resid_vec^2)), 6)
  
  # --- Assemble per-parameter table ---
  diag_df <- data.frame(
    Model        = model_label,
    Parameter    = param_names,
    Geweke_Z     = round(geweke_z, 4),
    IF           = round(IF_vals,  4),
    Bias         = bias,         # model-level, repeated per param
    Rel_Bias     = rel_bias,
    AcceptRate   = round(accept_vec, 4),
    RMSE         = rmse,
    Converged    = ifelse(abs(geweke_z) <= 2, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  return(diag_df)
}

# ==============================================================================
# SECTION 8: DIC Computation
# ==============================================================================
compute_DIC <- function(samples, Y, X_list, b0, k,
                        idx_rho, idx_beta, idx_omega,
                        idx_r, idx_phi, idx_b) {
  n_samples <- nrow(samples)
  
  rho_bar   <- mean(samples[, idx_rho])
  beta_bar  <- colMeans(samples[, idx_beta, drop = FALSE])
  r_bar     <- mean(samples[, idx_r])
  phi_bar   <- mean(samples[, idx_phi])
  omega_bar <- if (k > 0) colMeans(samples[, idx_omega, drop = FALSE]) else NULL
  b_bar     <- if (k > 0) sapply(1:k, function(i) mode_fn(samples[, idx_b[i]])) else NULL
  
  llik_bar <- log_likelihood(Y, rho_bar, beta_bar, omega_bar,
                             r_bar, phi_bar, X_list, b_bar, b0, k)
  D_bar    <- -2 * llik_bar
  
  llik_smp <- sapply(1:n_samples, function(s) {
    omega_s <- if (k > 0) samples[s, idx_omega] else NULL
    b_s     <- if (k > 0) sapply(1:k, function(i) samples[s, idx_b[i]]) else NULL
    log_likelihood(Y,
                   samples[s, idx_rho],
                   samples[s, idx_beta],
                   omega_s,
                   samples[s, idx_r],
                   samples[s, idx_phi],
                   X_list, b_s, b0, k)
  })
  
  E_D <- mean(-2 * llik_smp)
  pD  <- E_D - D_bar
  DIC <- D_bar + 2 * pD
  
  return(list(DIC   = round(DIC,   4),
              pD    = round(pD,    4),
              D_bar = round(D_bar, 4),
              E_D   = round(E_D,   4)))
}

# ==============================================================================
# SECTION 9: In-Sample Prediction and Residual Diagnostics
# ==============================================================================
compute_predicted <- function(res, Y, X_list, b0_mcmc = 3,
                              model_label, year_start = 2014) {
  samples <- res$samples
  k       <- res$k
  
  rho_bar   <- mean(samples[, res$idx_rho])
  beta_bar  <- colMeans(samples[, res$idx_beta, drop = FALSE])
  r_bar     <- mean(samples[, res$idx_r])
  phi_bar   <- mean(samples[, res$idx_phi])
  omega_bar <- if (k > 0) colMeans(samples[, res$idx_omega, drop = FALSE]) else NULL
  b_bar     <- if (k > 0) sapply(1:k, function(i) mode_fn(samples[, res$idx_b[i]])) else NULL
  
  n        <- length(Y)
  lambda_t <- rep(mean(Y[Y > 0]), n)
  mu_t     <- numeric(n)
  var_t    <- numeric(n)
  
  for (t in (b0_mcmc + 1):n) {
    
    # Conditional mean
    lambda_t[t] <- beta_bar[1] + beta_bar[2] * Y[t - 1] + beta_bar[3] * lambda_t[t - 1]
    if (k >= 1) {
      for (i in 1:k) {
        lambda_t[t] <- lambda_t[t] + omega_bar[i] * X_list[[i]][1, t - b_bar[i]]
      }
    }
    lambda_t[t] <- max(lambda_t[t], 1e-6)
    
    # ZIBNB marginal mean: E[Y_t] = (1 - rho) * lambda_t
    mu_t[t] <- (1 - rho_bar) * lambda_t[t]
    
    # ZIBNB marginal variance (requires phi > 2)
    if (phi_bar > 2) {
      gamma_t    <- (phi_bar - 1) / r_bar * lambda_t[t]
      EY2m1      <- r_bar * (r_bar + 1) * gamma_t * (gamma_t + 1) /
        ((phi_bar - 1) * (phi_bar - 2))
      EY_bnb     <- lambda_t[t]
      EY2_bnb    <- EY2m1 + EY_bnb
      Var_bnb    <- EY2_bnb - EY_bnb^2
      var_t[t]   <- (1 - rho_bar) * Var_bnb + rho_bar * (1 - rho_bar) * lambda_t[t]^2
    } else {
      var_t[t]   <- mu_t[t] * (1 + mu_t[t])
    }
    var_t[t] <- max(var_t[t], 1e-6)
  }
  
  residuals <- Y - mu_t
  std_resid <- residuals / sqrt(var_t)
  std_resid[1:b0_mcmc] <- NA
  
  time_seq    <- seq(from = as.Date(paste0(year_start, "-01-01")),
                     by = "week", length.out = n)
  year_breaks <- seq(from = as.Date(paste0(year_start, "-01-01")),
                     to = max(time_seq), by = "3 years")
  
  safe_label <- gsub("[: ]", "_", model_label)
  pdf(paste0("prediction_diagnostics_CenluzLepto_ZIBNB_", safe_label, ".pdf"),
      width = 12, height = 8)
  
  layout(matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = TRUE),
         heights = c(2, 1.2))
  
  par(mar = c(4, 4, 3, 2))
  plot(time_seq, Y, type = "l", lty = 2, col = "blue",
       xlab = "", ylab = "Cases",
       main = paste0("Cenluz Lepto ZIBNB-INGARCHX — ", model_label),
       cex.main = 1.1, xaxt = "n")
  axis.Date(1, at = year_breaks, labels = format(year_breaks, "%Y"),
            las = 1, cex.axis = 0.85)
  lines(time_seq, mu_t, col = "red", lwd = 1.5)
  legend("topleft", legend = c("Observed", "Predicted"),
         lty = c(2, 1), col = c("blue", "red"), bty = "n", cex = 0.85)
  
  par(mar = c(4, 4, 3, 1))
  plot(std_resid, type = "l", col = "black",
       main = "Standardized Residuals",
       xlab = "Time", ylab = "Residual", cex.main = 0.95)
  abline(h = 0, col = "red", lwd = 1.5)
  
  par(mar = c(4, 4, 3, 1))
  acf(na.omit(std_resid), main = "ACF of Residuals",
      col = "black", cex.main = 0.95)
  
  par(mar = c(4, 4, 3, 1))
  acf(na.omit(std_resid)^2, main = "ACF of Squared Residuals",
      col = "black", cex.main = 0.95)
  
  dev.off()
  cat("Prediction plot saved for Cenluz Lepto ZIBNB:", safe_label, "\n")
  
  return(list(Y_pred    = mu_t,
              std_resid = std_resid,
              lambda_t  = lambda_t,
              rho_bar   = rho_bar))
}

# ==============================================================================
# SECTION 10: Traceplots and ACF Plots
# ==============================================================================
plot_empirical <- function(res, param_names, model_label) {
  n_cont  <- length(param_names)
  samples <- res$samples
  nr      <- ceiling(n_cont / 3)
  
  pdf(paste0("traceplot_CenluzLepto_ZIBNB_", model_label, ".pdf"))
  par(mfrow = c(nr, 3))
  for (i in 1:n_cont) {
    plot(samples[, i], type = "l",
         main = param_names[i], ylab = "Value", xlab = "Iteration")
  }
  dev.off()
  
  pdf(paste0("acf_CenluzLepto_ZIBNB_", model_label, ".pdf"))
  par(mfrow = c(nr, 3))
  for (i in 1:n_cont) {
    acf(samples[, i], main = param_names[i])
  }
  dev.off()
  
  cat("Plots saved for Cenluz Lepto ZIBNB", model_label, "\n")
}

# ==============================================================================
# SECTION 11: Run All Summaries, Diagnostics, Predictions
# ==============================================================================
# Parameter names: rho, beta1, beta2, beta3, [omega1, omega2], r, phi
param_m1 <- c("rho", "beta1", "beta2", "beta3", "r", "phi")
param_m2 <- c("rho", "beta1", "beta2", "beta3", "omega1", "r", "phi")
param_m3 <- c("rho", "beta1", "beta2", "beta3", "omega1", "r", "phi")
param_m4 <- c("rho", "beta1", "beta2", "beta3", "omega1", "omega2", "r", "phi")

summary_m1 <- summarize_empirical(res_m1, param_m1)
summary_m2 <- summarize_empirical(res_m2, param_m2, b_names = "b1")
summary_m3 <- summarize_empirical(res_m3, param_m3, b_names = "b1")
summary_m4 <- summarize_empirical(res_m4, param_m4, b_names = c("b1", "b2"))

cat("\n========== Model 1 Summary (ZIBNB, No Exo) ==========\n")
print(summary_m1)
cat("\n========== Model 2 Summary (ZIBNB, Rainfall Only) ==========\n")
print(summary_m2)
cat("\n========== Model 3 Summary (ZIBNB, Max Temp Only) ==========\n")
print(summary_m3)
cat("\n========== Model 4 Summary (ZIBNB, Rain + Max Temp) ==========\n")
print(summary_m4)

diag_m1 <- compute_geweke_IF(res_m1, param_m1, "Model 1: No Exo")
diag_m2 <- compute_geweke_IF(res_m2, param_m2, "Model 2: Rainfall Only")
diag_m3 <- compute_geweke_IF(res_m3, param_m3, "Model 3: Max Temp Only")
diag_m4 <- compute_geweke_IF(res_m4, param_m4, "Model 4: Rain + Max Temp")

cat("\n========== DIC Comparison ==========\n")

dic_m1 <- compute_DIC(res_m1$samples, Y, NULL,             3, 0,
                      res_m1$idx_rho, res_m1$idx_beta, res_m1$idx_omega,
                      res_m1$idx_r,   res_m1$idx_phi,  res_m1$idx_b)

dic_m2 <- compute_DIC(res_m2$samples, Y, list(X1_mat),    3, 1,
                      res_m2$idx_rho, res_m2$idx_beta, res_m2$idx_omega,
                      res_m2$idx_r,   res_m2$idx_phi,  res_m2$idx_b)

dic_m3 <- compute_DIC(res_m3$samples, Y, list(X2_mat),    3, 1,
                      res_m3$idx_rho, res_m3$idx_beta, res_m3$idx_omega,
                      res_m3$idx_r,   res_m3$idx_phi,  res_m3$idx_b)

dic_m4 <- compute_DIC(res_m4$samples, Y, list(X1_mat, X2_mat), 3, 2,
                      res_m4$idx_rho, res_m4$idx_beta, res_m4$idx_omega,
                      res_m4$idx_r,   res_m4$idx_phi,  res_m4$idx_b)

dic_table <- data.frame(
  Model = c("Model 1: No Exo",
            "Model 2: Rainfall Only",
            "Model 3: Max Temp Only",
            "Model 4: Rain + Max Temp"),
  D_bar = c(dic_m1$D_bar, dic_m2$D_bar, dic_m3$D_bar, dic_m4$D_bar),
  E_D   = c(dic_m1$E_D,   dic_m2$E_D,   dic_m3$E_D,   dic_m4$E_D),
  pD    = c(dic_m1$pD,    dic_m2$pD,    dic_m3$pD,    dic_m4$pD),
  DIC   = c(dic_m1$DIC,   dic_m2$DIC,   dic_m3$DIC,   dic_m4$DIC)
)
print(dic_table)
cat("\nBest model (lowest DIC):",
    dic_table$Model[which.min(dic_table$DIC)], "\n")

pred_m1 <- compute_predicted(res_m1, Y, NULL,
                             model_label = "Model1_No_Exo",       year_start = 2014)
pred_m2 <- compute_predicted(res_m2, Y, list(X1_mat),
                             model_label = "Model2_Rainfall_Only", year_start = 2014)
pred_m3 <- compute_predicted(res_m3, Y, list(X2_mat),
                             model_label = "Model3_MaxTemp_Only",  year_start = 2014)
pred_m4 <- compute_predicted(res_m4, Y, list(X1_mat, X2_mat),
                             model_label = "Model4_Rain_MaxTemp",  year_start = 2014)

plot_empirical(res_m1, param_m1, "model1")
plot_empirical(res_m2, param_m2, "model2")
plot_empirical(res_m3, param_m3, "model3")
plot_empirical(res_m4, param_m4, "model4")

# ==============================================================================
# SECTION 11b: Comprehensive Diagnostics Table (Geweke, IF, Bias, RelBias, AR, RMSE)
# ==============================================================================
comp_diag_m1 <- compute_comprehensive_diag(
  res_m1, param_m1, pred_m1, Y, model_label = "Model1_No_Exo")
comp_diag_m2 <- compute_comprehensive_diag(
  res_m2, param_m2, pred_m2, Y, model_label = "Model2_Rainfall_Only")
comp_diag_m3 <- compute_comprehensive_diag(
  res_m3, param_m3, pred_m3, Y, model_label = "Model3_MaxTemp_Only")
comp_diag_m4 <- compute_comprehensive_diag(
  res_m4, param_m4, pred_m4, Y, model_label = "Model4_Rain_MaxTemp")

comprehensive_diag_all <- rbind(comp_diag_m1, comp_diag_m2,
                                comp_diag_m3, comp_diag_m4)

cat("\n========== Comprehensive Diagnostics (All Models) ==========\n")
print(comprehensive_diag_all)

# ==============================================================================
# SECTION 12: Save All Results
# ==============================================================================
write.csv(summary_m1, "CenluzLepto_ZIBNB_result_model1.csv",      row.names = FALSE)
write.csv(summary_m2, "CenluzLepto_ZIBNB_result_model2.csv",      row.names = FALSE)
write.csv(summary_m3, "CenluzLepto_ZIBNB_result_model3.csv",      row.names = FALSE)
write.csv(summary_m4, "CenluzLepto_ZIBNB_result_model4.csv",      row.names = FALSE)
write.csv(dic_table,  "CenluzLepto_ZIBNB_DIC_comparison.csv",     row.names = FALSE)
write.csv(diag_m1,    "CenluzLepto_ZIBNB_diagnostics_model1.csv", row.names = FALSE)
write.csv(diag_m2,    "CenluzLepto_ZIBNB_diagnostics_model2.csv", row.names = FALSE)
write.csv(diag_m3,    "CenluzLepto_ZIBNB_diagnostics_model3.csv", row.names = FALSE)
write.csv(diag_m4,    "CenluzLepto_ZIBNB_diagnostics_model4.csv", row.names = FALSE)

# New: comprehensive diagnostics CSV (Geweke, IF, Bias, Rel.Bias, AcceptRate, RMSE)
write.csv(comprehensive_diag_all,
          "CenluzLepto_ZIBNB_comprehensive_diagnostics.csv",
          row.names = FALSE)

cat("\nAll Cenluz Lepto ZIBNB-INGARCHX results saved.\n")
