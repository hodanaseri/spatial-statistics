# Load necessary libraries
library(coda)  # For diagnostic tools for MCMC
library(MCMCpack)  # For MALA function

# Set the parameters for the Beta distribution
alpha <- 7.2
beta <- 8.7

# Define the log-posterior of Beta distribution
log_posterior <- function(x) {
  (alpha - 1) * log(x) + (beta - 1) * log(1 - x)
}

# Define the gradient of the log-posterior
grad_log_posterior <- function(x) {
  (alpha - 1) / x - (beta - 1) / (1 - x)
}

# MALA function to generate samples
mala_sampler <- function(log_posterior, grad_log_posterior, n, step_size, initial) {
  samples <- numeric(n)
  samples[1] <- initial
  for (i in 2:n) {
    current <- samples[i - 1]
    proposal_mean <- current + step_size * grad_log_posterior(current) / 2
    proposal <- rnorm(1, mean = proposal_mean, sd = sqrt(step_size))
    log_accept_ratio <- log_posterior(proposal) - log_posterior(current)
    log_accept_ratio <- log_accept_ratio + dnorm(current, mean = proposal_mean, sd = sqrt(step_size), log = TRUE)
    log_accept_ratio <- log_accept_ratio - dnorm(proposal, mean = current + step_size * grad_log_posterior(proposal) / 2, sd = sqrt(step_size), log = TRUE)
    if (log(runif(1)) < log_accept_ratio) {
      samples[i] <- proposal
    } else {
      samples[i] <- current
    }
  }
  return(samples)
}

# Generate samples using MALA
set.seed(123)
n_samples <- 10000
step_size <- 0.01
initial_value <- 0.5  # Starting point for MALA

samples_mala <- mala_sampler(log_posterior, grad_log_posterior, n_samples, step_size, initial_value)

# Diagnostics
mcmc_samples <- as.mcmc(samples_mala)
plot(mcmc_samples, main = "Trace and Density Plot of MALA Samples")

# Mixing diagnostics
autocorr.plot(mcmc_samples, main = "Autocorrelation Plot")
effectiveSize(mcmc_samples)



alpha <- 7.2
beta <- 8.7
log_posterior <- function(x) {
  if (x <= 0 || x >= 1) return(-Inf)
 (alpha - 1) * log(x) + (beta - 1) * log(1 - x)}
grad_log_posterior <- function(x) {
  if (x <= 0 || x >= 1) return(0) 
  (alpha - 1) / x - (beta - 1) / (1 - x)
}
mala_sampler <- function(log_posterior, grad_log_posterior, n, step_size, initial) {
  samples <- numeric(n)
  samples[1] <- initial
  for (i in 2:n) {
    current <- samples[i - 1]
    proposal_mean <- current + step_size * grad_log_posterior(current) / 2 
    proposal <- rnorm(1, mean = proposal_mean, sd = sqrt(step_size)) 
    if (proposal <= 0) proposal <- 1e-10
    if (proposal >= 1) proposal <- 1 - 1e-10
    log_accept_ratio <- log_posterior(proposal) - log_posterior(current)
    if (!is.finite(log_accept_ratio)) log_accept_ratio <- -Inf
    log_accept_ratio <- log_accept_ratio + dnorm(current, mean = proposal_mean, sd = sqrt(step_size), log = TRUE)
    log_accept_ratio <- log_accept_ratio - dnorm(proposal, mean = current + step_size * grad_log_posterior(proposal) / 2, sd = sqrt(step_size), log = TRUE)
    
    
    
    if (log(runif(1)) < log_accept_ratio) {
      samples[i] <- proposal
    } else {
      samples[i] <- current
    }
  }
  return(samples)
}

set.seed(123)
n_samples <- 10000
step_size <- 0.01
initial_value <- 0.5
samples_mala <- mala_sampler(log_posterior, grad_log_posterior, n_samples, step_size, initial_value)

mcmc_samples <- as.mcmc(samples_mala)
plot(mcmc_samples, main = "Trace and Density Plot of MALA Samples")
autocorr.plot(mcmc_samples, main = "Autocorrelation Plot")
effectiveSize(mcmc_samples)

step_size <- 0.1
samples_mala <- mala_sampler(log_posterior, grad_log_posterior, n_samples, step_size, initial_value)



step_size <- 0.1
samples_mala <- mala_sampler(log_posterior, grad_log_posterior, n_samples, step_size, initial_value)
burn_in <- 5000
samples_after_burn_in <- samples_mala[(burn_in + 1):n_samples]
thinned_samples <- samples_after_burn_in[seq(1, length(samples_after_burn_in), by = 10)]
mcmc_samples <- as.mcmc(thinned_samples)
plot(mcmc_samples, main = "Trace and Density Plot of Thinned MALA Samples after Burn-in")
