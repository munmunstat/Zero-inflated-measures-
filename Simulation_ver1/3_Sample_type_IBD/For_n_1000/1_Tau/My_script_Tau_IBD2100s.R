############################################
# SLURM JOB SCRIPT: Kendall Tau Correlation for IBD sample
############################################

library(copula)
library(MASS)
library(Kendall)
library(cubature)
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(phyloseq)

# -----------------------------------------------
# Accept the job array index from Slurm
# -----------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
cat("Running Slurm array task:", task_id, "\n")

# ---------------------------------------------------------
# Define copula-based tau function
# ---------------------------------------------------------

tau <- function(mu, sigma, pi, R) {
  
  # Create a bivariate Gaussian copula from the covariance matrix R
  gaussCop <- normalCopula( param = R, dim = 2,dispstr = "un")
  
  # Define the Copula CDF and density functions
  copula_cdf <- function(u, v) {
    pCopula(c(u, v), copula =  gaussCop)
  }
  
  copula_pdf <- function(u, v) {
    dCopula(c(u, v), copula =  gaussCop)
  }
  
  # Function for the second term integral
  second_term_fun <- function(u) {
    2 * copula_cdf(u[1], pi[2]) * copula_pdf(u[1], u[2])
  }
  
  # Integrate over u from pi[1] to 1 for the second term
  second_term_integral <- hcubature(second_term_fun, lower = c(pi[1], 0), upper = c(1,pi[2]))$integral
  
  # Function for the third term integral
  third_term_fun <- function(u) {
    2 * copula_cdf(pi[1], u[2]) * copula_pdf(u[1], u[2])
  }
  
  # Integrate over u from pi[2] to 1 for the third term
  third_term_integral <- hcubature(third_term_fun, lower = c(0, pi[2]), upper = c(pi[1],1))$integral
  
  # Function for the fourth term (double integral)
  fourth_term_fun <- function(u) {
    4 * copula_cdf(u[1], u[2]) * copula_pdf(u[1], u[2])
  }
  
  # Double integration for the fourth term
  fourth_term_integral <- hcubature(fourth_term_fun, lower = c(pi[1], pi[2]), upper = c(1, 1))$integral
  
  # Calculate the total expected value (Kendall's Tau)
  total_expectation <- copula_cdf(pi[1], pi[2])^2 + second_term_integral + third_term_integral + fourth_term_integral
  
  # Return the result for Kendall's Tau
  return(total_expectation - 1)
}


# ---------------------------------------------------------
# Matrix calculation applying the tau() function
# ---------------------------------------------------------
tauMatrix <- function(mu, sigma, pi, r, f) {
  J <- length(mu)
  tau <- matrix(0, J, J)
  for (i in 1:(J-1)) {
    cat("i =", i, "\n")
    for (j in (i+1):J) {
      cat("\t j =", j, "\n")
      idx <- c(i,j)
      tau[i,j] <- f(mu=mu[idx], sigma=sigma[idx], pi=pi[idx], R=r[i,j])
    }
  }
  ut <- lower.tri(tau)
  tau[ut] <- t(tau)[ut]
  return(tau)
}

# ---------------------------------------------------------
# Empirical tau star calculation
# ---------------------------------------------------------

calculate_tau_star <- function(X, Y) {
  n <- length(X)
  
  p_00 <- sum(X == 0 & Y == 0) / n
  p_01 <- sum(X == 0 & Y > 0) / n
  p_10 <- sum(X > 0 & Y == 0) / n
  p_11 <- sum(X > 0 & Y > 0) / n
  
  X_pos <- X[X > 0 & Y > 0]
  Y_pos <- Y[X > 0 & Y > 0]
  
  if (length(X_pos) >= 3 && length(Y_pos) >= 3) {
    tau_11 <- Kendall(X_pos, Y_pos)$tau[1]
  } else {
    tau_11 <- 0
  }
  
  X10 <- X[Y == 0]
  X11 <- X[Y > 0]
  
  if (length(X10) == 0 || length(X11) == 0) {
    p1 <- 0
  } else {
    count_1 <- sum(sapply(X10, function(x) sum(X11 < x)))
    p1 <- count_1 / (length(X11) * length(X10))
  }
  
  Y01 <- Y[X == 0]
  Y11 <- Y[X > 0]
  
  if (length(Y01) == 0 || length(Y11) == 0) {
    p2 <- 0
  } else {
    count_2 <- sum(sapply(Y01, function(y) sum(Y11 < y)))
    p2 <- count_2 / (length(Y01) * length(Y11))
  }
  
  tau_star <- p_11^2 * tau_11 +
    2 * (p_00 * p_11 - p_01 * p_10) +
    2 * p_11 * (p_10 * (1 - 2 * p1) + p_01 * (1 - 2 * p2))
  
  return(tau_star)
}

# ---------------------------------------------------------
# Matrix calculation applying tau star
# ---------------------------------------------------------

kendallTaustar <- function(data) {
  J <- ncol(data)
  tau_star_matrix <- matrix(0, J, J)
  for (i in 1:(J-1)) {
    for (j in (i+1):J) {
      tau_star_matrix[i, j] <- calculate_tau_star(data[, i], data[, j])
    }
  }
  tau_star_matrix[lower.tri(tau_star_matrix)] <- 
    t(tau_star_matrix)[lower.tri(tau_star_matrix)]
  return(tau_star_matrix)
}

# -----------------------------------------------
# Load simulated data and compute rarefied tau
# -----------------------------------------------

sim_data_file <- paste0("IBD2_sim_data_", task_id, ".RData")
if (!file.exists(sim_data_file)) {
  stop(paste("Simulated data file", sim_data_file, "not found!"))
}
load(sim_data_file)

ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
rarefied_counts <- as(otu_table(ps_rare), "matrix")

taustar_matrix_result <- kendallTaustar(rarefied_counts)
print(taustar_matrix_result)

# -----------------------------------------------
# Load fitted data and compute copula-based tau
# -----------------------------------------------

fit_data_file <- paste0("fit_IBD2_data_", task_id, ".RData")
if (!file.exists(fit_data_file)) {
  stop(paste("Fitted data file", fit_data_file, "not found!"))
}
load(fit_data_file)

mu    <- fitted_simdata$EM_fit$fit$mu
pi    <- fitted_simdata$EM_fit$fit$pi0
R     <- fitted_simdata$EM_fit$fit$Sigma
sigma <- fitted_simdata$EM_fit$fit$sigma

tau_matrix_result <- tauMatrix(mu, sigma, pi, R, tau)
print(tau_matrix_result)

# -----------------------------------------------
# Load true tau matrix
# -----------------------------------------------

load("IBD_truetau_matrix_result.RData")

# -----------------------------------------------
# Compare metrics
# -----------------------------------------------

lower_tri_new <- tau_matrix_result[lower.tri(tau_matrix_result)]
lower_tri_taustar <- taustar_matrix_result[lower.tri(taustar_matrix_result)]
lower_tri_trutau <- IBD_truetau_matrix_result[lower.tri(IBD_truetau_matrix_result)]

mse_newdeveloped <- mean((lower_tri_new - lower_tri_trutau)^2)
mse_taustar <- mean((lower_tri_taustar - lower_tri_trutau)^2)

mad_new <- mean(abs(lower_tri_new - lower_tri_trutau))
mad_taustar <- mean(abs(lower_tri_taustar - lower_tri_trutau))

cat("Task", task_id, "MSE (copula-based tau):", mse_newdeveloped, "\n")
cat("Task", task_id, "MAD (copula-based tau):", mad_new, "\n")
cat("Task", task_id, "MSE (tau star):", mse_taustar, "\n")
cat("Task", task_id, "MAD (tau star):", mad_taustar, "\n")

save.image(file = paste0("TI2_simulation_task_", task_id, ".RData"))
cat("Task", task_id, "completed.\n")
