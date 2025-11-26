############################################
# SLURM JOB SCRIPT: Spearman's rho Correlation for Stool sample
############################################
library(VC2copula)
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
# Define copula-based spearman function
# ---------------------------------------------------------

spearman_rho <- function(mu, sigma, pi, R) {
  
  # Create a Gaussian copula from the covariance matrix R
  gaussCop <- normalCopula(param = R, dim = 2, dispstr = "un")
  
  # Define the Copula CDF and density functions
  copula_cdf <- function(u, v) {
    pCopula(c(u, v), copula = gaussCop)
  }
  
  # Step 4: Derivatives of Copula
  partial_duij <- function(u, pik) {
    dduCopula(c(u, pik), copula = gaussCop)
  }
  
  partial_duik <- function(pij, uik) {
    ddvCopula(c(pij, uik), copula = gaussCop)
  }
  
  # Step 5: Integrals for rho^* calculation
  first_part <- 0.25 * pi[1]^2 * pi[2]^2
  
  second_term_fun_rho <- function(u) {
    0.5 * pi[2] * u[1] * partial_duij(u[1], pi[2])
  }
  second_term <- hcubature(second_term_fun_rho, lower = pi[1], upper = 1)$integral
  
  third_term_fun_rho <- function(u) {
    0.5 * pi[1] * u[1] * partial_duik(pi[1], u[1])
  }
  third_term <- hcubature(third_term_fun_rho, lower = pi[2], upper = 1)$integral
  
  double_integral_fun_rho <- function(u) {
    copula_cdf(u[1], u[2])
  }
  double_integral_term <- hcubature(double_integral_fun_rho, lowerLimit = c(pi[1], pi[2]), upperLimit = c(1, 1))$integral
  
  # Step 6: Additional terms
  pi_term <- pi[1] * pi[2] * copula_cdf(pi[1], pi[2])
  
  constant_term <- -(pi[1]^2 + pi[2]^2) / 2
  
  int_C_pi_j_uik <- function(uik) {
    pi[1] * copula_cdf(pi[1], uik)
  }
  integral_pi_j_uik_val <- hcubature(int_C_pi_j_uik, lower = pi[2], upper = 1)$integral
  
  int_C_uij_pi_k <- function(uij) {
    pi[2] * copula_cdf(uij, pi[2])
  }
  integral_uij_pi_k_val <- hcubature(int_C_uij_pi_k, lower = pi[1], upper = 1)$integral
  
  # Step 7: Final rho^* calculation
  rho_star <- 12 * (first_part + second_term + third_term + double_integral_term  + pi_term + 
                      integral_pi_j_uik_val + integral_uij_pi_k_val + constant_term) - 3
  
  return(rho_star)
}



# ---------------------------------------------------------
# Matrix calculation applying the rho() function
# ---------------------------------------------------------
# f is the tau function to apply to each element
spearmanMatrix <- function(mu, sigma, pi, r, f) {
  J <- length(mu)
  spearman_rho <- matrix(0, J, J)
  for (i in 1:(J-1)) {
    cat("i =", i, "\n")
    for (j in (i+1):J) {
      cat("\t j =", j, "\n")
      idx <- c(i,j)
      spearman_rho [i,j] <- f(mu=mu[idx], sigma=sigma[idx], pi=pi[idx], R=r[i,j])
    }
  }
  ut <- lower.tri(spearman_rho)
  spearman_rho[ut] <- t( spearman_rho)[ut]
  return(spearman_rho)
}


# ---------------------------------------------------------
# Empirical rho star calculation
# ---------------------------------------------------------

# rho star for a pair of variables from simulated data
calculate_rho_star <- function(X, Y) {
  n <- length(X)
  
  # Calculate probabilities
  p_00 <- sum(X == 0 & Y == 0) / n
  p_01 <- sum(X == 0 & Y > 0) / n
  p_10 <- sum(X > 0 & Y == 0) / n
  p_11 <- sum(X > 0 & Y > 0) / n
  p_1_plus <- p_10 + p_11  # P(X > 0)
  p_plus_1 <- p_01 + p_11  # P(Y > 0)
  
  # Compute rho_11 on non-zero values
  X_pos <- X[X > 0 & Y > 0]
  Y_pos <- Y[X > 0 & Y > 0]
  if (length(X_pos) >= 3 && length(Y_pos) >= 3) {
    rho_11 <- cor(X_pos, Y_pos, method = "spearman")
  } else {
    rho_11 <- 0
  }
  
  # Final rho star formula
  rho_star <- p_11 * p_1_plus * p_plus_1 * rho_11 + 3 * (p_00 * p_11 - p_10 * p_01)
  return(rho_star)
}

# ---------------------------------------------------------
# Matrix calculation applying rho star
# ---------------------------------------------------------

rho_star_matrix <- function(data) {
  J <- ncol(data)
  rho_mat <- matrix(0, J, J)
  
  for (i in 1:(J-1)) {
    cat("i =", i, "\n")
    for (j in (i+1):J) {
      cat("\t j =", j, "\n")
      rho_mat[i, j] <- calculate_rho_star(data[, i], data[, j])
    }
  }
  
  rho_mat[lower.tri(rho_mat)] <- t(rho_mat)[lower.tri(rho_mat)]
  diag(rho_mat) <- 1
  return(rho_mat)
}

# -----------------------------------------------
# Load simulated data and compute rarefied rho
# -----------------------------------------------

sim_data_file <- paste0("sim_data_", task_id, ".RData")
if (!file.exists(sim_data_file)) {
  stop(paste("Simulated data file", sim_data_file, "not found!"))
}
load(sim_data_file)

ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
rarefied_counts <- as(otu_table(ps_rare), "matrix")

rho_star_result <- rho_star_matrix(rarefied_counts)
print(rho_star_result)

# -----------------------------------------------
# Load fitted data and compute copula-based rho
# -----------------------------------------------

fit_data_file <- paste0("fit_data_", task_id, ".RData")
if (!file.exists(fit_data_file)) {
  stop(paste("Fitted data file", fit_data_file, "not found!"))
}
load(fit_data_file)

mu    <- fitted_simdata$EM_fit$fit$mu
pi    <- fitted_simdata$EM_fit$fit$pi0
R     <- fitted_simdata$EM_fit$fit$Sigma
sigma <- fitted_simdata$EM_fit$fit$sigma

spearman_matrix_result <- spearmanMatrix(mu, sigma, pi, R, spearman_rho)
print(spearman_matrix_result)

# -----------------------------------------------
# Load true rho matrix
# -----------------------------------------------

load("truspearman_matrix_result.RData")

# -----------------------------------------------
# Compare metrics
# -----------------------------------------------

lower_tri_new <- spearman_matrix_result [lower.tri(spearman_matrix_result)]
lower_tri_rhostar <- rho_star_result[lower.tri(rho_star_result)]
lower_tri_trurho <- truspearman_matrix_result[lower.tri(truspearman_matrix_result)]


mse_newdeveloped <- mean((lower_tri_new - lower_tri_trurho)^2)
mse_rhostar <- mean((lower_tri_rhostar - lower_tri_trurho)^2)

mad_new <- mean(abs(lower_tri_new - lower_tri_trurho))
mad_rhostar <- mean(abs(lower_tri_rhostar - lower_tri_trurho))


cat("Task", task_id, "MSE (copula-based rho):", mse_newdeveloped, "\n")
cat("Task", task_id, "MAD (copula-based rho):", mad_new, "\n")
cat("Task", task_id, "MSE (tau star):", mse_rhostar, "\n")
cat("Task", task_id, "MAD (tau star):", mad_rhostar, "\n")

save.image(file = paste0("SS_simulation_task_", task_id, ".RData"))
cat("Task", task_id, "completed.\n")
