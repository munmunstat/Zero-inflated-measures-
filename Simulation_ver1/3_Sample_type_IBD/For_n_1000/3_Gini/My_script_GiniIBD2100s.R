##########################################################
# SLURM JOB SCRIPT: Gini Correlation IBD Sample n=1000
##########################################################

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

# -----------------------------------------------
# Gini Function Definitions
# -----------------------------------------------

GiniCorrelation <- function(X, Y) {
  F_X <- ecdf(X)
  F_Y <- ecdf(Y)
  cov_X_G_Y <- cov(X, F_Y(Y))
  cov_X_F_X <- cov(X, F_X(X))
  Gini_corr <- cov_X_G_Y / cov_X_F_X
  return(Gini_corr)
}

GiniMatrix_original_formula <- function(data) {
  J <- ncol(data)
  Gini_corr_new <- matrix(0, J, J)
  for (i in 1:(J-1)) {
    for (j in (i+1):J) {
      Gini_corr_new[i, j] <- GiniCorrelation(data[, i], data[, j])
    }
  }
  Gini_corr_new[lower.tri(Gini_corr_new)] <- t(Gini_corr_new)[lower.tri(Gini_corr_new)]
  return(Gini_corr_new)
}

# ---------------------------------------------------------
# Define the new Gini function integrating over Aij-space
# ---------------------------------------------------------

Gini_corr_new_Aspace <- function(mu, sigma, pi, R) {
  gaussCop <- normalCopula(param = R, dim = 2, dispstr = "un")
  
  A_lower <- qlnorm(pi[1], meanlog = mu[1], sdlog = sigma[1])
  
  numerator_part1_new <- function(Aij, mu, sigma) {
    Aij * dlnorm(Aij, meanlog = mu, sdlog = sigma)
  }
  
  numerator_part2_new <- function(x, mu, sigma, gaussCop) {
    Aij <- x[1]
    uik <- x[2]
    uij <- plnorm(Aij, meanlog = mu, sdlog = sigma)
    cop <- dCopula(c(uij, uik), copula = gaussCop)
    val <- Aij * uik * cop * dlnorm(Aij, meanlog = mu, sdlog = sigma)
    return(val)
  }
  
  numerator_part3_new <- function(x, mu, sigma, gaussCop) {
    Aij <- x[1]
    uik <- x[2]
    uij <- plnorm(Aij, meanlog = mu, sdlog = sigma)
    cop <- dCopula(c(uij, uik), copula = gaussCop)
    val <- Aij * cop * dlnorm(Aij, meanlog = mu, sdlog = sigma)
    return(val)
  }
  
  denominator_integral_new <- function(Aij, mu, sigma) {
    uij <- plnorm(Aij, meanlog = mu, sdlog = sigma)
    val <- Aij * uij * dlnorm(Aij, meanlog = mu, sdlog = sigma)
    return(val)
  }
  
  integral_numerator1 <- hcubature(
    f = function(a) numerator_part1_new(a, mu[1], sigma[1]),
    lowerLimit = A_lower,
    upperLimit = Inf
  )$integral
  
  integral_numerator2 <- hcubature(
    f = function(x) numerator_part2_new(x, mu[1], sigma[1], gaussCop),
    lowerLimit = c(A_lower, pi[2]),
    upperLimit = c(Inf, 1)
  )$integral
  
  integral_numerator3 <- hcubature(
    f = function(x) numerator_part3_new(x, mu[1], sigma[1], gaussCop),
    lowerLimit = c(A_lower, 0),
    upperLimit = c(Inf, pi[2])
  )$integral
  
  integral_denominator <- hcubature(
    f = function(a) denominator_integral_new(a, mu[1], sigma[1]),
    lowerLimit = A_lower,
    upperLimit = Inf
  )$integral
  
  numerator <- 0.5 * pi[2] * integral_numerator3 +
    integral_numerator2 -
    0.5 * integral_numerator1
  
  denominator <- integral_denominator -
    ((1 - pi[1]^2) / 2) * integral_numerator1
  
  Gini_newly_derived <- numerator / denominator
  return(Gini_newly_derived)
}

# ---------------------------------------------------------
# Matrix calculation applying the Gini_corr_new_Aspace
# ---------------------------------------------------------

GiniMatrix <- function(mu, sigma, pi, r, f) {
  J <- length(mu)
  Gini_corr_new_Aspace_mat <- matrix(0, J, J)
  for (i in 1:(J-1)) {
    cat("i =", i, "\n")
    for (j in (i+1):J) {
      cat("\t j =", j, "\n")
      idx <- c(i, j)
      Gini_corr_new_Aspace_mat[i, j] <- f(
        mu = mu[idx],
        sigma = sigma[idx],
        pi = pi[idx],
        R = r[i, j]
      )
    }
  }
  ut <- lower.tri(Gini_corr_new_Aspace_mat)
  Gini_corr_new_Aspace_mat[ut] <- t(Gini_corr_new_Aspace_mat)[ut]
  return(Gini_corr_new_Aspace_mat)
}

# -----------------------------------------------
# Load simulated data and compute rarefied Gini
# -----------------------------------------------

sim_data_file <- paste0("IBD2_sim_data_", task_id, ".RData")

if (!file.exists(sim_data_file)) {
  stop(paste("Simulated data file", sim_data_file, "not found!"))
}

load(sim_data_file)

ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
set.seed(123)
ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
rarefied_counts <- as(otu_table(ps_rare), "matrix")
Gini_matrix_original <- GiniMatrix_original_formula(rarefied_counts)

# -----------------------------------------------
# Load fitted data and compute new Gini
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

Gini_matrix_newly_developed_result <- GiniMatrix(mu = mu, sigma = sigma,pi = pi, r = R, f = Gini_corr_new_Aspace)

# -----------------------------------------------
# Load True Gini Matrix
# -----------------------------------------------

load("IBD_Gini_matrix_true_developed.RData")

# -----------------------------------------------
# Compare metrics
# -----------------------------------------------

lower_tri_new <- Gini_matrix_newly_developed_result[lower.tri(Gini_matrix_newly_developed_result)]
lower_tri_original <- Gini_matrix_original[lower.tri(Gini_matrix_original)]
lower_tri_truGini <- IBD_Gini_matrix_true_developed[lower.tri(IBD_Gini_matrix_true_developed)]

mse_Gini_develope <- mean((lower_tri_new - lower_tri_truGini)^2)
mad_Gini_develope <- mean(abs(lower_tri_new - lower_tri_truGini))

mse_Gini_original <- mean((lower_tri_original - lower_tri_truGini)^2)
mad_Gini_original <- mean(abs(lower_tri_original - lower_tri_truGini))

cat("Task", task_id, "MSE (developed Gini):", mse_Gini_develope, "\n")
cat("Task", task_id, "MAD (developed Gini):", mad_Gini_develope, "\n")
cat("Task", task_id, "MSE (original Gini):", mse_Gini_original, "\n")
cat("Task", task_id, "MAD (original Gini):", mad_Gini_original, "\n")

# Save all results
save.image(file = paste0("IBD2_Gini_simulation_task_", task_id, ".RData"))
cat("Task", task_id, "completed.\n")
