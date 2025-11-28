library(copula)
library(cubature)
library(MASS)
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Kendall)
library(phyloseq)
library(VC2copula)




# ---------------------------------------------------------
# Read Data from AGP file and fiting SParrsdossa
# ---------------------------------------------------------


# Load file
otu_small_df <- read.csv("otu_small.csv", row.names = 1)

# Subset first 50 taxa (rows) and first 500 samples (columns)
otu_small_subset <- otu_small_df[1:50, 1:500]

# Check dimensions
dim(otu_small_subset)
# Should print [1] 50 500

# Fit SparseDOSSA2
fitted_AGP <- fit_SparseDOSSA2(
  data = otu_small_subset,
  lambda = 0,
  control = list(verbose = TRUE)
)

# Save fitted object
save(fitted_AGP, file = "fitted_AGP.RData")



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

# ---------------------------------------------------------
# Newly developed Rho calculation
# ---------------------------------------------------------

# Define the function to compute Spearman's Rho including rho^* calculation without a separate function
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
# Matrix Rho calculation
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
# Empirical tau rho calculation
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
# Empirical Matrix rho star calculation
# ---------------------------------------------------------

# Matrix version: compute rho_star for all variable pairs in a data matrix
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

# -----------------------------------------------
# Gini Matrix
# -----------------------------------------------

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


### Extract parameter from fitted data
mu    <- fitted_AGP$EM_fit$fit$mu
pi    <- fitted_AGP$EM_fit$fit$pi0
R     <- fitted_AGP$EM_fit$fit$Sigma
sigma <- fitted_AGP$EM_fit$fit$sigma

# -------------------------------------------------------
# Compute AGP tau matrix newly developed
# -------------------------------------------------------

AGP_tau_matrix_result <- tauMatrix(mu, sigma, pi, R, tau)
save(AGP_tau_matrix_result, file = "AGP_tau_matrix_result.RData")
write.csv(AGP_tau_matrix_result, "AGP_tau_matrix_result.csv")

# -------------------------------------------------------
# Compute AGP empirical taustar
# -----------------------------------------------------

# Create phyloseq object
ps <- phyloseq(
  otu_table(otu_small_subset, taxa_are_rows = TRUE)
)

# Rarefy to even depth
ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)

# Convert rarefied counts to matrix
rarefied_counts <- as(otu_table(ps_rare), "matrix")

# Compute Kendall's tau*
AGP_taustar_matrix_result <- kendallTaustar(rarefied_counts)

save(AGP_taustar_matrix_result, file = "AGP_taustar_matrix_result.RData")
write.csv(AGP_taustar_matrix_result, "AGP_taustar_matrix_result.csv")

# -------------------------------------------------------
# Compute newly developed Spearman matrix
# -------------------------------------------------------

AGP_spearman_matrix_result <- spearmanMatrix(mu, sigma, pi, R, spearman_rho)

save(AGP_spearman_matrix_result, file = "AGP_spearman_matrix_result.RData")
write.csv(AGP_spearman_matrix_result, "AGP_spearman_matrix_result.csv")

# -------------------------------------------------------
# Compute rho_star matrix on rarefied counts
# -------------------------------------------------------

AGP_rho_star_result <- rho_star_matrix(rarefied_counts)

save(AGP_rho_star_result, file = "AGP_rho_star_result.RData")
write.csv(AGP_rho_star_result, "AGP_rho_star_result.csv")

# -------------------------------------------------------
# Compute newly developed Gini matrix
# -------------------------------------------------------

AGP_Gini_matrix_newly_developed_result <- GiniMatrix(
  mu = mu,
  sigma = sigma,
  pi = pi,
  r = R,
  f = Gini_corr_new_Aspace
)

save(AGP_Gini_matrix_newly_developed_result, file = "AGP_Gini_matrix_newly_developed_result.RData")
write.csv(AGP_Gini_matrix_newly_developed_result, "AGP_Gini_matrix_newly_developed_result.csv")

# -------------------------------------------------------
# Compute Gini matrix using original formula on rarefied data
# -------------------------------------------------------

AGP_Gini_matrix_original <- GiniMatrix_original_formula(rarefied_counts)

save(AGP_Gini_matrix_original, file = "AGP_Gini_matrix_original.RData")
write.csv(AGP_Gini_matrix_original, "AGP_Gini_matrix_original.csv")

