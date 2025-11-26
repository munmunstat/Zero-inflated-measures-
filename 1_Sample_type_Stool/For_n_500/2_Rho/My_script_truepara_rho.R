library(VC2copula)
library(copula)
library(MASS)
library(Kendall)
library(cubature)
# Function to simulate Zero-Inflated Log-Normal data using Gaussian copula
simZILN <- function(n, mu, sigma, pi, R) {
  J <- length(mu)
  
  # Check if input lengths are consistent
  if (J != length(pi) | J != length(sigma))
    stop("mu, sigma, and pi must all have the same length")
  
  # Generate copula data using Gaussian Copula
  gCop <- normalCopula(param = P2p(R), dim = J, dispstr = "un")
  z <- rCopula(n, gCop)
  
  # Initialize the matrix for Zero-Inflated Log-Normal data
  r <- matrix(0, nrow(z), J)
  
  # Iterate over each dimension
  for (j in 1:J) {
    # Apply Zero-Inflated Log-Normal transformation
    r[, j] <- qzilnorm(z[, j], mu[j], sigma[j], pi[j])
  }
  
  # Return simulated data and copula data
  return(list(dat = r, copula.dat = z))
}
# Quantile function of the Zero-Inflated Log-Normal with direct scaling logic
qzilnorm <- function(u, meanlog, sdlog, pi) {
  # Initialize a result vector
  r <- rep(0, length(u))
  
  # Logic to determine valid values where u > pi
  valid_u <- u > pi
  
  # Scale u values (only for valid_u) and assign 0 where u <= pi
  scaled_u <- (u - pi) / (1 - pi)
  scaled_u[!valid_u] <- 0  # For u <= pi, assign zero
  
  # Apply qlnorm for only valid scaled values, else assign zero
  r <- ifelse(valid_u, qlnorm(scaled_u, meanlog = meanlog, sdlog = sdlog), 0)
  
  return(r)  # Return the result
}



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

# New tauMatrix function for Kendall's correlation using cor()
spearmanMatrix_package <- function(data) {
  J <- ncol(data)  # Get the number of columns (features) in the dataset
  spearman_rho <- matrix(0, J, J)  # Initialize the matrix
  
  # Calculate Kendall's Tau for each pair of columns
  for (i in 1:(J-1)) {
    for (j in (i+1):J) {
      # Use Kendall correlation from the cor() function
      spearman_rho[i,j] <- cor(data[, i], data[, j], method = "spearman")
    }
  }
  
  # Fill the lower triangular part of the matrix to make it symmetric
  spearman_rho[lower.tri(spearman_rho)] <- t(spearman_rho)[lower.tri(spearman_rho)]
  
  return(spearman_rho)
}


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

#############
#####Function Finish here



##########
#Fitting true parameters value 
###############
# Using the pre-trained dataset for "Stool"
Stool_simulation_1 <- SparseDOSSA2(template = "Stool",  # choose from "Stool", "Vaginal" or "IBD"
                                   new_features = TRUE,  # should new features be simulated
                                   n_sample = 1000,  # number of samples to simulate
                                   n_feature = 100,  # number of features to simulate (when 'new_features = TRUE')
                                   verbose = FALSE)  # return detailed info



Sigma_full <- Stool_simulation_1$template$EM_fit$fit$Sigma
pi0_full <- Stool_simulation_1$template$EM_fit$fit$pi0
mu_full<-Stool_simulation_1$template$EM_fit$fit$mu
sigma_full<-Stool_simulation_1$template$EM_fit$fit$sigma
length(sigma_full)
# Adjust this vector as needed to select other features
selected_indices <- 1:50
R<- Sigma_full[selected_indices, selected_indices]
pi<-unname(pi0_full[selected_indices])

mu<-unname(mu_full[selected_indices])
sigma<-unname(sigma_full[selected_indices])

J <- length(mu)
dim(R)
n <- 1000
# Simulate Zero-Inflated Log-Normal data using Gaussian Copula
ziln_result <- simZILN(n, mu, sigma, pi, R)

# Extract simulated data
ziln_data<-ziln_result$dat
copula_data <- ziln_result$copula.dat



# Use the spearmanMatrix function to calculate Tau for all pairs
truspearman_matrix_result <- spearmanMatrix(mu, sigma, pi, R, spearman_rho)
print(spearman_matrix_result)
save(truspearman_matrix_result, file = "truspearman_matrix_result.RData")



# Use the new spearmanMatrix function using package
truspearman_rho_matrix_result_package <- spearmanMatrix_package(ziln_data)
print(spearman_rho_matrix_result_package)
save(truspearman_rho_matrix_result_package , file = "truspearman_rho_matrix_result_package.RData")


# Compute rho_star matrix
truerho_star_result <- rho_star_matrix(ziln_data)
print(rho_star_result)
save(truerho_star_result, file = "truerho_star_result.RData")





