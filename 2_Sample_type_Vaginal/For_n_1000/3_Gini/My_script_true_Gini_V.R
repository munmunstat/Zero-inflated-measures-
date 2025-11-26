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


#' Gini Correlation- $\Gamma(x, y) = \frac{\text{Cov}(x, F_{Y}(y))}{\text{Cov}(x, F_{X}(x))},$




GiniCorrelation <- function(X, Y) {
  #Empirical CDF functions
  F_X <- ecdf(X)
  F_Y <- ecdf(Y)
  
  #Covariances
  cov_X_G_Y <- cov(X, F_Y(Y))  # Covariance between X and the empirical CDF of Y
  cov_X_F_X <- cov(X, F_X(X))  # Covariance between X and the empirical CDF of X
  
  # Gini correlation calculation
  Gini_corr <- cov_X_G_Y / cov_X_F_X
  return(Gini_corr)
}


GiniMatrix_original_formula <- function(data) {
  J <- ncol(data)  # Get the number of columns (features) in the dataset
  Gini_corr_new <- matrix(0, J, J)  # Initialize the matrix
  
  # Calculate Gini correlation for each pair of columns
  for (i in 1:(J-1)) {
    for (j in (i+1):J) {
      # Use the custom GiniCorrelation function
      Gini_corr_new[i, j] <- GiniCorrelation(data[, i], data[, j])
    }
  }
  
  # Fill the lower triangular part of the matrix to make it symmetric
  Gini_corr_new[lower.tri(Gini_corr_new)] <- t(Gini_corr_new)[lower.tri(Gini_corr_new)]
  
  return(Gini_corr_new)
}


# Define the new Gini function integrating over Aij
Gini_corr_new_Aspace <- function(mu, sigma, pi, R) {
  # Create bivariate Gaussian copula
  gaussCop <- normalCopula(param = R, dim = 2, dispstr = "un")
  
  # Lower integration limit in Aij-space
  A_lower <- qlnorm(pi[1], meanlog = mu[1], sdlog = sigma[1])
  
  # Numerator part 1 integrand
  numerator_part1_new <- function(Aij, mu, sigma) {
    Aij * dlnorm(Aij, meanlog = mu, sdlog = sigma)
  }
  
  # Numerator part 2 integrand
  numerator_part2_new <- function(x, mu, sigma, gaussCop) {
    Aij <- x[1]
    uik <- x[2]
    uij <- plnorm(Aij, meanlog = mu, sdlog = sigma)
    cop <- dCopula(c(uij, uik), copula = gaussCop)
    val <- Aij * uik * cop * dlnorm(Aij, meanlog = mu, sdlog = sigma)
    return(val)
  }
  
  # Numerator part 3 integrand
  numerator_part3_new <- function(x, mu, sigma, gaussCop) {
    Aij <- x[1]
    uik <- x[2]
    uij <- plnorm(Aij, meanlog = mu, sdlog = sigma)
    cop <- dCopula(c(uij, uik), copula = gaussCop)
    val <- Aij * cop * dlnorm(Aij, meanlog = mu, sdlog = sigma)
    return(val)
  }
  
  # Denominator integrand
  denominator_integral_new <- function(Aij, mu, sigma) {
    uij <- plnorm(Aij, meanlog = mu, sdlog = sigma)
    val <- Aij * uij * dlnorm(Aij, meanlog = mu, sdlog = sigma)
    return(val)
  }
  
  # Perform all numerical integrations
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
  
  # Compute final numerator and denominator
  numerator <- 0.5 * pi[2] * integral_numerator3 +
    integral_numerator2 -
    0.5 * integral_numerator1
  
  denominator <- integral_denominator -
    ((1 - pi[1]^2) / 2) * integral_numerator1
  
  Gini_newly_derived <- numerator / denominator
  return(Gini_newly_derived)
}


# f is the Gini function to apply to each element
GiniMatrix <- function(mu, sigma, pi, r, f) {
  J <- length(mu)
  Gini_corr_new_Aspace<- matrix(0, J, J)
  for (i in 1:(J-1)) {
    cat("i =", i, "\n")
    for (j in (i+1):J) {
      cat("\t j =", j, "\n")
      idx <- c(i,j)
      Gini_corr_new_Aspace[i,j] <- f(mu=mu[idx], sigma=sigma[idx], pi=pi[idx], R=r[i,j])
    }
  }
  ut <- lower.tri(Gini_corr_new_Aspace)
  Gini_corr_new_Aspace[ut] <- t(Gini_corr_new_Aspace)[ut]
  return(Gini_corr_new_Aspace)
}


##########################
####Functions finish here 
############################

# Using the pre-trained dataset for "Stool"
Vaginal_simulation_1 <- SparseDOSSA2(template = "Vaginal",  # choose from "Stool", "Vaginal" or "IBD"
                                     new_features = TRUE,  # should new features be simulated
                                     n_sample = 1000,  # number of samples to simulate
                                     n_feature = 100,  # number of features to simulate (when 'new_features = TRUE')
                                     verbose = FALSE)  # return detailed info




Sigma_full <- Vaginal_simulation_1$template$EM_fit$fit$Sigma
pi0_full <- Vaginal_simulation_1$template$EM_fit$fit$pi0
mu_full <- Vaginal_simulation_1$template$EM_fit$fit$m
sigma_full <- Vaginal_simulation_1$template$EM_fit$fit$sigma
selected_indices <- 1:50
R <- Sigma_full[selected_indices, selected_indices]
pi <- unname(pi0_full[selected_indices])
mu <- unname(mu_full[selected_indices])
sigma <- unname(sigma_full[selected_indices])



n=1000

J <- length(mu)
# Simulate Zero-Inflated Log-Normal data using Gaussian Copula
ziln_result <- simZILN(n, mu, sigma, pi, R)

# Extract simulated data
ziln_data<-ziln_result$dat
copula_data <- ziln_result$copula.dat



V_Gini_matrix_true_developed<-GiniMatrix(mu, sigma, pi, R, Gini_corr_new_Aspace)



save(V_Gini_matrix_true_developed,file = "V_Gini_matrix_true_developed.RData")














