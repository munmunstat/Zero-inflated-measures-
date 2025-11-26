#for vaginal sample finding value of true parameter values
library(copula)
library(cubature)
library(MASS)
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Kendall)


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


# Define the function to compute Kendall's Tau
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
  fourth_term_integral <- hcubature(fourth_term_fun, lowerLimit = c(pi[1], pi[2]), upperLimit = c(1, 1))$integral
  
  # Calculate the total expected value (Kendall's Tau)
  total_expectation <- copula_cdf(pi[1], pi[2])^2 + second_term_integral + third_term_integral + fourth_term_integral
  
  # Return the result for Kendall's Tau
  return(total_expectation - 1)
}



# f is the tau function to apply to each element
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


# New tauMatrix function for Kendall's correlation using cor()
kendallTauMatrix_package <- function(data) {
  J <- ncol(data)  # Get the number of columns (features) in the dataset
  tau <- matrix(0, J, J)  # Initialize the matrix
  
  # Calculate Kendall's Tau for each pair of columns
  for (i in 1:(J-1)) {
    for (j in (i+1):J) {
      # Use Kendall correlation from the cor() function
      tau[i,j] <- cor(data[, i], data[, j], method = "kendall")
    }
  }
  
  # Fill the lower triangular part of the matrix to make it symmetric
  tau[lower.tri(tau)] <- t(tau)[lower.tri(tau)]
  
  return(tau)
}

################################
#Step 5:Using Empirical formula
################################
calculate_tau_star <- function(X, Y) {
  # Step 3: Calculate probabilities
  n <- length(X)  # Assuming X and Y are the same length
  
  p_00 <- sum(X == 0 & Y == 0) / n
  p_01 <- sum(X == 0 & Y > 0) / n
  p_10 <- sum(X > 0 & Y == 0) / n
  p_11 <- sum(X > 0 & Y > 0) / n
  
  # Step 4: Filter data for Kendall's tau
  X_pos <- X[X > 0 & Y > 0]
  Y_pos <- Y[X > 0 & Y > 0]
  
  
  
  # Check if there are at least 2 unique non-identical pairs
  if (length(X_pos) >= 3 && length(Y_pos) >= 3) {
    tau_11 <- Kendall(X_pos, Y_pos)$tau[1]
  } else {
    tau_11 <- 0  # Assign a default value if there are not enough pairs
    
  }
  # Calculate Kendall's tau for positive observations
  #tau_11 <- Kendall(X_pos, Y_pos)$tau[1]
  
  # Print tau_11 to track its value during execution
  cat("tau_11 for this pair:", tau_11, "\n")
  
  # Step 5: Calculate p1
  # X10: X where Y == 0
  # X11: X where Y > 0
  X10 <- X[Y == 0]
  X11 <- X[Y > 0]
  
  # Count pairs (X10, X11) such that X10 > X11
  
  count_1 <- 0
  for (i in seq_along(X10)) {
    count_1=count_1+sum(X11<X10[i])
    
  }
  
  # Probability p1
  p1 <- count_1 / (length(X11)*length(X10))
  
  # Step 6: Calculate p2
  #Y01: Y where X == 0
  #Y11: Y where X > 0
  Y01 <- Y[X == 0]
  Y11 <- Y[X > 0]
  
  # Count pairs (Y01, Y11) such that Y01 > Y11
  # Count pairs (Y01, Y11) such that Y01 > Y11
  count_2 <- 0
  for (i in seq_along(Y01)) {
    count_2=count_2+sum(Y11<Y01[i])
    
  }
  
  
  # Probability p2
  p2 <- count_2 / (length(Y01)*length(Y11))
  
  # Step 7: Calculate tau*
  tau_star <- p_11^2 * tau_11 + 
    2 * (p_00 * p_11 - p_01 * p_10) + 
    2 * p_11 * (p_10 * (1 - 2 * p1) + p_01 * (1 - 2 * p2))
  
  #Return results
  return(tau_star)
}



# New tauMatrix function for Kendall's correlation using cor()
kendallTaustar<- function(data) {
  J <- ncol(data)  # Get the number of columns (features) in the dataset
  tau_starmatrix  <- matrix(0, J, J)  # Initialize the matrix
  
  # Calculate Kendall's Tau for each pair of columns
  for (i in 1:(J-1)) {
    for (j in (i+1):J) {
      # Use Kendall correlation from the cor() function
      tau_starmatrix[i,j] <- calculate_tau_star(data[, i], data[, j])
    }
  }
  
  # Fill the lower triangular part of the matrix to make it symmetric
  tau_starmatrix [lower.tri(tau_starmatrix )] <- t(tau_starmatrix )[lower.tri(tau_starmatrix)]
  
  return(tau_starmatrix)
}

####Functions finish here 
####



##########
#Fitting true parameters value 
###############

load("Vaginal_simulation_1_sparseD.RData")

Sigma_full <- Vaginal_simulation_1$template$EM_fit$fit$Sigma
pi0_full <- Vaginal_simulation_1$template$EM_fit$fit$pi0
mu_full<-Vaginal_simulation_1$template$EM_fit$fit$mu
sigma_full<-Vaginal_simulation_1$template$EM_fit$fit$sigma
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



library(Kendall)
# Testing CDF of positive part
idx <- c(1,2)

# Regular Kendall's tau
cor(ziln_data[,idx], method="kendall")


########
# Use the tauMatrix function to calculate Tau for all pairs
########

V_truetau_matrix_result <- tauMatrix(mu, sigma, pi, R, tau)
print(V_truetau_matrix_result)
save(V_truetau_matrix_result, file = "V_truetau_matrix_result.RData")


V_truetau_star<-calculate_tau_star(X=ziln_data[, 1],Y=ziln_data[, 2])
print(V_truetau_star)
save(V_truetau_star, file = "V_truetau_star.RData")

V_truetaustar_matrix_result<-kendallTaustar(ziln_data)
print(V_truetaustar_matrix_result)
save(V_truetaustar_matrix_result, file = "V_truetaustar_matrix_result.RData")



