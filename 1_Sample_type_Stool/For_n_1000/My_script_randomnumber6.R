# 1000 Random Number generation 
# Load necessary libraries
library(copula)
library(MASS)
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)


# Accept the job array index from Slurm
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])  # Get the array index

# Set a unique seed for reproducibility
set.seed(123 + task_id)

# Define your simulation function
simulate_multinomial_data <- function(n, pi, mu, sigma, R, size_meanlog = 11.9, size_sdlog = 0.47) {
  p <- length(pi)
  if (!(length(mu) == p && length(sigma) == p && ncol(R) == p && nrow(R) == p)) {
    stop("Dimensions of pi, mu, sigma, and R must match.")
  }
  param <- P2p(R)
  gCop <- normalCopula(param = param, dim = p, dispstr = "un")
  u <- rCopula(n, gCop)
  ziln_data <- matrix(0, nrow = n, ncol = p)
  for (j in 1:p) {
    valid_u <- u[, j] > pi[j]
    scaled_u <- (u[, j] - pi[j]) / (1 - pi[j])
    scaled_u[!valid_u] <- 0
    ziln_data[, j] <- ifelse(valid_u, qlnorm(scaled_u, meanlog = mu[j], sdlog = sigma[j]), 0)
  }
  Relative_ziln <- t(apply(ziln_data, 1, function(i) {
    if (sum(i) == 0) {
      return(c(rep(0, (length(i) - 1)), 1))
    } else {
      return(i / sum(i))
    }
  }))
  random_numbers <- rlnorm(n = n, meanlog = size_meanlog, sdlog = size_sdlog)
  size_vector <- round(pmax(1, pmin(random_numbers, 100000)))
  multinomial_data <- t(sapply(1:n, function(i) {
    rmultinom(n = 1, size = size_vector[i], prob = Relative_ziln[i, ])
  }))
  return(multinomial_data)
}


library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)


# Using the pre-trained dataset for "Stool"
Stool_simulation_1 <- SparseDOSSA2(template = "Stool",  # choose from "Stool", "Vaginal" or "IBD"
                                   new_features = TRUE,  # should new features be simulated
                                   n_sample = 1000,  # number of samples to simulate
                                   n_feature = 100,  # number of features to simulate (when 'new_features = TRUE')
                                   verbose = FALSE)  # return detailed info



#save(Stool_simulation_1 , file = "Stool_simulation_1_sparseD.RData")

Sigma_full <- Stool_simulation_1$template$EM_fit$fit$Sigma
pi0_full <- Stool_simulation_1$template$EM_fit$fit$pi0
mu_full <- Stool_simulation_1$template$EM_fit$fit$mu
sigma_full <- Stool_simulation_1$template$EM_fit$fit$sigma
selected_indices <- 1:50
R <- Sigma_full[selected_indices, selected_indices]
pi <- unname(pi0_full[selected_indices])
mu <- unname(mu_full[selected_indices])
sigma <- unname(sigma_full[selected_indices])

# Run the simulation
sim_data <- simulate_multinomial_data(n = 1000, pi = pi, mu = mu, sigma = sigma, R = R)

# Save the output with a unique file name using the task_id
save(sim_data,file=paste0("sim_data6_", task_id, ".RData"))
