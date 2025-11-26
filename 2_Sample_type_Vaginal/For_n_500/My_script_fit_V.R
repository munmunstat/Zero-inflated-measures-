
library(copula)
library(cubature)
library(MASS)
library(SparseDOSSA2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Kendall)


# ---- Accept job array index from Slurm ----
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

if (is.na(task_id)) {
  stop("No task_id passed to the script!")
}

# ---- Load simulated data file ----
sim_data_file <- paste0("V_sim_data_", task_id, ".RData")

if (!file.exists(sim_data_file)) {
  stop(paste("Simulated data file", sim_data_file, "not found!"))
}

load(sim_data_file)  # This loads sim_data

if (!exists("sim_data")) {
  stop("sim_data variable not found in loaded file")
}

# ---- Fit SparseDOSSA2 ----
fitted_simdata <- fit_SparseDOSSA2(
  data = t(sim_data),
  lambda = 0,
  control = list(verbose = TRUE)
)

# ---- Save fitted result ----
fit_file <- paste0("fit_V_data_", task_id, ".RData")
save(fitted_simdata, file = fit_file)

cat("Done processing simulation", task_id, "\n")
cat("Saved result to", fit_file, "\n")
