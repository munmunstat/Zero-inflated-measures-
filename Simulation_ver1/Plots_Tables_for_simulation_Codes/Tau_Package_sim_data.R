library(Kendall)
library(phyloseq)


#tau Stoll sample n=500


# Load required libraries
library(phyloseq)
library(Kendall)

# Load the true Kendall tau matrix
load("truetau_matrix_result.RData")  # Should load truetau_matrix_result (50x50)

# Define your Kendall tau function using cor()
kendallTauMatrix_package <- function(data) {
  J <- ncol(data)
  tau <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      tau[i, j] <- cor(data[, i], data[, j], method = "kendall")
    }
  }
  tau[lower.tri(tau)] <- t(tau)[lower.tri(tau)]
  colnames(tau) <- colnames(data)
  rownames(tau) <- colnames(data)
  return(tau)
}

# Initialize vectors to hold MSE and MAD values
mse_list <- numeric(0)
mad_list <- numeric(0)

# Initialize skip counter and log
skipped_files <- 0
skipped_file_names <- c()

# List of all simulation files
sim_files <- list.files(pattern = "^sim_data_\\d+\\.RData$")

# Loop over each simulation file
for (file in sim_files) {
  cat("Processing:", file, "\n")
  load(file)  # Loads sim_data (500 samples × 50 taxa)
  
  # Assign correct column names to match the true matrix
  colnames(sim_data) <- colnames(truetau_matrix_result)
  
  # Rarefy using phyloseq
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix names again after rarefaction (they might be lost)
  colnames(rarefied_counts) <- colnames(truetau_matrix_result)
  
  # ✅ Check if number of taxa has been reduced
  if (ncol(rarefied_counts) < ncol(truetau_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Kendall's tau matrix
  tau_package <- kendallTauMatrix_package(rarefied_counts)
  
  # Extract lower triangle of both matrices
  lower_tri_package <- tau_package[lower.tri(tau_package)]
  lower_tri_true <- truetau_matrix_result[lower.tri(truetau_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# ✅ Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)

# Print summary
cat("Average MSE:", avg_mse, "\n")
cat("Overall RMSE:", overall_rmse, "\n")
cat("Average MAD:", avg_mad, "\n")
cat("Total skipped files:", skipped_files, "\n")



# Print summary with 5 digits after decimal
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")

# Save final summary
results_summary <- data.frame(
  Avg_MSE_PackageTau = avg_mse,
  Overall_RMSE_PackageTau = overall_rmse,
  Avg_MAD_PackageTau = avg_mad,
  Total_Skipped = skipped_files
)
write.csv(results_summary, "summary_package_kendall_tau_MSE_MAD.csv", row.names = FALSE)

# Optionally save skipped file names
writeLines(skipped_file_names, "skipped_sim_data_files.txt")



##Tau Stool n=1000

# Load required libraries
library(phyloseq)
library(Kendall)

# Load the true Kendall tau matrix
load("truetau_matrix_result.RData")  # Should load truetau_matrix_result (50x50)

# Define your Kendall tau function using cor()
kendallTauMatrix_package <- function(data) {
  J <- ncol(data)
  tau <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      tau[i, j] <- cor(data[, i], data[, j], method = "kendall")
    }
  }
  tau[lower.tri(tau)] <- t(tau)[lower.tri(tau)]
  colnames(tau) <- colnames(data)
  rownames(tau) <- colnames(data)
  return(tau)
}

# Initialize vectors to hold MSE and MAD values
mse_list <- numeric(0)
mad_list <- numeric(0)

# Initialize skip counter and log
skipped_files <- 0
skipped_file_names <- c()

# List of all simulation files
sim_files <- list.files(pattern = "^sim_data6_\\d+\\.RData$")

# Loop over each simulation file
for (file in sim_files) {
  cat("Processing:", file, "\n")
  load(file)  # Loads sim_data (500 samples × 50 taxa)
  
  # Assign correct column names to match the true matrix
  colnames(sim_data) <- colnames(truetau_matrix_result)
  
  # Rarefy using phyloseq
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix names again after rarefaction (they might be lost)
  colnames(rarefied_counts) <- colnames(truetau_matrix_result)
  
  # ✅ Check if number of taxa has been reduced
  if (ncol(rarefied_counts) < ncol(truetau_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Kendall's tau matrix
  tau_package <- kendallTauMatrix_package(rarefied_counts)
  
  # Extract lower triangle of both matrices
  lower_tri_package <- tau_package[lower.tri(tau_package)]
  lower_tri_true <- truetau_matrix_result[lower.tri(truetau_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

#Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)

# Print summary
cat("Average MSE:", avg_mse, "\n")
cat("Overall RMSE:", overall_rmse, "\n")
cat("Average MAD:", avg_mad, "\n")
cat("Total skipped files:", skipped_files, "\n")



# Print summary with 5 digits after decimal
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")
# Save final summary
results_summary <- data.frame(
  Avg_MSE_PackageTau = avg_mse,
  Overall_RMSE_PackageTau = overall_rmse,
  Avg_MAD_PackageTau = avg_mad,
  Total_Skipped = skipped_files
)
write.csv(results_summary, "summary_package_tau_MSE_MAD_S1000.csv", row.names = FALSE)


#-----------------------
#tau vaginal smaple n=500
#------------------------
# Load required libraries
library(phyloseq)
library(Kendall)

# Load the true Kendall tau matrix
load("V_truetau_matrix_result.RData")  # Should load truetau_matrix_result (50x50)

# Define your Kendall tau function using cor()
kendallTauMatrix_package <- function(data) {
  J <- ncol(data)
  tau <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      tau[i, j] <- cor(data[, i], data[, j], method = "kendall")
    }
  }
  tau[lower.tri(tau)] <- t(tau)[lower.tri(tau)]
  colnames(tau) <- colnames(data)
  rownames(tau) <- colnames(data)
  return(tau)
}

# Initialize vectors to hold MSE and MAD values
mse_list <- numeric(0)
mad_list <- numeric(0)

# Initialize skip counter and log
skipped_files <- 0
skipped_file_names <- c()

# List of all simulation files
sim_files <- list.files(pattern = "^V_sim_data_\\d+\\.RData$")

# Loop over each simulation file
for (file in sim_files) {
  cat("Processing:", file, "\n")
  load(file)  # Loads sim_data (500 samples × 50 taxa)
  
  # Assign correct column names to match the true matrix
  colnames(sim_data) <- colnames(V_truetau_matrix_result)
  
  # Rarefy using phyloseq
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix names again after rarefaction (they might be lost)
  colnames(rarefied_counts) <- colnames(V_truetau_matrix_result)
  
  # ✅ Check if number of taxa has been reduced
  if (ncol(rarefied_counts) < ncol(V_truetau_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Kendall's tau matrix
  tau_package <- kendallTauMatrix_package(rarefied_counts)
  
  # Extract lower triangle of both matrices
  lower_tri_package <- tau_package[lower.tri(tau_package)]
  lower_tri_true <- V_truetau_matrix_result[lower.tri(V_truetau_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# ✅ Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)



# Print summary with 5 digits after decimal
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")

# Save final summary
results_summary <- data.frame(
  Avg_MSE_PackageTau = avg_mse,
  Overall_RMSE_PackageTau = overall_rmse,
  Avg_MAD_PackageTau = avg_mad,
  Total_Skipped = skipped_files
)
write.csv(results_summary, "summary_package_kendall_tau_MSE_MAD.csv", row.names = FALSE)

# Optionally save skipped file names
writeLines(skipped_file_names, "skipped_sim_data_files.txt")






#-----------------------
#tau vaginal smaple n=1000
#------------------------
# Load required libraries
library(phyloseq)
library(Kendall)

# Load the true Kendall tau matrix
load("V_truetau_matrix_result.RData")  # Should load truetau_matrix_result (50x50)

# Define your Kendall tau function using cor()
kendallTauMatrix_package <- function(data) {
  J <- ncol(data)
  tau <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      tau[i, j] <- cor(data[, i], data[, j], method = "kendall")
    }
  }
  tau[lower.tri(tau)] <- t(tau)[lower.tri(tau)]
  colnames(tau) <- colnames(data)
  rownames(tau) <- colnames(data)
  return(tau)
}

# Initialize vectors to hold MSE and MAD values
mse_list <- numeric(0)
mad_list <- numeric(0)

# Initialize skip counter and log
skipped_files <- 0
skipped_file_names <- c()

# List of all simulation files
sim_files <- list.files(pattern = "^V2_sim_data_\\d+\\.RData$")

# Loop over each simulation file
for (file in sim_files) {
  cat("Processing:", file, "\n")
  load(file)  # Loads sim_data (500 samples × 50 taxa)
  
  # Assign correct column names to match the true matrix
  colnames(sim_data) <- colnames(V_truetau_matrix_result)
  
  # Rarefy using phyloseq
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix names again after rarefaction (they might be lost)
  colnames(rarefied_counts) <- colnames(V_truetau_matrix_result)
  
  # ✅ Check if number of taxa has been reduced
  if (ncol(rarefied_counts) < ncol(V_truetau_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Kendall's tau matrix
  tau_package <- kendallTauMatrix_package(rarefied_counts)
  
  # Extract lower triangle of both matrices
  lower_tri_package <- tau_package[lower.tri(tau_package)]
  lower_tri_true <- V_truetau_matrix_result[lower.tri(V_truetau_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# ✅ Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)



# Print summary with 5 digits after decimal
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")

# Save final summary
results_summary <- data.frame(
  Avg_MSE_PackageTau = avg_mse,
  Overall_RMSE_PackageTau = overall_rmse,
  Avg_MAD_PackageTau = avg_mad,
  Total_Skipped = skipped_files
)
write.csv(results_summary, "summary_package_tau_MSE_MAD_V1000.csv", row.names = FALSE)

# Optionally save skipped file names
writeLines(skipped_file_names, "skipped_sim_data_files.txt")




#-----------------------
#tau IBD smaple n=500
#------------------------
# Load required libraries
library(phyloseq)
library(Kendall)

# Load the true Kendall tau matrix
load("IBD_truetau_matrix_result.RData")  # Should load truetau_matrix_result (50x50)

# Define your Kendall tau function using cor()
kendallTauMatrix_package <- function(data) {
  J <- ncol(data)
  tau <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      tau[i, j] <- cor(data[, i], data[, j], method = "kendall")
    }
  }
  tau[lower.tri(tau)] <- t(tau)[lower.tri(tau)]
  colnames(tau) <- colnames(data)
  rownames(tau) <- colnames(data)
  return(tau)
}

# Initialize vectors to hold MSE and MAD values
mse_list <- numeric(0)
mad_list <- numeric(0)

# Initialize skip counter and log
skipped_files <- 0
skipped_file_names <- c()

# List of all simulation files
sim_files <- list.files(pattern = "^IBD_sim_data_\\d+\\.RData$")

# Loop over each simulation file
for (file in sim_files) {
  cat("Processing:", file, "\n")
  load(file)  # Loads sim_data (500 samples × 50 taxa)
  
  # Assign correct column names to match the true matrix
  colnames(sim_data) <- colnames(IBD_truetau_matrix_result)
  
  # Rarefy using phyloseq
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix names again after rarefaction (they might be lost)
  colnames(rarefied_counts) <- colnames(IBD_truetau_matrix_result)
  
  # ✅ Check if number of taxa has been reduced
  if (ncol(rarefied_counts) < ncol(IBD_truetau_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Kendall's tau matrix
  tau_package <- kendallTauMatrix_package(rarefied_counts)
  
  # Extract lower triangle of both matrices
  lower_tri_package <- tau_package[lower.tri(tau_package)]
  lower_tri_true <- IBD_truetau_matrix_result[lower.tri(IBD_truetau_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# ✅ Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)



# Print summary with 5 digits after decimal
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")

# Save final summary
results_summary <- data.frame(
  Avg_MSE_PackageTau = avg_mse,
  Overall_RMSE_PackageTau = overall_rmse,
  Avg_MAD_PackageTau = avg_mad,
  Total_Skipped = skipped_files
)
write.csv(results_summary, "summary_package_kendall_tau_MSE_MAD.csv", row.names = FALSE)

# Optionally save skipped file names
writeLines(skipped_file_names, "skipped_sim_data_files.txt")


#-----------------------
#tau IBD smaple n=1000
#------------------------
# Load required libraries
library(phyloseq)
library(Kendall)

# Load the true Kendall tau matrix
load("IBD_truetau_matrix_result.RData")  # Should load truetau_matrix_result (50x50)

# Define your Kendall tau function using cor()
kendallTauMatrix_package <- function(data) {
  J <- ncol(data)
  tau <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      tau[i, j] <- cor(data[, i], data[, j], method = "kendall")
    }
  }
  tau[lower.tri(tau)] <- t(tau)[lower.tri(tau)]
  colnames(tau) <- colnames(data)
  rownames(tau) <- colnames(data)
  return(tau)
}

# Initialize vectors to hold MSE and MAD values
mse_list <- numeric(0)
mad_list <- numeric(0)

# Initialize skip counter and log
skipped_files <- 0
skipped_file_names <- c()

# List of all simulation files
sim_files <- list.files(pattern = "^IBD2_sim_data_\\d+\\.RData$")

# Loop over each simulation file
for (file in sim_files) {
  cat("Processing:", file, "\n")
  load(file)  # Loads sim_data (500 samples × 50 taxa)
  
  # Assign correct column names to match the true matrix
  colnames(sim_data) <- colnames(IBD_truetau_matrix_result)
  
  # Rarefy using phyloseq
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix names again after rarefaction (they might be lost)
  colnames(rarefied_counts) <- colnames(IBD_truetau_matrix_result)
  
  # ✅ Check if number of taxa has been reduced
  if (ncol(rarefied_counts) < ncol(IBD_truetau_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Kendall's tau matrix
  tau_package <- kendallTauMatrix_package(rarefied_counts)
  
  # Extract lower triangle of both matrices
  lower_tri_package <- tau_package[lower.tri(tau_package)]
  lower_tri_true <- IBD_truetau_matrix_result[lower.tri(IBD_truetau_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# ✅ Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)



# Print summary with 5 digits after decimal
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")

# Save final summary
results_summary <- data.frame(
  Avg_MSE_PackageTau = avg_mse,
  Overall_RMSE_PackageTau = overall_rmse,
  Avg_MAD_PackageTau = avg_mad,
  Total_Skipped = skipped_files
)
write.csv(results_summary, "summary_package_tau_MSE_MAD_IBD1000.csv", row.names = FALSE)

# Optionally save skipped file names
writeLines(skipped_file_names, "skipped_sim_data_files.txt")



