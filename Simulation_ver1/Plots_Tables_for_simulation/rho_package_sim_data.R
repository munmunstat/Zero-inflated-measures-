# Load required libraries

#Rho stool sample n=500

library(phyloseq)
library(Kendall)
# Load the true Spearman rho matrix
load("truspearman_matrix_result.RData")  # Should load truespearman_matrix_result (50x50)
ls()
# Define Spearman rho function using cor()
spearmanMatrix_package <- function(data) {
  J <- ncol(data)
  rho <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      rho[i, j] <- cor(data[, i], data[, j], method = "spearman")
    }
  }
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  colnames(rho) <- colnames(data)
  rownames(rho) <- colnames(data)
  return(rho)
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
  load(file)  # Should load sim_data (500 × 50 matrix)
  
  # Assign correct column names
  colnames(sim_data) <- colnames(truspearman_matrix_result)
  
  # Rarefy
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix column names again (may be lost in rarefaction)
  colnames(rarefied_counts) <- colnames(truspearman_matrix_result)
  
  # Skip if taxa reduced
  if (ncol(rarefied_counts) < ncol(truspearman_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Spearman rho matrix
  rho_package <- spearmanMatrix_package(rarefied_counts)
  
  # Extract lower triangles
  lower_tri_package <- rho_package[lower.tri(rho_package)]
  lower_tri_true <- truspearman_matrix_result[lower.tri(truspearman_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)


# Print summary (5 decimal digits)
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")




# Prepare summary as data frame
summary_df <- data.frame(
  Method = "Spearman_Package",
  Average_MSE = round(avg_mse, 5),
  Overall_RMSE = round(overall_rmse, 5),
  Average_MAD = round(avg_mad, 5),
  Skipped_Files = skipped_files
)

# Save to CSV
write.csv(summary_df, file = "spearman_package_summary_S500.csv", row.names = FALSE)




##Rho Stool n=1000
# Load required libraries
library(phyloseq)
library(Kendall)
# Load the true Spearman rho matrix
load("truspearman_matrix_result.RData")  # Should load truespearman_matrix_result (50x50)

# Define Spearman rho function using cor()
spearmanMatrix_package <- function(data) {
  J <- ncol(data)
  rho <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      rho[i, j] <- cor(data[, i], data[, j], method = "spearman")
    }
  }
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  colnames(rho) <- colnames(data)
  rownames(rho) <- colnames(data)
  return(rho)
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
  load(file)  # Should load sim_data (500 × 50 matrix)
  
  # Assign correct column names
  colnames(sim_data) <- colnames(truspearman_matrix_result)
  
  # Rarefy
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix column names again (may be lost in rarefaction)
  colnames(rarefied_counts) <- colnames(truspearman_matrix_result)
  
  # Skip if taxa reduced
  if (ncol(rarefied_counts) < ncol(truspearman_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Spearman rho matrix
  rho_package <- spearmanMatrix_package(rarefied_counts)
  
  # Extract lower triangles
  lower_tri_package <- rho_package[lower.tri(rho_package)]
  lower_tri_true <- truspearman_matrix_result[lower.tri(truspearman_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)

# Print summary (5 decimal digits)
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")



# Prepare summary as data frame
summary_df <- data.frame(
  Method = "Spearman_Package",
  Average_MSE = round(avg_mse, 5),
  Overall_RMSE = round(overall_rmse, 5),
  Average_MAD = round(avg_mad, 5),
  Skipped_Files = skipped_files
)

# Save to CSV
write.csv(summary_df, file = "spearman_package_summary_S1000.csv", row.names = FALSE)




#spearman Vaginal n=500

# Load required libraries
library(phyloseq)
library(Kendall)
# Load the true Spearman rho matrix
load("V_truspearman_matrix_result.RData")  # Should load truespearman_matrix_result (50x50)
ls()
# Define Spearman rho function using cor()
spearmanMatrix_package <- function(data) {
  J <- ncol(data)
  rho <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      rho[i, j] <- cor(data[, i], data[, j], method = "spearman")
    }
  }
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  colnames(rho) <- colnames(data)
  rownames(rho) <- colnames(data)
  return(rho)
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
  load(file)  # Should load sim_data (500 × 50 matrix)
  
  # Assign correct column names
  colnames(sim_data) <- colnames(V_truspearman_matrix_result)
  
  # Rarefy
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix column names again (may be lost in rarefaction)
  colnames(rarefied_counts) <- colnames(V_truspearman_matrix_result)
  
  # Skip if taxa reduced
  if (ncol(rarefied_counts) < ncol(V_truspearman_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Spearman rho matrix
  rho_package <- spearmanMatrix_package(rarefied_counts)
  
  # Extract lower triangles
  lower_tri_package <- rho_package[lower.tri(rho_package)]
  lower_tri_true <- V_truspearman_matrix_result[lower.tri(V_truspearman_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)


# Print summary (5 decimal digits)
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")




# Prepare summary as data frame
summary_df <- data.frame(
  Method = "Spearman_Package",
  Average_MSE = round(avg_mse, 5),
  Overall_RMSE = round(overall_rmse, 5),
  Average_MAD = round(avg_mad, 5),
  Skipped_Files = skipped_files
)

# Save to CSV
write.csv(summary_df, file = "spearman_package_summary_V500.csv", row.names = FALSE)




## Vaginal n=1000
# Load required libraries
library(phyloseq)
library(Kendall)
# Load the true Spearman rho matrix
load("V_truspearman_matrix_result.RData")  # Should load truespearman_matrix_result (50x50)

# Define Spearman rho function using cor()
spearmanMatrix_package <- function(data) {
  J <- ncol(data)
  rho <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      rho[i, j] <- cor(data[, i], data[, j], method = "spearman")
    }
  }
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  colnames(rho) <- colnames(data)
  rownames(rho) <- colnames(data)
  return(rho)
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
  load(file)  # Should load sim_data (500 × 50 matrix)
  
  # Assign correct column names
  colnames(sim_data) <- colnames(V_truspearman_matrix_result)
  
  # Rarefy
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix column names again (may be lost in rarefaction)
  colnames(rarefied_counts) <- colnames(V_truspearman_matrix_result)
  
  # Skip if taxa reduced
  if (ncol(rarefied_counts) < ncol(V_truspearman_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Spearman rho matrix
  rho_package <- spearmanMatrix_package(rarefied_counts)
  
  # Extract lower triangles
  lower_tri_package <- rho_package[lower.tri(rho_package)]
  lower_tri_true <- V_truspearman_matrix_result[lower.tri(V_truspearman_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)

# Print summary (5 decimal digits)
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")



# Prepare summary as data frame
summary_df <- data.frame(
  Method = "Spearman_Package",
  Average_MSE = round(avg_mse, 5),
  Overall_RMSE = round(overall_rmse, 5),
  Average_MAD = round(avg_mad, 5),
  Skipped_Files = skipped_files
)

# Save to CSV
write.csv(summary_df, file = "spearman_package_summary_S1000.csv", row.names = FALSE)


#---------------
#Ibd n=500
#----------------

# Load required libraries
library(phyloseq)
library(Kendall)
# Load the true Spearman rho matrix
load("IBD_truspearman_matrix_result.RData")  # Should load truespearman_matrix_result (50x50)
ls()
# Define Spearman rho function using cor()
spearmanMatrix_package <- function(data) {
  J <- ncol(data)
  rho <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      rho[i, j] <- cor(data[, i], data[, j], method = "spearman")
    }
  }
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  colnames(rho) <- colnames(data)
  rownames(rho) <- colnames(data)
  return(rho)
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
  load(file)  # Should load sim_data (500 × 50 matrix)
  
  # Assign correct column names
  colnames(sim_data) <- colnames(IBD_truspearman_matrix_result)
  
  # Rarefy
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix column names again (may be lost in rarefaction)
  colnames(rarefied_counts) <- colnames(IBD_truspearman_matrix_result)
  
  # Skip if taxa reduced
  if (ncol(rarefied_counts) < ncol(IBD_truspearman_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Spearman rho matrix
  rho_package <- spearmanMatrix_package(rarefied_counts)
  
  # Extract lower triangles
  lower_tri_package <- rho_package[lower.tri(rho_package)]
  lower_tri_true <- IBD_truspearman_matrix_result[lower.tri(IBD_truspearman_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)

# Print summary (5 decimal digits)
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")



# Prepare summary as data frame
summary_df <- data.frame(
  Method = "Spearman_Package",
  Average_MSE = round(avg_mse, 5),
  Overall_RMSE = round(overall_rmse, 5),
  Average_MAD = round(avg_mad, 5),
  Skipped_Files = skipped_files
)

# Save to CSV
write.csv(summary_df, file = "spearman_package_summary_I500.csv", row.names = FALSE)

#---------------
#Ibd n=1000
#----------------

# Load required libraries
library(phyloseq)
library(Kendall)
# Load the true Spearman rho matrix
load("IBD_truspearman_matrix_result.RData")  # Should load truespearman_matrix_result (50x50)
ls()
# Define Spearman rho function using cor()
spearmanMatrix_package <- function(data) {
  J <- ncol(data)
  rho <- matrix(0, J, J)
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      rho[i, j] <- cor(data[, i], data[, j], method = "spearman")
    }
  }
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  colnames(rho) <- colnames(data)
  rownames(rho) <- colnames(data)
  return(rho)
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
  load(file)  # Should load sim_data (500 × 50 matrix)
  
  # Assign correct column names
  colnames(sim_data) <- colnames(IBD_truspearman_matrix_result)
  
  # Rarefy
  ps <- phyloseq(otu_table(sim_data, taxa_are_rows = FALSE))
  ps_rare <- rarefy_even_depth(ps, replace = FALSE, rngseed = 123, verbose = FALSE)
  rarefied_counts <- as(otu_table(ps_rare), "matrix")
  
  # Fix column names again (may be lost in rarefaction)
  colnames(rarefied_counts) <- colnames(IBD_truspearman_matrix_result)
  
  # Skip if taxa reduced
  if (ncol(rarefied_counts) < ncol(IBD_truspearman_matrix_result)) {
    warning(paste("Skipping", file, "- rarefied taxa reduced to", ncol(rarefied_counts)))
    skipped_files <- skipped_files + 1
    skipped_file_names <- c(skipped_file_names, file)
    next
  }
  
  # Compute Spearman rho matrix
  rho_package <- spearmanMatrix_package(rarefied_counts)
  
  # Extract lower triangles
  lower_tri_package <- rho_package[lower.tri(rho_package)]
  lower_tri_true <- IBD_truspearman_matrix_result[lower.tri(IBD_truspearman_matrix_result)]
  
  # Compute MSE and MAD
  mse <- mean((lower_tri_package - lower_tri_true)^2)
  mad <- median(abs(lower_tri_package - lower_tri_true))
  
  # Store results
  mse_list <- c(mse_list, mse)
  mad_list <- c(mad_list, mad)
}

# Final aggregation
avg_mse <- mean(mse_list)
avg_mad <- mean(mad_list)
overall_rmse <- sqrt(avg_mse)

# Print summary (5 decimal digits)
cat("Average MSE:", formatC(avg_mse, format = "f", digits = 5), "\n")
cat("Overall RMSE:", formatC(overall_rmse, format = "f", digits = 5), "\n")
cat("Average MAD:", formatC(avg_mad, format = "f", digits = 5), "\n")
cat("Total skipped files:", skipped_files, "\n")



# Prepare summary as data frame
summary_df <- data.frame(
  Method = "Spearman_Package",
  Average_MSE = round(avg_mse, 5),
  Overall_RMSE = round(overall_rmse, 5),
  Average_MAD = round(avg_mad, 5),
  Skipped_Files = skipped_files
)

# Save to CSV
write.csv(summary_df, file = "spearman_package_summary_IVd1000.csv", row.names = FALSE)

