##########Final simulationr result
#Tau Stool sample n=500

# Define the folder path where your .RData files are stored
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# Get the list of all .RData files
file_list <- list.files(path = folder_path, pattern = "^TS_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize a data frame to store the dimensions for each iteration
dimensions_check <- data.frame(
  File = character(),
  Truetau_Dim = character(),
  Tau_Dim = character(),
  stringsAsFactors = FALSE
)

# Loop through all the .RData files
for (file in file_list) {
  # Load the .RData file
  load(file)
  
  # Get the dimensions of the matrices
  truetau_dim <- dim(truetau_matrix_result)
  tau_dim <- dim(tau_matrix_result)
  
  # Add the dimensions to the data frame
  dimensions_check <- rbind(
    dimensions_check,
    data.frame(
      File = basename(file),
      Truetau_Dim = paste(truetau_dim, collapse = "x"),
      Tau_Dim = paste(tau_dim, collapse = "x"),
      stringsAsFactors = FALSE
    )
  )
}

# Print the dimensions check
print(dimensions_check)

# Save the dimensions check result to a CSV file (optional)
#write.csv(dimensions_check, file = "dimensions_check.csv", row.names = FALSE)



library(data.table)  # Load the data.table package

# Initialize an empty data table
results_table <- data.table(Task_ID = integer(), mse_newdeveloped = numeric(), mse_taustar = numeric())

# Loop through each file and extract the variables
for (file in file_list) {
  # Load the .RData file
  load(file)
  
  # Extract the task ID from the file name
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  # Append the values to the results table efficiently using rbindlist
  results_table <- rbindlist(
    list(results_table, 
         data.table(Task_ID = task_id, mse_newdeveloped = mse_newdeveloped, mse_taustar = mse_taustar))
  )
}

# Sort the results table by Task_ID using data.table syntax
# Format the MSE values to prevent scientific notation
results_table[, mse_newdeveloped := format(mse_newdeveloped, scientific = FALSE)]
results_table[, mse_taustar := format(mse_taustar, scientific = FALSE)]

# Print the updated results table
print(results_table)


#Convert the MSE columns back to numeric
results_table[, mse_newdeveloped := as.numeric(mse_newdeveloped)]
results_table[, mse_taustar := as.numeric(mse_taustar)]

# Compute RMSE for each method
rmse_newdeveloped <- sqrt(mean(results_table$mse_newdeveloped))
rmse_taustar <- sqrt(mean(results_table$mse_taustar))

# Print the RMSEs
cat("RMSE (Newly Developed):", rmse_newdeveloped, "\n")
cat("RMSE (Tau*):", rmse_taustar, "\n")

# Print the RMSEs with 5 decimal places
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_newdeveloped), "\n")
cat("RMSE (Tau*):", sprintf("%.5f", rmse_taustar), "\n")






library(data.table)

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^TS_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_newdeveloped = numeric(),
                            mse_taustar = numeric(),
                            mad_new = numeric(),
                            mad_taustar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mse_newdeveloped,
    mse_taustar = mse_taustar,
    mad_new = mad_new,
    mad_taustar = mad_taustar
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_newdeveloped))
rmse_tau <- sqrt(mean(results_table$mse_taustar))
mean_mad_new <- mean(results_table$mad_new)
mean_mad_tau <- mean(results_table$mad_taustar)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Tau*):", sprintf("%.5f", rmse_tau), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Tau*):", sprintf("%.5f", mean_mad_tau), "\n")
######

#############################
#Rho Stool sample n=500
############################

library(data.table)

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^SS_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_newdeveloped = numeric(),
                            mse_rhostar = numeric(),
                            mad_new = numeric(),
                            mad_rhostar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mse_newdeveloped,
    mse_rhostar = mse_rhostar,
    mad_new = mad_new,
    mad_rhostar = mad_rhostar
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_table$mse_rhostar))
mean_mad_new <- mean(results_table$mad_new)
mean_mad_rho <- mean(results_table$mad_rhostar)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")

#########
library(data.table)

folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^SS_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results summary table
results_summary <- data.table(Task_ID = integer(),
                              mse_newdeveloped = numeric(),
                              mse_rhostar = numeric(),
                              mad_new = numeric(),
                              mad_rhostar = numeric())

# Initialize detailed table
results_detailed <- data.table(Task_ID = integer(),
                               simulation_index = integer(),
                               mse_newdeveloped = numeric(),
                               mse_rhostar = numeric(),
                               mad_new = numeric(),
                               mad_rhostar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Should load mse_newdeveloped, mse_rhostar, mad_new, mad_rhostar
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  # --- Summary table (averaged per task)
  results_summary <- rbindlist(list(results_summary, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mean(mse_newdeveloped),
    mse_rhostar = mean(mse_rhostar),
    mad_new = mean(mad_new),
    mad_rhostar = mean(mad_rhostar)
  )))
  
  # --- Detailed table (simulation-wise)
  n_sim <- length(mse_newdeveloped)
  detailed_dt <- data.table(
    Task_ID = rep(task_id, n_sim),
    simulation_index = seq_len(n_sim),
    mse_newdeveloped = mse_newdeveloped,
    mse_rhostar = mse_rhostar,
    mad_new = mad_new,
    mad_rhostar = mad_rhostar
  )
  
  results_detailed <- rbindlist(list(results_detailed, detailed_dt))
}

# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar)

cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")

# Print detailed results to console
cat("\n==== DETAILED SIMULATION RESULTS ====\n")
print(results_detailed)

######


# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar, na.rm = TRUE))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar, na.rm = TRUE)

# Print summary
cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")



##############################
## Tau vaginal sample n=500
#############################

library(data.table)

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^TS_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_newdeveloped = numeric(),
                            mse_taustar = numeric(),
                            mad_new = numeric(),
                            mad_taustar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mse_newdeveloped,
    mse_taustar = mse_taustar,
    mad_new = mad_new,
    mad_taustar = mad_taustar
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_newdeveloped))
rmse_tau <- sqrt(mean(results_table$mse_taustar))
mean_mad_new <- mean(results_table$mad_new)
mean_mad_tau <- mean(results_table$mad_taustar)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Tau*):", sprintf("%.5f", rmse_tau), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Tau*):", sprintf("%.5f", mean_mad_tau), "\n")
######
library(data.table)
# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^TV_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_newdeveloped = numeric(),
                            mse_taustar = numeric(),
                            mad_new = numeric(),
                            mad_taustar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mse_newdeveloped,
    mse_taustar = mse_taustar,
    mad_new = mad_new,
    mad_taustar = mad_taustar
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_newdeveloped))
rmse_tau <- sqrt(mean(results_table$mse_taustar))
mean_mad_new <- mean(results_table$mad_new)
mean_mad_tau <- mean(results_table$mad_taustar)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Tau*):", sprintf("%.5f", rmse_tau), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Tau*):", sprintf("%.5f", mean_mad_tau), "\n")

#####################################
#spearman rho Vaginal sample n=500
######################################
# Define folder path
library(data.table)

folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^SV_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results summary table
results_summary <- data.table(Task_ID = integer(),
                              mse_newdeveloped = numeric(),
                              mse_rhostar = numeric(),
                              mad_new = numeric(),
                              mad_rhostar = numeric())

# Initialize detailed table
results_detailed <- data.table(Task_ID = integer(),
                               simulation_index = integer(),
                               mse_newdeveloped = numeric(),
                               mse_rhostar = numeric(),
                               mad_new = numeric(),
                               mad_rhostar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Should load mse_newdeveloped, mse_rhostar, mad_new, mad_rhostar
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  # --- Summary table (averaged per task)
  results_summary <- rbindlist(list(results_summary, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mean(mse_newdeveloped),
    mse_rhostar = mean(mse_rhostar),
    mad_new = mean(mad_new),
    mad_rhostar = mean(mad_rhostar)
  )))
  
  # --- Detailed table (simulation-wise)
  n_sim <- length(mse_newdeveloped)
  detailed_dt <- data.table(
    Task_ID = rep(task_id, n_sim),
    simulation_index = seq_len(n_sim),
    mse_newdeveloped = mse_newdeveloped,
    mse_rhostar = mse_rhostar,
    mad_new = mad_new,
    mad_rhostar = mad_rhostar
  )
  
  results_detailed <- rbindlist(list(results_detailed, detailed_dt))
}

# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar)

cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")

# Print detailed results to console
cat("\n==== DETAILED SIMULATION RESULTS ====\n")
print(results_detailed)

######


# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar, na.rm = TRUE))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar, na.rm = TRUE)

# Print summary
cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")



##########################
###Tau for IBD n=500
###########################

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^TI_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_newdeveloped = numeric(),
                            mse_taustar = numeric(),
                            mad_new = numeric(),
                            mad_taustar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mse_newdeveloped,
    mse_taustar = mse_taustar,
    mad_new = mad_new,
    mad_taustar = mad_taustar
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_newdeveloped))
rmse_tau <- sqrt(mean(results_table$mse_taustar))
mean_mad_new <- mean(results_table$mad_new)
mean_mad_tau <- mean(results_table$mad_taustar)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Tau*):", sprintf("%.5f", rmse_tau), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Tau*):", sprintf("%.5f", mean_mad_tau), "\n")

###############################
###Gini for stool sample n=500
################################

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^Gini_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_Gini_develope= numeric(),
                            mse_Gini_original = numeric(),
                            mad_Gini_develope = numeric(),
                            mad_Gini_original = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_Gini_develope = mse_Gini_develope,
    mse_Gini_original = mse_Gini_original,
    mad_Gini_develope = mad_Gini_develope,
    mad_Gini_original = mad_Gini_original
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_Gini_develope))
rmse_gini <- sqrt(mean(results_table$mse_Gini_original))
mean_mad_new <- mean(results_table$mad_Gini_develope)
mean_mad_gini <- mean(results_table$mad_Gini_original)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Gini_original):", sprintf("%.5f", rmse_gini), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Gini_original):", sprintf("%.5f", mean_mad_gini), "\n")


###########################
### spearman IBD n=500
##########################

library(data.table)

folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^SI_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results summary table
results_summary <- data.table(Task_ID = integer(),
                              mse_newdeveloped = numeric(),
                              mse_rhostar = numeric(),
                              mad_new = numeric(),
                              mad_rhostar = numeric())

# Initialize detailed table
results_detailed <- data.table(Task_ID = integer(),
                               simulation_index = integer(),
                               mse_newdeveloped = numeric(),
                               mse_rhostar = numeric(),
                               mad_new = numeric(),
                               mad_rhostar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Should load mse_newdeveloped, mse_rhostar, mad_new, mad_rhostar
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  # --- Summary table (averaged per task)
  results_summary <- rbindlist(list(results_summary, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mean(mse_newdeveloped),
    mse_rhostar = mean(mse_rhostar),
    mad_new = mean(mad_new),
    mad_rhostar = mean(mad_rhostar)
  )))
  
  # --- Detailed table (simulation-wise)
  n_sim <- length(mse_newdeveloped)
  detailed_dt <- data.table(
    Task_ID = rep(task_id, n_sim),
    simulation_index = seq_len(n_sim),
    mse_newdeveloped = mse_newdeveloped,
    mse_rhostar = mse_rhostar,
    mad_new = mad_new,
    mad_rhostar = mad_rhostar
  )
  
  results_detailed <- rbindlist(list(results_detailed, detailed_dt))
}

# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar)


cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")

# Print detailed results to console
cat("\n==== DETAILED SIMULATION RESULTS ====\n")
print(results_detailed)


# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar, na.rm = TRUE))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar, na.rm = TRUE)

# Print summary
cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")


##################################
#stool sample kendell tau n=1000
#################################

library(data.table)

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^TS2_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_newdeveloped = numeric(),
                            mse_taustar = numeric(),
                            mad_new = numeric(),
                            mad_taustar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mse_newdeveloped,
    mse_taustar = mse_taustar,
    mad_new = mad_new,
    mad_taustar = mad_taustar
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_newdeveloped))
rmse_tau <- sqrt(mean(results_table$mse_taustar))
mean_mad_new <- mean(results_table$mad_new)
mean_mad_tau <- mean(results_table$mad_taustar)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Tau*):", sprintf("%.5f", rmse_tau), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Tau*):", sprintf("%.5f", mean_mad_tau), "\n")


#############################
##Stool spearman n=1000
##############################
library(data.table)

folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^SS2_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results summary table
results_summary <- data.table(Task_ID = integer(),
                              mse_newdeveloped = numeric(),
                              mse_rhostar = numeric(),
                              mad_new = numeric(),
                              mad_rhostar = numeric())

# Initialize detailed table
results_detailed <- data.table(Task_ID = integer(),
                               simulation_index = integer(),
                               mse_newdeveloped = numeric(),
                               mse_rhostar = numeric(),
                               mad_new = numeric(),
                               mad_rhostar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Should load mse_newdeveloped, mse_rhostar, mad_new, mad_rhostar
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  # --- Summary table (averaged per task)
  results_summary <- rbindlist(list(results_summary, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mean(mse_newdeveloped),
    mse_rhostar = mean(mse_rhostar),
    mad_new = mean(mad_new),
    mad_rhostar = mean(mad_rhostar)
  )))
  
  # --- Detailed table (simulation-wise)
  n_sim <- length(mse_newdeveloped)
  detailed_dt <- data.table(
    Task_ID = rep(task_id, n_sim),
    simulation_index = seq_len(n_sim),
    mse_newdeveloped = mse_newdeveloped,
    mse_rhostar = mse_rhostar,
    mad_new = mad_new,
    mad_rhostar = mad_rhostar
  )
  
  results_detailed <- rbindlist(list(results_detailed, detailed_dt))
}

# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar)

cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")

# Print detailed results to console
cat("\n==== DETAILED SIMULATION RESULTS ====\n")
print(results_detailed)

######


# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar, na.rm = TRUE))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar, na.rm = TRUE)

# Print summary
cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")


#######################
#Gini Vaginal n=500
#####################
library(data.table)

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^V_Gini_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_Gini_develope= numeric(),
                            mse_Gini_original = numeric(),
                            mad_Gini_develope = numeric(),
                            mad_Gini_original = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_Gini_develope = mse_Gini_develope,
    mse_Gini_original = mse_Gini_original,
    mad_Gini_develope = mad_Gini_develope,
    mad_Gini_original = mad_Gini_original
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_Gini_develope))
rmse_gini <- sqrt(mean(results_table$mse_Gini_original))
mean_mad_new <- mean(results_table$mad_Gini_develope)
mean_mad_gini <- mean(results_table$mad_Gini_original)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Gini_original):", sprintf("%.5f", rmse_gini), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Gini_original):", sprintf("%.5f", mean_mad_gini), "\n")

#####################
###Gini stool n=1000
#######################

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^GiniS2_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_Gini_develope= numeric(),
                            mse_Gini_original = numeric(),
                            mad_Gini_develope = numeric(),
                            mad_Gini_original = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_Gini_develope = mse_Gini_develope,
    mse_Gini_original = mse_Gini_original,
    mad_Gini_develope = mad_Gini_develope,
    mad_Gini_original = mad_Gini_original
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_Gini_develope))
rmse_gini <- sqrt(mean(results_table$mse_Gini_original))
mean_mad_new <- mean(results_table$mad_Gini_develope)
mean_mad_gini <- mean(results_table$mad_Gini_original)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Gini_original):", sprintf("%.5f", rmse_gini), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Gini_original):", sprintf("%.5f", mean_mad_gini), "\n")

#############################
#vaginal spearman n=1000
###############################
library(data.table)

folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^SV2_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results summary table
results_summary <- data.table(Task_ID = integer(),
                              mse_newdeveloped = numeric(),
                              mse_rhostar = numeric(),
                              mad_new = numeric(),
                              mad_rhostar = numeric())

# Initialize detailed table
results_detailed <- data.table(Task_ID = integer(),
                               simulation_index = integer(),
                               mse_newdeveloped = numeric(),
                               mse_rhostar = numeric(),
                               mad_new = numeric(),
                               mad_rhostar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Should load mse_newdeveloped, mse_rhostar, mad_new, mad_rhostar
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  # --- Summary table (averaged per task)
  results_summary <- rbindlist(list(results_summary, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mean(mse_newdeveloped),
    mse_rhostar = mean(mse_rhostar),
    mad_new = mean(mad_new),
    mad_rhostar = mean(mad_rhostar)
  )))
  
  # --- Detailed table (simulation-wise)
  n_sim <- length(mse_newdeveloped)
  detailed_dt <- data.table(
    Task_ID = rep(task_id, n_sim),
    simulation_index = seq_len(n_sim),
    mse_newdeveloped = mse_newdeveloped,
    mse_rhostar = mse_rhostar,
    mad_new = mad_new,
    mad_rhostar = mad_rhostar
  )
  
  results_detailed <- rbindlist(list(results_detailed, detailed_dt))
}

# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar)

cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")

# Print detailed results to console
cat("\n==== DETAILED SIMULATION RESULTS ====\n")
print(results_detailed)

######


# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar, na.rm = TRUE))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar, na.rm = TRUE)

# Print summary
cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")


###############################
###Tau vaginal sample n=1000
#################################

library(data.table)
# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^TV2_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_newdeveloped = numeric(),
                            mse_taustar = numeric(),
                            mad_new = numeric(),
                            mad_taustar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mse_newdeveloped,
    mse_taustar = mse_taustar,
    mad_new = mad_new,
    mad_taustar = mad_taustar
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_newdeveloped))
rmse_tau <- sqrt(mean(results_table$mse_taustar))
mean_mad_new <- mean(results_table$mad_new)
mean_mad_tau <- mean(results_table$mad_taustar)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Tau*):", sprintf("%.5f", rmse_tau), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Tau*):", sprintf("%.5f", mean_mad_tau), "\n")


##########################
##GIni IBD sample n=500
#######################

library(data.table)

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^IBD_Gini_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_Gini_develope= numeric(),
                            mse_Gini_original = numeric(),
                            mad_Gini_develope = numeric(),
                            mad_Gini_original = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_Gini_develope = mse_Gini_develope,
    mse_Gini_original = mse_Gini_original,
    mad_Gini_develope = mad_Gini_develope,
    mad_Gini_original = mad_Gini_original
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_Gini_develope))
rmse_gini <- sqrt(mean(results_table$mse_Gini_original))
mean_mad_new <- mean(results_table$mad_Gini_develope)
mean_mad_gini <- mean(results_table$mad_Gini_original)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Gini_original):", sprintf("%.5f", rmse_gini), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Gini_original):", sprintf("%.5f", mean_mad_gini), "\n")



########################
#Gini Vaginal n=500
###############################
library(data.table)

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^V2_Gini_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_Gini_develope= numeric(),
                            mse_Gini_original = numeric(),
                            mad_Gini_develope = numeric(),
                            mad_Gini_original = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_Gini_develope = mse_Gini_develope,
    mse_Gini_original = mse_Gini_original,
    mad_Gini_develope = mad_Gini_develope,
    mad_Gini_original = mad_Gini_original
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_Gini_develope))
rmse_gini <- sqrt(mean(results_table$mse_Gini_original))
mean_mad_new <- mean(results_table$mad_Gini_develope)
mean_mad_gini <- mean(results_table$mad_Gini_original)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Gini_original):", sprintf("%.5f", rmse_gini), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Gini_original):", sprintf("%.5f", mean_mad_gini), "\n")


#########################
###Tau for IBD n=1000
######################
library(data.table)

# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^TI2_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_newdeveloped = numeric(),
                            mse_taustar = numeric(),
                            mad_new = numeric(),
                            mad_taustar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mse_newdeveloped,
    mse_taustar = mse_taustar,
    mad_new = mad_new,
    mad_taustar = mad_taustar
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_newdeveloped))
rmse_tau <- sqrt(mean(results_table$mse_taustar))
mean_mad_new <- mean(results_table$mad_new)
mean_mad_tau <- mean(results_table$mad_taustar)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Tau*):", sprintf("%.5f", rmse_tau), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Tau*):", sprintf("%.5f", mean_mad_tau), "\n")



#######################
#rho IBD n=1000
####################
library(data.table)

folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^SI2_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results summary table
results_summary <- data.table(Task_ID = integer(),
                              mse_newdeveloped = numeric(),
                              mse_rhostar = numeric(),
                              mad_new = numeric(),
                              mad_rhostar = numeric())

# Initialize detailed table
results_detailed <- data.table(Task_ID = integer(),
                               simulation_index = integer(),
                               mse_newdeveloped = numeric(),
                               mse_rhostar = numeric(),
                               mad_new = numeric(),
                               mad_rhostar = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Should load mse_newdeveloped, mse_rhostar, mad_new, mad_rhostar
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  # --- Summary table (averaged per task)
  results_summary <- rbindlist(list(results_summary, data.table(
    Task_ID = task_id,
    mse_newdeveloped = mean(mse_newdeveloped),
    mse_rhostar = mean(mse_rhostar),
    mad_new = mean(mad_new),
    mad_rhostar = mean(mad_rhostar)
  )))
  
  # --- Detailed table (simulation-wise)
  n_sim <- length(mse_newdeveloped)
  detailed_dt <- data.table(
    Task_ID = rep(task_id, n_sim),
    simulation_index = seq_len(n_sim),
    mse_newdeveloped = mse_newdeveloped,
    mse_rhostar = mse_rhostar,
    mad_new = mad_new,
    mad_rhostar = mad_rhostar
  )
  
  results_detailed <- rbindlist(list(results_detailed, detailed_dt))
}

# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar)


cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")

# Print detailed results to console
cat("\n==== DETAILED SIMULATION RESULTS ====\n")
print(results_detailed)


# Compute and print overall summary
rmse_new <- sqrt(mean(results_summary$mse_newdeveloped))
rmse_rho <- sqrt(mean(results_summary$mse_rhostar, na.rm = TRUE))
mean_mad_new <- mean(results_summary$mad_new)
mean_mad_rho <- mean(results_summary$mad_rhostar, na.rm = TRUE)

# Print summary
cat("==== SUMMARY METRICS ====\n")
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (rho*):", sprintf("%.5f", rmse_rho), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (rho*):", sprintf("%.5f", mean_mad_rho), "\n")





######################
###Gini IBD n=1000
#########################


library(data.table)
# Define folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# List all matching .RData files
file_list <- list.files(path = folder_path, pattern = "^IBD2_Gini_simulation_task_.*\\.RData$", full.names = TRUE)

# Initialize results table
results_table <- data.table(Task_ID = integer(),
                            mse_Gini_develope= numeric(),
                            mse_Gini_original = numeric(),
                            mad_Gini_develope = numeric(),
                            mad_Gini_original = numeric())

# Loop through each file
for (file in file_list) {
  load(file)  # Assumes mse_newdeveloped, mse_taustar, mad_new, mad_taustar are in the file
  
  task_id <- as.numeric(gsub("^.*_task_(\\d+)\\.RData$", "\\1", basename(file)))
  
  results_table <- rbindlist(list(results_table, data.table(
    Task_ID = task_id,
    mse_Gini_develope = mse_Gini_develope,
    mse_Gini_original = mse_Gini_original,
    mad_Gini_develope = mad_Gini_develope,
    mad_Gini_original = mad_Gini_original
  )))
}

# Compute RMSE and average MAD
rmse_new <- sqrt(mean(results_table$mse_Gini_develope))
rmse_gini <- sqrt(mean(results_table$mse_Gini_original))
mean_mad_new <- mean(results_table$mad_Gini_develope)
mean_mad_gini <- mean(results_table$mad_Gini_original)

# Print results (5 decimal places)
cat("RMSE (Newly Developed):", sprintf("%.5f", rmse_new), "\n")
cat("RMSE (Gini_original):", sprintf("%.5f", rmse_gini), "\n")
cat("MAD (Newly Developed):", sprintf("%.5f", mean_mad_new), "\n")
cat("MAD (Gini_original):", sprintf("%.5f", mean_mad_gini), "\n")




