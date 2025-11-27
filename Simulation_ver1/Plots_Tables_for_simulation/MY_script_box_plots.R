
############################################################
#Box plot for Sample type stool n=500 for comparison of Tau
############################################################

library(ggplot2)
library(dplyr)
library(phyloseq)

# Define the folder path where your .RData files are stored
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# Get the list of all .RData files
file_list <- list.files(path = folder_path, pattern = "^TS_simulation_task_.*\\.RData$", full.names = TRUE)

# Step 1: Check dimensions
dimensions_check <- data.frame(
  File = character(),
  Truetau_Dim = character(),
  Tau_Dim = character(),
  stringsAsFactors = FALSE
)

for (file in file_list) {
  load(file)
  truetau_dim <- dim(truetau_matrix_result)
  tau_dim <- dim(tau_matrix_result)
  dimensions_check <- rbind(dimensions_check, data.frame(
    File = basename(file),
    Truetau_Dim = paste(truetau_dim, collapse = "x"),
    Tau_Dim = paste(tau_dim, collapse = "x"),
    stringsAsFactors = FALSE
  ))
}

# Step 2: Exclude problematic files
exclude_files <- dimensions_check$File[dimensions_check$Tau_Dim %in% c("49x49", "48x48")]
filtered_file_list <- file_list[!basename(file_list) %in% exclude_files]

# Step 3: Define the selected matrix pairs
selected_pairs <- list(c(1, 2), c(12, 5), c(10, 20), c(15, 25), c(8, 30),
                       c(3, 18), c(22, 7), c(14, 19), c(4, 6), c(11, 13))
pair_names <- sapply(selected_pairs, function(pair) paste(pair, collapse = ","))

# Step 4: Initialize lists for values
selected_elements_newdev <- vector("list", length(selected_pairs))
selected_elements_taustar <- vector("list", length(selected_pairs))
names(selected_elements_newdev) <- pair_names
names(selected_elements_taustar) <- pair_names  # 

# Step 5: Loop over filtered files
for (file in filtered_file_list) {
  load(file)
  
  # Check if dimensions match to proceed
  if (!all(dim(truetau_matrix_result) == dim(tau_matrix_result))) {
    cat("Dimension mismatch in file:", basename(file), "\n")
    next
  }
  
  # Compute difference matrices
  diff_newdev <- truetau_matrix_result - tau_matrix_result
  diff_taustar <- truetau_matrix_result - taustar_matrix_result
  
  # Extract values for selected pairs
  for (i in seq_along(selected_pairs)) {
    pair <- selected_pairs[[i]]
    selected_elements_newdev[[i]] <- c(selected_elements_newdev[[i]], diff_newdev[pair[1], pair[2]])
    selected_elements_taustar[[i]] <- c(selected_elements_taustar[[i]], diff_taustar[pair[1], pair[2]])
  }
}

# Step 6: Prepare data for boxplot
boxplot_data_newdev <- do.call(rbind, lapply(names(selected_elements_newdev), function(name) {
  data.frame(Pair = name, Value = selected_elements_newdev[[name]], Difference_Type = "truetau - newdev")
}))

boxplot_data_taustar <- do.call(rbind, lapply(names(selected_elements_taustar), function(name) {
  data.frame(Pair = name, Value = selected_elements_taustar[[name]], Difference_Type = "truetau - rarified_taustar")
}))

boxplot_data <- rbind(boxplot_data_newdev, boxplot_data_taustar)

# Step 7: Plot
ggplot(boxplot_data, aes(x = Pair, y = Value, fill = Difference_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = c("skyblue", "lightcoral")) +
  labs(title = "Comparison of True and Estimated Differences for NewDev vs Rarefied Tau* (n = 500)",
       x = "Selected Pairs (Actual Indices)",
       y = "Difference Values",
       fill = "Difference Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 100, hjust = 1))




# Make the legend order explicit
boxplot_data$Difference_Type <- factor(
  boxplot_data$Difference_Type,
  levels = c("truetau - newdev", "truetau - rarified_taustar")
)

ggplot(boxplot_data, aes(x = Pair, y = Value, fill = Difference_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  scale_fill_manual(
    breaks = c("truetau - newdev", "truetau - rarified_taustar"),
    labels = c(
      expression("true " * tau^"*" - hat(tau^"*")[plain("ZILN")]),   # true τ − τ̂_ZILN
      expression("true " * tau^"*" - "Emp " *hat(tau^"*"))     # true τ − Empirical τ*
    ),
    values = c("skyblue", "lightcoral")
  ) +
  labs(
    title = bquote("Comparison of True and Estimated Differences for " ~ hat(tau^"*")[ZILN] ~ " vs Empirical " ~ hat(tau^"*") ~ " (n = 500, Template=Stool)"),
    x = "Selected Pairs (Actual Indices)",
    y = "Difference Values",
    fill = "Difference Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 100, hjust = 1))



library(ggplot2)
library(Cairo)



pdf("TS5500.pdf", width = 7, height = 5, useDingbats = FALSE)

# ensure legend order
boxplot_data$Difference_Type <- factor(
  boxplot_data$Difference_Type,
  levels = c("truetau - newdev", "truetau - rarified_taustar")
)

pretty_pairs <- paste0("(", boxplot_data$Pair, ")")
boxplot_data$Pair <- factor(pretty_pairs, levels = unique(pretty_pairs))

library(ggplot2)

ggplot(boxplot_data, aes(Pair, Value, fill = Difference_Type)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(
    breaks = c("truetau - newdev", "truetau - rarified_taustar"),
    labels = c(
      expression("true " * tau^"*" - hat(tau)["ZI"]^"*"),
      expression("true " * tau^"*" - hat(tau)["Emp"]^"*")
    ),
    values = c("#74add1", "lightcoral")
  ) +
  # no fixed range; let ggplot choose nice limits & breaks
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 9),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = bquote("Comparison of true " * tau^"*" * " differences for n=500"),
    x = "Selected pairs (indices)",
    y = "Difference values",
    fill = "Difference Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
    legend.position = c(0.98, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    legend.key.size = unit(10, "pt"),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", margin = margin(b = 6))
  )

dev.off()





library(ggplot2)
library(dplyr)

# Define the folder path
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# Step 1: Get list of relevant files
file_list <- list.files(path = folder_path, pattern = "^SI2_simulation_task_.*\\.RData$", full.names = TRUE)

# Step 2: Check matrix dimensions
dimensions_check <- data.frame(
  File = character(),
  Truerho_Dim = character(),
  rho_Dim = character(),
  stringsAsFactors = FALSE
)

for (file in file_list) {
  load(file)
  Truerho_Dim <- dim(IBD_truspearman_matrix_result)
  rho_Dim <- dim(spearman_matrix_result)
  dimensions_check <- rbind(dimensions_check, data.frame(
    File = basename(file),
    Truerho_Dim = paste(Truerho_Dim, collapse = "x"),
    rho_Dim = paste(rho_Dim, collapse = "x"),
    stringsAsFactors = FALSE
  ))
}

# Step 3: Exclude problematic files (e.g., small dimensions)
exclude_files <- dimensions_check$File[dimensions_check$rho_Dim %in% c("49x49", "48x48")]
filtered_file_list <- file_list[!basename(file_list) %in% exclude_files]

# Step 4: Define selected matrix index pairs (i,j)
selected_pairs <- list(c(2, 3), c(12, 5), c(10, 20), c(15, 25), c(8, 30),
                       c(26, 32), c(22, 7), c(45,47), c(4, 6), c(11, 13))
pair_names <- sapply(selected_pairs, function(pair) paste(pair, collapse = ","))

# Step 5: Initialize lists for pairwise difference values
selected_elements_newdev <- vector("list", length(selected_pairs))
selected_elements_rhostar <- vector("list", length(selected_pairs))  # <- Corrected name
names(selected_elements_newdev) <- pair_names
names(selected_elements_rhostar) <- pair_names

# Step 6: Loop over filtered files
for (file in filtered_file_list) {
  load(file)
  
  # Check dimension match before subtraction
  if (!all(dim(IBD_truspearman_matrix_result) == dim(spearman_matrix_result))) {
    cat("Dimension mismatch in file:", basename(file), "\n")
    next
  }
  
  # Difference matrices
  diff_newdev <- IBD_truspearman_matrix_result - spearman_matrix_result
  diff_rhostar <- IBD_truspearman_matrix_result - rho_star_result
  
  # Extract difference values at the selected pairs
  for (i in seq_along(selected_pairs)) {
    pair <- selected_pairs[[i]]
    selected_elements_newdev[[i]] <- c(selected_elements_newdev[[i]], diff_newdev[pair[1], pair[2]])
    selected_elements_rhostar[[i]] <- c(selected_elements_rhostar[[i]], diff_rhostar[pair[1], pair[2]])
  }
}

# Step 7: Prepare data for plotting
boxplot_data_newdev <- do.call(rbind, lapply(names(selected_elements_newdev), function(name) {
  data.frame(Pair = name, Value = selected_elements_newdev[[name]], Difference_Type = "truerho - newdev")
}))

boxplot_data_rhostar <- do.call(rbind, lapply(names(selected_elements_rhostar), function(name) {
  data.frame(Pair = name, Value = selected_elements_rhostar[[name]], Difference_Type = "truerho - rarefied_rhostar")
}))

boxplot_data <- rbind(boxplot_data_newdev, boxplot_data_rhostar)

# Step 8: Plot the boxplot
ggplot(boxplot_data, aes(x = Pair, y = Value, fill = Difference_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = c("skyblue", "lightcoral")) +
  labs(title = "Comparison of True and Estimated Spearman Differences (IBD Sample, n = 500)",
       x = "Selected Taxa Pairs (Index)",
       y = "Difference (True - Estimated)",
       fill = "Difference Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 100, hjust = 1))





# --- Save as PDF (vector, high quality) ---
pdf("SI2.pdf", width = 7, height = 5, useDingbats = FALSE)

# 1) Lock legend order (must match your data)
boxplot_data$Difference_Type <- factor(
  boxplot_data$Difference_Type,
  levels = c("truerho - newdev", "truerho - rarefied_rhostar")
)

# 2) Clean pair labels: ensure exactly one set of parentheses like "(i,j)"
pairs_clean <- gsub("[()]", "", boxplot_data$Pair)      # remove any existing ()
pretty_pairs <- sprintf("(%s)", pairs_clean)            # add one pair of ()
boxplot_data$Pair <- factor(pretty_pairs, levels = unique(pretty_pairs))

# 3) Plot
ggplot(boxplot_data, aes(Pair, Value, fill = Difference_Type)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(
    breaks = c("truerho - newdev", "truerho - rarefied_rhostar"),
    labels = c(
      # true ρ* − ρ̂*_{Emp}
      expression("true " * rho^"*" - hat(rho)[ZI]^"*"),
      # true ρ* − ρ̂*_{ZILN}
      
      expression("true " * rho^"*" - hat(rho)[Emp]^"*")
    ),
    values = c("#74add1", "lightcoral")
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 9),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = bquote("Comparison of true " * rho^"*" * " differences for n=1000"),
    x = "Selected pairs (indices)",
    y = "Difference values",
    fill = "Difference Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
    legend.position = c(0.98, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    legend.key.size = unit(10, "pt"),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", margin = margin(b = 6))
  )

dev.off()





########
#####Box PLot for GINI

library(ggplot2)
library(dplyr)

# Path to your Gini simulation .RData files
folder_path <- "C:/Users/HP/OneDrive/Documents/Research 1st year"

# -------------------------
# LIST ALL FILES
# -------------------------

file_list <- list.files(
  path = folder_path,
  pattern = "^Gini_simulation_task_.*\\.RData$",
  full.names = TRUE
)

# -------------------------
# DIMENSION CHECK BLOCK
# -------------------------

# Initialize data frame to store dimensions
dimensions_check <- data.frame(
  File = character(),
  TruGini_Dim = character(),
  NewGini_Dim = character(),
  OriginalGini_Dim = character(),
  stringsAsFactors = FALSE
)

for (file in file_list) {
  load(file)
  
  # Check if all required objects exist
  if (!exists("Gini_matrix_tru_developed") ||
      !exists("Gini_matrix_newly_developed_result") ||
      !exists("Gini_matrix_original")) {
    cat("Skipping file due to missing objects:", basename(file), "\n")
    next
  }
  
  # Get dimensions
  tru_dim <- dim(Gini_matrix_tru_developed)
  new_dim <- dim(Gini_matrix_newly_developed_result)
  original_dim <- dim(Gini_matrix_original)
  
  # Store dimensions
  dimensions_check <- rbind(
    dimensions_check,
    data.frame(
      File = basename(file),
      TruGini_Dim = paste(tru_dim, collapse = "x"),
      NewGini_Dim = paste(new_dim, collapse = "x"),
      OriginalGini_Dim = paste(original_dim, collapse = "x"),
      stringsAsFactors = FALSE
    )
  )
}

# -------------------------
# FILTER OUT BAD FILES
# -------------------------

# Exclude files where dimensions don't match true Gini matrix
exclude_files <- dimensions_check$File[
  !(dimensions_check$TruGini_Dim == dimensions_check$NewGini_Dim &
      dimensions_check$TruGini_Dim == dimensions_check$OriginalGini_Dim)
]

cat("Number of excluded files due to dimension mismatch:", length(exclude_files), "\n")

filtered_file_list <- file_list[!basename(file_list) %in% exclude_files]

# -------------------------
# SELECTED PAIRS
# -------------------------

selected_pairs <- list(c(2, 3), c(12, 5),  c(10, 20),c(15, 25),c(8, 30),c(26, 32)
                       ,c(22, 7),c(45, 47),c(4, 6),c(47, 49))

pair_names <- sapply(selected_pairs, function(pair) paste(pair, collapse = ","))

# Initialize result lists
selected_elements_newdev <- vector("list", length(selected_pairs))
selected_elements_original <- vector("list", length(selected_pairs))
names(selected_elements_newdev) <- pair_names
names(selected_elements_original) <- pair_names

# -------------------------
# LOOP OVER FILTERED FILES
# -------------------------

for (file in filtered_file_list) {
  load(file)
  
  diff_newdev <- Gini_matrix_tru_developed - Gini_matrix_newly_developed_result
  diff_original <- Gini_matrix_tru_developed - Gini_matrix_original
  
  for (i in seq_along(selected_pairs)) {
    pair <- selected_pairs[[i]]
    val_newdev <- diff_newdev[pair[1], pair[2]]
    val_original <- diff_original[pair[1], pair[2]]
    
    selected_elements_newdev[[i]] <- c(
      selected_elements_newdev[[i]],
      val_newdev
    )
    selected_elements_original[[i]] <- c(
      selected_elements_original[[i]],
      val_original
    )
  }
}

# -------------------------
# BUILD DATA FRAMES
# -------------------------

boxplot_data_newdev <- do.call(rbind, lapply(names(selected_elements_newdev), function(name) {
  data.frame(
    Pair = name,
    Value = selected_elements_newdev[[name]],
    Difference_Type = "truGini - newGini"
  )
}))

boxplot_data_original <- do.call(rbind, lapply(names(selected_elements_original), function(name) {
  data.frame(
    Pair = name,
    Value = selected_elements_original[[name]],
    Difference_Type = "truGini - originalGini"
  )
}))

# Merge datasets
boxplot_data <- rbind(boxplot_data_newdev, boxplot_data_original)

# -------------------------
# PLOT BOXPLOT
# -------------------------

ggplot(boxplot_data, aes(x = Pair, y = Value, fill = Difference_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = c("skyblue", "lightcoral")) +
  labs(
    title = "Differences Between True and Estimated Gini Correlations",
    x = "Selected Pairs (Variable Indices)",
    y = "Difference Values",
    fill = "Difference Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(boxplot_data, aes(x = Pair, y = Value, fill = Difference_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = c("skyblue", "lightcoral")) +
  labs(
    title = "Differences Between True and Estimated Gini Correlations",
    x = "Selected Pairs (Variable Indices)",
    y = "Difference Values",
    fill = "Difference Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(-0.5, 0.5))


# Create one combined dataset
plot_data <- boxplot_data %>%
  mutate(
    Method = case_when(
      Difference_Type == "truGini - newGini" ~ "TruGini_New - Newly Developed",
      Difference_Type == "truGini - originalGini" ~ "TruGini_New - Original",
      TRUE ~ Difference_Type
    )
  )




###########
#########
library(ggplot2)
library(dplyr)

# --- Save as PDF (vector, high quality) ---
pdf("Gini_boxplot.pdf", width = 7, height = 5, useDingbats = FALSE)

# 1) Lock legend order to match your data order
boxplot_data$Difference_Type <- factor(
  boxplot_data$Difference_Type,
  levels = c("truGini - newGini", "truGini - originalGini")
)

# 2) Clean pair labels: enforce consistent (i,j) format
pairs_clean <- gsub("[()]", "", boxplot_data$Pair)
pretty_pairs <- sprintf("(%s)", pairs_clean)
boxplot_data$Pair <- factor(pretty_pairs, levels = unique(pretty_pairs))

# 3) Plot
ggplot(boxplot_data, aes(Pair, Value, fill = Difference_Type)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(
    breaks = c("truGini - newGini", "truGini - originalGini"),
    labels = c(
      # true Γ* − Γ̂*_{ZILN}
      expression("true " * Gamma^"*" - hat(Gamma)[ZI]^"*"),
      # true Γ* − Γ̂*_{Emp}
      expression("true " * Gamma^"*" - hat(Gamma)[Org])
    ),
    values = c("#74add1", "lightcoral")
  ) +
  scale_y_continuous(
    #limits = c(-0.5, 0.5),
    breaks = seq(-0.5, 0.5, by = 0.1),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(-0.5, 0.5))+
  labs(
    title = bquote("Comparison of true " * Gamma^"*" * " differences for n=500"),
    x = "Selected pairs (indices)",
    y = "Difference values",
    fill = "Difference Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
    legend.position = c(0.98, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    legend.key.size = unit(10, "pt"),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", margin = margin(b = 6))
  )

dev.off()
