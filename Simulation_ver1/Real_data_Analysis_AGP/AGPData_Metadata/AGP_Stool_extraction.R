# ------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------
library(biomformat)
library(data.table)

# ------------------------------------------------------------------------
# 1. Read BIOM file
# ------------------------------------------------------------------------
# Replace with actual file path
biom_file <- "deblur_125nt_no_blooms.biom"

# Read HDF5 BIOM
x <- read_hdf5_biom(biom_file)

# ------------------------------------------------------------------------
# 2. Convert BIOM to OTU matrix
# ------------------------------------------------------------------------
otu_list <- x$data
otu <- do.call(rbind, otu_list)

# Add sample names as colnames
colnames(otu) <- sapply(x$columns, function(c) c$id)

# ------------------------------------------------------------------------
# 3. Extract taxonomy strings
# ------------------------------------------------------------------------
otu_taxonomies <- sapply(x$rows, function(r) {
  paste(r$metadata$taxonomy, collapse="|")
})

# ------------------------------------------------------------------------
# 4. Read metadata
# ------------------------------------------------------------------------
metadata <- fread("10317_20250630-161401.txt")

# Ensure sample_name column is character
metadata$sample_name <- as.character(metadata$sample_name)

# ------------------------------------------------------------------------
# 5. Identify stool samples
# ------------------------------------------------------------------------
sample_types <- metadata$sample_type
names(sample_types) <- metadata$sample_name

# Find stool sample IDs
stool_samples <- names(sample_types)[sample_types == "Stool"]

# Find overlap between stool samples and biom file columns
samples_in_both <- intersect(colnames(otu), stool_samples)

# Subset the OTU table
otu_stool <- otu[, samples_in_both]

# Check dimension (full stool dataset)
dim(otu_stool)
# [1] 32954 9507 (example full size)

# ------------------------------------------------------------------------
# 6. Subset smaller dataset
# ------------------------------------------------------------------------
# Grab first 100 taxa and first 1000 stool samples
otu_small <- otu_stool[1:100, 1:1000]

# ------------------------------------------------------------------------
# 7. Create a single tidy data frame with taxonomy
# ------------------------------------------------------------------------
# Convert to data frame
otu_small_df <- as.data.frame(otu_small)

# Add TaxonID (1 to 100)
otu_small_df$TaxonID <- 1:nrow(otu_small_df)

# Add Taxonomy column
otu_small_df$Taxonomy <- otu_taxonomies[1:nrow(otu_small_df)]

# Reorder columns: TaxonID, Taxonomy, then sample columns
otu_small_with_taxonomy <- otu_small_df[, c("TaxonID", "Taxonomy", colnames(otu_small_df)[1:ncol(otu_small)])]

# Check dimensions
dim(otu_small_with_taxonomy)
# Should be:
# 100 rows x (2 + 1000) columns = 100 rows x 1002 columns

# Preview first few rows
head(otu_small_with_taxonomy[, 1:6])

# ------------------------------------------------------------------------
# 8. Save as .csv and .RData
# ------------------------------------------------------------------------

# Save as CSV
write.csv(
  otu_small_with_taxonomy,
  file = "otu_small_with_taxonomy.csv",
  row.names = FALSE
)

# Save as RData
save(
  otu_small_with_taxonomy,
  file = "otu_small_with_taxonomy.RData"
)
