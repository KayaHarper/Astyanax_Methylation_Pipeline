#Calculating the perc change in methylation per CpG site across the genome in all 
#feeding conditions, tissues, and pop
# Load necessary libraries
library(data.table)
library(pheatmap)
library(doParallel)
library(parallelly) # For availableCores()
library(wesanderson)
library(circlize)

# Define paths to cov files
cov_files <- list.files("path/to/cov/files",
                        pattern = "\\.cov(\\.gz)?$", full.names = TRUE)

# Step 1: Function to read and clean a single .cov file
read_cov_file <- function(file) {
  df <- fread(file, header = FALSE, col.names = c("Chromosome", "Start", "End",
                                                  "Methylation_Percent",
                                                  "Methylated_Reads", "Unmethylated_Reads"))
  df$CpG_ID <- paste(df$Chromosome, df$Start, sep = ":")
  df$Coverage <- df$Methylated_Reads + df$Unmethylated_Reads
  return(df[, .(CpG_ID, Methylation_Percent, Coverage)])
}

# Step 2: Group files based on their names
mol_Files_Fed <- grep("Mol-.*F", cov_files, value = TRUE)
mol_Files_Starved <- grep("Mol-.*S", cov_files, value = TRUE)
pa_Files_Fed <- grep("PA-.*F", cov_files, value = TRUE)
pa_Files_Starved <- grep("PA-.*S", cov_files, value = TRUE)
sf_Files_Fed <- grep("SF-.*F", cov_files, value = TRUE)
sf_Files_Starved <- grep("SF-.*S", cov_files, value = TRUE)

# Combine groups into a list
file_groups <- list(
  "Mol_Fed" = mol_Files_Fed,
  "Mol_Starved" = mol_Files_Starved,
  "PA_Fed" = pa_Files_Fed,
  "PA_Starved" = pa_Files_Starved,
  "SF_Fed" = sf_Files_Fed,
  "SF_Starved" = sf_Files_Starved
)

# Step 3: Set up parallel backend
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Step 4: Read files in parallel and process
read_group <- function(group_files) {
  foreach(file = group_files,
          .packages = "data.table",
          .export = c("read_cov_file")) %dopar% {
            read_cov_file(file)
          }
}

# Read and process all groups
processed_data <- lapply(file_groups, read_group)

# Stop the cluster
stopCluster(cl)

# Step 5: Process conditions
process_condition <- function(data_list) {
  if (length(data_list) == 0) return(NULL)
  
  # Merge replicate data
  # Rename columns before merging to avoid duplicates
  data_list <- lapply(seq_along(data_list), function(i) {
    setnames(data_list[[i]],
             old = c("Methylation_Percent", "Coverage"),
             new = c(paste0("Methylation_Percent_", i), paste0("Coverage_", i)))
    data_list[[i]]
  })
  
  # Merge replicate data
  merged_data <- Reduce(function(x, y) merge(x, y, by = "CpG_ID", all = TRUE), data_list)
  
  
  # Calculate pooled coverage and average methylation
  coverage_cols <- grep("^Coverage", colnames(merged_data), value = TRUE)
  methylation_cols <- grep("^Methylation_Percent", colnames(merged_data), value = TRUE)
  
  merged_data$Pooled_Coverage <- rowSums(merged_data[, ..coverage_cols], na.rm = TRUE)
  merged_data$Average_Methylation <- rowMeans(merged_data[, ..methylation_cols], na.rm = TRUE)
  
  # Filter rows
  filtered_data <- merged_data[!is.na(merged_data$Average_Methylation) &
                                 merged_data$Pooled_Coverage >= 10, ]
  return(filtered_data[, .(CpG_ID, Average_Methylation)])
}

# Apply filtering to all groups
final_data <- lapply(processed_data, process_condition)

# Remove NULL results
final_data <- Filter(Negate(is.null), final_data)

# Rename columns to include condition names
condition_names <- names(processed_data)  # Assuming processed_data has named conditions
final_data <- lapply(seq_along(final_data), function(i) {
  setnames(final_data[[i]], old = "Average_Methylation", new = paste0("Average_Methylation_", condition_names[i]))
  final_data[[i]]
})

# Combine results into a heatmap-ready matrix
heatmap_data <- Reduce(function(x, y) {
  merge(x, y, by = "CpG_ID", all = TRUE)
}, final_data)



# Convert data to matrix
heatmap_matrix <- as.matrix(heatmap_data[, -1, with = FALSE])
rownames(heatmap_matrix) <- heatmap_data$CpG_ID


# Step 7: Create heatmap with custom labels and angle for both rows and columns
# Custom row names
custom_row_labels <- c("Molino Fed", "Molino Starved", "Pachon Fed",
                       "Pachon Starved", "Surface Fed", "Surface Starved")

# Create matrices for correlation coefficients and p-values
corvals <- pvals <- matrix(NA, 6, 6)
colnames(corvals) <- colnames(pvals) <- custom_row_labels
row.names(corvals) <- row.names(pvals) <- custom_row_labels


options(digits = 10)
# Compute correlations and p-values for the upper triangle only
for (i in 1:6) { # Rows
  for (j in i:6) { # Upper triangle (including diagonal)
    x <- heatmap_matrix[, i]
    y <- heatmap_matrix[, j]
    yg <- y[!is.na(y) & !is.na(x)]
    xg <- x[!is.na(y) & !is.na(x)]
    res <- cor.test(yg, xg)
    pvals[i, j] <- res$p.value
    corvals[i, j] <- res$estimate
  }
}


# Mask the lower triangle of the correlation matrix
corvals[lower.tri(corvals)] <- NA
pvals[lower.tri(pvals)] <- NA
# Create a combined matrix for p-values and correlations
masked_correlations <- corvals
masked_correlations[pvals > 0.05] <- NA  # Mask correlations with p > 0.05


# Define a custom three-color gradient
col_fun <- colorRamp2(c(min(masked_correlations, na.rm = TRUE), 0, max(masked_correlations, na.rm = TRUE)),
                      c("#E6A0C4", "#C6CDF7", "#7294D4"))

# Generate the heatmap with the new color scale
pheatmap(
  mat = masked_correlations,  # Only significant correlations are shown in color
  color = col_fun(seq(min(masked_correlations, na.rm = TRUE), max(masked_correlations, na.rm = TRUE), length.out = 50)),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = round(corvals, 2),  # Show all correlation coefficients (rounded to 2 decimals)
  number_color = "black",
  na_col = "grey",  # Non-significant correlations are grayed out
  fontsize_number = 10,
  main = "Correlations of % CpG Sites in Liver Tissue",
  angle_col = 45  # Rotate column labels by 45 degrees
) 