#Determining the Diff Mentholated Cytosine & Diff Mentholated Regions from 
#the .cov files generated from our Bismark pipeline 
#Load required libraries 

library(methylKit)
library(genomation)
library(UpSetR)
library(ComplexHeatmap)
library(tcltk)
library(GenomicRanges)
library(rtracklayer)  # for reading GFF files
library(rstudioapi)
library(IRanges)
library(circlize)
library(magick)#  methylation analysis

# select the directory containing .bismark.cov files (must be one file, with all test files inside)
selected_folder <- selectDirectory(
  caption = "Select Directory", 
)
# lists all files in selected folder (double check file list)
file.list <- list.files(path = selected_folder, pattern = "*.bismark.cov", full.names = TRUE)

extracted_names <- sub(".*/", "", file.list)  # Remove path
extracted_names <- sub("^([^_]*)_[^_]*_([^_]*)_.*$", "\\1_\\2", extracted_names)
sample.ids <- c(extracted_names)
name_prefixes <- sub("-.*$", "", extracted_names)  # Extract everything before the first dash
print(name_prefixes)
# Create a treatment vector with 1 for unique names and 0 otherwise
unique_names <- unique(name_prefixes)
if (length(unique_names) == 1) {
  # Use the second character of the section after the dash
  section_after_dash <- sub("^[^-]*-", "", extracted_names)
  second_char <- substr(section_after_dash, 2, 2)
  treatment <- ifelse(second_char == "F", 0, 1)
  print("option1")
} else if (unique_names[1] == 'SF'){
  treatment <- as.numeric(name_prefixes == unique_names[2])
} else if (unique_names[2] == 'SF'){
  treatment <- as.numeric(name_prefixes == unique_names[1])
} else {
  # Create a treatment vector with 1 for the first unique name and 0 otherwise
  treatment <- as.numeric(name_prefixes == unique_names[1])
  print("catchall")
}
print(sample.ids)
print(treatment)

# read and combine all coverage files into a list (methylRaw objects) 
methylRawList.obj <- lapply(seq_along(file.list), function(i) {
  obj <- methRead(file.list[i],              # File path (file.list)
                  sample.id = sample.ids[i], # reads in sample ID (above)
                  assembly = "GCA_023375975.1", # Genome assembly used (Housekeeping, but it is required)
                  treatment = treatment[i],  # Treatment condition (treatment condition above)
                  context = "CpG",           # Methylation context (CpG),  housekeeping 
                  pipeline = "bismarkCoverage",
                  mincov = 6)
  # Pipeline used to generate the data,also housekeeping 
  # Remove unwanted contigs from the data
  obj <- obj[!grepl("^(JALAHU010000|NW_0260400)", obj$chr), ] 
  return(obj)  # Return as a filtered object
})

# turns my  list of all  methylRaw objects to methylRawList objects, required for calculateDiffMeth
methylRawList.obj <- methylRawList(methylRawList.obj, treatment = treatment)

# filters out extreme coverage values (our min low count is 10 hi count is 99.9)
filtered.myobj <- filterByCoverage(methylRawList.obj, lo.count = 6, hi.perc = 99.9)

# scales the read coverage values for each sample to a common scale across the samples
normalized.myobj <- normalizeCoverage(filtered.myobj)

# combine all samples into a single object for differential methylation analysis
methylBase.obj <- unite(normalized.myobj, destrand = FALSE, min.per.group = 2L)

# calculate differential methylation at cytosine level (DMCs) MN adjusts for overdispersion (higher varibility then the model assumes) 
#across biological replicates by using a Negative Binomial model.
#fast.fisher recommended by methylKit because it is statistical robust in small samples and not comp. intensive.

# calculate differential methylation at regions  (DMRs), with window size 1000 and step size 500
dmr <- tileMethylCounts(methylBase.obj, win.size = 1000,  # Window size (1kb)
                        step.size = 500,                 # Step size (500bp)
                        cov.bases = 10)                 # minimum coverage threshold (also set before as 10)

# calculate differential methylation for the tiled regions
dmr_diff <- calculateDiffMeth(dmr, test = "fast.fisher", overdispersion = "MN")

# filter significant DMRs based on methylation difference and p-value thresholds (15 - 30% dif is standard, so is p = 0.05 or q =0.01)
dmr_filtered <- getMethylDiff(dmr_diff, difference = 15, qvalue = 0.05)

# Save the filter results to a CSV file
write.csv(dmr_filtered, "/Users/kayaharper/Desktop/Liver/DMR/dmr_filltered_SFL_SSL_0.05.csv", row.names = FALSE)

# Save the unfiltered results to a CSV file
write.csv(dmr_diff, "/Users/kayaharper/Desktop/Liver/DMR/dmr_diff_SFL_SSL_0.05.csv", row.names = FALSE)



# # Convert DMR filtered results to BED format and save
# dmr_gr <- GRanges(seqnames = dmr_filtered$chr,
#                   ranges = IRanges(start = dmr_filtered$start, end = dmr_filtered$end))
# export(dmr_gr, "dmr_pachon_liver_fed_starved.bed")


#_________________________________________________________________________________

# 
# #DMR HEATMAP
# ## Generate a per-sample methylation matrix
# meth_matrix_samples <- percMethylation(methylBase.obj)
# 
# # Subset the matrix for DMRs
# dmr_rows <- as.integer(rownames(dmr_filtered))  # Get row indices of DMRs in the original data
# subset_dmr_matrix <- meth_matrix_samples[dmr_rows, ]
# 
# # Add annotations for treatment
# # Correcting the HeatmapAnnotation for DMRs
# treatment_annotation_DMR <- HeatmapAnnotation(
#   Treatment = factor(treatment, levels = c(0, 1), labels = c("fed", "starved")),
#   col = list(Treatment = c("fed" = "skyblue", "starved" = "pink"))
# )
# 
# 
# # Create the heatmap for DMRs
# heatmap_DMR <- Heatmap(
#   subset_dmr_matrix,
#   name = "Molino Fed Liver vs Brain - DMRs",
#   top_annotation = treatment_annotation_DMR,
#   show_row_names = FALSE,
#   show_column_names = TRUE,
#   col = colorRamp2(c(0, 50, 100), c("purple", "white", "darkblue")),
#   cluster_rows = TRUE,
#   cluster_columns = TRUE,
#   heatmap_legend_param = list(title = "Methylation %", legend_direction = "horizontal",
#                               main = "Differential Methylation Across Liver and Brain Tissues")
# )
# 
# # Draw the heatmap
# draw(heatmap_DMR)
# 
# #_________________________________________________________________________________
# 
# 
# 
# 
# 







dmc_diff <- calculateDiffMeth(methylBase.obj, overdispersion = "MN", test = "fast.fisher")

# filter significant DMCs based on methylation difference and p-value thresholds (15 - 30% dif is standard, so is p = 0.05 or q =0.01)
dmc_filtered <- getMethylDiff(dmc_diff, difference = 15, qvalue = 0.05)

# saves the results to a CSV file
write.csv(dmc_filtered, "dmc_results_starved_fed_surface_Liver.csv", row.names = FALSE)

# Convert DMC filtered results to BED format and save
dmc_gr <- GRanges(seqnames = dmc_filtered$chr,
                  ranges = IRanges(start = dmc_filtered$start, end = dmc_filtered$end))
export(dmc_gr, "dmc_results_starved_fed_surface_Liver.bed")

# 
# 
# #DMC HEATMAP
# # generate a matrix of percentage methylation values for further analysis (not sure if i will use)
# meth_matrix <- percMethylation(methylBase.obj) 
# 
# # Subset the matrix if needed (e.g., top DMCs/DMRs based on significance or other criteria)
# # For example, use rows from dmc_filtered
# dmc_rows <- as.integer(rownames(dmc_filtered))  # Get row indices of DMCs in the original data
# subset_meth_matrix <- meth_matrix[dmc_rows, ]
# 
# # Add annotations (optional)
# # Create annotation bars for treatment conditions
# # Correct the HeatmapAnnotation to reflect Surface Fish (SF) as control and Mol as test
# treatment_annotation_DMC <- HeatmapAnnotation(
#   Treatment = factor(treatment, levels = c(0, 1), labels = c("fed", "starved")),
#   col = list(Treatment = c("fed" = "seagreen", "starved" = "pink"))
# )
# # Create the heatmap
# heatmap_DMC <- Heatmap(
#   subset_meth_matrix,
#   name = "Molino Fed Liver vs Brain",
#   top_annotation = treatment_annotation_DMC,
#   show_row_names = FALSE,  # Set to TRUE if you want to see the genomic regions/cytosines
#   show_column_names = TRUE, # Set to TRUE if you want to see sample IDs
#   col = colorRamp2(c(0, 50, 100), c("white", "skyblue", "darkblue")), # Color scale: blue = low, red = high methylation
#   cluster_rows = TRUE, # Cluster rows based on similarity
#   cluster_columns = TRUE, # Cluster columns (samples)
#   
#   heatmap_legend_param = list(title = "Methylation %", legend_direction = "horizontal")
# )
# 
# # Draw the heatmap
# draw(heatmap_DMC)
# 




# Histogram of the density of the methylation windows  

library(ggplot2)
# 0) Make sure you're using the *unfiltered* differential object
#    (your 'dmr_diff' from calculateDiffMeth on the tiled object)
dd <- as.data.frame(getData(dmr_diff))  # should have 'meth.diff'

# 1) Diagnostics
summary(dd$meth.diff)
table(sign(dd$meth.diff))      # expect some -1 and +1 (and 0s)
range(dd$meth.diff, na.rm=TRUE)

# 2) If you used getMethylDiff() earlier, ensure type="all":
# dmr_all <- getMethylDiff(dmr_diff, difference = 0, qvalue = 1, type = "all")
# dd <- as.data.frame(getData(dmr_all))

# 3) Plot the *signed* differences (NOT abs)
library(ggplot2)
ggplot(dd, aes(x = meth.diff)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60) +
  labs(x = "Window methylation difference (%)  [treatment=1 minus treatment=0]",
       y = "Density",
       title = "Distribution of per-window methylation differences") +
  coord_cartesian(xlim = c(-60, 60)) +   # adjust as needed
  theme_bw()


# 
# 
# 
