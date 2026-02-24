#Tissue specific diff methylation, taken from DMR intersected with a annotated 
#gene or a promoter region
# Load DMRs and genomic features
dmr <- read.csv("/path/to/file")
gf <- read.table("/path/to/file", 
                 sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)

# Initialize empty data frames for upregulated and downregulated DMRs
results_up <- data.frame()   # For upregulated DMRs (+ % methylation)
results_down <- data.frame() # For downregulated DMRs (- % methylation)

# Loop through each DMR
for (i in 1:nrow(dmr)) {
  # Extract the chromosome and define flanking regions
  cur <- gf[gf$V1 == dmr$chr[i], ]
  cur <- cur[(dmr$start[i] >= (cur$V4 - 1000)) & (dmr$end[i] <= (cur$V5 + 1000)), ]
  
  if (dmr$meth.diff[i] > 0) {
    results_down <- rbind(results_down, cur)  # methylation up → downregulation
  } else if (dmr$meth.diff[i] < 0) {
    results_up <- rbind(results_up, cur)      # methylation down → upregulation
  }
}

# Filter for only gene rows *and* drop NA
gene_lines_up <- results_up[results_up$V3 == "gene" & !is.na(results_up$V9), ]
gene_ids_up <- sub("ID=gene-([^;]+);.*", "\\1", gene_lines_up$V9)

gene_lines_down <- results_down[results_down$V3 == "gene" & !is.na(results_down$V9), ]
gene_ids_down <- sub("ID=gene-([^;]+);.*", "\\1", gene_lines_down$V9)
 
gene_ids_up_unique <- unique(gene_ids_up)
gene_ids_down_unique <- unique(gene_ids_down)


write.table(gene_ids_up_unique, file = "/path/to/file",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(gene_ids_down_unique, file ="/path/to/file",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)