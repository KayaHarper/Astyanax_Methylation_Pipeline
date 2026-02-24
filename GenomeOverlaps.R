#Tissue specific diff methylation, taken from DMR intersected with a annotated 
#gene or a promoter region
#load necessary libraries
library(methylKit)
library(genomation)
library(GenomicRanges)
library(rtracklayer)  # for reading GFF files
library(rstudioapi)
library(IRanges)


selected_folder <- selectDirectory(caption = "Name of your data input folder") #data folder select 



file.list <- list.files(path = selected_folder, pattern = "*.bismark.cov", full.names = TRUE) # makes list of .bismark cov files in folder
print(file.list)


if (length(file.list) > 0) {
  sample_ids <- as.list(paste0("SF_Liver", seq_along(file.list)))  
  
  treatment_vector <- c(1, 1, 1, 1, 1, 1)  
  print(treatment_vector)
  
  myobj <- methRead( 
    location = as.list(file.list), 
    sample.id = sample_ids,
    pipeline = "bismarkCoverage",
    assembly = "GCF_023375975.1",  
    treatment = treatment_vector,  
    mincov = 10  
  )
  
  for (i in seq_along(myobj)) {  # corrected loop
    myobj[[i]] <- myobj[[i]][!grepl("^(JALAHU010000|NW_0260400)", myobj[[i]]$chr), ]
  }
  
  print(myobj)
} else {
  print("No .bismark.cov files found in the selected folder.")
}


# filter methylation data based on read coverage
edited.cov.filt <- filterByCoverage(myobj, lo.count=1, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

#  merge data into a single object
meth <- unite(edited.cov.filt, destrand=FALSE)

coverage_sum <- meth$coverage1+meth$coverage2+meth$coverage3+meth$coverage4+meth$coverage5+meth$coverage6

# Summing all numCs columns
numCs_sum <- meth$numCs1 + meth$numCs2 + meth$numCs3 + meth$numCs4 + meth$numCs5 + meth$numCs6

# Summing all numTs columns
numTs_sum <- meth$numTs1 + meth$numTs2 + meth$numTs3 + meth$numTs4 + meth$numTs5 + meth$numTs6

methylDF <- meth[ ,1:4]

# Combining these columns into the meth data frame
meth$coverage_sum <- coverage_sum
meth$numCs_sum <- numCs_sum
meth$numTs_sum <- numTs_sum

# View the updated data frame
head(meth)


methGrange <- as(meth, "GRanges")

# load and update bed file (CpG island features) with correct chromosome names

# load the first bed file with CpG islands
bed_CpG <- read.table(file.choose(), 
                   header = FALSE, sep = "\t", 
                   col.names = c("chr", "start", "end", "name"),
                   stringsAsFactors = FALSE)

# load the second BED file (e.g., gene annotations)

bed_gene <- read.table(file.choose(), 
                       header = FALSE, sep = "\t", 
                       col.names = c("chr", "start", "end", "name"),
                       stringsAsFactors = FALSE)


#gff_file <- file.choose()  # Choose GFF file with gene annotations
# create GRanges for the first BED file
bedGrange1 <- GRanges(seqnames = bed_CpG$chr, 
                      ranges = IRanges(start = bed_CpG$start, end = bed_CpG$end), 
                      name = bed_CpG$name)

# create GRanges for the second BED file
bedGrange2 <- GRanges(seqnames = bed_gene$chr, 
                      ranges = IRanges(start = bed_gene$start, end = bed_gene$end), 
                      name = bed_gene$name)

# ecxtend CpG island regions by 4000 bp flanks on both sides
extended_bedGrange <- resize(bedGrange1, width = width(bedGrange1) + 4000, fix = "center")  # Corrected variable

# find overlaps between methylation data 
overlaps <- findOverlaps(methGrange, extended_bedGrange, minoverlap = 1, maxgap = 40000)

# extract methylated regions that overlap with extended CpG islands
methylated_cpg_regions <- methGrange[queryHits(overlaps)]

#  which genes overlap with these methylated CpG regions using the second BED file
gene_overlaps <- findOverlaps(extended_bedGrange, bedGrange2, maxgap = 4000) 

# extract the overlapping genes
overlapping_genes <- bedGrange2[subjectHits(gene_overlaps)]

#methylated_gene_list <- unique(overlapping_genes$seqnames)  # Adjust column name as needed

# print or save the list of methylated genes
print(overlapping_genes)
print(methylated_gene_list)
dataFrame <- as.data.frame(overlapping_genes)

#save as csv file 
write.csv(overlapping_genes, file = "Methylated_Gene_List.csv", row.names = FALSE)


