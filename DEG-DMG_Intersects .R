#DEG-DMG-inersects aka the code that makes the bar graph of all 3 tissues, 
#feeding condtions and populations comparisons on one graph 
# ---- LOAD LIBRARIES ----

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(apeglm)
library(clusterProfiler)
library(org.Dr.eg.db)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(patchwork)




# ---- SETTINGS ----
morph_labels <- c("mol" = "Molino", "pa" = "Pachón", "sf" = "Rio Choy")
morph_colors <- c("Molino" = "#ce9e91", "Pachón" = "#E57E81", "Rio Choy" = "#BEC5DF")
anno_colors <- list(
  Morph = morph_colors,
  Condition = c("Fed" = "#6cbfc2", "Starved" = "#C43A1A")
)
shape_map <- c("Rio Choy" = 16, "Pachón" = 17, "Molino" = 15)
output_dir <- "/Users/kayaharper/Desktop/GSEA_Results"
dir.create(output_dir, showWarnings = FALSE)

# ---- LOAD & CLEAN COUNT DATA ----
counts <- read.table("/path/to/gene_counts.txt",
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1, check.names = FALSE)
counts <- counts[, -(1:5)]
colnames(counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", basename(colnames(counts)))
# Define custom shapes for each tissue
tissue_shape_map <- c("liver" = 23, "brain" = 25, "muscle" = 4)

parse_sample <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  morph <- tolower(parts[1])
  fish_id <- parts[2]
  tissue <- tolower(parts[3])
  condition <- ifelse(grepl("F$", fish_id), "Fed", "Starved")
  return(c(morph = morph, fish_id = fish_id, tissue = tissue, condition = condition))
}

sample_info <- t(sapply(colnames(counts), parse_sample))
sampleTable <- as.data.frame(sample_info)
rownames(sampleTable) <- colnames(counts)




run_DE_by_tissue_morph <- function(dds, tissue_name, morph_name) {
  dds_sub <- dds[, dds$tissue == tissue_name & dds$morph == morph_name]
  dds_sub$tissue <- droplevels(dds_sub$tissue)
  dds_sub$morph <- droplevels(dds_sub$morph)
  dds_sub$condition <- droplevels(dds_sub$condition)
  
  if (length(unique(dds_sub$condition)) < 2) {
    stop(paste("Not enough conditions in", morph_name, tissue_name))
  }
  
  design(dds_sub) <- ~ condition
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("condition", "Fed", "Starved"))
  return(res)
}

# ---- DESEQ2 ANALYSIS ----
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleTable,
                              design = ~ morph + tissue + condition)
dds <- DESeq(dds)
# ==== DEG–DMG OVERLAP ANALYSIS: LIVER ====

# -- Prep data subsets --
bg_genes <- rownames(counts)

# ---- Helper Function: Apply padj and log2FC cutoffs ----
get_sig_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  res <- na.omit(res)
  res[which(res$padj < padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff), ]
}

# ---- Fisher's Exact Test Helper ----
run_fisher <- function(deg, dmg, background) {
  overlap <- length(intersect(deg, dmg))
  only_deg <- length(setdiff(deg, dmg))
  only_dmg <- length(setdiff(dmg, deg))
  neither <- length(setdiff(background, union(deg, dmg)))
  matrix <- matrix(c(overlap, only_deg, only_dmg, neither), nrow = 2)
  fisher.test(matrix)
}

# ---- Define consistent colors for all DEG–DMG bar plots ----
segment_colors <- c(
  "DEG only" = "#6cbfc2",    # teal-blue
  "Both"     = "#8888AA",    # muted purple
  "DMG only" = "#E57E81"     # pink-red
)


liver_samples <- rownames(sampleTable)[sampleTable$tissue == "liver"]
liver_counts <- counts[, liver_samples]
liver_sampleTable <- sampleTable[liver_samples, ]

# -- Fed vs Starved within each morph --
liver_dds <- DESeqDataSetFromMatrix(countData = liver_counts, colData = liver_sampleTable, design = ~ morph + condition)
liver_dds <- DESeq(liver_dds)

get_DEG_metabolic <- function(dds, morph_name) {
  dds_sub <- dds[, dds$morph == morph_name]
  dds_sub$condition <- droplevels(dds_sub$condition)
  design(dds_sub) <- ~ condition
  dds_sub <- DESeq(dds_sub)
  results(dds_sub, contrast = c("condition", "Fed", "Starved"))
}

deg_liver_metabolic <- list(
  Molino = get_sig_genes(get_DEG_metabolic(liver_dds, "mol")),
  Pachon = get_sig_genes(get_DEG_metabolic(liver_dds, "pa")),
  Surface = get_sig_genes(get_DEG_metabolic(liver_dds, "sf"))
)

dmg_liver_metabolic <- list(
  Molino = read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1,
  Pachon = read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1,
  Surface = read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
)

fisher_liver_metabolic <- list(
  Molino = run_fisher(rownames(deg_liver_metabolic$Molino), dmg_liver_metabolic$Molino, bg_genes),
  Pachon = run_fisher(rownames(deg_liver_metabolic$Pachon), dmg_liver_metabolic$Pachon, bg_genes),
  Surface = run_fisher(rownames(deg_liver_metabolic$Surface), dmg_liver_metabolic$Surface, bg_genes)
)

# Build overlap_counts_met
overlap_counts_met <- data.frame(
  Comparison = rep(c("Molino (FvS)", "Pachón (FvS)", "Rio Choy (FvS)"), each = 3),
  Category = rep(c("DEG only", "DMG only", "Both"), 3),
  Count = c(
    length(setdiff(rownames(deg_liver_metabolic$Molino), dmg_liver_metabolic$Molino)),
    length(setdiff(dmg_liver_metabolic$Molino, rownames(deg_liver_metabolic$Molino))),
    length(intersect(rownames(deg_liver_metabolic$Molino), dmg_liver_metabolic$Molino)),
    
    length(setdiff(rownames(deg_liver_metabolic$Pachon), dmg_liver_metabolic$Pachon)),
    length(setdiff(dmg_liver_metabolic$Pachon, rownames(deg_liver_metabolic$Pachon))),
    length(intersect(rownames(deg_liver_metabolic$Pachon), dmg_liver_metabolic$Pachon)),
    
    length(setdiff(rownames(deg_liver_metabolic$Surface), dmg_liver_metabolic$Surface)),
    length(setdiff(dmg_liver_metabolic$Surface, rownames(deg_liver_metabolic$Surface))),
    length(intersect(rownames(deg_liver_metabolic$Surface), dmg_liver_metabolic$Surface))
  )
)

# Add FvS sig labels
sig_labels_met <- data.frame(
  Comparison = c("Molino (FvS)", "Pachón (FvS)", "Rio Choy (FvS)"),
  p_value = c(fisher_liver_metabolic$Molino$p.value,
              fisher_liver_metabolic$Pachon$p.value,
              fisher_liver_metabolic$Surface$p.value)
)
sig_labels_met$label <- cut(sig_labels_met$p_value,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            labels = c("***", "**", "*", ""),
                            right = FALSE)
total_counts_met <- overlap_counts_met %>%
  group_by(Comparison) %>%
  summarise(total = sum(Count))
sig_labels_met <- left_join(sig_labels_met, total_counts_met, by = "Comparison")
sig_labels_met$ypos <- sig_labels_met$total + 20

# -- Fed-only cross-morph comparisons --
fed_liver_samples <- rownames(sampleTable)[sampleTable$condition == "Fed" & sampleTable$tissue == "liver"]
fed_liver_counts <- counts[, fed_liver_samples]
fed_liver_sampleTable <- sampleTable[fed_liver_samples, ]

fed_liver_dds <- DESeqDataSetFromMatrix(countData = fed_liver_counts,
                                        colData = fed_liver_sampleTable,
                                        design = ~ morph)
fed_liver_dds <- DESeq(fed_liver_dds)

deg_liver_morph <- list(
  MvP = get_sig_genes(results(fed_liver_dds, contrast = c("morph", "mol", "pa"))),
  MvS = get_sig_genes(results(fed_liver_dds, contrast = c("morph", "mol", "sf"))),
  PvS = get_sig_genes(results(fed_liver_dds, contrast = c("morph", "pa", "sf")))
)

dmg_liver_morph <- list(
  MvP = read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1,
  MvS = read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1,
  PvS = read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
)

fisher_liver_morph <- list(
  MvP = run_fisher(rownames(deg_liver_morph$MvP), dmg_liver_morph$MvP, rownames(fed_liver_dds)),
  MvS = run_fisher(rownames(deg_liver_morph$MvS), dmg_liver_morph$MvS, rownames(fed_liver_dds)),
  PvS = run_fisher(rownames(deg_liver_morph$PvS), dmg_liver_morph$PvS, rownames(fed_liver_dds))
)

# Bar segments for morph comparisons
overlap_counts_morph <- data.frame(
  Comparison = rep(c("Molino vs Pachón", "Molino vs Surface", "Pachón vs Surface"), each = 3),
  Category = rep(c("DEG only", "DMG only", "Both"), 3),
  Count = c(
    length(setdiff(rownames(deg_liver_morph$MvP), dmg_liver_morph$MvP)),
    length(setdiff(dmg_liver_morph$MvP, rownames(deg_liver_morph$MvP))),
    length(intersect(rownames(deg_liver_morph$MvP), dmg_liver_morph$MvP)),
    
    length(setdiff(rownames(deg_liver_morph$MvS), dmg_liver_morph$MvS)),
    length(setdiff(dmg_liver_morph$MvS, rownames(deg_liver_morph$MvS))),
    length(intersect(rownames(deg_liver_morph$MvS), dmg_liver_morph$MvS)),
    
    length(setdiff(rownames(deg_liver_morph$PvS), dmg_liver_morph$PvS)),
    length(setdiff(dmg_liver_morph$PvS, rownames(deg_liver_morph$PvS))),
    length(intersect(rownames(deg_liver_morph$PvS), dmg_liver_morph$PvS))
  )
)

# Morph p-values
sig_labels_morph <- data.frame(
  Comparison = c("Molino vs Pachón", "Molino vs Surface", "Pachón vs Surface"),
  p_value = c(fisher_liver_morph$MvP$p.value, fisher_liver_morph$MvS$p.value, fisher_liver_morph$PvS$p.value)
)
sig_labels_morph$label <- cut(sig_labels_morph$p_value,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                              labels = c("***", "**", "*", ""),
                              right = FALSE)
total_counts_morph <- overlap_counts_morph %>%
  group_by(Comparison) %>%
  summarise(total = sum(Count))
sig_labels_morph <- left_join(sig_labels_morph, total_counts_morph, by = "Comparison")
sig_labels_morph$ypos <- sig_labels_morph$total + 20

# -- Combine all liver data for final plot --
overlap_counts_liver_all <- rbind(overlap_counts_met, overlap_counts_morph)
sig_labels_liver_all <- rbind(sig_labels_met, sig_labels_morph)

# Rename & reorder
overlap_counts_liver_all$Comparison <- recode(overlap_counts_liver_all$Comparison,
                                              "Surface (FvS)" = "Rio Choy (FvS)",
                                              "Molino vs Surface" = "Rio Choy vs Molino",
                                              "Pachón vs Surface" = "Rio Choy vs Pachón")
sig_labels_liver_all$Comparison <- recode(sig_labels_liver_all$Comparison,
                                          "Molino vs Surface" = "Rio Choy vs Molino",
                                          "Pachón vs Surface" = "Rio Choy vs Pachón")

ordered_comparisons <- c("Rio Choy (FvS)", "Pachón (FvS)", "Molino (FvS)",
                         "Rio Choy vs Pachón", "Rio Choy vs Molino", "Molino vs Pachón")

overlap_counts_liver_all$Comparison <- factor(overlap_counts_liver_all$Comparison, levels = ordered_comparisons)
sig_labels_liver_all$Comparison <- factor(sig_labels_liver_all$Comparison, levels = ordered_comparisons)
overlap_counts_liver_all$Category <- factor(overlap_counts_liver_all$Category, levels = c("DEG only", "Both", "DMG only"))

# Done! Now you can use `overlap_counts_liver_all` and `sig_labels_liver_all` to plot.





#Brain 



# ----- Brain Tissue DEG–DMG Overlap -----

# ---- Fed vs Starved ----
brain_samples <- rownames(sampleTable)[sampleTable$tissue == "brain"]
brain_counts <- counts[, brain_samples]
brain_sampleTable <- sampleTable[brain_samples, ]
brain_dds <- DESeqDataSetFromMatrix(countData = brain_counts, colData = brain_sampleTable, design = ~ morph + condition)
brain_dds <- DESeq(brain_dds)

res_Molino_brain <- get_DEG_metabolic(brain_dds, "mol")
res_Pachon_brain <- get_DEG_metabolic(brain_dds, "pa")
res_Surface_brain <- get_DEG_metabolic(brain_dds, "sf")

sig_Molino_brain <- get_sig_genes(res_Molino_brain)
sig_Pachon_brain <- get_sig_genes(res_Pachon_brain)
sig_Surface_brain <- get_sig_genes(res_Surface_brain)

deg_Molino_brain <- rownames(sig_Molino_brain)
deg_Pachon_brain <- rownames(sig_Pachon_brain)
deg_Surface_brain <- rownames(sig_Surface_brain)

dmg_Molino_brain <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
dmg_Pachon_brain <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
dmg_Surface_brain <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1

fisher_Molino_brain <- run_fisher(deg_Molino_brain, dmg_Molino_brain, bg_genes)
fisher_Pachon_brain <- run_fisher(deg_Pachon_brain, dmg_Pachon_brain, bg_genes)
fisher_Surface_brain <- run_fisher(deg_Surface_brain, dmg_Surface_brain, bg_genes)

overlap_counts_brain <- data.frame(
  Comparison = rep(c("Molino (FvS)", "Pachón (FvS)", "Surface (FvS)"), each = 3),
  Category = rep(c("DEG only", "DMG only", "Both"), 3),
  Count = c(
    length(setdiff(deg_Molino_brain, dmg_Molino_brain)),
    length(setdiff(dmg_Molino_brain, deg_Molino_brain)),
    length(intersect(deg_Molino_brain, dmg_Molino_brain)),
    
    length(setdiff(deg_Pachon_brain, dmg_Pachon_brain)),
    length(setdiff(dmg_Pachon_brain, deg_Pachon_brain)),
    length(intersect(deg_Pachon_brain, dmg_Pachon_brain)),
    
    length(setdiff(deg_Surface_brain, dmg_Surface_brain)),
    length(setdiff(dmg_Surface_brain, deg_Surface_brain)),
    length(intersect(deg_Surface_brain, dmg_Surface_brain))
  )
)

sig_labels_brain <- data.frame(
  Comparison = c("Molino (FvS)", "Pachón (FvS)", "Surface (FvS)"),
  p_value = c(fisher_Molino_brain$p.value, fisher_Pachon_brain$p.value, fisher_Surface_brain$p.value)
)
sig_labels_brain$label <- cut(sig_labels_brain$p_value,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                              labels = c("***", "**", "*", ""),
                              right = FALSE)
total_counts_brain <- overlap_counts_brain %>%
  group_by(Comparison) %>%
  summarise(total = sum(Count))
sig_labels_brain <- left_join(sig_labels_brain, total_counts_brain, by = "Comparison")
sig_labels_brain$ypos <- sig_labels_brain$total + 20

# ---- Cross-Morph Comparisons ----
fed_brain_samples <- rownames(sampleTable)[sampleTable$condition == "Fed" & sampleTable$tissue == "brain"]
fed_brain_counts <- counts[, fed_brain_samples]
fed_brain_sampleTable <- sampleTable[fed_brain_samples, ]

fed_brain_dds <- DESeqDataSetFromMatrix(countData = fed_brain_counts, colData = fed_brain_sampleTable, design = ~ morph)
fed_brain_dds <- DESeq(fed_brain_dds)

res_MvP_brain <- results(fed_brain_dds, contrast = c("morph", "mol", "pa"))
res_MvS_brain <- results(fed_brain_dds, contrast = c("morph", "mol", "sf"))
res_PvS_brain <- results(fed_brain_dds, contrast = c("morph", "pa", "sf"))

sig_MvP_brain <- get_sig_genes(res_MvP_brain)
sig_MvS_brain <- get_sig_genes(res_MvS_brain)
sig_PvS_brain <- get_sig_genes(res_PvS_brain)

deg_MvP_brain <- rownames(sig_MvP_brain)
deg_MvS_brain <- rownames(sig_MvS_brain)
deg_PvS_brain <- rownames(sig_PvS_brain)

dmg_MvP_brain <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
dmg_MvS_brain <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
dmg_PvS_brain <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1

fisher_MvP_brain <- run_fisher(deg_MvP_brain, dmg_MvP_brain, bg_genes)
fisher_MvS_brain <- run_fisher(deg_MvS_brain, dmg_MvS_brain, bg_genes)
fisher_PvS_brain <- run_fisher(deg_PvS_brain, dmg_PvS_brain, bg_genes)

overlap_counts_brain_morph <- data.frame(
  Comparison = rep(c("Molino vs Pachón", "Molino vs Surface", "Pachón vs Surface"), each = 3),
  Category = rep(c("DEG only", "DMG only", "Both"), 3),
  Count = c(
    length(setdiff(deg_MvP_brain, dmg_MvP_brain)),
    length(setdiff(dmg_MvP_brain, deg_MvP_brain)),
    length(intersect(deg_MvP_brain, dmg_MvP_brain)),
    
    length(setdiff(deg_MvS_brain, dmg_MvS_brain)),
    length(setdiff(dmg_MvS_brain, deg_MvS_brain)),
    length(intersect(deg_MvS_brain, dmg_MvS_brain)),
    
    length(setdiff(deg_PvS_brain, dmg_PvS_brain)),
    length(setdiff(dmg_PvS_brain, deg_PvS_brain)),
    length(intersect(deg_PvS_brain, dmg_PvS_brain))
  )
)

sig_labels_brain_morph <- data.frame(
  Comparison = c("Molino vs Pachón", "Molino vs Surface", "Pachón vs Surface"),
  p_value = c(fisher_MvP_brain$p.value, fisher_MvS_brain$p.value, fisher_PvS_brain$p.value)
)
sig_labels_brain_morph$label <- cut(sig_labels_brain_morph$p_value,
                                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                    labels = c("***", "**", "*", ""),
                                    right = FALSE)
total_counts_brain_morph <- overlap_counts_brain_morph %>%
  group_by(Comparison) %>%
  summarise(total = sum(Count))
sig_labels_brain_morph <- left_join(sig_labels_brain_morph, total_counts_brain_morph, by = "Comparison")
sig_labels_brain_morph$ypos <- sig_labels_brain_morph$total + 20

# Combine FvS and morph comparisons
overlap_counts_brain_all <- rbind(overlap_counts_brain, overlap_counts_brain_morph)
sig_labels_brain_all <- rbind(sig_labels_brain, sig_labels_brain_morph)

# Recode and reorder
overlap_counts_brain_all$Comparison <- recode(overlap_counts_brain_all$Comparison,
                                              "Surface (FvS)" = "Rio Choy (FvS)",
                                              "Molino vs Surface" = "Rio Choy vs Molino",
                                              "Pachón vs Surface" = "Rio Choy vs Pachón"
)
sig_labels_brain_all$Comparison <- recode(sig_labels_brain_all$Comparison,
                                          "Molino vs Surface" = "Rio Choy vs Molino",
                                          "Pachón vs Surface" = "Rio Choy vs Pachón"
)
overlap_counts_brain_all$Comparison <- factor(overlap_counts_brain_all$Comparison,
                                              levels = c("Rio Choy (FvS)", "Pachón (FvS)", "Molino (FvS)",
                                                         "Rio Choy vs Pachón", "Rio Choy vs Molino", "Molino vs Pachón"))
sig_labels_brain_all$Comparison <- factor(sig_labels_brain_all$Comparison, levels = levels(overlap_counts_brain_all$Comparison))


overlap_counts_brain_all$Category <- factor(overlap_counts_brain_all$Category, 
                                             levels = c("DEG only", "Both", "DMG only"))

# PLOT
ggplot(overlap_counts_brain_all, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = segment_colors) +
  geom_text(data = sig_labels_brain_all,
            aes(x = Comparison, y = ypos, label = label),
            inherit.aes = FALSE,
            size = 6, fontface = "bold") +
  labs(
    y = "Gene Count",
    title = "Shared and Unique DEGs & DMGs in Brain (per Comparison)",
    fill = "Gene Set"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))















#muscle 

# ----- Muscle Tissue DEG–DMG Overlap -----

# ---- Fed vs Starved ----
muscle_samples <- rownames(sampleTable)[sampleTable$tissue == "muscle"]
muscle_counts <- counts[, muscle_samples]
muscle_sampleTable <- sampleTable[muscle_samples, ]
muscle_dds <- DESeqDataSetFromMatrix(countData = muscle_counts, colData = muscle_sampleTable, design = ~ morph + condition)
muscle_dds <- DESeq(muscle_dds)

res_Molino_muscle <- get_DEG_metabolic(muscle_dds, "mol")
res_Pachon_muscle <- get_DEG_metabolic(muscle_dds, "pa")
res_Surface_muscle <- get_DEG_metabolic(muscle_dds, "sf")

sig_Molino_muscle <- get_sig_genes(res_Molino_muscle)
sig_Pachon_muscle <- get_sig_genes(res_Pachon_muscle)
sig_Surface_muscle <- get_sig_genes(res_Surface_muscle)

deg_Molino_muscle <- rownames(sig_Molino_muscle)
deg_Pachon_muscle <- rownames(sig_Pachon_muscle)
deg_Surface_muscle <- rownames(sig_Surface_muscle)

dmg_Molino_muscle <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
dmg_Pachon_muscle <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
dmg_Surface_muscle <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1

bg_genes <- rownames(counts)


fisher_Molino_muscle <- run_fisher(deg_Molino_muscle, dmg_Molino_muscle, bg_genes)
fisher_Pachon_muscle <- run_fisher(deg_Pachon_muscle, dmg_Pachon_muscle, bg_genes)
fisher_Surface_muscle <- run_fisher(deg_Surface_muscle, dmg_Surface_muscle, bg_genes)


overlap_counts_muscle <- data.frame(
  Comparison = rep(c("Molino (FvS)", "Pachón (FvS)", "Surface (FvS)"), each = 3),
  Category = rep(c("DEG only", "DMG only", "Both"), 3),
  Count = c(
    length(setdiff(deg_Molino_muscle, dmg_Molino_muscle)),
    length(setdiff(dmg_Molino_muscle, deg_Molino_muscle)),
    length(intersect(deg_Molino_muscle, dmg_Molino_muscle)),
    
    length(setdiff(deg_Pachon_muscle, dmg_Pachon_muscle)),
    length(setdiff(dmg_Pachon_muscle, deg_Pachon_muscle)),
    length(intersect(deg_Pachon_muscle, dmg_Pachon_muscle)),
    
    length(setdiff(deg_Surface_muscle, dmg_Surface_muscle)),
    length(setdiff(dmg_Surface_muscle, deg_Surface_muscle)),
    length(intersect(deg_Surface_muscle, dmg_Surface_muscle))
  )
)

sig_labels_muscle <- data.frame(
  Comparison = c("Molino (FvS)", "Pachón (FvS)", "Surface (FvS)"),
  p_value = c(fisher_Molino_muscle$p.value, fisher_Pachon_muscle$p.value, fisher_Surface_muscle$p.value)
)
sig_labels_muscle$label <- cut(sig_labels_muscle$p_value,
                               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                               labels = c("***", "**", "*", ""),
                               right = FALSE)
total_counts_muscle <- overlap_counts_muscle %>%
  group_by(Comparison) %>%
  summarise(total = sum(Count))
sig_labels_muscle <- left_join(sig_labels_muscle, total_counts_muscle, by = "Comparison")
sig_labels_muscle$ypos <- sig_labels_muscle$total + 20

# ---- Cross-Morph Comparisons (Fed only) ----
fed_muscle_samples <- rownames(sampleTable)[sampleTable$condition == "Fed" & sampleTable$tissue == "muscle"]
fed_muscle_counts <- counts[, fed_muscle_samples]
fed_muscle_sampleTable <- sampleTable[fed_muscle_samples, ]

fed_muscle_dds <- DESeqDataSetFromMatrix(countData = fed_muscle_counts, colData = fed_muscle_sampleTable, design = ~ morph)
fed_muscle_dds <- DESeq(fed_muscle_dds)

res_MvP_muscle <- results(fed_muscle_dds, contrast = c("morph", "mol", "pa"))
res_MvS_muscle <- results(fed_muscle_dds, contrast = c("morph", "mol", "sf"))
res_PvS_muscle <- results(fed_muscle_dds, contrast = c("morph", "pa", "sf"))

sig_MvP_muscle <- get_sig_genes(res_MvP_muscle)
sig_MvS_muscle <- get_sig_genes(res_MvS_muscle)
sig_PvS_muscle <- get_sig_genes(res_PvS_muscle)

deg_MvP_muscle <- rownames(sig_MvP_muscle)
deg_MvS_muscle <- rownames(sig_MvS_muscle)
deg_PvS_muscle <- rownames(sig_PvS_muscle)

dmg_MvP_muscle <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
dmg_MvS_muscle <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1
dmg_PvS_muscle <- read.table("/Path/to/intersects/Intersects.txt", header = FALSE)$V1

fisher_MvP_muscle <- run_fisher(deg_MvP_muscle, dmg_MvP_muscle, bg_genes)
fisher_MvS_muscle <- run_fisher(deg_MvS_muscle, dmg_MvS_muscle, bg_genes)
fisher_PvS_muscle <- run_fisher(deg_PvS_muscle, dmg_PvS_muscle, bg_genes)

overlap_counts_muscle_morph <- data.frame(
  Comparison = rep(c("Molino vs Pachón", "Molino vs Surface", "Pachón vs Surface"), each = 3),
  Category = rep(c("DEG only", "DMG only", "Both"), 3),
  Count = c(
    length(setdiff(deg_MvP_muscle, dmg_MvP_muscle)),
    length(setdiff(dmg_MvP_muscle, deg_MvP_muscle)),
    length(intersect(deg_MvP_muscle, dmg_MvP_muscle)),
    
    length(setdiff(deg_MvS_muscle, dmg_MvS_muscle)),
    length(setdiff(dmg_MvS_muscle, deg_MvS_muscle)),
    length(intersect(deg_MvS_muscle, dmg_MvS_muscle)),
    
    length(setdiff(deg_PvS_muscle, dmg_PvS_muscle)),
    length(setdiff(dmg_PvS_muscle, deg_PvS_muscle)),
    length(intersect(deg_PvS_muscle, dmg_PvS_muscle))
  )
)

sig_labels_muscle_morph <- data.frame(
  Comparison = c("Molino vs Pachón", "Molino vs Surface", "Pachón vs Surface"),
  p_value = c(fisher_MvP_muscle$p.value, fisher_MvS_muscle$p.value, fisher_PvS_muscle$p.value)
)
sig_labels_muscle_morph$label <- cut(sig_labels_muscle_morph$p_value,
                                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                     labels = c("***", "**", "*", ""),
                                     right = FALSE)
total_counts_muscle_morph <- overlap_counts_muscle_morph %>%
  group_by(Comparison) %>%
  summarise(total = sum(Count))
sig_labels_muscle_morph <- left_join(sig_labels_muscle_morph, total_counts_muscle_morph, by = "Comparison")
sig_labels_muscle_morph$ypos <- sig_labels_muscle_morph$total + 20

# Combine all
overlap_counts_muscle_all <- rbind(overlap_counts_muscle, overlap_counts_muscle_morph)
sig_labels_muscle_all <- rbind(sig_labels_muscle, sig_labels_muscle_morph)

# Rename and reorder
overlap_counts_muscle_all$Comparison <- recode(overlap_counts_muscle_all$Comparison,
                                               "Surface (FvS)" = "Rio Choy (FvS)",
                                               "Molino vs Surface" = "Rio Choy vs Molino",
                                               "Pachón vs Surface" = "Rio Choy vs Pachón"
)
sig_labels_muscle_all$Comparison <- recode(sig_labels_muscle_all$Comparison,
                                           "Molino vs Surface" = "Rio Choy vs Molino",
                                           "Pachón vs Surface" = "Rio Choy vs Pachón"
)

overlap_counts_muscle_all$Comparison <- factor(overlap_counts_muscle_all$Comparison,
                                               levels = c("Rio Choy (FvS)", "Pachón (FvS)", "Molino (FvS)",
                                                          "Rio Choy vs Pachón", "Rio Choy vs Molino", "Molino vs Pachón"))
sig_labels_muscle_all$Comparison <- factor(sig_labels_muscle_all$Comparison,
                                           levels = levels(overlap_counts_muscle_all$Comparison))

overlap_counts_muscle_all$Category <- factor(overlap_counts_muscle_all$Category,
                                             levels = c("DEG only", "Both", "DMG only"))

# PLOT
ggplot(overlap_counts_muscle_all, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = segment_colors) +
  geom_text(data = sig_labels_muscle_all,
            aes(x = Comparison, y = ypos, label = label),
            inherit.aes = FALSE,
            size = 6, fontface = "bold") +
  labs(
    y = "Gene Count",
    title = "Shared and Unique DEGs & DMGs in Muscle (per Comparison)",
    fill = "Gene Set"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))



# Intersect genes: DEG ∩ DMG
muscle_intersects <- list(
  Molino_FvS = intersect(deg_Molino_muscle, dmg_Molino_muscle),
  Pachon_FvS = intersect(deg_Pachon_muscle, dmg_Pachon_muscle),
  Surface_FvS = intersect(deg_Surface_muscle, dmg_Surface_muscle),
  Molino_vs_Pachon = intersect(deg_MvP_muscle, dmg_MvP_muscle),
  Molino_vs_Surface = intersect(deg_MvS_muscle, dmg_MvS_muscle),
  Pachon_vs_Surface = intersect(deg_PvS_muscle, dmg_PvS_muscle)
)

# Write to CSV
write.csv(stack(muscle_intersects), 
          "/Users/Tissue.csv",
          row.names = FALSE)


# Intersect genes: DEG ∩ DMG
brain_intersects <- list(
  Molino_FvS = intersect(deg_Molino_brain, dmg_Molino_brain),
  Pachon_FvS = intersect(deg_Pachon_brain, dmg_Pachon_brain),
  Surface_FvS = intersect(deg_Surface_brain, dmg_Surface_brain),
  Molino_vs_Pachon = intersect(deg_MvP_brain, dmg_MvP_brain),
  Molino_vs_Surface = intersect(deg_MvS_brain, dmg_MvS_brain),
  Pachon_vs_Surface = intersect(deg_PvS_brain, dmg_PvS_brain)
)

# Write to CSV
write.csv(stack(brain_intersects), 
          "/Users/Tissue.csv",
          row.names = FALSE)










































#NEW _________________


# ====================== DEG–DMG (Direction-Aware) — REWRITE ======================

# ---- LOAD LIBRARIES ----
library(DESeq2)
library(tidyverse)
library(ggplot2)

# ---- SETTINGS ----
morph_labels <- c("mol" = "Molino", "pa" = "Pachón", "sf" = "Rio Choy")
output_dir <- "/Users/kayaharper/Desktop/GSEA_Results"
dir.create(output_dir, showWarnings = FALSE)

# ---- LOAD & CLEAN COUNT DATA ----
counts <- read.table("/Users/kayaharper/Desktop/Projects/Methyl_Seq/RNA-Seq/gene_counts.txt",
                     header = TRUE, sep = "\t", comment.char = "#", row.names = 1, check.names = FALSE)
counts <- counts[, -(1:5)]
colnames(counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", basename(colnames(counts)))

parse_sample <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  morph <- tolower(parts[1])
  fish_id <- parts[2]
  tissue <- tolower(parts[3])
  condition <- ifelse(grepl("F$", fish_id), "Fed", "Starved")
  return(c(morph = morph, fish_id = fish_id, tissue = tissue, condition = condition))
}
sample_info <- t(sapply(colnames(counts), parse_sample))
sampleTable <- as.data.frame(sample_info)
rownames(sampleTable) <- colnames(counts)

# ---- Build global DESeq object (used, then subset) ----
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleTable,
                              design = ~ morph + tissue + condition)
dds <- DESeq(dds)

# ====================== HELPERS ======================

get_deg_direction <- function(res, alpha = 0.05, lfc = 1) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::drop_na(padj, log2FoldChange) %>%
    dplyr::filter(padj < alpha, abs(log2FoldChange) >= lfc) %>%
    dplyr::mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
    dplyr::select(gene, Direction)
}

run_fisher <- function(setA, setB, background) {
  a_and_b <- length(intersect(setA, setB))
  a_only  <- length(setdiff(setA, setB))
  b_only  <- length(setdiff(setB, setA))
  neither <- length(setdiff(background, union(setA, setB)))
  fisher.test(matrix(c(a_and_b, a_only, b_only, neither), nrow = 2))
}

# Summarize directional overlaps for one comparison
directional_overlap_summary <- function(deg_df, dmg_hyper, dmg_hypo, background,
                                        label_comparison, label_tissue) {
  if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
  deg_up   <- deg_df$gene[deg_df$Direction == "Up"]
  deg_down <- deg_df$gene[deg_df$Direction == "Down"]

  tibble::tibble(
    Tissue = label_tissue,
    Comparison = label_comparison,
    Category = factor(c("Hyper ∩ Down","Hypo ∩ Up","Hyper ∩ Up","Hypo ∩ Down"),
                      levels = c("Hyper ∩ Down","Hypo ∩ Up","Hyper ∩ Up","Hypo ∩ Down")),
    Count = c(
      length(intersect(dmg_hyper, deg_down)),  # repression-consistent
      length(intersect(dmg_hypo,  deg_up)),    # de-repression-consistent
      length(intersect(dmg_hyper, deg_up)),    # opposite
      length(intersect(dmg_hypo,  deg_down))   # opposite
    ),
    p_value = c(
      run_fisher(dmg_hyper, deg_down, background)$p.value,
      run_fisher(dmg_hypo,  deg_up,   background)$p.value,
      NA_real_, NA_real_
    )
  ) %>%
    dplyr::mutate(sig = cut(p_value,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            labels = c("***","**","*",""),
                            right = FALSE))
}

# read DMG U/D files if they exist; return NULL if missing so we can skip gracefully
read_dmg_ud <- function(path_U, path_D) {
  if (!file.exists(path_U) || !file.exists(path_D)) {
    message("Skipping (missing U or D files):\n  U: ", path_U, "\n  D: ", path_D)
    return(NULL)
  }
  list(hyper = scan(path_U, what = "character", quiet = TRUE),
       hypo  = scan(path_D, what = "character", quiet = TRUE))
}

# plot helper
plot_directional_overlap <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    warning("No data to plot.")
    return(invisible(NULL))
  }
  fill_cols <- c("Hyper ∩ Down"="#6cbfc2","Hypo ∩ Up"="#6cbfc2",
                 "Hyper ∩ Up"="#E57E81","Hypo ∩ Down"="#E57E81")

  df2 <- df %>%
    group_by(Tissue, Comparison) %>%
    mutate(ypos = Count + max(1, round(max(Count, na.rm = TRUE)*0.05))) %>%
    ungroup()

  ggplot(df2, aes(x = Comparison, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fill_cols) +
    geom_text(
      data = subset(df2, !is.na(sig) & Category %in% c("Hyper ∩ Down","Hypo ∩ Up")),
      aes(x = Comparison, y = ypos, label = sig),
      inherit.aes = FALSE, size = 5, fontface = "bold"
    ) +
    facet_wrap(~ Tissue, nrow = 1, scales = "free_x") +
    labs(y = "Gene Count", title = "Directional DEG–DMG Overlap", fill = "Overlap") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          strip.text = element_text(size = 14, face = "bold"))
}

# background universe
bg_genes <- rownames(counts)

# ====================== COMPARISON BUILDERS ======================
# 1) Brain: within-morph Fed vs Starved (you provided U/D files)

brain_dds <- {
  keep <- sampleTable$tissue == "brain"
  DESeqDataSetFromMatrix(counts[, keep], sampleTable[keep, ], design = ~ morph + condition) %>% DESeq()
}

get_deg_fvs_brain <- function(morph_code) {
  dds_sub <- brain_dds[, brain_dds$morph == morph_code]
  dds_sub$condition <- droplevels(dds_sub$condition)
  design(dds_sub) <- ~ condition
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("condition", "Fed", "Starved"))
  get_deg_direction(res, alpha = 0.05, lfc = 1)
}

# DMG paths you listed for Brain (Metabolic) — U (hyper), D (hypo)
brain_paths <- list(
  "Molino (FvS)" = list(
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_U_FMB_SMB.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_D_FMB_SMB.txt",
    morph_code = "mol"
  ),
  "Pachón (FvS)" = list(
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_U_FPB_SPB.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_D_FPB_SPB.txt",
    morph_code = "pa"
  ),
  "Rio Choy (FvS)" = list(
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_U_FSB_SSB.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_D_FSB_SSB.txt",
    morph_code = "sf"
  )
)

brain_rows <- list()
for (label in names(brain_paths)) {
  info <- brain_paths[[label]]
  dmg_sets <- read_dmg_ud(info$U, info$D)
  if (is.null(dmg_sets)) next
  deg_df <- get_deg_fvs_brain(info$morph_code)
  row <- directional_overlap_summary(deg_df, dmg_sets$hyper, dmg_sets$hypo,
                                     background = bg_genes,
                                     label_comparison = label, label_tissue = "Brain")
  if (!is.null(row)) brain_rows[[label]] <- row
}
brain_df <- dplyr::bind_rows(brain_rows)

# 2) Liver: cross-morph (Fed only) — you provided U/D files under Liver/Morph/Fed
# build a DESeq for Fed-only liver, design ~ morph
liver_dds_fed <- {
  keep <- sampleTable$tissue == "liver" & sampleTable$condition == "Fed"
  DESeqDataSetFromMatrix(counts[, keep], sampleTable[keep, ], design = ~ morph) %>% DESeq()
}

get_deg_morph_liver <- function(level_a, level_b) {
  res <- results(liver_dds_fed, contrast = c("morph", level_a, level_b))
  get_deg_direction(res, alpha = 0.05, lfc = 1)
}
# map liver DMG U/D files to correct comparison labels
# map liver DMG U/D files to correct comparison labels
liver_paths <- list(
  "Molino vs Pachón" = list(
    # FPL vs FML  → Pachón vs Molino
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_U_FPL_FML.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_D_FPL_FML.txt",
    # DESeq contrast: Molino vs Pachón
    level_a = "mol",
    level_b = "pa"
  ),
  
  "Rio Choy vs Molino" = list(
    # FSL vs FML → Surface vs Molino
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_U_FML_FSL.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_D_FML_FSL.txt",
    level_a = "sf",
    level_b = "mol"
  ),
  
  "Rio Choy vs Pachón" = list(
    # FSL vs FPL → Surface vs Pachón
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_U_FPL_FSL.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_D_FPL_FSL.txt",
    level_a = "sf",
    level_b = "pa"
  )
)

liver_rows <- list()
for (label in names(liver_paths)) {
  info <- liver_paths[[label]]
  dmg_sets <- read_dmg_ud(info$U, info$D)
  if (is.null(dmg_sets)) next
  deg_df <- get_deg_morph_liver(info$level_a, info$level_b)
  row <- directional_overlap_summary(deg_df, dmg_sets$hyper, dmg_sets$hypo,
                                     background = bg_genes,
                                     label_comparison = label, label_tissue = "Liver")
  if (!is.null(row)) liver_rows[[label]] <- row
}
liver_df <- dplyr::bind_rows(liver_rows)

# 3) Muscle: cross-morph (Fed only) — you provided U/D files under Muscle/Morph/Fed
muscle_dds_fed <- {
  keep <- sampleTable$tissue == "muscle" & sampleTable$condition == "Fed"
  DESeqDataSetFromMatrix(counts[, keep], sampleTable[keep, ], design = ~ morph) %>% DESeq()
}

get_deg_morph_muscle <- function(level_a, level_b) {
  res <- results(muscle_dds_fed, contrast = c("morph", level_a, level_b))
  get_deg_direction(res, alpha = 0.05, lfc = 1)
}

muscle_paths <- list(
  "Molino vs Pachón" = list(
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_U_FMM_FPM.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_D_FMM_FPM.txt",
    level_a = "mol", level_b = "pa"
  ),
  "Rio Choy vs Molino" = list(
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_U_FSM_FMM.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_D_FSM_FMM.txt",
    level_a = "sf", level_b = "mol"
  ),
  "Rio Choy vs Pachón" = list(
    U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_U_FSM_FPM.txt",
    D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_D_FSM_FPM.txt",
    level_a = "sf", level_b = "pa"
  )
)

muscle_rows <- list()
for (label in names(muscle_paths)) {
  info <- muscle_paths[[label]]
  dmg_sets <- read_dmg_ud(info$U, info$D)
  if (is.null(dmg_sets)) next
  deg_df <- get_deg_morph_muscle(info$level_a, info$level_b)
  row <- directional_overlap_summary(deg_df, dmg_sets$hyper, dmg_sets$hypo,
                                     background = bg_genes,
                                     label_comparison = label, label_tissue = "Muscle")
  if (!is.null(row)) muscle_rows[[label]] <- row
}
muscle_df <- dplyr::bind_rows(muscle_rows)

# ====================== COMBINE & PLOT ======================
dir_all <- dplyr::bind_rows(brain_df, liver_df, muscle_df)

# Consistent ordering in facets and x-axis
dir_all$Tissue <- factor(dir_all$Tissue, levels = c("Liver","Brain","Muscle"))

# For prettier x-order within each tissue
cmp_levels <- c("Rio Choy (FvS)","Pachón (FvS)","Molino (FvS)",
                "Rio Choy vs Pachón","Rio Choy vs Molino","Molino vs Pachón")
dir_all$Comparison <- factor(dir_all$Comparison, levels = cmp_levels)

p <- plot_directional_overlap(dir_all)
print(p)

ggsave(filename = file.path(output_dir, "Directional_DEG_DMG_Overlap_All_Tissues.pdf"),
       plot = p, width = 12, height = 5)

































# ====================== ALL 18 COMPARISONS: TOTAL DEG / TOTAL DMG / TOTAL INTERSECT ======================

library(DESeq2)
library(tidyverse)

# ---- helpers ----
get_deg_total <- function(res, alpha = 0.05, lfc = 1) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::drop_na(padj, log2FoldChange) %>%
    dplyr::filter(padj < alpha, abs(log2FoldChange) >= lfc) %>%
    dplyr::pull(gene) %>%
    unique()
}

read_gene_list <- function(path) {
  # scan is fast + avoids factor weirdness
  unique(trimws(scan(path, what = "character", quiet = TRUE)))
}

read_dmg_ud <- function(path_U, path_D) {
  if (!file.exists(path_U) || !file.exists(path_D)) {
    message("Missing DMG files:\n  U: ", path_U, "\n  D: ", path_D)
    return(NULL)
  }
  U <- read_gene_list(path_U)
  D <- read_gene_list(path_D)
  list(U = U, D = D, all = unique(c(U, D)))
}

count_row <- function(tissue, comparison, deg_genes, dmg_sets) {
  if (is.null(dmg_sets)) return(NULL)
  deg_genes <- unique(trimws(deg_genes))
  dmg_all   <- unique(trimws(dmg_sets$all))
  
  tibble(
    Tissue = tissue,
    Comparison = comparison,
    n_DEG_total = length(deg_genes),
    n_DMG_total = length(dmg_all),
    n_intersect_total = length(intersect(deg_genes, dmg_all))
  )
}

# ---- DMG PATH MAPS (your real files) ----

# FvS within morph: 3 per tissue
dmg_paths_fvs <- list(
  Liver = list(
    "Molino (FvS)"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Metabolic/Intersects_U_FML_SML.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Metabolic/Intersects_D_FML_SML.txt",
                            morph_code = "mol"),
    "Pachón (FvS)"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Metabolic/Intersects_U_FPL_SPL.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Metabolic/Intersects_D_FPL_SPL.txt",
                            morph_code = "pa"),
    "Rio Choy (FvS)" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Metabolic/Intersects_U_FSL_SSL.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Metabolic/Intersects_D_FSL_SSL.txt",
                            morph_code = "sf")
  ),
  Brain = list(
    "Molino (FvS)"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_U_FMB_SMB.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_D_FMB_SMB.txt",
                            morph_code = "mol"),
    "Pachón (FvS)"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_U_FPB_SPB.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_D_FPB_SPB.txt",
                            morph_code = "pa"),
    "Rio Choy (FvS)" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_U_FSB_SSB.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Metabolic/Intersects_D_FSB_SSB.txt",
                            morph_code = "sf")
  ),
  Muscle = list(
    "Molino (FvS)"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Metabolic/Intersects_U_FMM_SMM.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Metabolic/Intersects_D_FMM_SMM.txt",
                            morph_code = "mol"),
    "Pachón (FvS)"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Metabolic/Intersects_U_SPM_FPM.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Metabolic/Intersects_D_SPM_FPM.txt",
                            morph_code = "pa"),
    "Rio Choy (FvS)" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Metabolic/Intersects_U_FSM_SSM.txt",
                            D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Metabolic/Intersects_D_FSM_SSM.txt",
                            morph_code = "sf")
  )
)

# Fed cross-morph: 3 per tissue
dmg_paths_morph <- list(
  Liver = list(
    "Molino vs Pachón"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_U_FPL_FML.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_D_FPL_FML.txt",
                                level_a = "mol", level_b = "pa"),
    "Rio Choy vs Molino" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_U_FML_FSL.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_D_FML_FSL.txt",
                                level_a = "sf", level_b = "mol"),
    "Rio Choy vs Pachón" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_U_FPL_FSL.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Liver/Morph/Fed/Intersects_D_FPL_FSL.txt",
                                level_a = "sf", level_b = "pa")
  ),
  Brain = list(
    "Molino vs Pachón"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Morph/Fed/Intersects_U_FPB_FMB.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Morph/Fed/Intersects_D_FPB_FMB.txt",
                                level_a = "mol", level_b = "pa"),
    "Rio Choy vs Molino" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Morph/Fed/Intersects_U_FSB_FMB.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Morph/Fed/Intersects_D_FSB_FMB.txt",
                                level_a = "sf", level_b = "mol"),
    "Rio Choy vs Pachón" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Morph/Fed/Intersects_U_FSB_FPB.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Brain/Morph/Fed/Intersects_D_FSB_FPB.txt",
                                level_a = "sf", level_b = "pa")
  ),
  Muscle = list(
    "Molino vs Pachón"   = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_U_FMM_FPM.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_D_FMM_FPM.txt",
                                level_a = "mol", level_b = "pa"),
    "Rio Choy vs Molino" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_U_FSM_FMM.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_D_FSM_FMM.txt",
                                level_a = "sf", level_b = "mol"),
    "Rio Choy vs Pachón" = list(U = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_U_FSM_FPM.txt",
                                D = "/Users/kayaharper/Desktop/Projects/Methyl_Seq/Gene_Intersects/Muscle/Morph/Fed/Intersects_D_FSM_FPM.txt",
                                level_a = "sf", level_b = "pa")
  )
)

# ---- Build DESeq objects per tissue ----
make_dds_tissue <- function(tissue) {
  keep <- sampleTable$tissue == tissue
  DESeqDataSetFromMatrix(counts[, keep], sampleTable[keep, ], design = ~ morph + condition) %>% DESeq()
}

brain_dds  <- make_dds_tissue("brain")
liver_dds  <- make_dds_tissue("liver")
muscle_dds <- make_dds_tissue("muscle")

# Fed-only DESeq for cross-morph per tissue
make_dds_fed_tissue <- function(tissue) {
  keep <- sampleTable$tissue == tissue & sampleTable$condition == "Fed"
  DESeqDataSetFromMatrix(counts[, keep], sampleTable[keep, ], design = ~ morph) %>% DESeq()
}

brain_dds_fed  <- make_dds_fed_tissue("brain")
liver_dds_fed  <- make_dds_fed_tissue("liver")
muscle_dds_fed <- make_dds_fed_tissue("muscle")

dds_by_tissue <- list(brain = brain_dds, liver = liver_dds, muscle = muscle_dds)
dds_fed_by_tissue <- list(brain = brain_dds_fed, liver = liver_dds_fed, muscle = muscle_dds_fed)

# ---- Functions to get DEGs for each comparison type ----
get_deg_fvs <- function(dds_tissue, morph_code, alpha=0.05, lfc=1) {
  dds_sub <- dds_tissue[, dds_tissue$morph == morph_code]
  dds_sub$condition <- droplevels(dds_sub$condition)
  design(dds_sub) <- ~ condition
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("condition", "Fed", "Starved"))
  get_deg_total(res, alpha=alpha, lfc=lfc)
}

get_deg_morph <- function(dds_fed_tissue, level_a, level_b, alpha=0.05, lfc=1) {
  res <- results(dds_fed_tissue, contrast = c("morph", level_a, level_b))
  get_deg_total(res, alpha=alpha, lfc=lfc)
}

# ---- Build the 18-row table ----
rows <- list()

# 9 FvS rows
for (tissue_name in names(dmg_paths_fvs)) {
  tissue_code <- tolower(tissue_name)  # "Liver" -> "liver"
  dds_tissue <- dds_by_tissue[[tissue_code]]
  
  for (label in names(dmg_paths_fvs[[tissue_name]])) {
    info <- dmg_paths_fvs[[tissue_name]][[label]]
    dmg_sets <- read_dmg_ud(info$U, info$D)
    deg_genes <- get_deg_fvs(dds_tissue, info$morph_code)
    rows[[paste(tissue_name, label, sep=" | ")]] <- count_row(tissue_name, label, deg_genes, dmg_sets)
  }
}

# 9 cross-morph Fed rows
for (tissue_name in names(dmg_paths_morph)) {
  tissue_code <- tolower(tissue_name)
  dds_fed_tissue <- dds_fed_by_tissue[[tissue_code]]
  
  for (label in names(dmg_paths_morph[[tissue_name]])) {
    info <- dmg_paths_morph[[tissue_name]][[label]]
    dmg_sets <- read_dmg_ud(info$U, info$D)
    deg_genes <- get_deg_morph(dds_fed_tissue, info$level_a, info$level_b)
    rows[[paste(tissue_name, label, sep=" | ")]] <- count_row(tissue_name, label, deg_genes, dmg_sets)
  }
}

deg_dmg_totals_18 <- bind_rows(rows)

# nice ordering (optional)
cmp_levels <- c("Rio Choy (FvS)","Pachón (FvS)","Molino (FvS)",
                "Rio Choy vs Pachón","Rio Choy vs Molino","Molino vs Pachón")

deg_dmg_totals_18 <- deg_dmg_totals_18 %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Liver","Brain","Muscle")),
    Comparison = recode(Comparison, "Surface (FvS)" = "Rio Choy (FvS)"),
    Comparison = factor(Comparison, levels = cmp_levels)
  ) %>%
  arrange(Tissue, Comparison)

print(deg_dmg_totals_18, n = Inf)

write.csv(deg_dmg_totals_18,
          file = file.path(output_dir, "DEG_DMG_Totals_All18_ByTissue_Comparison.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")





# ====================== GENE LISTS: ALL DEG∩DMG INTERSECTS (ALL 18) ======================

# Returns a tidy table of intersect genes for one comparison
make_intersect_gene_table <- function(tissue, comparison, deg_genes, dmg_sets) {
  if (is.null(dmg_sets)) return(NULL)
  dmg_all <- unique(c(dmg_sets$U, dmg_sets$D))
  inter <- sort(unique(intersect(deg_genes, dmg_all)))
  if (length(inter) == 0) return(NULL)
  
  tibble(
    Tissue = tissue,
    Comparison = comparison,
    GeneID = inter
  )
}

intersect_rows <- list()

# ---- 9 FvS comparisons (within morph) ----
for (tissue_name in names(dmg_paths_fvs)) {
  tissue_code <- tolower(tissue_name)
  dds_tissue <- dds_by_tissue[[tissue_code]]
  
  for (label in names(dmg_paths_fvs[[tissue_name]])) {
    info <- dmg_paths_fvs[[tissue_name]][[label]]
    dmg_sets <- read_dmg_ud(info$U, info$D)
    if (is.null(dmg_sets)) next
    
    deg_genes <- get_deg_fvs(dds_tissue, info$morph_code)
    
    out <- make_intersect_gene_table(tissue_name, label, deg_genes, dmg_sets)
    if (!is.null(out)) intersect_rows[[paste(tissue_name, label, sep=" | ")]] <- out
  }
}

# ---- 9 Fed cross-morph comparisons ----
for (tissue_name in names(dmg_paths_morph)) {
  tissue_code <- tolower(tissue_name)
  dds_fed_tissue <- dds_fed_by_tissue[[tissue_code]]
  
  for (label in names(dmg_paths_morph[[tissue_name]])) {
    info <- dmg_paths_morph[[tissue_name]][[label]]
    dmg_sets <- read_dmg_ud(info$U, info$D)
    if (is.null(dmg_sets)) next
    
    deg_genes <- get_deg_morph(dds_fed_tissue, info$level_a, info$level_b)
    
    out <- make_intersect_gene_table(tissue_name, label, deg_genes, dmg_sets)
    if (!is.null(out)) intersect_rows[[paste(tissue_name, label, sep=" | ")]] <- out
  }
}

deg_dmg_intersect_genes_all18 <- bind_rows(intersect_rows)

# Save ONE combined file (recommended)
write.csv(
  deg_dmg_intersect_genes_all18,
  file = file.path(output_dir, "DEG_DMG_Intersect_Genes_All18_ByTissue_Comparison.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Optional: also save a wide-ish summary of counts per comparison (handy quick check)
deg_dmg_intersect_counts_check <- deg_dmg_intersect_genes_all18 %>%
  count(Tissue, Comparison, name = "n_intersect_genes") %>%
  arrange(factor(Tissue, levels = c("Liver","Brain","Muscle")), Comparison)

write.csv(
  deg_dmg_intersect_counts_check,
  file = file.path(output_dir, "DEG_DMG_Intersect_Genes_Counts_Check_All18.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

