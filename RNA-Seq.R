
# @kayaharper@tamu.edu
# RNA-Seq Analysis in Astyanax mexicanus 
# ---- LIBRARIES ----

  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(apeglm)
  library(clusterProfiler)
  library(org.Dr.eg.db)
  library(tidyverse)
  library(RColorBrewer)
  library(matrixStats)

# Make sure dplyr verbs win conflicts (optional but helpful)
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange

# ---- PATHS ----
counts_path <- "/path/to/file"
output_dir  <- "/path/to/file"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- SETTINGS ----
alpha <- 0.05
lfc_thresh <- 1

morph_labels <- c("mol" = "Molino", "pa" = "Pachón", "sf" = "Rio Choy")
morph_colors <- c("Molino" = "#ce9e91", "Pachón" = "#E57E81", "Rio Choy" = "#BEC5DF")

anno_colors <- list(
  Morph = morph_colors,
  Condition = c("Fed" = "#6cbfc2", "Starved" = "#C43A1A")
)

shape_map <- c("Rio Choy" = 16, "Pachón" = 17, "Molino" = 15)

group_order <- c(
  "Rio Choy_Fed", "Rio Choy_Starved",
  "Pachón_Fed",  "Pachón_Starved",
  "Molino_Fed",  "Molino_Starved"
)

heatmap_palette <- colorRampPalette(c("#65999e", "#FFFFFF", "#9467bd"))(100)

# ---- LOAD & CLEAN COUNT DATA ----
counts <- read.table(
  counts_path,
  header = TRUE, sep = "\t", comment.char = "#",
  row.names = 1, check.names = FALSE
)

# featureCounts often has annotation columns; you removed 1:5 in your script.
counts <- counts[, -(1:5)]
colnames(counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", basename(colnames(counts)))

# ---- PARSE SAMPLE INFO ----
parse_sample <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  morph <- tolower(parts[1])
  fish_id <- parts[2]
  tissue <- tolower(parts[3])
  condition <- ifelse(grepl("F$", fish_id), "Fed", "Starved")
  c(morph = morph, fish_id = fish_id, tissue = tissue, condition = condition)
}

sample_info <- t(sapply(colnames(counts), parse_sample))
sampleTable <- as.data.frame(sample_info, stringsAsFactors = FALSE)
rownames(sampleTable) <- colnames(counts)

# enforce factors + consistent levels
sampleTable$morph <- factor(sampleTable$morph, levels = names(morph_labels))
sampleTable$tissue <- factor(sampleTable$tissue, levels = c("liver", "brain", "muscle"))
sampleTable$condition <- factor(sampleTable$condition, levels = c("Fed", "Starved"))

# ---- BUILD MASTER DDS (for convenience; we subset from this) ----
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sampleTable,
  design = ~ morph + tissue + condition
)

# ---- HELPERS ----
run_DE_by_tissue_morph <- function(dds, tissue_name, morph_name) {
  dds_sub <- dds[, dds$tissue == tissue_name & dds$morph == morph_name]
  dds_sub$tissue <- droplevels(dds_sub$tissue)
  dds_sub$morph <- droplevels(dds_sub$morph)
  dds_sub$condition <- droplevels(dds_sub$condition)
  
  if (length(unique(dds_sub$condition)) < 2) {
    stop(paste("Not enough conditions in", morph_name, tissue_name))
  }
  
  design(dds_sub) <- ~ condition
  dds_sub <- DESeq(dds_sub, quiet = TRUE)
  results(dds_sub, contrast = c("condition", "Fed", "Starved"))
}

# ---- GSEA FUNCTION (plot + store top terms) ----
run_GSEA_plot <- function(de_result, tissue_label, comparison_name, morph_label) {
  res_df <- as.data.frame(de_result) %>%
    rownames_to_column("gene") %>%
    drop_na(pvalue, log2FoldChange) %>%
    mutate(rank_metric = -log10(pvalue) * sign(log2FoldChange))
  
  gene_list <- sort(setNames(res_df$rank_metric, res_df$gene), decreasing = TRUE)
  
  gsea_res <- gseGO(
    geneList = gene_list,
    OrgDb = org.Dr.eg.db,
    ont = "BP",
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE,
    eps = 0
  )
  
  # empty safe return
  if (nrow(gsea_res@result) == 0) {
    return(list(
      gsea_result = gsea_res,
      top_terms_df = tibble(),
      gsea_plot = ggplot() +
        theme_void() +
        ggtitle(paste("No enriched GO terms in", comparison_name)) +
        theme(plot.title = element_text(size = 14, face = "italic", hjust = 0.5))
    ))
  }
  
  gsea_res@result <- gsea_res@result %>%
    mutate(direction = ifelse(NES > 0, "Down in Starved", "Up in Starved"))
  
  top_terms <- gsea_res@result %>%
    group_by(direction) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  gsea_df <- gsea_res@result %>%
    filter(ID %in% top_terms$ID) %>%
    mutate(
      direction = ifelse(NES > 0, "Down in Starved", "Up in Starved"),
      GeneRatio = (str_count(core_enrichment, ",") + 1) / setSize,
      p.adjust = as.numeric(p.adjust),
      Description = str_to_title(Description)
    )
  
  # plot
  gsea_plot <- ggplot(gsea_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
    geom_point(aes(size = setSize, color = p.adjust)) +
    scale_color_gradientn(colors = c("#65999e", "#9467bd"), name = "Adjusted p-value", trans = "reverse") +
    scale_size_continuous(name = "Gene Count") +
    facet_wrap(~direction, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(title = paste("GSEA -", tissue_label, "-", comparison_name), x = "GeneRatio", y = NULL) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 12),
      strip.text = element_text(size = 13, face = "bold")
    )
  
  top_terms_df <- gsea_df %>%
    transmute(
      Population = morph_label,
      Tissue = tissue_label,
      Direction = direction,
      GO_Term = ID,
      Description,
      p.adjust,
      NES,
      GeneRatio,
      setSize,
      leadingEdge = core_enrichment
    )
  
  list(gsea_result = gsea_res, top_terms_df = top_terms_df, gsea_plot = gsea_plot)
}

# ---- HEATMAP (per tissue) ----
make_tissue_heatmap <- function(dds, tissue_name) {
  dds_t <- dds[, dds$tissue == tissue_name]
  dds_t$tissue <- droplevels(dds_t$tissue)
  dds_t$morph <- droplevels(dds_t$morph)
  dds_t$condition <- droplevels(dds_t$condition)
  
  # model condition while controlling for morph
  design(dds_t) <- ~ morph + condition
  dds_t <- DESeq(dds_t, quiet = TRUE)
  
  res <- results(dds_t, contrast = c("condition", "Starved", "Fed"))
  deg <- res[which(res$padj < alpha & !is.na(res$padj)), ]
  deg_genes <- rownames(deg)
  
  if (length(deg_genes) == 0) {
    message("No DEGs for heatmap in tissue: ", tissue_name)
    return(list(n_degs = 0))
  }
  
  vsd <- vst(dds_t, blind = FALSE)
  mat <- assay(vsd)[deg_genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  
  # annotation
  meta <- as.data.frame(colData(vsd)) %>%
    mutate(
      morph = recode(as.character(morph), sf = "Rio Choy", pa = "Pachón", mol = "Molino"),
      condition = str_to_title(as.character(condition)),
      group = paste(morph, condition, sep = "_")
    )
  
  # group means
  mat_grouped <- as.data.frame(t(mat_z)) %>%
    mutate(group = meta$group) %>%
    group_by(group) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    column_to_rownames("group") %>%
    t()
  
  # enforce consistent order (keep only those present)
  keep_groups <- intersect(group_order, colnames(mat_grouped))
  mat_grouped <- mat_grouped[, keep_groups, drop = FALSE]
  
  anno_group <- data.frame(
    Morph = gsub("_.*", "", keep_groups),
    Condition = gsub(".*_", "", keep_groups),
    row.names = keep_groups
  )
  
  colnames(mat_grouped) <- gsub("_", " ", colnames(mat_grouped))
  
  pheatmap(
    mat_grouped,
    annotation_col = anno_group,
    annotation_colors = anno_colors,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    gaps_col = c(2, 4),
    scale = "none",
    color = heatmap_palette,
    fontsize_col = 18,
    fontsize = 14,
    angle_col = 45,
    border_color = NA,
    breaks = seq(-2, 2, length.out = 101),
    main = paste0("DEGs in ", str_to_title(tissue_name), ": Starved vs Fed (controlling for morph)")
  )
  
  list(n_degs = length(deg_genes))
}

# ---- DEG UNION (for PCA per tissue) ----
get_deg_union_by_tissue <- function(dds, tissue_name) {
  unique(unlist(lapply(names(morph_labels), function(morph_name) {
    res <- run_DE_by_tissue_morph(dds, tissue_name = tissue_name, morph_name = morph_name)
    rownames(res)[which(res$padj < alpha & !is.na(res$padj))]
  })))
}

# ---- PCA + SCREE (per tissue, DEG-union) ----
run_pca_and_scree <- function(dds, tissue_name, plot_all_pcs = TRUE) {
  deg_union <- get_deg_union_by_tissue(dds, tissue_name)
  
  dds_t <- dds[, dds$tissue == tissue_name]
  dds_t$tissue <- droplevels(dds_t$tissue)
  dds_t$morph <- droplevels(dds_t$morph)
  dds_t$condition <- droplevels(dds_t$condition)
  
  design(dds_t) <- ~ morph + condition
  dds_t <- DESeq(dds_t, quiet = TRUE)
  
  vsd <- vst(dds_t, blind = FALSE)
  
  if (length(deg_union) < 2) {
    message("Too few DEG-union genes for PCA in tissue: ", tissue_name)
    return(invisible(NULL))
  }
  
  mat <- assay(vsd)[deg_union, , drop = FALSE]
  pca <- prcomp(t(mat), scale. = TRUE)
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  
  meta <- as.data.frame(colData(vsd)) %>%
    mutate(
      Morph = factor(morph, levels = names(morph_labels), labels = morph_labels),
      Condition = factor(str_to_title(as.character(condition)), levels = c("Fed", "Starved"))
    )
  
  # scree
  scree_df <- data.frame(PC = seq_along(percentVar), Variance = percentVar)
  p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
    geom_line(group = 1, linewidth = 1.2) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = 1:min(10, nrow(scree_df))) +
    labs(
      title = paste("Scree Plot -", str_to_title(tissue_name), "(DEG Union)"),
      x = "Principal Component",
      y = "% Variance Explained"
    ) +
    theme_minimal(base_size = 14)
  print(p_scree)
  
  plot_custom_pca <- function(pc_x, pc_y) {
    ggplot(data.frame(pca$x, meta), aes(
      x = .data[[paste0("PC", pc_x)]],
      y = .data[[paste0("PC", pc_y)]],
      color = Condition,
      shape = Morph
    )) +
      geom_point(size = 8.5, alpha = 0.75, stroke = 1.2) +
      scale_color_manual(values = anno_colors$Condition) +
      scale_shape_manual(values = shape_map) +
      labs(
        title = paste0("PCA - ", str_to_title(tissue_name), " (DEG Union): PC", pc_x, " vs PC", pc_y),
        x = paste0("PC", pc_x, ": ", percentVar[pc_x], "% variance"),
        y = paste0("PC", pc_y, ": ", percentVar[pc_y], "% variance")
      ) +
      theme_minimal(base_size = 16)
  }
  
  # always PC1 vs PC2
  print(plot_custom_pca(1, 2))
  
  if (plot_all_pcs) {
    pcs <- combn(1:4, 2)
    for (i in 1:ncol(pcs)) print(plot_custom_pca(pcs[1, i], pcs[2, i]))
  }
  
  invisible(list(pca = pca, percentVar = percentVar))
}

# ---- MAIN LOOPS: DEGs CSV + DEG COUNTS + GSEA ----
tissues <- levels(sampleTable$tissue)
morphs  <- levels(sampleTable$morph)

all_gsea_top_terms <- list()
deg_counts <- list()

for (tissue in tissues) {
  for (morph in morphs) {
    
    message("DE + GSEA: ", morph, " / ", tissue)
    
    # DE
    de_res <- run_DE_by_tissue_morph(dds, tissue_name = tissue, morph_name = morph)
    
    # significant DEGs
    sig_df <- as.data.frame(de_res) %>%
      rownames_to_column("Gene ID") %>%
      drop_na(padj, log2FoldChange) %>%
      filter(padj < alpha & abs(log2FoldChange) >= lfc_thresh) %>%
      mutate(
        Population = morph_labels[[morph]],
        Tissue = str_to_title(tissue),
        `Regulation in starved (up/down)` = ifelse(log2FoldChange > 0, "Up", "Down"),
        log2FoldChange = round(log2FoldChange, 3),
        `P-Value` = signif(pvalue, 4),
        `Adjusted P-Value` = signif(padj, 4)
      ) %>%
      dplyr::select(Population, Tissue, `Regulation in starved (up/down)`,
                    `Gene ID`, log2FoldChange, `P-Value`, `Adjusted P-Value`)
    
    
    out_csv <- file.path(output_dir, paste0(
      "DEGs_", morph_labels[[morph]], "_", str_to_title(tissue), "_Starved_vs_Fed.csv"
    ))
    write.csv(sig_df, out_csv, row.names = FALSE)
    
    # store DEG counts for bar chart
    deg_counts[[paste(morph, tissue, sep = "_")]] <- tibble(
      Morph = morph_labels[[morph]],
      Tissue = str_to_title(tissue),
      Num_DEGs = nrow(sig_df)
    )
    
    # GSEA
    comp_name <- paste0(morph_labels[[morph]], "_Starved_vs_Fed")
    gsea_out <- run_GSEA_plot(de_res, str_to_title(tissue), comp_name, morph_labels[[morph]])
    print(gsea_out$gsea_plot)
    
    all_gsea_top_terms[[paste(morph, tissue, sep = "_")]] <- gsea_out$top_terms_df
  }
}

# write combined GSEA top terms
all_top_terms_df <- bind_rows(all_gsea_top_terms)
write.csv(all_top_terms_df, file.path(output_dir, "GSEA_Top_Terms_All.csv"), row.names = FALSE)

# ---- BAR CHART: # DEGs per Morph × Tissue ----
deg_df <- bind_rows(deg_counts)

p_deg_bar <- ggplot(deg_df, aes(x = Tissue, y = Num_DEGs, fill = Morph)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = morph_colors) +
  labs(y = "# of DEGs (Fed vs Starved within morph)", x = NULL, fill = "Morph") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  )
print(p_deg_bar)

ggsave(
  filename = file.path(output_dir, "DEG_Counts_BarChart.pdf"),
  plot = p_deg_bar,
  width = 8, height = 6
)

# ---- HEATMAPS: one per tissue (blue-white-purple) + print DEG counts used ----
heatmap_deg_counts <- list()

for (tissue in tissues) {
  message("Heatmap: ", tissue)
  hm_out <- make_tissue_heatmap(dds, tissue)
  heatmap_deg_counts[[tissue]] <- hm_out$n_degs
}

cat("\nNumber of DEGs used in each tissue heatmap (condition effect controlling for morph):\n")
cat("Liver:  ", heatmap_deg_counts[["liver"]],  " genes\n")
cat("Brain:  ", heatmap_deg_counts[["brain"]],  " genes\n")
cat("Muscle: ", heatmap_deg_counts[["muscle"]], " genes\n")

# ---- PCA + SCREE (per tissue; DEG union across morphs) ----
liver_results  <- run_pca_and_scree(dds, tissue_name = "liver",  plot_all_pcs = TRUE)
brain_results  <- run_pca_and_scree(dds, tissue_name = "brain",  plot_all_pcs = TRUE)
muscle_results <- run_pca_and_scree(dds, tissue_name = "muscle", plot_all_pcs = TRUE)

# ---- OPTIONAL: PCA all tissues (top variable genes) ----
vsd_all <- vst(dds, blind = FALSE)
top_var_genes <- head(order(rowVars(assay(vsd_all)), decreasing = TRUE), 1000)
mat_all <- assay(vsd_all)[top_var_genes, ]
pca_all <- prcomp(t(mat_all), scale. = TRUE)
percentVar_all <- round(100 * (pca_all$sdev^2 / sum(pca_all$sdev^2)))

meta_all <- as.data.frame(colData(vsd_all)) %>%
  mutate(
    Morph = factor(morph, levels = names(morph_labels), labels = morph_labels),
    Condition = factor(str_to_title(as.character(condition)), levels = c("Fed", "Starved")),
    Tissue = factor(tissue)
  )

p_all <- ggplot(data.frame(pca_all$x, meta_all), aes(PC1, PC2, color = Condition, shape = Tissue)) +
  geom_point(size = 6, alpha = 0.8, stroke = 1.1) +
  scale_color_manual(values = anno_colors$Condition) +
  labs(
    title = "PCA: All Tissues (Top 1000 Variable Genes)",
    x = paste0("PC1: ", percentVar_all[1], "% variance"),
    y = paste0("PC2: ", percentVar_all[2], "% variance")
  ) +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_all)

ggsave(
  filename = file.path(output_dir, "PCA_AllTissues_Top1000.pdf"),
  plot = p_all,
  width = 8, height = 6
)
