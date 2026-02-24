# DMR–DEG_DMR_Bar_Graph
#this graph determiens within what genomic feature our DMRs
#fall in within our DMG-DEG intersects -----------------------------------------------
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)

# --- Load GFF and define gene features ---
gff_path <- "/path/to/genomic/featurs"
gff <- import(gff_path)
genes_gr  <- gff[gff$type == "gene"]
exons_gr  <- gff[gff$type == "exon"]
cpg_gr    <- gff[gff$type == "CpG_island"]

# Define promoters (2kb around TSS)
promoters_gr   <- GenomicRanges::promoters(genes_gr, upstream = 1000, downstream = 1000)
gene_bodies_gr <- genes_gr

# --- Function to classify DMRs by feature --------------------------------
classify_dmrs <- function(dmr_path, comparison_name,
                          promoters_gr, exons_gr, gene_bodies_gr) {
  dmr_df <- read.csv(dmr_path)
  colnames(dmr_df)[1:3] <- c("chr", "start", "end")
  
  dmr_gr <- GRanges(
    seqnames = dmr_df$chr,
    ranges   = IRanges(start = dmr_df$start, end = dmr_df$end)
  )
  
  dmr_df$Feature <- "Intergenic"  # default
  
  # Overlaps
  hits_promoter <- queryHits(findOverlaps(dmr_gr, promoters_gr))
  hits_exon     <- queryHits(findOverlaps(dmr_gr, exons_gr))
  hits_gene     <- queryHits(findOverlaps(dmr_gr, gene_bodies_gr))
  
  # Hierarchical assignment
  dmr_df$Feature[hits_gene]     <- "GeneBody"
  dmr_df$Feature[hits_exon]     <- "Exon"
  dmr_df$Feature[hits_promoter] <- "Promoter"
  
  dmr_df$Comparison <- comparison_name
  dmr_df
}

# --- Helper: do full pipeline for one tissue -----------------------------
make_tissue_summary <- function(tissue_label,
                                dmr_paths,        # named list: Comparison -> path
                                intersect_csv) {  # CSV from stack(...), with 'values' column
  
  # 1) classify DMRs for all comparisons for this tissue
  dmr_list <- mapply(
    FUN  = function(path, comp) {
      classify_dmrs(path, comp, promoters_gr, exons_gr, gene_bodies_gr)
    },
    path = unlist(dmr_paths),
    comp = names(dmr_paths),
    SIMPLIFY = FALSE
  )
  
  all_dmrs_annotated <- dplyr::bind_rows(dmr_list)
  
  # 2) read DEG–DMG intersect genes for this tissue
  shared_df <- read.csv(intersect_csv)
  # from stack(), gene IDs are in 'values'
  shared_genes <- unique(shared_df$values)
  
  # 3) match these genes to the GFF genes
  gene_ids <- mcols(genes_gr)$Name   # change to $ID / $gene_id if that's what matches your LOC IDs
  shared_gene_gr <- genes_gr[gene_ids %in% shared_genes]
  
  # 4) mark DMRs that overlap those shared DEG–DMG genes
  dmr_gr_all <- GRanges(
    seqnames = all_dmrs_annotated$chr,
    ranges   = IRanges(start = all_dmrs_annotated$start,
                       end   = all_dmrs_annotated$end)
  )
  
  shared_overlaps <- findOverlaps(dmr_gr_all, shared_gene_gr)
  
  all_dmrs_annotated$Shared_DEG_DMG <- FALSE
  all_dmrs_annotated$Shared_DEG_DMG[queryHits(shared_overlaps)] <- TRUE
  
  filtered_dmrs <- all_dmrs_annotated %>% filter(Shared_DEG_DMG)
  
  if (nrow(filtered_dmrs) == 0) {
    warning("No DMRs overlapping DEG–DMG genes for tissue: ", tissue_label)
    return(dplyr::tibble())
  }
  
  # 5) CpG island annotation
  dmr_gr <- GRanges(
    seqnames = filtered_dmrs$chr,
    ranges   = IRanges(filtered_dmrs$start, filtered_dmrs$end)
  )
  
  overlap_cpg <- findOverlaps(dmr_gr, cpg_gr)
  
  filtered_dmrs$CpG_status <- "no CpG island"
  filtered_dmrs$CpG_status[queryHits(overlap_cpg)] <- "CpG island"
  
  # 6) summarise for plotting
  summary_cpg <- filtered_dmrs %>%
    filter(Feature != "Intergenic") %>%
    group_by(Feature, CpG_status) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Category = paste(Feature, CpG_status, sep = "\n"),
      Category = factor(Category, levels = c(
        "Promoter\nCpG island", "Promoter\nno CpG island",
        "Exon\nCpG island", "Exon\nno CpG island",
        "GeneBody\nCpG island", "GeneBody\nno CpG island"
      ))
    )
  
  total_N <- sum(summary_cpg$Count)
  
  summary_cpg %>%
    mutate(
      Percent = Count / total_N * 100,
      label   = paste0("n=", Count, ", ", round(Percent, 1), "%"),
      Tissue  = tissue_label,
      Total_N = total_N
    )
}


# --- MUSCLE --------------------------------------------------------------
muscle_summary <- make_tissue_summary(
  tissue_label = "Muscle",
  dmr_paths = list(
    "Molino vs Surface" = "/path/to/file",
    "Pachón vs Molino"  = "/path/to/file",
    "Pachón vs Surface" = "/path/to/file"
  ),
  intersect_csv = "/path/to/file"
)

# --- LIVER (update paths) -----------------------------------------------
liver_summary <- make_tissue_summary(
  tissue_label = "Liver",
  dmr_paths = list(
    # TODO: update these to your actual liver DMR CSVs
    "Molino vs Surface" = "/path/to/file",
    "Pachón vs Molino"  = "/path/to/file",
    "Pachón vs Surface" = "/path/to/file"
  ),
  intersect_csv = "/path/to/file"
)

# --- BRAIN (update paths) -----------------------------------------------
brain_summary <- make_tissue_summary(
  tissue_label = "Brain",
  dmr_paths = list(
    # TODO: update these to your actual brain DMR CSVs
    "Molino vs Surface" = "/path/to/file",
    "Pachón vs Molino"  = "/path/to/file",
    "Pachón vs Surface" = "/path/to/file"
  ),
  intersect_csv = "/path/to/file"
)


# Combine (drop any empty ones if a tissue had 0 overlaps)
summary_all <- dplyr::bind_rows(muscle_summary, liver_summary, brain_summary) %>%
  filter(!is.na(Category))

# Make Tissue a nice factor order
summary_all$Tissue <- factor(summary_all$Tissue, levels = c("Liver", "Brain", "Muscle"))

# Data frame for N annotations per tissue
N_by_tissue <- summary_all %>%
  group_by(Tissue) %>%
  summarise(Total_N = unique(Total_N)[1]) %>%
  mutate(x = 1, y = 2)   # put label near the bottom of each bar


# Contingency table: Tissue × Feature Category
ct <- xtabs(Count ~ Tissue + Category, data = summary_all)
ct
chisq.test(ct)

ct_Liver_Muscle <- (ct[c("Liver", "Muscle"), ])
ct_Muscle_brain <- (ct[c("Brain", "Muscle"), ])
ct_liver_brain <- ct[c("Liver", "Brain"), ]
chisq.test(ct_liver_brain)
chisq.test(ct_liver_brain)
chisq.test(ct_Muscle_brain)

# --- Plot: all three tissues in one graph -------------------------------
ggplot(summary_all, aes(x = "", y = Percent, fill = Category)) +
  geom_col(width = 0.2) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4, color = "black"
  ) +
  geom_text(
    data = N_by_tissue,
    aes(x = x, y = y, label = paste0("N = ", Total_N)),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "Promoter\nCpG island"     = "#B9F2E1",  # pastel mint-aqua
    "Promoter\nno CpG island"  = "#6BC1B5",  # muted aqua-teal
    "Exon\nCpG island"         = "#E4D4FA",  # very light purple
    "Exon\nno CpG island"      = "#AF9ED6",  # lavender
    "GeneBody\nCpG island"     = "#C6D7FA",  # cool periwinkle
    "GeneBody\nno CpG island"  = "#9FA8C6"   # cool grey-blue-purple
  )) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  facet_wrap(~ Tissue, nrow = 1) +
  labs(
    x = NULL,
    y = "% of DMR–DEG Intersections",
    fill = NULL,
    title = "DMR–DEG Localization by CpG Island and Region (Liver, Brain, Muscle)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  )
