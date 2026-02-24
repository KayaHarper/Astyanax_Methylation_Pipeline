#Average percentage methylation at each CpG site in the genome 
# ===== Clean % methylation dot plot (faceted) =====
library(methylKit)
library(rstudioapi)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# 1) Pick folder
selected_folder <- selectDirectory(caption = "Select Directory Containing Bismark .cov Files")
if (is.null(selected_folder) || selected_folder == "") stop("No directory selected.")

# Grab both .cov and .cov.gz
file.list <- list.files(
  path = selected_folder,
  pattern = "bismark\\.cov(\\.gz)?$",
  full.names = TRUE
)

if (length(file.list) == 0) stop("No .bismark.cov files found.")
sample.ids <- basename(file.list)

cat("Files detected:\n"); print(file.list)

# 2) Parse condition robustly from the token after the dash (handles -1F_, -2S_, -10F_, etc.)
section_after_dash <- sub("^[^-]*-", "", sample.ids)      # e.g., "1F_Liver_S..."
id_token           <- sub("_.*", "", section_after_dash)  # e.g., "1F"
cond_char          <- substr(id_token, nchar(id_token), nchar(id_token))  # "F" or "S"

# treatment for methylKit (0 = Fed, 1 = Starved)
treatment <- ifelse(cond_char == "F", 0L, 1L)

# 3) Read coverage
methylRawList.obj <- lapply(seq_along(file.list), function(i) {
  methRead(
    location  = file.list[i],             # use 'location=' for clarity
    sample.id = sample.ids[i],
    assembly  = "GCF_023375975.1",        # use your consistent assembly tag
    treatment = treatment[i],
    context   = "CpG",
    pipeline  = "bismarkCoverage"
  )
})
cat("Successfully loaded all methylation files!\n")

# 4) Unite and compute % methylation per sample
methylRawList.obj <- new("methylRawList", methylRawList.obj)
methylBase.obj <- unite(methylRawList.obj)

percent_meth_matrix <- percMethylation(methylBase.obj)
sample_means <- colMeans(percent_meth_matrix, na.rm = TRUE)

# 5) Build plotting data (note the comma after Condition!)
methylation_data <- data.frame(
  Sample = sample.ids,
  Morph   = sub("-.*", "", sample.ids),
  Tissue  = sub("^[^-]+-[^_]+_([^_]+)_.*", "\\1", sample.ids),
  Condition = ifelse(cond_char == "F", "Fed", "Starved"),
  Methylation_Percent = sample_means,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Morph = recode(Morph, "SF" = "Rio Choy", "PA" = "Pachón", "Mol" = "Molino", .default = Morph),
    Morph = factor(Morph, levels = c("Rio Choy", "Pachón", "Molino")),
    Tissue = factor(Tissue, levels = c("Brain", "Liver", "Muscle")),
    Condition = factor(Condition, levels = c("Fed", "Starved"))
  )

# (Optional) quick sanity check
# with(methylation_data, table(Condition, Tissue, Morph))

# 6) Per-group summary (mean ± SE)
summary_methylation <- summarise(
  methylation_data,
  n = dplyr::n(),
  Mean = mean(Methylation_Percent, na.rm = TRUE),
  SD   = stats::sd(Methylation_Percent,  na.rm = TRUE),
  SE   = if_else(n > 1, SD / sqrt(n), 0),
  .by  = c(Morph, Tissue, Condition)
)

# ---- Stars when CI doesn't overlap Fed·Rio Choy of the same tissue ----
z <- 1.96  # 95% CI

# CI for every group
summary_ci <- summary_methylation %>%
  dplyr::mutate(
    ci_low  = Mean - z * SE,
    ci_high = Mean + z * SE
  )

# Reference per tissue = Fed · Rio Choy
ref_ci <- summary_ci %>%
  dplyr::filter(Morph == "Rio Choy", Condition == "Fed") %>%
  dplyr::select(Tissue, ref_low = ci_low, ref_high = ci_high, ref_Mean = Mean, ref_SE = SE)

# Join ref to all rows in same tissue and flag non-overlap
summary_with_ref <- summary_ci %>%
  dplyr::left_join(ref_ci, by = "Tissue") %>%
  dplyr::mutate(
    non_overlap =
      !is.na(ref_low) &                                  # we have a reference in this tissue
      !(Morph == "Rio Choy" & Condition == "Fed") &      # don't star the reference itself
      (ci_high < ref_low | ci_low > ref_high)            # no CI overlap with Fed·Rio Choy
  )

# Star positions (a bit above the higher CI between the two)
sig_df <- summary_with_ref %>%
  dplyr::filter(non_overlap) %>%
  dplyr::transmute(
    Morph, Tissue, Condition,
    y = pmax(Mean + SE, ref_high, na.rm = TRUE) + 0.5,   # clear the taller interval
    label = "*"
  )

# (optional) see who got starred
sig_df %>% dplyr::arrange(Tissue, Condition, Morph) %>%
  dplyr::select(Tissue, Condition, Morph) %>% print()


# 8) Plot
pos_dodge <- position_dodge(width = 0.5)



# Draw summaries first, replicates last so dots are visible
pos_jit <- position_jitter(width = 0.18, height = 0)


ggplot(methylation_data, aes(x = Morph, y = Methylation_Percent, color = Condition)) +
  # summaries behind
  geom_errorbar(
    data = summary_methylation,
    inherit.aes = FALSE,
    aes(x = Morph, ymin = Mean - SE, ymax = Mean + SE, color = Condition),
    width = 0.18, size = 0.9, show.legend = FALSE
  ) +
  geom_point(
    data = summary_methylation,
    inherit.aes = FALSE,
    aes(x = Morph, y = Mean, color = Condition),
    size = 4.2, stroke = 1.1, shape = 21, fill = "white", alpha = 1,
    show.legend = FALSE
  ) +
  # replicates on top (beeswarm-like that works with singletons)
  geom_quasirandom(
    data = methylation_data,
    aes(x = Morph, y = Methylation_Percent),
    groupOnX = TRUE, width = 0.22, varwidth = FALSE,
    size = 2.6, alpha = 0.75, shape = 16, show.legend = FALSE
  ) +
  # stars
  geom_text(
    data = sig_df,
    inherit.aes = FALSE,
    aes(x = Morph, y = y, label = label),
    size = 5
  ) +
  facet_grid(rows = vars(Tissue), cols = vars(Condition), switch = "y", drop = FALSE) +
  scale_color_manual(values = c("Fed" = "#6cbfc2", "Starved" = "#C43A1A")) +
  labs(x = "Population", y = "% Methylation") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = NA, color = NA),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  )














md <- getData(methylBase.obj)

numCs_cols <- grep("^numCs", names(md), value = TRUE)
numTs_cols <- grep("^numTs", names(md), value = TRUE)

global_weighted <- sapply(seq_along(numCs_cols), function(i) {
  C <- sum(md[[numCs_cols[i]]], na.rm = TRUE)
  T <- sum(md[[numTs_cols[i]]], na.rm = TRUE)
  100 * C / (C + T)
})

# make sure it lines up with your samples
names(global_weighted) <- sample.ids


# ===== Compare unweighted vs weighted CpG methylation =====

compare_methylation <- data.frame(
  Sample = sample.ids,
  Morph   = sub("-.*", "", sample.ids),
  Tissue  = sub("^[^-]+-[^_]+_([^_]+)_.*", "\\1", sample.ids),
  Condition = ifelse(cond_char == "F", "Fed", "Starved"),
  Unweighted_Mean_CpG_Methylation = as.numeric(sample_means),
  Weighted_Global_CpG_Methylation = as.numeric(global_weighted[sample.ids]),
  Difference_Weighted_minus_Unweighted =
    as.numeric(global_weighted[sample.ids]) - as.numeric(sample_means),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    Morph = dplyr::recode(Morph, "SF" = "Rio Choy", "PA" = "Pachón", "Mol" = "Molino", .default = Morph),
    Tissue = factor(Tissue, levels = c("Brain", "Liver", "Muscle")),
    Condition = factor(Condition, levels = c("Fed", "Starved"))
  )

print(compare_methylation)

