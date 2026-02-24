#__________________________Cross_Population_Venn_______________________________-
library(tidyverse)
library(VennDiagram)
library(grid)
library(ggplot2)

# === Step 1: Define file paths ===
folder <- "/path/to/folder"
files <- list.files(folder, pattern = "^Intersects_.*\\.txt$", full.names = TRUE)

# === Step 2: Morph name mapping ===
morph_from_code <- function(code) {
  if (grepl("S", code)) return("Rio Choy")
  if (grepl("P", code)) return("Pachón")
  if (grepl("M", code)) return("Molino")
  return("Unknown")
  
}

# === Step 3: Extract morph pair from filename ===
get_morph_pair <- function(fname) {
  short <- gsub("Intersects_|\\.txt", "", basename(fname))
  parts <- unlist(strsplit(short, "_"))
  morph1 <- morph_from_code(parts[1])
  morph2 <- morph_from_code(parts[2])
  paste(sort(c(morph1, morph2)), collapse = " vs ")
}

# === Step 4–5: Load gene sets and force correct order and color mapping ===
grouped_files <- split(files, sapply(files, get_morph_pair))
gene_sets <- lapply(grouped_files, function(file_group) {
  genes <- unlist(lapply(file_group, function(f) {
    read.table(f, stringsAsFactors = FALSE, header = FALSE)[[1]]
  }))
  unique(trimws(genes[!is.na(genes) & genes != ""]))
})

# === Manually set desired Venn order and color mapping ===
venn_order <- c("Molino vs Rio Choy", "Molino vs Pachón", "Pachón vs Rio Choy")
venn_colors <- c("#8E9B4D", "#8C4F76", "#4DA3B6")  # green, purple, blue

# Check that all required comparisons exist
missing <- setdiff(venn_order, names(gene_sets))
if (length(missing) > 0) stop(paste("Missing expected comparisons:", paste(missing, collapse = ", ")))

# Reorder gene sets to match color and position
gene_sets <- gene_sets[venn_order]
names(gene_sets) <- venn_order

# === Step 6: Overlap calculations ===
area1 <- length(setdiff(gene_sets[[1]], union(gene_sets[[2]], gene_sets[[3]])))
area2 <- length(setdiff(gene_sets[[2]], union(gene_sets[[1]], gene_sets[[3]])))
area3 <- length(setdiff(gene_sets[[3]], union(gene_sets[[1]], gene_sets[[2]])))
a12   <- length(setdiff(intersect(gene_sets[[1]], gene_sets[[2]]), gene_sets[[3]]))
a13   <- length(setdiff(intersect(gene_sets[[1]], gene_sets[[3]]), gene_sets[[2]]))
a23   <- length(setdiff(intersect(gene_sets[[2]], gene_sets[[3]]), gene_sets[[1]]))
a123  <- length(Reduce(intersect, gene_sets))
total <- area1 + area2 + area3 + a12 + a13 + a23 + a123

make_label <- function(count) {
  pct <- round(100 * count / total, 1)
  paste0(count, "\n(", pct, "%)")
}
labels <- lapply(c(area1, area2, a12, area3, a13, a23, a123), make_label)

# === Step 7: Draw Venn ===
venn.plot <- draw.triple.venn(
  area1 = length(gene_sets[[1]]),
  area2 = length(gene_sets[[2]]),
  area3 = length(gene_sets[[3]]),
  n12   = length(intersect(gene_sets[[1]], gene_sets[[2]])),
  n13   = length(intersect(gene_sets[[1]], gene_sets[[3]])),
  n23   = length(intersect(gene_sets[[2]], gene_sets[[3]])),
  n123  = a123,
  category = names(gene_sets),
  fill = venn_colors,
  alpha = 0.7,
  cex = 0.01,
  cat.cex = 1.5,
  cat.fontface = "bold",
  label.col = "white",
  print.mode = "raw",
  ind = TRUE,
  scaled = FALSE
)

draw_label <- function(label, x, y, fontsize = 14) {
  grid.text(label, x = x, y = y, gp = gpar(fontsize = fontsize, fontface = "bold"))
}

# === Step 8: Add count labels manually ===
grid.text("Gut", y = 1, gp = gpar(fontsize = 18, fontface = "bold"))
draw_label(labels[[1]], x = 0.25, y = 0.75, fontsize = 20)  # Molino vs Rio
draw_label(labels[[2]], x = 0.75, y = 0.75, fontsize = 20)  # Molino vs Pachón
draw_label(labels[[3]], x = 0.5,  y = 0.8,  fontsize = 20)  # a12
draw_label(labels[[4]], x = 0.5,  y = 0.25, fontsize = 20)  # Pachón vs Rio
draw_label(labels[[5]], x = 0.3,  y = 0.45, fontsize = 20)  # a13
draw_label(labels[[6]], x = 0.7,  y = 0.45, fontsize = 20)  # a23
draw_label(labels[[7]], x = 0.5,  y = 0.55, fontsize = 20)  # a123

#__________________________metabolic_Venn_______________________________

library(tidyverse)
library(VennDiagram)
library(svDialogs)

# === Step 1: Select folder ===
folder <- dlg_dir(title = "Select folder with Intersects_ files")$res
files <- list.files(folder, pattern = "^Intersects_[UD]_.*\\.txt$", full.names = TRUE)

# === Step 2: Extract morph from filename ===
extract_morph <- function(filepath) {
  parts <- str_match(basename(filepath), "^Intersects_[UD]_(\\w{3,4})_(\\w{3,4})\\.txt")
  triplet1 <- parts[2]
  triplet2 <- parts[3]
  
  # Extract morph letters
  morph1 <- substr(triplet1, 2, 2)
  morph2 <- substr(triplet2, 2, 2)
  
  # Ensure both triplets refer to the same morph
  if (morph1 != morph2) return(NA_character_)
  
  case_when(
    morph1 == "M" ~ "Molino",
    morph1 == "P" ~ "Pachón",
    morph1 == "S" ~ "Rio Choy",
    TRUE ~ NA_character_
  )
}

# === Step 3: Read and merge U + D gene sets per morph ===
gene_sets <- list()

for (file in files) {
  morph <- extract_morph(file)
  if (!is.na(morph)) {
    genes <- read.table(file, stringsAsFactors = FALSE)[[1]]
    gene_sets[[morph]] <- unique(c(gene_sets[[morph]], genes))
  }
}

# === Step 4: Check ===
if (length(gene_sets) != 3) {
  stop("You need exactly 3 morphs represented. Check file naming or morph assignment.")
}

# === Step 5: Set final morph order and colors ===
gene_sets <- gene_sets[c("Pachón", "Rio Choy", "Molino")]
pretty_labels <- names(gene_sets)

venn_colors <- c(
  "Pachón" = "#E57E81",   # rose pink
  "Rio Choy" = "#9CA7CE", # blue
  "Molino" = "#ce9e91"    # beige pink
)[pretty_labels]

# === Step 6: Overlap counts ===
area1 <- length(setdiff(gene_sets[[1]], union(gene_sets[[2]], gene_sets[[3]])))
area2 <- length(setdiff(gene_sets[[2]], union(gene_sets[[1]], gene_sets[[3]])))
area3 <- length(setdiff(gene_sets[[3]], union(gene_sets[[1]], gene_sets[[2]])))
a12   <- length(setdiff(intersect(gene_sets[[1]], gene_sets[[2]]), gene_sets[[3]]))
a13   <- length(setdiff(intersect(gene_sets[[1]], gene_sets[[3]]), gene_sets[[2]]))
a23   <- length(setdiff(intersect(gene_sets[[2]], gene_sets[[3]]), gene_sets[[1]]))
a123  <- length(Reduce(intersect, gene_sets))
total <- area1 + area2 + area3 + a12 + a13 + a23 + a123

# === Step 7: Labels ===
make_label <- function(count) {
  pct <- round(100 * count / total, 1)
  paste0(count, "\n(", pct, "%)")
}
labels <- lapply(c(area1, area2, a12, area3, a13, a23, a123), make_label)

# === Step 8: Plot ===
venn.plot <- draw.triple.venn(
  area1 = length(gene_sets[[1]]),
  area2 = length(gene_sets[[2]]),
  area3 = length(gene_sets[[3]]),
  n12 = length(intersect(gene_sets[[1]], gene_sets[[2]])),
  n13 = length(intersect(gene_sets[[1]], gene_sets[[3]])),
  n23 = length(intersect(gene_sets[[2]], gene_sets[[3]])),
  n123 = a123,
  category = pretty_labels,
  fill = venn_colors,
  alpha = 0.75,
  cex = 0.01,
  cat.cex = 1.5,
  cat.fontface = "bold",
  label.col = "white",
  print.mode = "raw",
  ind = TRUE,
  scaled = FALSE
)

# === Step 9: Add labels ===
draw_label <- function(label, x, y, fontsize = 16) {
  grid.text(label, x = x, y = y,
            gp = gpar(fontsize = fontsize, fontface = "bold"))
}

grid.text("", y = 0.99, gp = gpar(fontsize = 22, fontface = "bold"))


draw_label(labels[[1]], x = 0.25, y = 0.75,  fontsize = 25)  # Pachón
draw_label(labels[[2]], x = 0.75, y = 0.75,  fontsize = 25)  # Rio Choy
 draw_label(labels[[4]], x = 0.5,  y = 0.25,  fontsize = 25)  # Molino
 draw_label(labels[[3]], x = 0.51,  y = 0.65,  fontsize = 12)   # P ∩ R
 draw_label(labels[[5]], x = 0.35, y = 0.52,  fontsize = 15)  # P ∩ M
 draw_label(labels[[6]], x = 0.63, y = 0.52,  fontsize = 15)  # R ∩ M draw_label(labels[[7]], x = 0.5,  y = 0.55,  fontsize = 20)  # All 3
