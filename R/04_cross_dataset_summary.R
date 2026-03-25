# =============================================================================
# 04_cross_dataset_summary.R
# Cross-dataset DEG overlap and summary
#
# After running scripts 01–03, this script:
#   - Loads significant DEGs from all three datasets
#   - Identifies any genes that appear across multiple datasets
#   - Produces a clean summary table
#
# Note: Probe IDs differ across platforms; gene symbol matching is used.
# =============================================================================

source("R/00_setup.R")

# =============================================================================
# 1. Load Results
# =============================================================================

# Check all result files exist before proceeding
result_files <- c(
  "results/GSE30540_significant_DEGs.csv",
  "results/GSE4107_significant_DEGs.csv",
  "results/GSE32323_tissue_significant_DEGs.csv"
)

missing <- result_files[!file.exists(result_files)]
if (length(missing) > 0) {
  stop("Missing result files. Run scripts 01–03 first:\n", paste(missing, collapse = "\n"))
}

degs_30540 <- read.csv("results/GSE30540_significant_DEGs.csv", row.names = 1)
degs_4107  <- read.csv("results/GSE4107_significant_DEGs.csv",  row.names = 1)
degs_32323 <- read.csv("results/GSE32323_tissue_significant_DEGs.csv", row.names = 1)

# =============================================================================
# 2. Summary Table
# =============================================================================

summary_df <- data.frame(
  Dataset     = c("GSE30540", "GSE4107", "GSE32323"),
  Comparison  = c(
    "CIN-high vs CIN-low",
    "CRC mucosa vs Healthy controls",
    "Tumour vs Normal tissue (paired)"
  ),
  Total_DEGs  = c(nrow(degs_30540), nrow(degs_4107), nrow(degs_32323)),
  Upregulated = c(
    sum(degs_30540$logFC > 0),
    sum(degs_4107$logFC  > 0),
    sum(degs_32323$logFC > 0)
  ),
  Downregulated = c(
    sum(degs_30540$logFC < 0),
    sum(degs_4107$logFC  < 0),
    sum(degs_32323$logFC < 0)
  )
)

print(summary_df)
write.csv(summary_df, "results/cross_dataset_summary.csv", row.names = FALSE)
message("Saved: results/cross_dataset_summary.csv")

# =============================================================================
# 3. Summary Bar Chart
# =============================================================================

plot_df <- summary_df %>%
  select(Dataset, Upregulated, Downregulated) %>%
  tidyr::pivot_longer(cols = c(Upregulated, Downregulated),
                      names_to = "Direction", values_to = "Count") %>%
  mutate(Count_plot = ifelse(Direction == "Downregulated", -Count, Count))

summary_plot <- ggplot(plot_df, aes(x = Dataset, y = Count_plot, fill = Direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Upregulated" = "#d7191c", "Downregulated" = "#2c7bb6")) +
  scale_y_continuous(labels = abs) +
  labs(
    title    = "Differentially Expressed Genes Across Three CRC Datasets",
    subtitle = "|log2FC| > 1 and adjusted p-value < 0.05 (BH correction)",
    x        = NULL,
    y        = "Number of significant DEGs",
    fill     = NULL,
    caption  = "Up = upregulated in the condition of interest; Down = downregulated"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("figures/cross_dataset_DEG_summary.png", summary_plot, width = 8, height = 5.5, dpi = 300)
message("Saved: figures/cross_dataset_DEG_summary.png")

# =============================================================================
# 4. Gene Overlap (if gene symbols are available as rownames)
# =============================================================================
# Note: GEO probe IDs are platform-specific. If rownames are gene symbols
# (e.g. after annotation), we can find overlaps. Otherwise skip this section.

genes_30540 <- rownames(degs_30540)
genes_4107  <- rownames(degs_4107)
genes_32323 <- rownames(degs_32323)

overlap_all   <- Reduce(intersect, list(genes_30540, genes_4107, genes_32323))
overlap_30_41 <- intersect(genes_30540, genes_4107)
overlap_30_32 <- intersect(genes_30540, genes_32323)
overlap_41_32 <- intersect(genes_4107,  genes_32323)

cat("\n--- Gene/Probe Overlaps ---\n")
cat("GSE30540 ∩ GSE4107:            ", length(overlap_30_41), "\n")
cat("GSE30540 ∩ GSE32323:           ", length(overlap_30_32), "\n")
cat("GSE4107  ∩ GSE32323:           ", length(overlap_41_32), "\n")
cat("All three datasets:            ", length(overlap_all),   "\n\n")

if (length(overlap_all) > 0) {
  cat("Genes present in all three datasets:\n")
  print(overlap_all)
  write.csv(data.frame(Gene = overlap_all),
            "results/cross_dataset_shared_genes.csv", row.names = FALSE)
}

message("=== Cross-dataset summary complete ===")
message("All outputs saved to results/ and figures/")
