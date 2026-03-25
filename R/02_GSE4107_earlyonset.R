# =============================================================================
# 02_GSE4107_earlyonset.R
# Differential expression: CRC mucosa vs healthy controls
# Dataset: GSE4107 (Sabates-Bellver et al., 2007)
#
# Question: Which genes are upregulated in CRC mucosa vs healthy tissue?
# These may serve as early detection biomarkers.
# =============================================================================

source("R/00_setup.R")

# =============================================================================
# 1. Data Retrieval
# =============================================================================

message("Downloading GSE4107 from GEO...")
gse <- getGEO("GSE4107", GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]

expr  <- exprs(eset)
pheno <- pData(eset)

message(sprintf("Dataset dimensions: %d probes x %d samples", nrow(expr), ncol(expr)))

# =============================================================================
# 2. Normalization
# =============================================================================
# GSE4107 expression values are on a linear scale (range ~1 to 360,000).
# Log2 transformation is required before DEA.

if (max(expr, na.rm = TRUE) > 100) {
  message("Data appears to be on linear scale — applying log2(x+1) transformation.")
  expr <- log2(expr + 1)
}

# =============================================================================
# 2. Phenotype Parsing
# =============================================================================

# Inspect available columns
# print(colnames(pheno))

sample_group <- ifelse(grepl("healthy control", pheno$title, ignore.case = TRUE),
                       "Normal", "CRC")
sample_group <- factor(sample_group, levels = c("Normal", "CRC"))

message(table(sample_group))

# =============================================================================
# 3. Quality Control — PCA
# =============================================================================

pca_result <- prcomp(t(expr), scale. = TRUE)
pca_df <- data.frame(
  PC1   = pca_result$x[, 1],
  PC2   = pca_result$x[, 2],
  Group = sample_group
)

pca_var <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = c("Normal" = "#1a9641", "CRC" = "#d7191c")) +
  labs(
    title    = "PCA — GSE4107: CRC Mucosa vs Healthy Controls",
    subtitle = "Microarray gene expression, colon mucosa biopsies",
    x        = paste0("PC1 (", pca_var[1], "% variance)"),
    y        = paste0("PC2 (", pca_var[2], "% variance)"),
    colour   = "Group"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/GSE4107_PCA.png", pca_plot, width = 7, height = 5, dpi = 300)
message("Saved: figures/GSE4107_PCA.png")

# =============================================================================
# 4. Differential Expression Analysis — limma
# =============================================================================

keep         <- !is.na(sample_group)
expr         <- expr[, keep]
sample_group <- droplevels(sample_group[keep])

design <- model.matrix(~ sample_group)
colnames(design) <- c("Intercept", "CRC_vs_Normal")

fit  <- lmFit(expr, design)
fit  <- eBayes(fit)

all_results <- topTable(fit, coef = "CRC_vs_Normal", number = Inf, sort.by = "P")
write.csv(all_results, "results/GSE4107_all_DEGs.csv", row.names = TRUE)

# Threshold: BH-adjusted p < 0.05 and |log2FC| > 1
# Data is log2 transformed — fold changes are now meaningful.
sig_degs <- all_results %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05)

message(sprintf("Significant DEGs: %d total (%d up in CRC, %d down in CRC)",
                nrow(sig_degs),
                sum(sig_degs$logFC > 0),
                sum(sig_degs$logFC < 0)))

write.csv(sig_degs, "results/GSE4107_significant_DEGs.csv", row.names = TRUE)

# =============================================================================
# 5. Volcano Plot
# =============================================================================

all_results <- all_results %>%
  mutate(
    significance = case_when(
      logFC > 1  & adj.P.Val < 0.05 ~ "Upregulated in CRC",
      logFC < -1 & adj.P.Val < 0.05 ~ "Downregulated in CRC",
      TRUE                           ~ "Not significant"
    )
  )

top_labels <- all_results %>%
  filter(significance != "Not significant") %>%
  arrange(adj.P.Val) %>%
  head(15)

volcano_plot <- ggplot(all_results, aes(x = logFC, y = -log10(adj.P.Val), colour = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "Upregulated in CRC"   = "#d7191c",
    "Downregulated in CRC" = "#2c7bb6",
    "Not significant"      = "grey70"
  )) +
  labs(
    title  = "Volcano Plot - GSE4107: CRC Mucosa vs Healthy Controls",
    x      = "log2 Fold Change",
    y      = "-log10(Adjusted P-value)",
    colour = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

if (nrow(top_labels) > 0) {
  volcano_plot <- volcano_plot +
    geom_text_repel(data = top_labels, aes(label = rownames(top_labels)),
                    size = 3, max.overlaps = 20)
}

ggsave("figures/GSE4107_volcano.png", volcano_plot, width = 8, height = 6, dpi = 300)
message("Saved: figures/GSE4107_volcano.png")

# =============================================================================
# 6. Heatmap of Top 40 DEGs
# =============================================================================

top40_probes <- sig_degs %>%
  slice_min(adj.P.Val, n = 40) %>%
  rownames()

if (length(top40_probes) >= 2) {
  heatmap_mat        <- expr[top40_probes, ]
  heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

  annotation_col <- data.frame(Group = sample_group)
  rownames(annotation_col) <- colnames(heatmap_mat)
  ann_colors <- list(Group = c(Normal = "#1a9641", CRC = "#d7191c"))

  png("figures/GSE4107_heatmap_top40.png", width = 10, height = 10, units = "in", res = 300)
  pheatmap(
    heatmap_mat_scaled,
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    show_colnames     = FALSE,
    fontsize_row      = 7,
    cluster_cols      = TRUE,
    cluster_rows      = TRUE,
    color             = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100),
    main              = "Top 40 DEGs — GSE4107 (CRC vs Normal Mucosa)"
  )
  dev.off()
  message("Saved: figures/GSE4107_heatmap_top40.png")
} else {
  message("Fewer than 2 significant DEGs — skipping heatmap.")
}

message("=== GSE4107 analysis complete ===")
