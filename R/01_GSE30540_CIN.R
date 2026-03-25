# =============================================================================
# 01_GSE30540_CIN.R
# Differential expression: CIN-high vs CIN-low in colorectal cancer
# Dataset: GSE30540 (Habermann et al., 2012)
#
# Question: Which genes discriminate CIN-high from CIN-low CRC phenotypes?
# CIN status is stored in pheno$description for this dataset.
# =============================================================================

source("R/00_setup.R")

# =============================================================================
# 1. Data Retrieval
# =============================================================================

message("Downloading GSE30540 from GEO...")
gse   <- getGEO("GSE30540", GSEMatrix = TRUE, getGPL = FALSE)
eset  <- gse[[1]]
expr  <- exprs(eset)
pheno <- pData(eset)

message(sprintf("Dataset dimensions: %d probes x %d samples", nrow(expr), ncol(expr)))

# =============================================================================
# 2. Phenotype Parsing
# =============================================================================
# CIN status lives in pheno$description ("CIN-high" / "CIN-low")

cin_group <- factor(pheno$description, levels = c("CIN-low", "CIN-high"))
message("Group counts:"); print(table(cin_group))

# Keep only samples with a valid CIN label
keep      <- !is.na(cin_group)
expr      <- expr[, keep]
cin_group <- droplevels(cin_group[keep])

# =============================================================================
# 3. Quality Control - PCA
# =============================================================================

pca_result <- prcomp(t(expr), scale. = TRUE)
pca_df <- data.frame(
  PC1   = pca_result$x[, 1],
  PC2   = pca_result$x[, 2],
  Group = cin_group
)
pca_var <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = c("CIN-low" = "#2c7bb6", "CIN-high" = "#d7191c")) +
  labs(
    title    = "PCA - GSE30540: CIN-high vs CIN-low",
    subtitle = "Microarray gene expression, CRC patient samples",
    x        = paste0("PC1 (", pca_var[1], "% variance)"),
    y        = paste0("PC2 (", pca_var[2], "% variance)"),
    colour   = "CIN Status"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/GSE30540_PCA.png", pca_plot, width = 7, height = 5, dpi = 300)
message("Saved: figures/GSE30540_PCA.png")

# =============================================================================
# 4. Differential Expression Analysis - limma
# =============================================================================

design <- model.matrix(~ cin_group)
colnames(design) <- c("Intercept", "CIN_high_vs_low")

fit <- lmFit(expr, design)
fit <- eBayes(fit)

all_results <- topTable(fit, coef = "CIN_high_vs_low", number = Inf, sort.by = "P")
write.csv(all_results, "results/GSE30540_all_DEGs.csv", row.names = TRUE)

# Threshold: nominal p < 0.05 and |log2FC| > 1
# Note: With 54,675 probes and only 10 CIN-low samples, BH-adjusted p-values
# do not reach 0.05 (minimum adj.P ~ 0.40). Nominal p-value threshold is used
# here as an exploratory analysis, consistent with the original publication.
# Results should be interpreted as hypothesis-generating, not confirmatory.
sig_degs <- all_results %>%
  filter(abs(logFC) > 1, P.Value < 0.05)

message(sprintf("Significant DEGs: %d total (%d up in CIN-high, %d down in CIN-high)",
                nrow(sig_degs),
                sum(sig_degs$logFC > 0),
                sum(sig_degs$logFC < 0)))

write.csv(sig_degs, "results/GSE30540_significant_DEGs.csv", row.names = TRUE)

# =============================================================================
# 5. Volcano Plot
# =============================================================================

all_results <- all_results %>%
  mutate(significance = case_when(
    logFC > 1  & adj.P.Val < 0.05 ~ "Upregulated",
    logFC < -1 & adj.P.Val < 0.05 ~ "Downregulated",
    TRUE                           ~ "Not significant"
  ))

top_labels <- all_results %>%
  filter(significance != "Not significant") %>%
  arrange(adj.P.Val) %>%
  head(15)

volcano_plot <- ggplot(all_results, aes(x = logFC, y = -log10(adj.P.Val), colour = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "Upregulated"     = "#d7191c",
    "Downregulated"   = "#2c7bb6",
    "Not significant" = "grey70"
  )) +
  labs(
    title  = "Volcano Plot - GSE30540: CIN-high vs CIN-low",
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

ggsave("figures/GSE30540_volcano.png", volcano_plot, width = 8, height = 6, dpi = 300)
message("Saved: figures/GSE30540_volcano.png")

# =============================================================================
# 6. DEG Count Bar Chart
# =============================================================================

up_count   <- sum(sig_degs$logFC > 0)
down_count <- sum(sig_degs$logFC < 0)

if (up_count + down_count > 0) {
  dist_df <- data.frame(
    Expression = c("Upregulated in CIN-high", "Downregulated in CIN-high"),
    Count      = c(up_count, down_count)
  )
  bar_plot <- ggplot(dist_df, aes(x = Expression, y = Count, fill = Expression)) +
    geom_col(width = 0.55, show.legend = FALSE) +
    geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
    scale_fill_manual(values = c(
      "Upregulated in CIN-high"   = "#d7191c",
      "Downregulated in CIN-high" = "#2c7bb6"
    )) +
    labs(
      title    = "Differentially Expressed Genes - GSE30540",
      subtitle = paste0("CIN-high vs CIN-low | Total significant DEGs: ", nrow(sig_degs)),
      x = NULL, y = "Number of genes"
    ) +
    theme_classic(base_size = 13) +
    ylim(0, max(dist_df$Count) * 1.15)

  ggsave("figures/GSE30540_DEG_counts.png", bar_plot, width = 6, height = 5, dpi = 300)
  message("Saved: figures/GSE30540_DEG_counts.png")
} else {
  message("No significant DEGs — skipping bar chart.")
}

# =============================================================================
# 7. Heatmap of Top DEGs
# =============================================================================

n_heatmap <- min(40, nrow(sig_degs))

if (n_heatmap >= 2) {
  top_probes         <- sig_degs %>% arrange(adj.P.Val) %>% head(n_heatmap) %>% rownames()
  heatmap_mat        <- expr[top_probes, ]
  heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

  annotation_col <- data.frame(CIN = cin_group)
  rownames(annotation_col) <- colnames(heatmap_mat)
  ann_colors <- list(CIN = c("CIN-low" = "#2c7bb6", "CIN-high" = "#d7191c"))

  png("figures/GSE30540_heatmap_top40.png", width = 10, height = 10, units = "in", res = 300)
  pheatmap(
    heatmap_mat_scaled,
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    show_colnames     = FALSE,
    fontsize_row      = 7,
    cluster_cols      = TRUE,
    cluster_rows      = TRUE,
    color             = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100),
    main              = paste0("Top ", n_heatmap, " DEGs - GSE30540 (CIN-high vs CIN-low)")
  )
  dev.off()
  message("Saved: figures/GSE30540_heatmap_top40.png")
} else {
  message("Fewer than 2 significant DEGs - skipping heatmap.")
}

message("=== GSE30540 analysis complete ===")
