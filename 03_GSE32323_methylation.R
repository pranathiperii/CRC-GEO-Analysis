# =============================================================================
# 03_GSE32323_methylation.R
# Dataset: GSE32323 (Khamas et al., 2012)
#
# This dataset contains two sub-experiments:
#   A) Matched tumour vs normal tissue (17 patient pairs)
#   B) 5-aza-2'-deoxycytidine treatment in 5 CRC cell lines
#
# Analysis A identifies tumour vs normal DEGs in patient tissue.
# Analysis B identifies epigenetically silenced genes reactivated by 
# demethylation — genes up after 5-aza are candidates for methylation silencing.
# =============================================================================

source("R/00_setup.R")

# =============================================================================
# 1. Data Retrieval
# =============================================================================

message("Downloading GSE32323 from GEO...")
gse   <- getGEO("GSE32323", GSEMatrix = TRUE, getGPL = FALSE)
eset  <- gse[[1]]
expr  <- exprs(eset)
pheno <- pData(eset)

message(sprintf("Dataset dimensions: %d probes x %d samples", nrow(expr), ncol(expr)))

# =============================================================================
# 2. Split into two sub-experiments based on title
# =============================================================================

titles <- pheno$title

# Sub-experiment A: matched patient tissue (normal vs cancer)
tissue_idx    <- grep("patient", titles, ignore.case = TRUE)
tissue_group  <- ifelse(grepl("normal", titles[tissue_idx], ignore.case = TRUE),
                        "Normal", "Cancer")
tissue_group  <- factor(tissue_group, levels = c("Normal", "Cancer"))
expr_tissue   <- expr[, tissue_idx]

message("--- Sub-experiment A: Patient tissue ---")
print(table(tissue_group))

# Sub-experiment B: cell lines, 5-aza vs control
cellline_idx   <- grep("cell line", titles, ignore.case = TRUE)
cellline_group <- ifelse(grepl("5aza|5-aza", titles[cellline_idx], ignore.case = TRUE),
                         "AzaTreated", "Control")
cellline_group <- factor(cellline_group, levels = c("Control", "AzaTreated"))
expr_cellline  <- expr[, cellline_idx]

message("--- Sub-experiment B: Cell lines ---")
print(table(cellline_group))

# =============================================================================
# 3. PCA — Patient Tissue
# =============================================================================

pca_result <- prcomp(t(expr_tissue), scale. = TRUE)
pca_df <- data.frame(
  PC1   = pca_result$x[, 1],
  PC2   = pca_result$x[, 2],
  Group = tissue_group
)
pca_var <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = c("Normal" = "#1a9641", "Cancer" = "#d7191c")) +
  labs(
    title    = "PCA - GSE32323: Tumour vs Normal Tissue",
    subtitle = "17 matched patient pairs",
    x        = paste0("PC1 (", pca_var[1], "% variance)"),
    y        = paste0("PC2 (", pca_var[2], "% variance)"),
    colour   = "Tissue"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/GSE32323_tissue_PCA.png", pca_plot, width = 7, height = 5, dpi = 300)
message("Saved: figures/GSE32323_tissue_PCA.png")

# =============================================================================
# 4. DEA — Patient Tissue (paired design using patient as blocking factor)
# =============================================================================

# Extract patient IDs for paired design
patient_ids <- gsub("patient (\\d+b),.*", "\\1", titles[tissue_idx])
patient_ids <- factor(patient_ids)

design_tissue <- model.matrix(~ patient_ids + tissue_group)

fit_tissue  <- lmFit(expr_tissue, design_tissue)
fit_tissue  <- eBayes(fit_tissue)

results_tissue <- topTable(fit_tissue, coef = "tissue_groupCancer",
                            number = Inf, sort.by = "P")
write.csv(results_tissue, "results/GSE32323_tissue_all_DEGs.csv", row.names = TRUE)

sig_tissue <- results_tissue %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05)

message(sprintf("Tissue DEGs (adj.P<0.05, |logFC|>1): %d total (%d up in tumour, %d down in tumour)",
                nrow(sig_tissue),
                sum(sig_tissue$logFC > 0),
                sum(sig_tissue$logFC < 0)))

write.csv(sig_tissue, "results/GSE32323_tissue_significant_DEGs.csv", row.names = TRUE)

# =============================================================================
# 5. Volcano Plot — Patient Tissue
# =============================================================================

results_tissue <- results_tissue %>%
  mutate(significance = case_when(
    logFC > 1  & adj.P.Val < 0.05 ~ "Upregulated in tumour",
    logFC < -1 & adj.P.Val < 0.05 ~ "Downregulated in tumour",
    TRUE                           ~ "Not significant"
  ))

top_labels <- results_tissue %>%
  filter(significance != "Not significant") %>%
  arrange(adj.P.Val) %>%
  head(15)

volcano_tissue <- ggplot(results_tissue, aes(x = logFC, y = -log10(adj.P.Val), colour = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "Upregulated in tumour"   = "#d7191c",
    "Downregulated in tumour" = "#1a9641",
    "Not significant"         = "grey70"
  )) +
  labs(
    title  = "Volcano Plot - GSE32323: Tumour vs Normal (paired)",
    x      = "log2 Fold Change",
    y      = "-log10(Adjusted P-value)",
    colour = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

if (nrow(top_labels) > 0) {
  volcano_tissue <- volcano_tissue +
    geom_text_repel(data = top_labels, aes(label = rownames(top_labels)),
                    size = 3, max.overlaps = 20)
}

ggsave("figures/GSE32323_tissue_volcano.png", volcano_tissue, width = 8, height = 6, dpi = 300)
message("Saved: figures/GSE32323_tissue_volcano.png")

# =============================================================================
# 6. DEA — Cell Lines, 5-aza vs Control
# =============================================================================

pca_cell <- prcomp(t(expr_cellline), scale. = TRUE)
pca_cell_df <- data.frame(
  PC1   = pca_cell$x[, 1],
  PC2   = pca_cell$x[, 2],
  Group = cellline_group
)
pca_var_cell <- round(summary(pca_cell)$importance[2, 1:2] * 100, 1)

pca_cell_plot <- ggplot(pca_cell_df, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = c("Control" = "#756bb1", "AzaTreated" = "#fd8d3c")) +
  labs(
    title    = "PCA - GSE32323: 5-aza-dC Treated vs Control (Cell Lines)",
    subtitle = "5 CRC cell lines",
    x        = paste0("PC1 (", pca_var_cell[1], "% variance)"),
    y        = paste0("PC2 (", pca_var_cell[2], "% variance)"),
    colour   = "Treatment"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/GSE32323_cellline_PCA.png", pca_cell_plot, width = 7, height = 5, dpi = 300)
message("Saved: figures/GSE32323_cellline_PCA.png")

design_cell <- model.matrix(~ cellline_group)
colnames(design_cell) <- c("Intercept", "AzaTreated_vs_Control")

fit_cell  <- lmFit(expr_cellline, design_cell)
fit_cell  <- eBayes(fit_cell)

results_cell <- topTable(fit_cell, coef = "AzaTreated_vs_Control",
                          number = Inf, sort.by = "P")
write.csv(results_cell, "results/GSE32323_cellline_all_DEGs.csv", row.names = TRUE)

# Use nominal p here — only 5 cell lines per group
sig_cell <- results_cell %>%
  filter(abs(logFC) > 1, P.Value < 0.05)

epi_candidates <- sig_cell %>% filter(logFC > 0) %>% arrange(P.Value)

message(sprintf("Cell line DEGs (nominal p<0.05, |logFC|>1): %d (%d up after 5-aza, %d down)",
                nrow(sig_cell),
                sum(sig_cell$logFC > 0),
                sum(sig_cell$logFC < 0)))
message(sprintf("Epigenetically silenced candidates (up after demethylation): %d", nrow(epi_candidates)))

write.csv(sig_cell,       "results/GSE32323_cellline_significant_DEGs.csv",   row.names = TRUE)
write.csv(epi_candidates, "results/GSE32323_epigenetic_candidates.csv",        row.names = TRUE)

# =============================================================================
# 7. Volcano Plot — Cell Lines
# =============================================================================

results_cell <- results_cell %>%
  mutate(significance = case_when(
    logFC > 1  & P.Value < 0.05 ~ "Reactivated by demethylation",
    logFC < -1 & P.Value < 0.05 ~ "Suppressed by demethylation",
    TRUE                         ~ "Not significant"
  ))

top_cell_labels <- results_cell %>%
  filter(significance != "Not significant") %>%
  arrange(P.Value) %>%
  head(15)

volcano_cell <- ggplot(results_cell, aes(x = logFC, y = -log10(P.Value), colour = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "Reactivated by demethylation" = "#fd8d3c",
    "Suppressed by demethylation"  = "#756bb1",
    "Not significant"              = "grey70"
  )) +
  labs(
    title    = "Volcano Plot - GSE32323: 5-aza-dC vs Control (Cell Lines)",
    subtitle = "Genes reactivated = candidate epigenetically silenced genes",
    x        = "log2 Fold Change",
    y        = "-log10(P-value)",
    colour   = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

if (nrow(top_cell_labels) > 0) {
  volcano_cell <- volcano_cell +
    geom_text_repel(data = top_cell_labels, aes(label = rownames(top_cell_labels)),
                    size = 3, max.overlaps = 20)
}

ggsave("figures/GSE32323_cellline_volcano.png", volcano_cell, width = 8, height = 6, dpi = 300)
message("Saved: figures/GSE32323_cellline_volcano.png")

# =============================================================================
# 8. Heatmap — Patient Tissue Top DEGs
# =============================================================================

n_heatmap <- min(40, nrow(sig_tissue))

if (n_heatmap >= 2) {
  top_probes         <- sig_tissue %>% arrange(adj.P.Val) %>% head(n_heatmap) %>% rownames()
  heatmap_mat        <- expr_tissue[top_probes, ]
  heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

  annotation_col <- data.frame(Tissue = tissue_group)
  rownames(annotation_col) <- colnames(heatmap_mat)
  ann_colors <- list(Tissue = c("Normal" = "#1a9641", "Cancer" = "#d7191c"))

  png("figures/GSE32323_tissue_heatmap.png", width = 10, height = 10, units = "in", res = 300)
  pheatmap(
    heatmap_mat_scaled,
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    show_colnames     = FALSE,
    fontsize_row      = 7,
    cluster_cols      = TRUE,
    cluster_rows      = TRUE,
    color             = colorRampPalette(c("#1a9641", "white", "#d7191c"))(100),
    main              = paste0("Top ", n_heatmap, " DEGs - GSE32323 Tumour vs Normal")
  )
  dev.off()
  message("Saved: figures/GSE32323_tissue_heatmap.png")
}

message("=== GSE32323 analysis complete ===")
