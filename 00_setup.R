# =============================================================================
# 00_setup.R
# Package installation and loading for CRC GEO Analysis
# =============================================================================

# --- Bioconductor packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("GEOquery", "limma", "Biobase")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# --- CRAN packages ---
cran_packages <- c("ggplot2", "pheatmap", "dplyr", "tidyr", "ggrepel", "RColorBrewer")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# --- Load all ---
library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)

message("All packages loaded successfully.")

# --- Create output directories if they don't exist ---
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
