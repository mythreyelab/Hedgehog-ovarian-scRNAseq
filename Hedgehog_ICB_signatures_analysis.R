# =============================================================================
# Hedgehog Pathway vs Published ICB Signatures (TISCH-style)
# Hypothesis: HH activation → immune cold / immunosuppressive phenotype
# Dataset: OV_GSE130000
# =============================================================================

library(hdf5r)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Set working directory to the repo root (folder containing OV_GSE130000_expression.h5
# and OV_GSE130000_CellMetainfo_table.tsv). Edit the path below or use setwd() interactively.
# Example:
#   setwd("/path/to/GSE130000")
# If the data files live one directory above this script (default layout in our deposit),
# uncomment the next line:
# setwd("..")
stopifnot(
  file.exists("OV_GSE130000_expression.h5"),
  file.exists("OV_GSE130000_CellMetainfo_table.tsv")
)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("=== Loading data ===\n")
h5 <- H5File$new("OV_GSE130000_expression.h5", mode = "r")
barcodes <- h5[["matrix/barcodes"]]$read()
genes    <- h5[["matrix/features/name"]]$read()
data_vals <- h5[["matrix/data"]]$read()
indices  <- h5[["matrix/indices"]]$read()
indptr   <- h5[["matrix/indptr"]]$read()
shape    <- h5[["matrix/shape"]]$read()
expr_mat <- sparseMatrix(i = indices + 1L, p = indptr, x = data_vals, dims = shape)
rownames(expr_mat) <- genes
colnames(expr_mat) <- barcodes
h5$close_all()

meta <- read.delim("OV_GSE130000_CellMetainfo_table.tsv", stringsAsFactors = FALSE)
rownames(meta) <- meta$Cell
common_cells <- intersect(colnames(expr_mat), rownames(meta))
expr_mat <- expr_mat[, common_cells]
meta <- meta[common_cells, ]
col_sums <- colSums(expr_mat)

cat("Total cells:", length(common_cells), "\n")

# =============================================================================
# 2. PUBLISHED ICB SIGNATURES
# =============================================================================

# References:
# Checkpoint: Auslander 2018, Ayers 2017
# CTL: Rooney 2015 (cytotoxic lymphocyte)
# IFNG: Ayers 2017 (IFN-gamma 6-gene)
# IMPRES: Auslander 2018 (15 immune checkpoint pairs)
# Inflammatory: Spranger 2015
# IPRES: Hugo 2016 (innate anti-PD1 resistance)
# T cell-inflamed: Ayers 2017 (18-gene T cell-inflamed GEP)
# TLS: Cabrita 2020 (tertiary lymphoid structures)
# T-quiescent: Peng 2019

icb_signatures <- list(

  Checkpoint = c(
    "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2",
    "TIGIT", "BTLA", "TNFRSF9", "ICOS", "IDO1", "CD276", "VTCN1"
  ),

  CTL = c(
    "GZMA", "GZMB", "GZMH", "GZMK", "PRF1",
    "GNLY", "NKG7", "KLRK1", "KLRD1", "CTSW", "CST7"
  ),

  IFNG = c(
    "IFNG", "STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA"
  ),

  IMPRES = c(
    "PDCD1", "CD274", "CTLA4", "CD28", "CD80", "CD86",
    "CD27", "CD40", "CD40LG", "ICOS", "ICOSLG",
    "TNFRSF4", "TNFRSF9", "TNFRSF18", "CD276"
  ),

  Inflammatory = c(
    "CCL2", "CCL3", "CCL4", "CCL5", "CXCL9", "CXCL10", "CXCL11",
    "IL1A", "IL1B", "IL6", "TNF", "IFNG",
    "IRF1", "STAT1", "NFKB1"
  ),

  IPRES = c(
    "AXL", "ROR2", "WNT5A", "LOXL2", "TWIST2", "TAGLN",
    "FAP", "VIM", "ACTA2", "COL4A1", "COL5A1", "COL5A2",
    "MMP1", "MMP2", "MMP3", "MMP9", "MMP10",
    "IL10", "VEGFA", "VEGFC", "TGFB1", "TGFB2"
  ),

  T_cell_inflamed = c(
    "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9",
    "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1",
    "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1",
    "TIGIT", "CCL5"
  ),

  TLS = c(
    "CD79B", "EIF1AY", "RBP5", "PTGDS", "CETP",
    "SFTPB", "GNG4", "DMBT1", "TNFRSF17",
    "LAMP3", "CCR7", "CCL19", "CCL21", "CXCL13"
  ),

  T_quiescent = c(
    "TCF7", "LEF1", "CCR7", "SELL", "IL7R",
    "CD28", "CD27", "BCL2", "FOXP1"
  )
)

# Check gene availability
cat("\nICB signature gene availability:\n")
for (sig in names(icb_signatures)) {
  found <- sum(icb_signatures[[sig]] %in% rownames(expr_mat))
  total <- length(icb_signatures[[sig]])
  missing <- icb_signatures[[sig]][!icb_signatures[[sig]] %in% rownames(expr_mat)]
  cat(sprintf("  %s: %d/%d found", sig, found, total))
  if (length(missing) > 0) cat(sprintf("  [missing: %s]", paste(missing, collapse=", ")))
  cat("\n")
}

# =============================================================================
# 3. SCORING
# =============================================================================

score_geneset <- function(mat, genes_in_set, col_sums) {
  genes_use <- intersect(genes_in_set, rownames(mat))
  if (length(genes_use) < 2) return(NULL)
  sub <- mat[genes_use, , drop = FALSE]
  norm_sub <- sweep(sub, 2, col_sums, "/") * 10000
  norm_sub <- log1p(norm_sub)
  score <- colMeans(as.matrix(norm_sub))
  return(list(score = score, n_genes = length(genes_use), genes_used = genes_use))
}

# Hedgehog genes
hh_all <- c("IHH", "SHH", "DHH", "PTCH1", "PTCH2", "HHIP", "GLI1", "GLI2")
hh_extended <- unique(c(hh_all, "SMO", "SUFU", "GLI3", "STK36", "KIF7",
                        "CDON", "BOC", "GAS1", "DISP1", "DISP2"))

# Score HH
hh_score_res <- score_geneset(expr_mat, hh_all, col_sums)
hh_ext_res   <- score_geneset(expr_mat, hh_extended, col_sums)

# Score all ICB signatures
icb_scores <- list()
for (sig in names(icb_signatures)) {
  res <- score_geneset(expr_mat, icb_signatures[[sig]], col_sums)
  if (!is.null(res)) icb_scores[[sig]] <- res
}

cat("\nAll signatures scored successfully:", length(icb_scores), "\n")

# Build master data frame
cell_df <- data.frame(
  Cell = common_cells,
  HH_core = hh_score_res$score[common_cells],
  HH_extended = if (!is.null(hh_ext_res)) hh_ext_res$score[common_cells] else NA,
  Celltype = meta$Celltype..major.lineage.,
  Tissue = meta$Tissue,
  Patient = meta$Patient,
  stringsAsFactors = FALSE
)
for (sig in names(icb_scores)) {
  cell_df[[sig]] <- icb_scores[[sig]]$score[common_cells]
}
cell_df$Tissue <- factor(cell_df$Tissue, levels = c("Tumor", "Metastatic", "Relapsed"))

# =============================================================================
# 4. CORRELATIONS: HH vs EACH ICB SIGNATURE
# =============================================================================

cat("\n=== HH Core vs ICB Signatures ===\n")

# --- 4A: All cells ---
cat("\n--- All cells ---\n")
cor_all <- data.frame()
for (sig in names(icb_scores)) {
  ct <- cor.test(cell_df$HH_core, cell_df[[sig]], method = "spearman")
  cor_all <- rbind(cor_all, data.frame(
    Signature = sig, Scope = "All_cells",
    N = nrow(cell_df),
    Spearman_rho = ct$estimate, pval = ct$p.value,
    stringsAsFactors = FALSE
  ))
}
cor_all$padj <- p.adjust(cor_all$pval, method = "BH")
cor_all$Sig <- ifelse(cor_all$padj < 0.05, "*", "")
print(cor_all %>% mutate(across(where(is.numeric), ~round(., 4))) %>%
        select(-pval) %>% arrange(Spearman_rho))

# --- 4B: Malignant cells only ---
cat("\n--- Malignant cells ---\n")
mal_df <- cell_df %>% filter(Celltype == "Malignant")
cor_mal <- data.frame()
for (sig in names(icb_scores)) {
  ct <- cor.test(mal_df$HH_core, mal_df[[sig]], method = "spearman")
  cor_mal <- rbind(cor_mal, data.frame(
    Signature = sig, Scope = "Malignant",
    N = nrow(mal_df),
    Spearman_rho = ct$estimate, pval = ct$p.value,
    stringsAsFactors = FALSE
  ))
}
cor_mal$padj <- p.adjust(cor_mal$pval, method = "BH")
cor_mal$Sig <- ifelse(cor_mal$padj < 0.05, "*", "")
print(cor_mal %>% mutate(across(where(is.numeric), ~round(., 4))) %>%
        select(-pval) %>% arrange(Spearman_rho))

# --- 4C: Malignant cells by tissue ---
cat("\n--- Malignant cells by tissue ---\n")
cor_mal_tis <- data.frame()
for (sig in names(icb_scores)) {
  for (tis in c("Tumor", "Metastatic", "Relapsed")) {
    sub <- mal_df %>% filter(Tissue == tis)
    if (nrow(sub) > 30) {
      ct <- cor.test(sub$HH_core, sub[[sig]], method = "spearman")
      cor_mal_tis <- rbind(cor_mal_tis, data.frame(
        Signature = sig, Tissue = tis,
        N = nrow(sub),
        Spearman_rho = ct$estimate, pval = ct$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}
cor_mal_tis$padj <- p.adjust(cor_mal_tis$pval, method = "BH")
cor_mal_tis$Sig <- ifelse(cor_mal_tis$padj < 0.05, "*", "")
cat("\nSignificant results:\n")
print(cor_mal_tis %>% filter(padj < 0.05) %>%
        mutate(across(where(is.numeric), ~round(., 4))) %>%
        select(-pval) %>% arrange(Spearman_rho))
cat("\nAll results:\n")
print(cor_mal_tis %>% mutate(across(where(is.numeric), ~round(., 4))) %>%
        select(-pval) %>% arrange(Signature, Tissue))

# --- 4D: Individual HH gene vs each ICB signature (all cells) ---
cat("\n--- Individual HH gene vs ICB signatures ---\n")
hh_genes_present <- hh_all[hh_all %in% rownames(expr_mat)]
gene_sig_cor <- data.frame()
for (g in hh_genes_present) {
  g_expr <- as.numeric(log1p(expr_mat[g, common_cells] / col_sums * 10000))
  for (sig in names(icb_scores)) {
    ct <- cor.test(g_expr, cell_df[[sig]], method = "spearman")
    gene_sig_cor <- rbind(gene_sig_cor, data.frame(
      Gene = g, Signature = sig,
      Spearman_rho = ct$estimate, pval = ct$p.value,
      stringsAsFactors = FALSE
    ))
  }
}
gene_sig_cor$padj <- p.adjust(gene_sig_cor$pval, method = "BH")
gene_sig_cor$Sig <- ifelse(gene_sig_cor$padj < 0.05, "*", "")

cat("\nSignificant gene-signature correlations:\n")
print(gene_sig_cor %>% filter(padj < 0.05) %>%
        mutate(across(where(is.numeric), ~round(., 4))) %>%
        select(-pval) %>% arrange(Spearman_rho))

# --- 4E: Per-patient pseudobulk ---
cat("\n--- Patient-level pseudobulk (malignant) ---\n")
sig_names <- names(icb_scores)
patient_df <- cell_df %>%
  filter(Celltype == "Malignant") %>%
  group_by(Patient, Tissue) %>%
  summarise(across(c(HH_core, all_of(sig_names)), mean), N = n(), .groups = "drop")

cor_patient <- data.frame()
for (sig in sig_names) {
  if (nrow(patient_df) >= 5) {
    ct <- cor.test(patient_df$HH_core, patient_df[[sig]], method = "spearman")
    cor_patient <- rbind(cor_patient, data.frame(
      Signature = sig,
      Spearman_rho = ct$estimate, pval = ct$p.value,
      stringsAsFactors = FALSE
    ))
  }
}
print(cor_patient %>% mutate(across(where(is.numeric), ~round(., 4))) %>% arrange(Spearman_rho))

# =============================================================================
# 5. VISUALIZATIONS
# =============================================================================

cat("\n=== Generating figures ===\n")
outdir <- "Hedgehog"

# --- Fig 1: Correlation heatmap — HH genes vs ICB signatures ---

heat_df <- gene_sig_cor %>%
  select(Gene, Signature, Spearman_rho, Sig) %>%
  mutate(Gene = factor(Gene, levels = c("IHH", "PTCH1", "PTCH2", "HHIP", "GLI1", "GLI2")),
         Signature = gsub("_", " ", Signature))

p1 <- ggplot(heat_df, aes(x = Signature, y = Gene, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Sig), size = 7, fontface = "bold") +
  scale_fill_gradient2(low = "#7B2D8E", mid = "white", high = "#E69F00",
                       midpoint = 0, limits = c(-0.15, 0.15),
                       name = "Spearman\nrho") +
  theme_minimal() +
  ggtitle("Hedgehog Genes vs ICB Signatures\n(All Cells, * = padj < 0.05)") +
  xlab("") + ylab("") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))

ggsave(file.path(outdir, "Fig9_HH_genes_vs_ICB_signatures_correlation_heatmap.pdf"), p1, width = 10, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig9_HH_genes_vs_ICB_signatures_correlation_heatmap.png"), p1, width = 10, height = 5, dpi = 300)
cat("Saved: Fig9_HH_genes_vs_ICB_signatures_correlation_heatmap\n")

# --- Fig 2: Bar plot — HH core pathway vs each ICB signature (malignant cells) ---

cor_mal$Signature <- factor(cor_mal$Signature,
                            levels = cor_mal$Signature[order(cor_mal$Spearman_rho)])

p2 <- ggplot(cor_mal, aes(x = Signature, y = Spearman_rho,
                           fill = ifelse(Spearman_rho > 0, "Positive", "Negative"))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = Sig,
                y = ifelse(Spearman_rho > 0,
                           Spearman_rho + 0.001,
                           Spearman_rho - 0.001)),
            hjust = 0.5, vjust = 0.5, size = 7, fontface = "bold") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("Positive" = "#E69F00", "Negative" = "#7B2D8E"), name = "") +
  coord_flip() +
  theme_bw() +
  ggtitle("HH Pathway Score vs ICB Signatures\n(Malignant Cells)") +
  ylab("Spearman rho") + xlab("") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "none")

ggsave(file.path(outdir, "Fig10_HH_pathway_vs_ICB_signatures_barplot_malignant.pdf"), p2, width = 8, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig10_HH_pathway_vs_ICB_signatures_barplot_malignant.png"), p2, width = 8, height = 5, dpi = 300)
cat("Saved: Fig10_HH_pathway_vs_ICB_signatures_barplot_malignant\n")

# --- Fig 3: Tissue-stratified heatmap (malignant cells) ---

heat_tis <- cor_mal_tis %>%
  mutate(label = paste0(round(Spearman_rho, 3), Sig),
         Signature = gsub("_", " ", Signature),
         Tissue = factor(Tissue, levels = c("Tumor", "Metastatic", "Relapsed")))

p3 <- ggplot(heat_tis, aes(x = Tissue, y = Signature, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "#7B2D8E", mid = "white", high = "#E69F00",
                       midpoint = 0, name = "Spearman\nrho") +
  theme_minimal() +
  ggtitle("HH Pathway vs ICB Signatures by Disease State\n(Malignant Cells, * = padj < 0.05)") +
  xlab("") + ylab("") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"))

ggsave(file.path(outdir, "Fig11_HH_pathway_vs_ICB_signatures_by_tissue.pdf"), p3, width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "Fig11_HH_pathway_vs_ICB_signatures_by_tissue.png"), p3, width = 8, height = 6, dpi = 300)
cat("Saved: Fig11_HH_pathway_vs_ICB_signatures_by_tissue\n")

# --- Fig 4: Scatter plots for strongest correlations (malignant, tumor) ---

top_sigs <- cor_mal_tis %>%
  filter(Tissue == "Tumor") %>%
  arrange(Spearman_rho) %>%
  head(4) %>%
  pull(Signature)

tumor_mal <- mal_df %>% filter(Tissue == "Tumor")

scatter_list <- list()
for (i in seq_along(top_sigs)) {
  sig <- top_sigs[i]
  rho <- cor_mal_tis %>% filter(Signature == sig, Tissue == "Tumor") %>% pull(Spearman_rho)
  pv  <- cor_mal_tis %>% filter(Signature == sig, Tissue == "Tumor") %>% pull(padj)
  scatter_list[[i]] <- ggplot(tumor_mal, aes(x = HH_core, y = .data[[sig]])) +
    geom_point(size = 0.3, alpha = 0.2, color = "#0072B2") +
    geom_smooth(method = "lm", color = "#D55E00", se = TRUE) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho=%.3f\npadj=%.2e", rho, pv),
             hjust = 1.1, vjust = 1.3, size = 3.5) +
    theme_bw() +
    ggtitle(gsub("_", " ", sig)) +
    xlab("HH Pathway Score") + ylab("ICB Score") +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
}

p4 <- wrap_plots(scatter_list, nrow = 1) +
  plot_annotation(
    title = "Strongest Negative HH–ICB Correlations (Primary Tumor, Malignant)",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(outdir, "Fig12_HH_vs_ICB_top_correlations_scatter_tumor.pdf"), p4, width = 16, height = 4.5, dpi = 300)
ggsave(file.path(outdir, "Fig12_HH_vs_ICB_top_correlations_scatter_tumor.png"), p4, width = 16, height = 4.5, dpi = 300)
cat("Saved: Fig12_HH_vs_ICB_top_correlations_scatter_tumor\n")

# --- Fig 5: GLI2 specifically vs all ICB signatures (hypothesis gene) ---

gli2_cors <- gene_sig_cor %>% filter(Gene == "GLI2")
gli2_cors$Signature <- factor(gsub("_", " ", gli2_cors$Signature),
                               levels = gsub("_", " ", gli2_cors$Signature[order(gli2_cors$Spearman_rho)]))

p5 <- ggplot(gli2_cors, aes(x = Signature, y = Spearman_rho,
                              fill = ifelse(Spearman_rho > 0, "Positive", "Negative"))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = Sig,
                y = ifelse(Spearman_rho > 0,
                           Spearman_rho + 0.001,
                           Spearman_rho - 0.001)),
            hjust = 0.5, vjust = 0.5, size = 7, fontface = "bold") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("Positive" = "#E69F00", "Negative" = "#7B2D8E"), name = "") +
  coord_flip() +
  theme_bw() +
  ggtitle("GLI2 Expression vs ICB Signatures\n(All Cells)") +
  ylab("Spearman rho") + xlab("") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "none")

ggsave(file.path(outdir, "Fig13_GLI2_vs_ICB_signatures_barplot.pdf"), p5, width = 8, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig13_GLI2_vs_ICB_signatures_barplot.png"), p5, width = 8, height = 5, dpi = 300)
cat("Saved: Fig13_GLI2_vs_ICB_signatures_barplot\n")

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

all_cors <- bind_rows(
  cor_all %>% mutate(Tissue = "All"),
  cor_mal %>% mutate(Tissue = "All_Malignant"),
  cor_mal_tis
)

write.csv(all_cors %>% mutate(across(where(is.numeric), ~round(., 6))),
          file.path(outdir, "HH_vs_ICB_signatures_correlations.csv"), row.names = FALSE)
write.csv(gene_sig_cor %>% mutate(across(where(is.numeric), ~round(., 6))),
          file.path(outdir, "HH_genes_vs_ICB_signatures_correlations.csv"), row.names = FALSE)
write.csv(cor_patient %>% mutate(across(where(is.numeric), ~round(., 4))),
          file.path(outdir, "HH_vs_ICB_patient_pseudobulk.csv"), row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
