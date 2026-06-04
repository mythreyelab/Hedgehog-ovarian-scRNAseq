# =============================================================================
# Hedgehog Pathway Activity & Correlation with ICB Signature
# Dataset: OV_GSE130000 (ovarian cancer scRNA-seq)
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

expr_mat <- sparseMatrix(
  i = indices + 1L, p = indptr, x = data_vals, dims = shape
)
rownames(expr_mat) <- genes
colnames(expr_mat) <- barcodes
h5$close_all()

meta <- read.delim("OV_GSE130000_CellMetainfo_table.tsv", stringsAsFactors = FALSE)
rownames(meta) <- meta$Cell
common_cells <- intersect(colnames(expr_mat), rownames(meta))
expr_mat <- expr_mat[, common_cells]
meta <- meta[common_cells, ]

cat("Total cells:", length(common_cells), "\n")

# =============================================================================
# 2. DEFINE HEDGEHOG PATHWAY GENES
# =============================================================================

hh_ligands <- c("SHH", "IHH", "DHH")
hh_receptors_effectors <- c("PTCH1", "PTCH2", "HHIP", "GLI1", "GLI2")
hh_all <- c(hh_ligands, hh_receptors_effectors)

# Extended Hedgehog pathway (for pathway scoring)
hh_extended <- c(
  hh_all,
  "SMO", "SUFU", "GLI3",        # additional core components
  "STK36", "KIF7", "CDON",      # signal transduction
  "BOC", "GAS1",                 # co-receptors
  "DISP1", "DISP2",             # ligand release
  "HHIP", "PTCH2"               # negative regulators
)
hh_extended <- unique(hh_extended)

cat("\nHedgehog genes in dataset:\n")
for (g in hh_all) {
  present <- g %in% rownames(expr_mat)
  cat(sprintf("  %s: %s\n", g, ifelse(present, "FOUND", "NOT FOUND")))
}
cat("\nExtended set found:", sum(hh_extended %in% rownames(expr_mat)), "/", length(hh_extended), "\n")

# =============================================================================
# 3. ICB (IMMUNE CHECKPOINT BLOCKADE) SIGNATURE GENES
# =============================================================================

# Immune checkpoint receptors and ligands
icb_genes <- c(
  # Inhibitory checkpoints (targets of ICB therapy)
  "PDCD1",     # PD-1
  "CD274",     # PD-L1
  "PDCD1LG2",  # PD-L2
  "CTLA4",
  "LAG3",
  "HAVCR2",    # TIM-3
  "TIGIT",
  "VSIR",      # VISTA
  "CD276",     # B7-H3
  "VTCN1",     # B7-H4
  "IDO1",
  "SIGLEC15",
  # Co-stimulatory / immune activation
  "CD80", "CD86",
  "TNFRSF9",   # 4-1BB
  "TNFRSF4",   # OX40
  "ICOS", "ICOSLG",
  # Cytotoxic / effector markers
  "GZMA", "GZMB", "PRF1", "IFNG",
  "CXCL9", "CXCL10", "CXCL11",
  # T cell markers
  "CD8A", "CD8B", "CD3E",
  # Antigen presentation
  "HLA-A", "HLA-B", "HLA-C",
  "B2M", "TAP1", "TAP2"
)

cat("\nICB signature genes in dataset:\n")
cat("  Found:", sum(icb_genes %in% rownames(expr_mat)), "/", length(icb_genes), "\n")
missing_icb <- icb_genes[!icb_genes %in% rownames(expr_mat)]
if (length(missing_icb) > 0) cat("  Missing:", paste(missing_icb, collapse = ", "), "\n")

# =============================================================================
# 4. SCORE FUNCTIONS
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

col_sums <- colSums(expr_mat)

# =============================================================================
# 5. PART 1 — HEDGEHOG PATHWAY ACTIVITY
# =============================================================================

cat("\n=== PART 1: Hedgehog Pathway Activity ===\n")

# --- 5A: Score across ALL cell types ---

hh_score_all <- score_geneset(expr_mat, hh_all, col_sums)
hh_score_ext <- score_geneset(expr_mat, hh_extended, col_sums)

cat("\nHedgehog core genes used:", paste(hh_score_all$genes_used, collapse = ", "), "\n")

# Build per-cell data
cell_df <- data.frame(
  Cell = common_cells,
  HH_score = hh_score_all$score[common_cells],
  HH_ext_score = if (!is.null(hh_score_ext)) hh_score_ext$score[common_cells] else NA,
  Celltype = meta$Celltype..major.lineage.,
  Tissue = meta$Tissue,
  Patient = meta$Patient,
  Cluster = as.character(meta$Cluster),
  stringsAsFactors = FALSE
)
cell_df$Tissue <- factor(cell_df$Tissue, levels = c("Tumor", "Metastatic", "Relapsed"))

# --- Mean HH score by cell type and tissue ---
cat("\nHedgehog pathway score by cell type and tissue:\n")
summary_ct <- cell_df %>%
  group_by(Celltype, Tissue) %>%
  summarise(
    N = n(),
    Mean_HH = round(mean(HH_score), 4),
    SD_HH = round(sd(HH_score), 4),
    .groups = "drop"
  ) %>%
  arrange(Celltype, Tissue)
print(as.data.frame(summary_ct))

# --- 5B: Individual gene expression (all cells) ---

hh_genes_present <- hh_all[hh_all %in% rownames(expr_mat)]
hh_norm <- log1p(sweep(as.matrix(expr_mat[hh_genes_present, ]), 2, col_sums, "/") * 10000)

gene_by_celltype <- data.frame()
for (g in hh_genes_present) {
  for (ct in unique(meta$Celltype..major.lineage.)) {
    cells_ct <- rownames(meta)[meta$Celltype..major.lineage. == ct]
    for (tis in c("Tumor", "Metastatic", "Relapsed")) {
      cells_tis <- intersect(cells_ct, rownames(meta)[meta$Tissue == tis])
      if (length(cells_tis) > 0) {
        gene_by_celltype <- rbind(gene_by_celltype, data.frame(
          Gene = g,
          Celltype = ct,
          Tissue = tis,
          Mean_Expr = mean(hh_norm[g, cells_tis]),
          Pct_Expressing = mean(hh_norm[g, cells_tis] > 0) * 100,
          N_cells = length(cells_tis),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}
gene_by_celltype$Tissue <- factor(gene_by_celltype$Tissue, levels = c("Tumor", "Metastatic", "Relapsed"))

cat("\nHedgehog gene expression summary (top expressing):\n")
print(gene_by_celltype %>%
        filter(Pct_Expressing > 1) %>%
        arrange(desc(Mean_Expr)) %>%
        head(30) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# --- 5C: Malignant cell focus ---

mal_cells <- rownames(meta)[meta$Celltype..major.lineage. == "Malignant"]
mal_mat <- expr_mat[, mal_cells]
mal_meta <- meta[mal_cells, ]
mal_col_sums <- colSums(mal_mat)

cat("\n--- Malignant cells only ---\n")
cat("Malignant cells:", length(mal_cells), "\n")

tumor_cells    <- mal_cells[mal_meta$Tissue == "Tumor"]
met_cells      <- mal_cells[mal_meta$Tissue == "Metastatic"]
relapsed_cells <- mal_cells[mal_meta$Tissue == "Relapsed"]

# Individual gene stats in malignant cells
mal_hh_norm <- log1p(sweep(as.matrix(mal_mat[hh_genes_present, ]), 2, mal_col_sums, "/") * 10000)

mal_gene_df <- data.frame(
  Gene = hh_genes_present,
  Tumor_Mean = sapply(hh_genes_present, function(g) mean(mal_hh_norm[g, tumor_cells])),
  Met_Mean = sapply(hh_genes_present, function(g) mean(mal_hh_norm[g, met_cells])),
  Relapsed_Mean = sapply(hh_genes_present, function(g) mean(mal_hh_norm[g, relapsed_cells])),
  Tumor_Pct = sapply(hh_genes_present, function(g) mean(mal_hh_norm[g, tumor_cells] > 0) * 100),
  Met_Pct = sapply(hh_genes_present, function(g) mean(mal_hh_norm[g, met_cells] > 0) * 100),
  Relapsed_Pct = sapply(hh_genes_present, function(g) mean(mal_hh_norm[g, relapsed_cells] > 0) * 100),
  stringsAsFactors = FALSE
)
mal_gene_df$Delta_Rel_Tum <- mal_gene_df$Relapsed_Mean - mal_gene_df$Tumor_Mean

# Wilcoxon tests
mal_gene_df$pval <- sapply(hh_genes_present, function(g) {
  wilcox.test(mal_hh_norm[g, tumor_cells], mal_hh_norm[g, relapsed_cells])$p.value
})
mal_gene_df$padj <- p.adjust(mal_gene_df$pval, method = "BH")
mal_gene_df$Sig <- ifelse(mal_gene_df$padj < 0.05, "*", "")

cat("\nHedgehog genes in malignant cells (Tumor vs Relapsed):\n")
print(mal_gene_df %>%
        mutate(across(where(is.numeric), ~round(., 4))) %>%
        select(-pval) %>%
        arrange(desc(abs(Delta_Rel_Tum))))

# =============================================================================
# 6. PART 2 — CORRELATION WITH ICB SIGNATURE
# =============================================================================

cat("\n=== PART 2: Hedgehog vs ICB Correlation ===\n")

# Score ICB signature
icb_score_res <- score_geneset(expr_mat, icb_genes, col_sums)
cat("ICB signature genes used:", icb_score_res$n_genes, "\n")

cell_df$ICB_score <- icb_score_res$score[common_cells]

# --- 6A: Per-cell correlation (all cells) ---

cor_all <- cor.test(cell_df$HH_score, cell_df$ICB_score, method = "spearman")
cat(sprintf("\nAll cells — Spearman rho: %.4f, p-value: %.2e\n", cor_all$estimate, cor_all$p.value))

# --- 6B: Correlation by cell type ---

cat("\nCorrelation by cell type:\n")
cor_by_ct <- cell_df %>%
  group_by(Celltype) %>%
  summarise(
    N = n(),
    Spearman_rho = cor(HH_score, ICB_score, method = "spearman"),
    pval = cor.test(HH_score, ICB_score, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(pval, method = "BH"),
         Sig = ifelse(padj < 0.05, "*", ""))
print(as.data.frame(cor_by_ct %>% mutate(across(where(is.numeric), ~round(., 4)))))

# --- 6C: Correlation by tissue in malignant cells ---

mal_df <- cell_df %>% filter(Celltype == "Malignant")

cat("\nCorrelation in malignant cells by tissue:\n")
cor_by_tis <- mal_df %>%
  group_by(Tissue) %>%
  summarise(
    N = n(),
    Spearman_rho = cor(HH_score, ICB_score, method = "spearman"),
    pval = cor.test(HH_score, ICB_score, method = "spearman")$p.value,
    .groups = "drop"
  )
print(as.data.frame(cor_by_tis %>% mutate(across(where(is.numeric), ~round(., 4)))))

# --- 6D: Per-patient pseudobulk correlation ---

patient_df <- cell_df %>%
  filter(Celltype == "Malignant") %>%
  group_by(Patient, Tissue) %>%
  summarise(
    HH_mean = mean(HH_score),
    ICB_mean = mean(ICB_score),
    N = n(),
    .groups = "drop"
  )

cat("\nPer-patient pseudobulk (malignant cells):\n")
print(as.data.frame(patient_df %>% mutate(across(where(is.numeric), ~round(., 4)))))

if (nrow(patient_df) >= 5) {
  cor_patient <- cor.test(patient_df$HH_mean, patient_df$ICB_mean, method = "spearman")
  cat(sprintf("\nPatient-level Spearman rho: %.4f, p-value: %.4f\n",
              cor_patient$estimate, cor_patient$p.value))
}

# --- 6E: Individual HH gene vs ICB score correlations ---

cat("\nIndividual HH gene vs ICB score (all cells):\n")
gene_icb_cor <- data.frame()
for (g in hh_genes_present) {
  g_expr <- log1p(expr_mat[g, common_cells] / col_sums * 10000)
  ct <- cor.test(as.numeric(g_expr), cell_df$ICB_score, method = "spearman")
  gene_icb_cor <- rbind(gene_icb_cor, data.frame(
    Gene = g,
    Spearman_rho = ct$estimate,
    pval = ct$p.value,
    stringsAsFactors = FALSE
  ))
}
gene_icb_cor$padj <- p.adjust(gene_icb_cor$pval, method = "BH")
gene_icb_cor$Sig <- ifelse(gene_icb_cor$padj < 0.05, "*", "")
print(gene_icb_cor %>% mutate(across(where(is.numeric), ~round(., 4))) %>% arrange(desc(abs(Spearman_rho))))

# =============================================================================
# 7. VISUALIZATIONS
# =============================================================================

cat("\n=== Generating figures ===\n")

tissue_colors <- c("Tumor" = "#009E73", "Metastatic" = "#7B2D8E", "Relapsed" = "#E69F00")
celltype_colors <- c(
  "CD8Tex"         = "#D55E00",
  "Fibroblasts"    = "#0072B2",
  "Malignant"      = "#009E73",
  "Mono/Macro"     = "#CC79A7",
  "Myofibroblasts" = "#56B4E9"
)

outdir <- "Hedgehog"

# --- Fig 1: HH pathway score violin by cell type and tissue ---

p1 <- ggplot(cell_df, aes(x = Tissue, y = HH_score, fill = Tissue)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  facet_wrap(~ Celltype, nrow = 1) +
  scale_fill_manual(values = tissue_colors) +
  theme_bw() +
  ggtitle("Hedgehog Pathway Activity by Cell Type and Disease State") +
  ylab("HH Pathway Score") + xlab("") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "Fig1_HH_pathway_violin_by_celltype_and_tissue.pdf"), p1, width = 16, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig1_HH_pathway_violin_by_celltype_and_tissue.png"), p1, width = 16, height = 5, dpi = 300)
cat("Saved: Fig1_HH_pathway_violin_by_celltype_and_tissue\n")

# --- Fig 2: Individual gene dot plot (malignant cells) ---

gene_plot_df <- mal_gene_df %>%
  select(Gene, Tumor_Mean, Met_Mean, Relapsed_Mean) %>%
  pivot_longer(cols = c(Tumor_Mean, Met_Mean, Relapsed_Mean),
               names_to = "Tissue", values_to = "Expr") %>%
  mutate(Tissue = gsub("_Mean", "", Tissue),
         Tissue = recode(Tissue, "Met" = "Metastatic"),
         Tissue = factor(Tissue, levels = c("Tumor", "Metastatic", "Relapsed")))

pct_plot_df <- mal_gene_df %>%
  select(Gene, Tumor_Pct, Met_Pct, Relapsed_Pct) %>%
  pivot_longer(cols = c(Tumor_Pct, Met_Pct, Relapsed_Pct),
               names_to = "Tissue", values_to = "Pct") %>%
  mutate(Tissue = gsub("_Pct", "", Tissue),
         Tissue = recode(Tissue, "Met" = "Metastatic"),
         Tissue = factor(Tissue, levels = c("Tumor", "Metastatic", "Relapsed")))

dot_df <- left_join(gene_plot_df, pct_plot_df, by = c("Gene", "Tissue"))

# Categorize genes
dot_df$Category <- ifelse(dot_df$Gene %in% hh_ligands, "Ligand", "Receptor/Effector")
dot_df$Gene <- factor(dot_df$Gene, levels = c(hh_ligands, hh_receptors_effectors))

p2 <- ggplot(dot_df, aes(x = Tissue, y = Gene, size = Pct, color = Expr)) +
  geom_point() +
  scale_color_gradient2(low = "#7B2D8E", mid = "#F7F7F7", high = "#E69F00",
                        midpoint = median(dot_df$Expr), name = "Mean Expr\n(log-norm)") +
  scale_size_continuous(range = c(1, 8), name = "% Expressing") +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  ggtitle("Hedgehog Gene Expression in Malignant Cells") +
  xlab("") + ylab("") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        strip.text.y = element_text(face = "bold", angle = 0),
        axis.text.x = element_text(face = "bold", size = 11))

ggsave(file.path(outdir, "Fig2_HH_gene_dotplot_malignant_cells.pdf"), p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "Fig2_HH_gene_dotplot_malignant_cells.png"), p2, width = 8, height = 6, dpi = 300)
cat("Saved: Fig2_HH_gene_dotplot_malignant_cells\n")

# --- Fig 3: HH gene heatmap across all cell types ---

ct_heatmap <- gene_by_celltype %>%
  mutate(Gene = factor(Gene, levels = c(hh_ligands, hh_receptors_effectors)))

p3 <- ggplot(ct_heatmap, aes(x = Celltype, y = Gene, fill = Mean_Expr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Pct_Expressing, 1)), size = 2.2) +
  scale_fill_gradient(low = "white", high = "#D55E00", name = "Mean Expr") +
  facet_wrap(~ Tissue, nrow = 1) +
  theme_minimal() +
  ggtitle("Hedgehog Genes Across Cell Types\n(numbers = % cells expressing)") +
  xlab("") + ylab("") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 11))

ggsave(file.path(outdir, "Fig3_HH_gene_heatmap_all_celltypes_all_tissues.pdf"), p3, width = 14, height = 6, dpi = 300)
ggsave(file.path(outdir, "Fig3_HH_gene_heatmap_all_celltypes_all_tissues.png"), p3, width = 14, height = 6, dpi = 300)
cat("Saved: Fig3_HH_gene_heatmap_all_celltypes_all_tissues\n")

# --- Fig 4: HH gene heatmap — primary tumor only ---

ct_heatmap_tumor <- gene_by_celltype %>%
  filter(Tissue == "Tumor") %>%
  mutate(Gene = factor(Gene, levels = c(hh_ligands, hh_receptors_effectors)))

p3b <- ggplot(ct_heatmap_tumor, aes(x = Celltype, y = Gene, fill = Mean_Expr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Pct_Expressing)), size = 3) +
  scale_fill_gradient(low = "white", high = "#D55E00", name = "Mean Expr\n(log-norm)") +
  theme_minimal() +
  ggtitle("Hedgehog Gene Expression Across Cell Types\n(Primary Tumor Only)") +
  xlab("") + ylab("") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"))

ggsave(file.path(outdir, "Fig4_HH_gene_heatmap_all_celltypes_tumor_only.pdf"), p3b, width = 8, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig4_HH_gene_heatmap_all_celltypes_tumor_only.png"), p3b, width = 8, height = 5, dpi = 300)
cat("Saved: Fig4_HH_gene_heatmap_all_celltypes_tumor_only\n")

# --- Fig 6: HH vs ICB scatter (per cell, malignant) ---

p4 <- ggplot(mal_df, aes(x = HH_score, y = ICB_score, color = Tissue)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = tissue_colors) +
  facet_wrap(~ Tissue, nrow = 1) +
  theme_bw() +
  ggtitle("Hedgehog vs ICB Signature — Malignant Cells") +
  xlab("Hedgehog Pathway Score") + ylab("ICB Signature Score") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "none")

# Add correlation text
cor_labels <- cor_by_tis %>%
  mutate(label = sprintf("rho=%.3f\np=%.2e", Spearman_rho, pval))
p4 <- p4 + geom_text(data = cor_labels,
                      aes(x = Inf, y = Inf, label = label),
                      hjust = 1.1, vjust = 1.3, size = 3.5, color = "black",
                      inherit.aes = FALSE)

ggsave(file.path(outdir, "Fig6_HH_vs_compositeICB_scatter_malignant.pdf"), p4, width = 14, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig6_HH_vs_compositeICB_scatter_malignant.png"), p4, width = 14, height = 5, dpi = 300)
cat("Saved: Fig6_HH_vs_compositeICB_scatter_malignant\n")

# --- Fig 7: HH vs ICB per-patient pseudobulk scatter ---

p5 <- ggplot(patient_df, aes(x = HH_mean, y = ICB_mean, color = Tissue, label = Patient)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  scale_color_manual(values = tissue_colors) +
  theme_bw() +
  ggtitle("Patient-Level: Hedgehog vs ICB Score (Malignant Cells)") +
  xlab("Mean HH Pathway Score") + ylab("Mean ICB Signature Score") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = "bottom")

ggsave(file.path(outdir, "Fig7_HH_vs_compositeICB_patient_pseudobulk.pdf"), p5, width = 8, height = 7, dpi = 300)
ggsave(file.path(outdir, "Fig7_HH_vs_compositeICB_patient_pseudobulk.png"), p5, width = 8, height = 7, dpi = 300)
cat("Saved: Fig7_HH_vs_compositeICB_patient_pseudobulk\n")

# --- Fig 8: Individual HH gene vs ICB correlation bar plot ---

gene_icb_cor$Gene <- factor(gene_icb_cor$Gene, levels = gene_icb_cor$Gene[order(gene_icb_cor$Spearman_rho)])

p6 <- ggplot(gene_icb_cor, aes(x = Gene, y = Spearman_rho,
                                fill = ifelse(Spearman_rho > 0, "Positive", "Negative"))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = Sig,
                y = ifelse(Spearman_rho > 0,
                           Spearman_rho + 0.001,
                           Spearman_rho - 0.001)),
            hjust = 0.5, vjust = 0.5, size = 7, fontface = "bold") +
  scale_fill_manual(values = c("Positive" = "#E69F00", "Negative" = "#7B2D8E"), name = "") +
  coord_flip() +
  theme_bw() +
  ggtitle("Hedgehog Gene Correlation with ICB Signature") +
  ylab("Spearman rho") + xlab("") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = "bottom")

ggsave(file.path(outdir, "Fig8_HH_genes_vs_compositeICB_barplot.pdf"), p6, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig8_HH_genes_vs_compositeICB_barplot.png"), p6, width = 7, height = 5, dpi = 300)
cat("Saved: Fig8_HH_genes_vs_compositeICB_barplot\n")

# --- Fig 5: HH pathway by patient ---

p7 <- ggplot(cell_df %>% filter(Celltype == "Malignant"),
             aes(x = Patient, y = HH_score, fill = Tissue)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", alpha = 0.5) +
  scale_fill_manual(values = tissue_colors) +
  theme_bw() +
  ggtitle("Hedgehog Pathway Score by Patient (Malignant Cells)") +
  ylab("HH Pathway Score") + xlab("") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(outdir, "Fig5_HH_pathway_violin_by_patient.pdf"), p7, width = 10, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig5_HH_pathway_violin_by_patient.png"), p7, width = 10, height = 5, dpi = 300)
cat("Saved: Fig5_HH_pathway_violin_by_patient\n")

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

write.csv(summary_ct, file.path(outdir, "HH_pathway_scores_by_celltype_tissue.csv"), row.names = FALSE)
write.csv(mal_gene_df %>% mutate(across(where(is.numeric), ~round(., 4))),
          file.path(outdir, "HH_gene_expression_malignant.csv"), row.names = FALSE)
write.csv(gene_by_celltype %>% mutate(across(where(is.numeric), ~round(., 4))),
          file.path(outdir, "HH_gene_expression_all_celltypes.csv"), row.names = FALSE)
write.csv(cor_by_ct %>% mutate(across(where(is.numeric), ~round(., 4))),
          file.path(outdir, "HH_ICB_correlation_by_celltype.csv"), row.names = FALSE)
write.csv(cor_by_tis %>% mutate(across(where(is.numeric), ~round(., 4))),
          file.path(outdir, "HH_ICB_correlation_malignant_by_tissue.csv"), row.names = FALSE)
write.csv(patient_df %>% mutate(across(where(is.numeric), ~round(., 4))),
          file.path(outdir, "HH_ICB_patient_pseudobulk.csv"), row.names = FALSE)
write.csv(gene_icb_cor %>% mutate(across(where(is.numeric), ~round(., 4))),
          file.path(outdir, "HH_gene_vs_ICB_correlation.csv"), row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All outputs saved to:", file.path(getwd(), outdir), "\n")
