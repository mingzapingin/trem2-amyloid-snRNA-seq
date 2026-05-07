# ============================================================
# export_shiny_data.R
#
# Extracts and copies all files needed for the R Shiny app
# into a single folder: ./shiny_data/
#
# Outputs:
#   shiny_data/sample_info.csv         — sample metadata
#   shiny_data/normalized_counts.csv   — global pseudobulk normalized counts
#   shiny_data/de_results.csv          — Q1 global DE full results
#   shiny_data/gsea_results.csv        — Q1 broad excitatory GSEA
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
})

select <- dplyr::select
filter <- dplyr::filter

OUT_DIR <- "./shiny_data"
dir.create(OUT_DIR, showWarnings = FALSE)

# --- Paths (edit if your layout differs) ---
SEURAT_PATH  <- "./processed/seurat_7m_final_annotated.rds"
DE_PATH      <- "./result/de_results/Q1_global_WT_vs_WT5XFAD/Q1_global_WT_vs_WT5XFAD_Global_full_results.csv"
GSEA_PATH    <- "./result/gsea_results/broad/Q1_broad_WT_vs_WT5XFAD/Excitatory_neuron_GSEA_ALL_databases.csv"

# =============================================================
# 1. Sample info (with per-sample QC metrics)
# =============================================================
message("Loading Seurat object...")
obj <- readRDS(SEURAT_PATH)

sample_info <- obj@meta.data %>%
  group_by(sample_id) %>%
  summarise(
    condition      = condition[1],
    genotype       = genotype[1],
    trem2_status   = trem2_status[1],
    amyloid_status = amyloid_status[1],
    n_cells        = n(),
    mean_nCount    = round(mean(nCount_RNA, na.rm = TRUE), 1),
    median_nCount  = round(median(nCount_RNA, na.rm = TRUE), 1),
    mean_nFeature  = round(mean(nFeature_RNA, na.rm = TRUE), 1),
    median_nFeature = round(median(nFeature_RNA, na.rm = TRUE), 1),
    mean_pct_mt    = round(mean(percent.mt, na.rm = TRUE), 3),
    median_pct_mt  = round(median(percent.mt, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(condition, sample_id)

write.csv(sample_info, file.path(OUT_DIR, "sample_info.csv"), row.names = FALSE)
message("  sample_info.csv: ", nrow(sample_info), " samples")

# =============================================================
# 2. Normalized counts (global pseudobulk)
# =============================================================
message("Aggregating pseudobulk counts...")

# Sum raw counts per sample across all cells
raw_counts <- AggregateExpression(
  obj,
  group.by = "sample_id",
  slot     = "counts",
  return.seurat = FALSE
)$RNA

raw_counts <- as.matrix(raw_counts)

# Normalize using DESeq2 size factors
message("  Normalizing with DESeq2 size factors...")
coldata <- data.frame(
  sample_id = colnames(raw_counts),
  condition = sample_info$condition[match(colnames(raw_counts), sample_info$sample_id)],
  row.names = colnames(raw_counts)
)

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData   = coldata,
  design    = ~ 1
)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Filter low-expression genes: keep genes with ≥5 counts in ≥2 samples
keep <- rowSums(raw_counts >= 5) >= 2
norm_counts <- norm_counts[keep, ]

norm_df <- as.data.frame(norm_counts) %>%
  rownames_to_column("gene")

write.csv(norm_df, file.path(OUT_DIR, "normalized_counts.csv"), row.names = FALSE)
message("  normalized_counts.csv: ", nrow(norm_df), " genes × ", ncol(norm_df) - 1, " samples")

# =============================================================
# 3. DE results (Q1 global)
# =============================================================
if (file.exists(DE_PATH)) {
  file.copy(DE_PATH, file.path(OUT_DIR, "de_results.csv"), overwrite = TRUE)
  de <- read_csv(DE_PATH, show_col_types = FALSE)
  message("  de_results.csv: ", nrow(de), " genes (Q1 global)")
} else {
  message("  WARNING: DE file not found at ", DE_PATH)
  message("  Run 09_07_7m_global.R first, then re-run this script")
}

# =============================================================
# 4. GSEA results (Q1 broad excitatory)
# =============================================================
if (file.exists(GSEA_PATH)) {
  file.copy(GSEA_PATH, file.path(OUT_DIR, "gsea_results.csv"), overwrite = TRUE)
  gsea <- read_csv(GSEA_PATH, show_col_types = FALSE)
  message("  gsea_results.csv: ", nrow(gsea), " pathways (Q1 excitatory)")
} else {
  message("  WARNING: GSEA file not found at ", GSEA_PATH)
  message("  Run 12_02_gsea_broad.R first, then re-run this script")
}

# =============================================================
# Done
# =============================================================
message("\n=== All files exported to: ", OUT_DIR, " ===")
message("  sample_info.csv         — ", nrow(sample_info), " samples")
message("  normalized_counts.csv   — ", nrow(norm_df), " genes")
if (exists("de"))   message("  de_results.csv          — ", nrow(de), " genes")
if (exists("gsea")) message("  gsea_results.csv        — ", nrow(gsea), " pathways")
message("\nCopy this folder to your Shiny project directory.")