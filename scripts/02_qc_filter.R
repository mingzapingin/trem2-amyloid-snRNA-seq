# ============================================================
# 02_qc_filter.R
# Compute QC metrics, visualise, filter cells, save.
#
# Inputs:  ./processed/seurat_7m_raw.rds
#          ./processed/seurat_15m_raw.rds
# Outputs: ./processed/seurat_7m_qc.rds
#          ./processed/seurat_15m_qc.rds
#          ./result/qc_summary_before_filter.csv
#          ./result/qc_summary_after_filter.csv
#          ./result/qc_compare_before_after.csv
#          ./result/filter_plot/*.png
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

dir.create("./result/filter_plot", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Load raw objects
# -----------------------------
seurat_7m  <- readRDS("./processed/seurat_7m_raw.rds")
seurat_15m <- readRDS("./processed/seurat_15m_raw.rds")

# -----------------------------
# 2. Add QC metrics
# -----------------------------
seurat_7m[["percent.mt"]]   <- PercentageFeatureSet(seurat_7m,  pattern = "^mt-")
seurat_15m[["percent.mt"]]  <- PercentageFeatureSet(seurat_15m, pattern = "^mt-")
seurat_7m[["percent.ribo"]] <- PercentageFeatureSet(seurat_7m,  pattern = "^Rpl|^Rps")
seurat_15m[["percent.ribo"]]<- PercentageFeatureSet(seurat_15m, pattern = "^Rpl|^Rps")

# -----------------------------
# 3. QC summary helper
# -----------------------------
qc_summary <- function(obj, stage_label) {
  obj@meta.data %>%
    group_by(sample_id) %>%
    summarise(
      cohort           = first(cohort),
      age_months       = first(age_months),
      condition        = first(condition),
      region           = first(region),
      n_cells          = n(),
      median_nCount    = median(nCount_RNA),
      median_nFeature  = median(nFeature_RNA),
      median_percent_mt = median(percent.mt),
      mean_percent_mt  = mean(percent.mt),
      .groups = "drop"
    ) %>%
    mutate(qc_stage = stage_label)
}

qc_before <- bind_rows(
  qc_summary(seurat_7m,  "before_filter_7m"),
  qc_summary(seurat_15m, "before_filter_15m")
)
write_csv(qc_before, "./result/qc_summary_before_filter.csv")

# -----------------------------
# 4. QC violin + scatter plots
# -----------------------------
plot_qc <- function(obj, prefix) {
  p1 <- VlnPlot(obj,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "sample_id", ncol = 3, pt.size = 0.1) + NoLegend()
  ggsave(file.path("./result/filter_plot", paste0(prefix, "_violin.png")),
         plot = p1, width = 14, height = 5, dpi = 300)

  p2 <- FeatureScatter(obj, "nCount_RNA", "percent.mt") +
    ggtitle(paste0(prefix, ": nCount vs percent.mt"))
  p3 <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA") +
    ggtitle(paste0(prefix, ": nCount vs nFeature"))
  ggsave(file.path("./result/filter_plot", paste0(prefix, "_scatter.png")),
         plot = p2 + p3, width = 12, height = 5, dpi = 300)
}

plot_qc(seurat_7m,  "7m_before")
plot_qc(seurat_15m, "15m_before")

# -----------------------------
# 5. Filter cells
# -----------------------------
# 7m:  percent.mt < 5
# 15m: percent.mt < 7.5  (relaxed — see project notes)

seurat_7m_filt <- subset(seurat_7m,
                          subset = nFeature_RNA > 200 &
                                   nFeature_RNA < 6000 &
                                   nCount_RNA   < 25000 &
                                   percent.mt   < 5)

seurat_15m_filt <- subset(seurat_15m,
                           subset = nFeature_RNA > 200 &
                                    nFeature_RNA < 6000 &
                                    nCount_RNA   < 25000 &
                                    percent.mt   < 7.5)

# -----------------------------
# 6. Post-filter summaries
# -----------------------------
qc_after <- bind_rows(
  qc_summary(seurat_7m_filt,  "after_filter_7m"),
  qc_summary(seurat_15m_filt, "after_filter_15m")
)
write_csv(qc_after, "./result/qc_summary_after_filter.csv")

qc_compare <- qc_before %>%
  select(sample_id, n_cells_before = n_cells, cohort, condition, region) %>%
  left_join(qc_after %>% select(sample_id, n_cells_after = n_cells),
            by = "sample_id") %>%
  mutate(cells_removed = n_cells_before - n_cells_after,
         pct_removed   = round(100 * cells_removed / n_cells_before, 2))

write_csv(qc_compare, "./result/qc_compare_before_after.csv")

plot_qc(seurat_7m_filt,  "7m_after")
plot_qc(seurat_15m_filt, "15m_after")

# -----------------------------
# 7. Save
# -----------------------------
saveRDS(seurat_7m_filt,  "./processed/seurat_7m_qc.rds")
saveRDS(seurat_15m_filt, "./processed/seurat_15m_qc.rds")

cat("\nCells kept by sample:\n")
print(qc_compare)
