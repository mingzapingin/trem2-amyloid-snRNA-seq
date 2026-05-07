# ============================================================
# 03_normalize_cluster.R
# Join layers → normalize → HVGs → scale → PCA → UMAP → cluster.
# Runs separately for 7m and 15m.
#
# Inputs:  ./processed/seurat_7m_qc.rds
#          ./processed/seurat_15m_qc.rds
# Outputs: ./processed/seurat_7m_clustered.rds
#          ./processed/seurat_15m_clustered.rds
#          ./result/cluster_plots/*.png
#          ./result/cluster_tables/*.csv
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)
dir.create("./result/cluster_tables", showWarnings = FALSE, recursive = TRUE)
dir.create("./result/cluster_plots",  showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Load QC-filtered objects
# -----------------------------
seurat_7m  <- readRDS("./processed/seurat_7m_qc.rds")
seurat_15m <- readRDS("./processed/seurat_15m_qc.rds")

# -----------------------------
# 2. Processing function
# -----------------------------
process_cohort <- function(obj, cohort_name,
                           nfeatures = 3000, npcs = 30,
                           dims_use = 1:20, resolution = 0.5,
                           regress_mt = FALSE) {
  message("Processing cohort: ", cohort_name)
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj, assay = "RNA")

  # Normalize → HVGs → Scale → PCA

  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                               nfeatures = nfeatures, verbose = FALSE)

  write.csv(
    data.frame(variable_gene = VariableFeatures(obj)),
    file.path("./result/cluster_tables",
              paste0(cohort_name, "_variable_features.csv")),
    row.names = FALSE
  )

  if (regress_mt) {
    obj <- ScaleData(obj, features = VariableFeatures(obj),
                     vars.to.regress = "percent.mt", verbose = FALSE)
  } else {
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  }

  obj <- RunPCA(obj, features = VariableFeatures(obj),
                npcs = npcs, verbose = FALSE)

  # --- Elbow + PCA plots ---
  p_elbow <- ElbowPlot(obj, ndims = npcs) +
    ggtitle(paste0(cohort_name, " elbow plot"))
  ggsave(file.path("./result/cluster_plots",
                    paste0(cohort_name, "_elbow.png")),
         plot = p_elbow, width = 6, height = 4, dpi = 300)

  p_pca_s <- DimPlot(obj, reduction = "pca", group.by = "sample_id") +
    ggtitle(paste0(cohort_name, " PCA by sample"))
  p_pca_c <- DimPlot(obj, reduction = "pca", group.by = "condition") +
    ggtitle(paste0(cohort_name, " PCA by condition"))
  ggsave(file.path("./result/cluster_plots",
                    paste0(cohort_name, "_pca_sample_condition.png")),
         plot = p_pca_s + p_pca_c, width = 12, height = 5, dpi = 300)

  # --- Neighbors → Cluster → UMAP ---
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)

  # --- UMAP plots ---
  p_cl <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                  label = TRUE) +
    ggtitle(paste0(cohort_name, " UMAP clusters"))
  ggsave(file.path("./result/cluster_plots",
                    paste0(cohort_name, "_umap_clusters.png")),
         plot = p_cl, width = 6, height = 5, dpi = 300)

  p_us <- DimPlot(obj, reduction = "umap", group.by = "sample_id") +
    ggtitle(paste0(cohort_name, " UMAP by sample"))
  p_uc <- DimPlot(obj, reduction = "umap", group.by = "condition") +
    ggtitle(paste0(cohort_name, " UMAP by condition"))
  ggsave(file.path("./result/cluster_plots",
                    paste0(cohort_name, "_umap_sample_condition.png")),
         plot = p_us + p_uc, width = 12, height = 5, dpi = 300)

  p_split <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                     split.by = "sample_id", label = FALSE, ncol = 4) +
    ggtitle(paste0(cohort_name, " UMAP split by sample"))
  ggsave(file.path("./result/cluster_plots",
                    paste0(cohort_name, "_umap_split_by_sample.png")),
         plot = p_split, width = 16, height = 10, dpi = 300)

  # --- Cluster composition tables ---
  cluster_tab <- table(obj$seurat_clusters, obj$sample_id)
  write.csv(as.data.frame.matrix(cluster_tab),
            file.path("./result/cluster_tables",
                      paste0(cohort_name, "_cluster_by_sample.csv")))
  write.csv(as.data.frame.matrix(prop.table(cluster_tab, margin = 2)),
            file.path("./result/cluster_tables",
                      paste0(cohort_name, "_cluster_prop_by_sample.csv")))

  # --- Per-sample UMAP + counts ---
  for (sid in unique(obj$sample_id)) {
    obj_sub <- subset(obj, subset = sample_id == sid)

    p_s <- DimPlot(obj_sub, reduction = "umap",
                   group.by = "seurat_clusters", label = TRUE) +
      ggtitle(paste0(cohort_name, " - ", sid, " clusters"))
    ggsave(file.path("./result/cluster_plots",
                      paste0(cohort_name, "_", sid, "_umap_clusters.png")),
           plot = p_s, width = 6, height = 5, dpi = 300)

    sample_counts <- obj_sub@meta.data %>%
      dplyr::count(seurat_clusters, name = "n_cells") %>%
      dplyr::mutate(sample_id  = sid,
                    proportion = n_cells / sum(n_cells)) %>%
      dplyr::arrange(as.numeric(as.character(seurat_clusters)))
    write.csv(sample_counts,
              file.path("./result/cluster_tables",
                        paste0(cohort_name, "_", sid, "_cluster_counts.csv")),
              row.names = FALSE)
  }

  # --- Save ---
  saveRDS(obj, file.path("./processed",
                          paste0("seurat_", cohort_name, "_clustered.rds")))
  obj
}

# -----------------------------
# 3. Run both cohorts
# -----------------------------
# Tuned settings from grid search (see project notes):
#   7m:  dims 1:12, resolution 0.1
#   15m: dims 1:12, resolution 0.1

seurat_7m_cl  <- process_cohort(seurat_7m,  "7m",
                                 npcs = 30, dims_use = 1:12, resolution = 0.1)
seurat_15m_cl <- process_cohort(seurat_15m, "15m",
                                 npcs = 40, dims_use = 1:12, resolution = 0.1)

cat("\n7m clusters:\n");  print(table(seurat_7m_cl$seurat_clusters))
cat("\n15m clusters:\n"); print(table(seurat_15m_cl$seurat_clusters))
