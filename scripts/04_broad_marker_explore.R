# ============================================================
# 04_broad_marker_explore.R
# Broad marker visualisation (dotplots, matrixplots, module
# scores, signature heatmaps) + FindAllMarkers for each cohort.
#
# This script is exploratory — its outputs help you decide the
# cluster-to-broad-type mapping applied in 05_assign_broad_types.R.
#
# Inputs:  ./processed/seurat_7m_clustered.rds
#          ./processed/seurat_15m_clustered.rds
# Outputs: ./processed/seurat_7m_broad_annot_scores.rds
#          ./processed/seurat_15m_broad_annot_scores.rds
#          ./result/annotation_plots/*.png
#          ./result/annotation_tables/*.csv
#          ./result/markers/*_all_markers.csv
#          ./result/markers/*_top10_markers_per_cluster.csv
# ============================================================

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(readr)

source("./scripts/00_utils.R")

dir.create("./result/annotation_plots",  showWarnings = FALSE, recursive = TRUE)
dir.create("./result/annotation_tables", showWarnings = FALSE, recursive = TRUE)
dir.create("./result/markers",           showWarnings = FALSE, recursive = TRUE)

# =============================================================
# A. Marker panel definitions
# =============================================================

# Panel 1: paper-like / simpler
marker_dict1_broad <- list(
  Neuron_general     = c("Grin1", "Syt1", "Rbfox3", "Snap25"),
  Excitatory_neurons = c("Grin2a", "Slc17a7", "Satb2", "Camk2a", "Nptx1"),
  Interneurons       = c("Gad1", "Gad2", "Sst", "Npy", "Tac1", "Penk", "Pde10a"),
  Oligodendrocytes   = c("Plp1", "Mbp", "Cldn11", "Mog"),
  Astrocytes         = c("Slc1a2", "Gja1", "Aqp4"),
  Microglia          = c("Hexb", "Csf1r", "C1qa", "P2ry12", "Trem2", "Tmem119"),
  OPCs               = c("Pdgfra", "Vcan", "Cspg4", "Olig1"),
  Endothelium        = c("Flt1", "Cldn5", "Vtn")
)

# Panel 2: compiled / tier-0 broad panel
marker_dict2_broad <- list(
  Pan_Neuronal    = c("Snap25", "Syt1", "Rbfox3"),
  Pan_Excitatory  = c("Slc17a7", "Neurod6", "Camk2a"),
  Pan_Inhibitory  = c("Gad1", "Gad2", "Slc32a1"),
  Astrocyte       = c("Slc1a2", "Gja1", "Aqp4", "Aldh1l1"),
  Microglia       = c("Hexb", "Csf1r", "C1qa", "P2ry12"),
  Oligodendrocyte = c("Plp1", "Mbp", "Mog", "Cldn11"),
  OPC             = c("Pdgfra", "Cspg4", "Vcan"),
  Endothelial     = c("Cldn5", "Flt1", "Pecam1"),
  Pericyte        = c("Vtn", "Pdgfrb")
)

# =============================================================
# B. Plotting helpers (use filter_marker_list from 00_utils.R)
# =============================================================

plot_marker_dotplot <- function(obj, marker_list, cohort_name, dict_name) {
  genes_use <- unique(unlist(marker_list))
  genes_use <- genes_use[genes_use %in% rownames(obj)]
  
  p <- DotPlot(obj, features = genes_use, group.by = "seurat_clusters") +
    RotatedAxis() +
    ggtitle(paste0(cohort_name, " ", dict_name, " broad marker dotplot")) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))
  
  ggsave(file.path("./result/annotation_plots",
                   paste0(cohort_name, "_", dict_name, "_dotplot.png")),
         plot = p, width = 14, height = 7, dpi = 300)
  p
}

plot_marker_matrixplot <- function(obj, marker_list, cohort_name, dict_name,
                                   group.by = "seurat_clusters") {
  marker_list <- lapply(marker_list, function(x) intersect(x, rownames(obj)))
  marker_list <- marker_list[lengths(marker_list) > 0]
  if (length(marker_list) == 0) return(NULL)
  
  p_grouped <- matrixplot(obj, features = marker_list, group.by = group.by)
  ggsave(file.path("./result/annotation_plots",
                   paste0(cohort_name, "_", dict_name, "_matrixplot_grouped.png")),
         plot = p_grouped, width = 12, height = 8, dpi = 300)
  
  p_flat <- matrixplot(obj, features = unique(unlist(marker_list)),
                       group.by = group.by)
  ggsave(file.path("./result/annotation_plots",
                   paste0(cohort_name, "_", dict_name, "_matrixplot_flat.png")),
         plot = p_flat, width = 12, height = 8, dpi = 300)
  
  list(grouped = p_grouped, flat = p_flat)
}

# --- Module scores ---

add_signature_scores <- function(obj, marker_list, prefix) {
  for (nm in names(marker_list)) {
    genes_use <- marker_list[[nm]]
    if (length(genes_use) < 2) {
      message("Skipping ", prefix, "_", nm, " (<2 genes)")
      next
    }
    score_name <- paste0(prefix, "_", nm)
    obj <- AddModuleScore(obj, features = list(genes_use),
                          name = paste0(score_name, "_tmp"),
                          assay = "RNA", search = FALSE)
    old_col <- paste0(score_name, "_tmp1")
    obj[[score_name]] <- obj[[old_col]]
    obj[[old_col]] <- NULL
  }
  obj
}

# --- Signature summary + plots ---

make_signature_summary <- function(obj, cohort_name, prefix) {
  score_cols <- grep(paste0("^", prefix, "_"),
                     colnames(obj@meta.data), value = TRUE)
  
  summary_tbl <- obj@meta.data %>%
    rownames_to_column("cell") %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_cols), mean), .groups = "drop")
  
  write.csv(summary_tbl,
            file.path("./result/annotation_tables",
                      paste0(cohort_name, "_", prefix,
                             "_cluster_signature_means.csv")),
            row.names = FALSE)
  summary_tbl
}

plot_signature_boxplots <- function(obj, cohort_name, prefix) {
  score_cols <- grep(paste0("^", prefix, "_"),
                     colnames(obj@meta.data), value = TRUE)
  df <- FetchData(obj, vars = c("seurat_clusters", score_cols)) %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = all_of(score_cols),
                 names_to = "signature", values_to = "score")
  
  p <- ggplot(df, aes(x = seurat_clusters, y = score)) +
    geom_boxplot(outlier.size = 0.2) +
    facet_wrap(~ signature, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(paste0(cohort_name, " ", prefix, " signature boxplots"))
  
  ggsave(file.path("./result/annotation_plots",
                   paste0(cohort_name, "_", prefix, "_signature_boxplots.png")),
         plot = p, width = 14, height = 10, dpi = 300)
  p
}

plot_signature_heatmap <- function(summary_tbl, cohort_name, prefix) {
  long_tbl <- summary_tbl %>%
    pivot_longer(cols = -seurat_clusters,
                 names_to = "signature", values_to = "mean_score")
  
  p <- ggplot(long_tbl, aes(x = seurat_clusters, y = signature,
                            fill = mean_score)) +
    geom_tile() + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(paste0(cohort_name, " ", prefix, " signature heatmap"))
  
  ggsave(file.path("./result/annotation_plots",
                   paste0(cohort_name, "_", prefix, "_signature_heatmap.png")),
         plot = p, width = 10, height = 6, dpi = 300)
  p
}

# =============================================================
# C. Main annotation function for one cohort
# =============================================================
annotate_broad_markers <- function(obj, cohort_name) {
  DefaultAssay(obj) <- "RNA"
  
  dict1_use <- filter_marker_list(marker_dict1_broad, obj)
  dict2_use <- filter_marker_list(marker_dict2_broad, obj)
  
  write.csv(stack(dict1_use),
            file.path("./result/annotation_tables",
                      paste0(cohort_name, "_dict1_genes_found.csv")),
            row.names = FALSE)
  write.csv(stack(dict2_use),
            file.path("./result/annotation_tables",
                      paste0(cohort_name, "_dict2_genes_found.csv")),
            row.names = FALSE)
  
  # Dotplots
  plot_marker_dotplot(obj, dict1_use, cohort_name, "dict1")
  plot_marker_dotplot(obj, dict2_use, cohort_name, "dict2")
  
  # Matrixplots
  plot_marker_matrixplot(obj, dict1_use, cohort_name, "dict1")
  plot_marker_matrixplot(obj, dict2_use, cohort_name, "dict2")
  
  # Module scores
  obj <- add_signature_scores(obj, dict1_use, prefix = "d1")
  obj <- add_signature_scores(obj, dict2_use, prefix = "d2")
  
  # Boxplots
  plot_signature_boxplots(obj, cohort_name, "d1")
  plot_signature_boxplots(obj, cohort_name, "d2")
  
  # Summary heatmaps
  d1_sum <- make_signature_summary(obj, cohort_name, "d1")
  d2_sum <- make_signature_summary(obj, cohort_name, "d2")
  plot_signature_heatmap(d1_sum, cohort_name, "d1")
  plot_signature_heatmap(d2_sum, cohort_name, "d2")
  
  obj
}

# =============================================================
# D. FindAllMarkers for one cohort
# =============================================================
find_markers_for_cohort <- function(obj, cohort_name,
                                    only.pos = TRUE,
                                    min.pct = 0.1,
                                    logfc.threshold = 0.1,
                                    test.use = "wilcox") {
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj, assay = "RNA")
  Idents(obj) <- "seurat_clusters"
  
  if (!"data" %in% Layers(obj[["RNA"]])) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  
  message("Layers in ", cohort_name, ": ",
          paste(Layers(obj[["RNA"]]), collapse = ", "))
  message("Cluster sizes:")
  print(table(Idents(obj)))
  
  markers <- FindAllMarkers(obj, assay = "RNA", slot = "data",
                            only.pos = only.pos,
                            min.pct = min.pct,
                            logfc.threshold = logfc.threshold,
                            test.use = test.use,
                            verbose = TRUE)
  
  write_csv(markers,
            file.path("./result/markers",
                      paste0(cohort_name, "_all_markers.csv")))
  
  if (nrow(markers) == 0) {
    warning("No markers found for ", cohort_name)
    return(markers)
  }
  
  top10 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()
  write_csv(top10,
            file.path("./result/markers",
                      paste0(cohort_name, "_top10_markers_per_cluster.csv")))
  markers
}

# =============================================================
# E. Run
# =============================================================
seurat_7m  <- readRDS("./processed/seurat_7m_clustered.rds")
seurat_15m <- readRDS("./processed/seurat_15m_clustered.rds")

# Broad marker visualisation + module scores
seurat_7m_annot  <- annotate_broad_markers(seurat_7m,  "7m")
seurat_15m_annot <- annotate_broad_markers(seurat_15m, "15m")

saveRDS(seurat_7m_annot,  "./processed/seurat_7m_broad_annot_scores.rds")
saveRDS(seurat_15m_annot, "./processed/seurat_15m_broad_annot_scores.rds")

# FindAllMarkers (uses clustered objects, not the annotated ones)
markers_7m  <- find_markers_for_cohort(seurat_7m,  "7m")
markers_15m <- find_markers_for_cohort(seurat_15m, "15m")