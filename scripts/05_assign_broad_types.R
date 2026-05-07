# ============================================================
# 05_assign_broad_types.R
# Apply the cluster → broad cell-type label mapping determined
# from the marker exploration in step 04.
#
# Edit the cluster_to_label_* vectors below when you update
# your assignments.
#
# Inputs:  ./processed/seurat_7m_clustered.rds
#          ./processed/seurat_15m_clustered.rds
# Outputs: ./processed/seurat_7m_broad_typed.rds
#          ./processed/seurat_15m_broad_typed.rds
#          ./result/annotation_tables/*_broad_annotation_table.csv
#          ./result/annotation_plots/*_umap_broad_type.png
# ============================================================

library(Seurat)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)

dir.create("./result/annotation_tables", showWarnings = FALSE, recursive = TRUE)
dir.create("./result/annotation_plots",  showWarnings = FALSE, recursive = TRUE)

# =============================================================
# 1. Cluster-to-label mappings   *** EDIT THESE ***
# =============================================================

cluster_to_label_7m <- c(
  "0"  = "Excitatory_neuron",
  "1"  = "Inhibitory_neuron",
  "2"  = "Oligodendrocyte",
  "3"  = "Inhibitory_neuron",
  "4"  = "Excitatory_neuron",
  "5"  = "Astrocyte",
  "6"  = "Microglia",
  "7"  = "OPC",
  "8"  = "Endothelial_Pericyte",
  "9"  = "Unknown",
  "10" = "Unknown",
  "11" = "Unknown"
)

cluster_to_label_15m <- c(
  "0"  = "Excitatory_neuron",
  "1"  = "Inhibitory_neuron",
  "2"  = "Excitatory_neuron",
  "3"  = "Excitatory_neuron",
  "4"  = "Astrocyte",
  "5"  = "Oligodendrocyte",
  "6"  = "Endothelial_Pericyte",
  "7"  = "OPC",
  "8"  = "Microglia"
)

# Canonical ordering for the factor levels
# "Unknown" is last — these clusters are parked here until resolved.
# Downstream scripts (06, 06.1) ignore Unknown cells automatically
# because they subset by specific broad_type values.
broad_type_levels <- c(
  "Excitatory_neuron", "Inhibitory_neuron",
  "Astrocyte", "Oligodendrocyte", "OPC",
  "Microglia", "Endothelial_Pericyte",
  "Unknown"
)

# =============================================================
# 2. Assignment function
# =============================================================

assign_broad_type <- function(obj, cluster_to_label) {
  clusters <- as.character(obj$seurat_clusters)
  
  missing <- setdiff(unique(clusters), names(cluster_to_label))
  if (length(missing) > 0) {
    stop("Missing broad labels for clusters: ",
         paste(missing, collapse = ", "))
  }
  
  obj$original_cluster <- obj$seurat_clusters
  obj$broad_type       <- factor(
    unname(cluster_to_label[clusters]),
    levels = broad_type_levels
  )
  obj
}

# =============================================================
# 3. Apply to both cohorts
# =============================================================

seurat_7m  <- readRDS("./processed/seurat_7m_clustered.rds")
seurat_15m <- readRDS("./processed/seurat_15m_clustered.rds")

seurat_7m  <- assign_broad_type(seurat_7m,  cluster_to_label_7m)
seurat_15m <- assign_broad_type(seurat_15m, cluster_to_label_15m)

# =============================================================
# 4. Save annotation tables
# =============================================================

for (info in list(
  list(obj = seurat_7m,  name = "7m"),
  list(obj = seurat_15m, name = "15m")
)) {
  write_csv(
    info$obj@meta.data %>%
      rownames_to_column("cell_id") %>%
      select(cell_id, sample_id, original_cluster, broad_type),
    file.path("./result/annotation_tables",
              paste0(info$name, "_broad_annotation_table.csv"))
  )
  
  # UMAP coloured by broad type
  p <- DimPlot(info$obj, group.by = "broad_type", label = TRUE, repel = TRUE) +
    ggtitle(paste0(info$name, " broad cell types"))
  ggsave(file.path("./result/annotation_plots",
                   paste0(info$name, "_umap_broad_type.png")),
         plot = p, width = 8, height = 6, dpi = 300)
  
  # Summary counts
  cat("\n", info$name, " broad type counts:\n")
  print(table(info$obj$broad_type))
  
  n_unk <- sum(info$obj$broad_type == "Unknown")
  if (n_unk > 0) {
    cat("  *** ", n_unk, " cells labelled Unknown (clusters: ",
        paste(unique(info$obj$original_cluster[info$obj$broad_type == "Unknown"]),
              collapse = ", "),
        ") — these will be excluded from subtype reclustering in 06 ***\n")
  }
}

# =============================================================
# 5. Save annotated objects
# =============================================================

saveRDS(seurat_7m,  "./processed/seurat_7m_broad_typed.rds")
saveRDS(seurat_15m, "./processed/seurat_15m_broad_typed.rds")