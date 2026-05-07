# ============================================================
# 07_transfer_labels.R
#
# For each cell type group, subset query cells from the
# broad-typed objects, find transfer anchors against the
# matching Allen reference built in 06, and transfer
# subclass + cluster labels onto query cells.
#
# Inputs:
#   ./processed/seurat_7m_broad_typed.rds
#   ./processed/seurat_15m_broad_typed.rds
#   ./processed/allen_references/allen_<group>_reference.rds
#
# Outputs:
#   ./processed/transfer_objects/seurat_<cohort>_<group>_transferred.rds
#   ./result/transfer_labels/<cohort>_<group>_*.png
#   ./result/transfer_labels/<cohort>_<group>_*.csv
# ============================================================

library(Seurat)
library(SeuratObject)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)

set.seed(1234)

# MapQuery uses future/parallel — default 500MB limit is too small
# for large reference objects. Set to 50 GB.
options(future.globals.maxSize = 50 * 1024^3)

dir.create("./processed/transfer_objects", recursive = TRUE, showWarnings = FALSE)
dir.create("./result/transfer_labels",     recursive = TRUE, showWarnings = FALSE)

RDS_DIR  <- "./processed/transfer_objects"
PLOT_DIR <- "./result/transfer_labels"
REF_DIR  <- "./processed/allen_references"

# =============================================================
# 1. Subgroup config: broad_type → reference file
# =============================================================

subgroup_config <- list(
  excitatory = list(
    broad_types = c("Excitatory_neuron"),
    ref_file    = file.path(REF_DIR, "allen_excitatory_reference.rds")
  ),
  inhibitory = list(
    broad_types = c("Inhibitory_neuron"),
    ref_file    = file.path(REF_DIR, "allen_inhibitory_reference.rds")
  ),
  microglia = list(
    broad_types = c("Microglia"),
    ref_file    = file.path(REF_DIR, "allen_microglia_reference.rds")
  ),
  astrocyte = list(
    broad_types = c("Astrocyte"),
    ref_file    = file.path(REF_DIR, "allen_astrocyte_reference.rds")
  ),
  oligodendrocyte = list(
    broad_types = c("Oligodendrocyte", "OPC"),
    ref_file    = file.path(REF_DIR, "allen_oligodendrocyte_reference.rds")
  ),
  endothelial_pericyte = list(
    broad_types = c("Endothelial_Pericyte"),
    ref_file    = file.path(REF_DIR, "allen_endothelial_pericyte_reference.rds")
  )
)

# =============================================================
# 2. Transfer settings (adjust if needed)
# =============================================================

TRANSFER_DIMS   <- 1:30   # PCA dims for anchor finding
N_FEATURES      <- 3000   # variable features for integration
K_ANCHOR        <- 5      # number of anchors per cell
K_WEIGHT        <- 50     # neighbors for weight transfer

# Labels to transfer from the Allen reference
# These must be column names in the reference metadata
LABELS_TO_TRANSFER <- c("subclass_label", "cluster_label")

# =============================================================
# 3. Cohort registry
# =============================================================

cohort_inputs <- list(
  "7m"  = "./processed/seurat_7m_broad_typed.rds",
  "15m" = "./processed/seurat_15m_broad_typed.rds"
)

# Which cohorts and groups to run?
cohorts_to_run <- c("7m", "15m")
groups_to_run  <- names(subgroup_config)
# groups_to_run <- c("microglia", "astrocyte")  # or just a subset

# =============================================================
# 4. Core transfer function
# =============================================================

transfer_labels_one <- function(query_obj, ref_obj,
                                cohort_name, group_name,
                                labels_to_transfer,
                                transfer_dims, n_features,
                                k_anchor, k_weight) {
  
  tag <- paste0(cohort_name, "_", group_name)
  message("\n--- ", tag, ": ", ncol(query_obj), " query cells, ",
          ncol(ref_obj), " ref cells ---")
  
  # --- Prepare query: normalize + find variable features ---
  DefaultAssay(query_obj) <- "RNA"
  query_obj <- tryCatch(
    JoinLayers(query_obj, assay = "RNA"),
    error = function(e) query_obj
  )
  query_obj <- NormalizeData(query_obj, verbose = FALSE)
  query_obj <- FindVariableFeatures(query_obj, nfeatures = n_features,
                                    verbose = FALSE)
  
  # --- Ensure reference is normalized ---
  DefaultAssay(ref_obj) <- "RNA"
  if (!"data" %in% Layers(ref_obj[["RNA"]])) {
    ref_obj <- NormalizeData(ref_obj, verbose = FALSE)
  }
  if (length(VariableFeatures(ref_obj)) == 0) {
    ref_obj <- FindVariableFeatures(ref_obj, nfeatures = n_features,
                                    verbose = FALSE)
  }
  
  # --- Gene name diagnostic ---
  query_genes <- rownames(query_obj)
  ref_genes   <- rownames(ref_obj)
  shared_genes <- intersect(query_genes, ref_genes)
  
  cat("  Query genes:  ", length(query_genes), "  (examples: ",
      paste(head(query_genes, 5), collapse = ", "), ")\n")
  cat("  Ref genes:    ", length(ref_genes), "  (examples: ",
      paste(head(ref_genes, 5), collapse = ", "), ")\n")
  cat("  Shared genes: ", length(shared_genes), "\n")
  
  if (length(shared_genes) < 100) {
    # Try case-insensitive match to diagnose casing issue
    query_lower <- tolower(query_genes)
    ref_lower   <- tolower(ref_genes)
    n_case_match <- length(intersect(query_lower, ref_lower))
    cat("  Case-insensitive overlap: ", n_case_match, "\n")
    
    if (n_case_match > 1000) {
      message("  Gene name CASING mismatch detected — converting reference ",
              "genes to match query format (first-letter-uppercase)")
      # Mouse gene convention: first letter uppercase, rest lowercase (e.g. Snap25)
      # Allen may use all-caps (SNAP25) or all-lowercase (snap25)
      new_ref_genes <- tools::toTitleCase(tolower(ref_genes))
      # But mouse genes are actually: first letter upper, rest lower
      # toTitleCase capitalizes each word — we want only first char upper
      new_ref_genes <- paste0(toupper(substr(ref_genes, 1, 1)),
                              tolower(substr(ref_genes, 2, nchar(ref_genes))))
      # Handle make.unique suffixes (e.g., "Gm12345.1")
      new_ref_genes <- make.unique(new_ref_genes, sep = ".")
      
      rownames(ref_obj[["RNA"]]) <- new_ref_genes
      ref_obj <- FindVariableFeatures(ref_obj, nfeatures = n_features,
                                      verbose = FALSE)
      ref_obj <- ScaleData(ref_obj, verbose = FALSE)
      ref_obj <- RunPCA(ref_obj, npcs = max(transfer_dims), verbose = FALSE)
      
      shared_genes <- intersect(rownames(query_obj), new_ref_genes)
      cat("  After case fix — shared genes: ", length(shared_genes), "\n")
    }
    
    if (length(shared_genes) < 100) {
      stop("Only ", length(shared_genes), " shared genes between query and ",
           "reference. Check gene name format.")
    }
  }
  
  # --- Select features: intersection of variable features ---
  query_var <- VariableFeatures(query_obj)
  ref_var   <- VariableFeatures(ref_obj)
  shared_var <- intersect(query_var, ref_var)
  
  if (length(shared_var) < 200) {
    # Fall back: use top shared genes by variance in both
    message("  Few shared variable features (", length(shared_var),
            ") — using top shared genes instead")
    shared_var <- shared_genes[seq_len(min(n_features, length(shared_genes)))]
  }
  
  cat("  Using ", length(shared_var), " features for anchor finding\n")
  
  # --- Find transfer anchors ---
  message("  Finding transfer anchors...")
  anchors <- FindTransferAnchors(
    reference  = ref_obj,
    query      = query_obj,
    dims       = transfer_dims,
    features   = shared_var,
    k.anchor   = k_anchor,
    normalization.method = "LogNormalize",
    reference.reduction  = "pca",
    verbose    = TRUE
  )
  
  message("  Found ", nrow(anchors@anchors), " anchors")
  
  # --- Transfer each label ---
  for (label_col in labels_to_transfer) {
    if (!label_col %in% colnames(ref_obj@meta.data)) {
      warning("  Label '", label_col, "' not in reference — skipping")
      next
    }
    
    message("  Transferring: ", label_col)
    
    predictions <- TransferData(
      anchorset = anchors,
      refdata   = ref_obj@meta.data[[label_col]],
      dims      = transfer_dims,
      k.weight  = k_weight,
      verbose   = FALSE
    )
    
    # Add predictions to query
    # TransferData returns: predicted.id, prediction.score.*, prediction.score.max
    pred_col     <- paste0("allen_", label_col)
    score_col    <- paste0("allen_", label_col, "_score")
    query_obj[[pred_col]]  <- predictions$predicted.id
    query_obj[[score_col]] <- predictions$prediction.score.max
  }
  
  # --- Map reference UMAP onto query ---
  message("  Mapping query onto reference UMAP...")
  ref_obj <- RunUMAP(ref_obj, dims = transfer_dims, return.model = TRUE,
                     verbose = FALSE)
  
  query_obj <- MapQuery(
    anchorset = anchors,
    reference = ref_obj,
    query     = query_obj,
    refdata   = list(subclass = "subclass_label"),
    reference.reduction = "pca",
    reduction.model     = "umap",
    verbose = FALSE
  )
  
  # --- Save transferred object ---
  out_path <- file.path(RDS_DIR, paste0("seurat_", tag, "_transferred.rds"))
  saveRDS(query_obj, out_path)
  message("  Saved: ", out_path)
  
  # --- Plots ---
  # UMAP colored by predicted subclass
  if ("allen_subclass_label" %in% colnames(query_obj@meta.data)) {
    p_sub <- DimPlot(query_obj, group.by = "allen_subclass_label",
                     reduction = "ref.umap",
                     label = TRUE, repel = TRUE) +
      ggtitle(paste(tag, "— predicted subclass")) +
      theme(legend.position = "right")
    ggsave(file.path(PLOT_DIR, paste0(tag, "_umap_predicted_subclass.png")),
           plot = p_sub, width = 10, height = 7, dpi = 200)
    
    # Prediction confidence
    p_score <- FeaturePlot(query_obj, features = "allen_subclass_label_score",
                           reduction = "ref.umap") +
      scale_color_viridis_c() +
      ggtitle(paste(tag, "— subclass prediction score"))
    ggsave(file.path(PLOT_DIR, paste0(tag, "_umap_subclass_score.png")),
           plot = p_score, width = 8, height = 7, dpi = 200)
  }
  
  # UMAP colored by predicted cluster
  if ("allen_cluster_label" %in% colnames(query_obj@meta.data)) {
    p_cl <- DimPlot(query_obj, group.by = "allen_cluster_label",
                    reduction = "ref.umap",
                    label = FALSE) +
      ggtitle(paste(tag, "— predicted cluster")) +
      NoLegend()
    ggsave(file.path(PLOT_DIR, paste0(tag, "_umap_predicted_cluster.png")),
           plot = p_cl, width = 8, height = 7, dpi = 200)
  }
  
  # UMAP by condition
  if ("condition" %in% colnames(query_obj@meta.data)) {
    p_cond <- DimPlot(query_obj, group.by = "condition",
                      reduction = "ref.umap", shuffle = TRUE) +
      ggtitle(paste(tag, "— condition"))
    ggsave(file.path(PLOT_DIR, paste0(tag, "_umap_condition.png")),
           plot = p_cond, width = 8, height = 7, dpi = 200)
  }
  
  # UMAP by sample
  p_samp <- DimPlot(query_obj, group.by = "orig.ident",
                    reduction = "ref.umap", shuffle = TRUE) +
    ggtitle(paste(tag, "— sample"))
  ggsave(file.path(PLOT_DIR, paste0(tag, "_umap_sample.png")),
         plot = p_samp, width = 8, height = 7, dpi = 200)
  
  # --- Tables ---
  # Prediction summary
  if ("allen_subclass_label" %in% colnames(query_obj@meta.data)) {
    sub_summary <- query_obj@meta.data %>%
      count(allen_subclass_label, name = "n_cells") %>%
      mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
      arrange(desc(n_cells))
    write_csv(sub_summary,
              file.path(PLOT_DIR, paste0(tag, "_predicted_subclass_counts.csv")))
    
    # Confidence summary per predicted subclass
    score_summary <- query_obj@meta.data %>%
      group_by(allen_subclass_label) %>%
      summarise(
        n_cells      = n(),
        mean_score   = round(mean(allen_subclass_label_score), 3),
        median_score = round(median(allen_subclass_label_score), 3),
        min_score    = round(min(allen_subclass_label_score), 3),
        pct_above_0.5 = round(100 * mean(allen_subclass_label_score > 0.5), 1),
        .groups = "drop"
      ) %>%
      arrange(desc(n_cells))
    write_csv(score_summary,
              file.path(PLOT_DIR, paste0(tag, "_subclass_confidence_summary.csv")))
  }
  
  if ("allen_cluster_label" %in% colnames(query_obj@meta.data)) {
    cl_summary <- query_obj@meta.data %>%
      count(allen_cluster_label, name = "n_cells") %>%
      mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
      arrange(desc(n_cells))
    write_csv(cl_summary,
              file.path(PLOT_DIR, paste0(tag, "_predicted_cluster_counts.csv")))
  }
  
  # Subclass × condition cross-tab (useful for DE planning)
  if ("allen_subclass_label" %in% colnames(query_obj@meta.data) &&
      "condition" %in% colnames(query_obj@meta.data)) {
    cross_tab <- query_obj@meta.data %>%
      count(allen_subclass_label, condition, name = "n_cells") %>%
      arrange(allen_subclass_label, condition)
    write_csv(cross_tab,
              file.path(PLOT_DIR, paste0(tag, "_subclass_by_condition.csv")))
  }
  
  # Subclass × sample cross-tab
  if ("allen_subclass_label" %in% colnames(query_obj@meta.data)) {
    samp_tab <- query_obj@meta.data %>%
      count(allen_subclass_label, orig.ident, name = "n_cells") %>%
      arrange(allen_subclass_label, orig.ident)
    write_csv(samp_tab,
              file.path(PLOT_DIR, paste0(tag, "_subclass_by_sample.csv")))
  }
  
  query_obj
}

# =============================================================
# 5. Run all cohorts × all groups
# =============================================================

all_results <- list()

for (cohort_name in cohorts_to_run) {
  rds_path <- cohort_inputs[[cohort_name]]
  message("\n\n========== Loading ", cohort_name, " ==========")
  obj_full <- readRDS(rds_path)
  cat(cohort_name, " broad_type distribution:\n")
  print(table(obj_full$broad_type))
  
  cohort_results <- list()
  
  for (group_name in groups_to_run) {
    cfg <- subgroup_config[[group_name]]
    
    # --- Check reference exists ---
    if (!file.exists(cfg$ref_file)) {
      warning("Reference not found: ", cfg$ref_file, " — skipping ", group_name)
      next
    }
    
    # --- Subset query cells ---
    cells_use <- rownames(obj_full@meta.data)[
      obj_full$broad_type %in% cfg$broad_types
    ]
    n_cells <- length(cells_use)
    message("\n>>> ", cohort_name, " — ", group_name,
            ": ", n_cells, " query cells")
    
    if (n_cells < 30) {
      warning("  Too few cells (", n_cells, ") — skipping")
      next
    }
    
    query_obj <- subset(obj_full, cells = cells_use)
    
    # --- Load reference ---
    message("  Loading reference: ", cfg$ref_file)
    ref_obj <- readRDS(cfg$ref_file)
    cat("  Reference: ", ncol(ref_obj), " cells, ",
        n_distinct(ref_obj$subclass_label), " subclasses\n")
    
    # --- Transfer ---
    result <- tryCatch({
      transfer_labels_one(
        query_obj  = query_obj,
        ref_obj    = ref_obj,
        cohort_name = cohort_name,
        group_name  = group_name,
        labels_to_transfer = LABELS_TO_TRANSFER,
        transfer_dims = TRANSFER_DIMS,
        n_features    = N_FEATURES,
        k_anchor      = K_ANCHOR,
        k_weight      = K_WEIGHT
      )
    }, error = function(e) {
      message("  !! Transfer FAILED for ", cohort_name, " ", group_name,
              ": ", e$message)
      NULL
    })
    
    if (!is.null(result)) {
      cohort_results[[group_name]] <- result
    }
    
    # Free reference memory
    rm(ref_obj); gc(verbose = FALSE)
  }
  
  all_results[[cohort_name]] <- cohort_results
}

# =============================================================
# 6. Summary
# =============================================================

cat("\n\n==================== SUMMARY ====================\n\n")

for (cohort_name in names(all_results)) {
  res <- all_results[[cohort_name]]
  cat(cohort_name, ":\n")
  
  if (length(res) == 0) {
    cat("  No transfers completed\n")
    next
  }
  
  for (group_name in names(res)) {
    obj <- res[[group_name]]
    n <- ncol(obj)
    
    if ("allen_subclass_label" %in% colnames(obj@meta.data)) {
      n_sub <- n_distinct(obj$allen_subclass_label)
      med_score <- round(median(obj$allen_subclass_label_score), 3)
      cat("  ", group_name, " — ", n, " cells → ",
          n_sub, " predicted subclasses",
          " (median score: ", med_score, ")\n")
    } else {
      cat("  ", group_name, " — ", n, " cells (no subclass label)\n")
    }
  }
}

cat("\nTransferred objects saved in: ", RDS_DIR, "\n")
cat("Plots and tables saved in:   ", PLOT_DIR, "\n")
