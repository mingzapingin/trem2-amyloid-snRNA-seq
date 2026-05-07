# ============================================================
# 06_build_allen_references.R
#
# Build subset reference Seurat objects from the Allen Brain
# Cell Atlas for reference-based label transfer in step 07.
#
# The full Allen dataset is 31,053 genes × 1,169,320 cells.
# We never load it all. Instead:
#
#   1. Read metadata CSV (lightweight)
#   2. Filter to cortex + hippocampus regions
#   3. Map Allen class/subclass → our subgroups
#   4. Downsample each subclass proportionally
#   5. Read only selected cells from HDF5
#   6. Build a Seurat reference object per subgroup
#   7. Normalize + PCA + UMAP each reference
#
# Inputs:
#   ./metadata/reference_metadata.csv
#   ./metadata/reference_dendrogram.RData
#   ./raw_data/reference/expression_matrix.hdf5
#
# Outputs:
#   ./processed/allen_references/allen_<group>_reference.rds
#   ./result/reference_subset/*.png
#   ./result/reference_subset/*.csv
#
# Requires: rhdf5  (BiocManager::install("rhdf5"))
# ============================================================

library(Seurat)
library(SeuratObject)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(rhdf5)
library(Matrix)

set.seed(1234)

dir.create("./processed/allen_references", recursive = TRUE, showWarnings = FALSE)
dir.create("./result/reference_subset",  recursive = TRUE, showWarnings = FALSE)

REF_DIR <- "./processed/allen_references"
SUBSET_DIR  <- "./result/reference_subset"

# =============================================================
# 1. File paths
# =============================================================

METADATA_CSV  <- "./metadata/reference_metadata.csv"
DENDRO_RDATA  <- "./metadata/reference_dendrogram.RData"
EXPRESSION_H5 <- "./raw_data/reference/expression_matrix.hdf5"

# =============================================================
# 2. Read metadata
# =============================================================

message("Reading metadata...")
ref_meta <- read_csv(METADATA_CSV, show_col_types = FALSE)
cat("Full reference: ", nrow(ref_meta), " cells\n\n")

# Quick summary of key columns
cat("=== class_label ===\n")
print(sort(table(ref_meta$class_label), decreasing = TRUE))

cat("\n=== subclass_label ===\n")
print(sort(table(ref_meta$subclass_label), decreasing = TRUE))

cat("\n=== region_label ===\n")
print(sort(table(ref_meta$region_label), decreasing = TRUE))

cat("\n=== neighborhood_label ===\n")
print(sort(table(ref_meta$neighborhood_label), decreasing = TRUE))

# Save full summary for reference
subclass_by_class <- ref_meta %>%
  count(class_label, subclass_label, neighborhood_label, name = "n_cells") %>%
  arrange(class_label, desc(n_cells))
write_csv(subclass_by_class,
          file.path(SUBSET_DIR, "allen_subclass_by_class_full.csv"))

# =============================================================
# 3. Region filter
# =============================================================
# Keep cortex + hippocampus regions relevant to your mouse brain
# snRNA-seq (5XFAD / TREM2 study).
#
# Allen region_label values for cortex + hippocampus include
# areas like: VISp, VISl, SSp, SSs, MOp, MOs, ACA, PL, ILA,
# RSP, TEa-PERI-ECT, AUD, CA, DG, SUB, ProS, HPF, etc.
#
# Strategy: EXCLUDE regions that are clearly NOT cortex or
# hippocampus (e.g., TH = thalamus, MB = midbrain, HY = hypo-
# thalamus, MY = medulla, CB = cerebellum, etc.).
# This is safer than trying to enumerate all cortex regions.

# First, see what's there
cat("\n\nAll unique region_label values:\n")
print(sort(unique(ref_meta$region_label)))

# Regions to EXCLUDE (subcortical, brainstem, etc.)
# Adjust this list after inspecting the output above.
regions_exclude <- c(
  # Add region_label values to exclude here.
  # Leave empty on first run to see what regions exist.
  # Example excludes (uncomment/edit based on your data):
  # "TH", "HY", "MB", "MY", "P", "CB"
)

if (length(regions_exclude) > 0) {
  ref_meta_filtered <- ref_meta %>%
    filter(!region_label %in% regions_exclude)
} else {
  # If no excludes set, keep everything for now
  # (you'll rely on class/subclass filtering below)
  ref_meta_filtered <- ref_meta
}

cat("\nAfter region filter: ", nrow(ref_meta_filtered), " cells\n")

# =============================================================
# 4. Subgroup mapping: Allen class/subclass → our groups
# =============================================================
# Map each Allen subclass_label to the subgroup we want to build
# a reference for. Cells not matching any group are dropped.
#
# This is the most important section to get right.
# Check the subclass_by_class table saved above.

# --- Mapping table ---
# Format: Allen subclass_label → our group name
# Allen class_label is used as a sanity check / fallback.

subgroup_mapping <- tribble(
  ~class_label,     ~subclass_label,      ~our_group,
  
  # --- Excitatory neurons (Glutamatergic) ---
  # Cortical layer types
  "Glutamatergic",  "L2/3 IT CTX",        "excitatory",
  "Glutamatergic",  "L4/5 IT CTX",        "excitatory",
  "Glutamatergic",  "L5 IT CTX",          "excitatory",
  "Glutamatergic",  "L5 PT CTX",          "excitatory",
  "Glutamatergic",  "L5/6 NP CTX",        "excitatory",
  "Glutamatergic",  "L6 IT CTX",          "excitatory",
  "Glutamatergic",  "L6 CT CTX",          "excitatory",
  "Glutamatergic",  "L6b CTX",            "excitatory",
  # Entorhinal / retrosplenial / perirhinal
  "Glutamatergic",  "L6b/CT ENT",         "excitatory",
  "Glutamatergic",  "L3 IT ENT",          "excitatory",
  "Glutamatergic",  "L2/3 IT ENTl",       "excitatory",
  "Glutamatergic",  "L2 IT ENTl",         "excitatory",
  "Glutamatergic",  "L2 IT ENTm",         "excitatory",
  "Glutamatergic",  "L6 IT ENTl",         "excitatory",
  "Glutamatergic",  "L5/6 IT TPE-ENT",    "excitatory",
  "Glutamatergic",  "L4 RSP-ACA",         "excitatory",
  "Glutamatergic",  "L2/3 IT RHP",        "excitatory",
  # Parahippocampal / other
  "Glutamatergic",  "L2/3 IT PPP",        "excitatory",
  "Glutamatergic",  "L5 PPP",             "excitatory",
  "Glutamatergic",  "NP PPP",             "excitatory",
  "Glutamatergic",  "NP SUB",             "excitatory",
  "Glutamatergic",  "Car3",               "excitatory",
  "Glutamatergic",  "CR",                 "excitatory",
  # Hippocampal
  "Glutamatergic",  "CA1-ProS",           "excitatory",
  "Glutamatergic",  "CA2-IG-FC",          "excitatory",
  "Glutamatergic",  "CA3",                "excitatory",
  "Glutamatergic",  "DG",                 "excitatory",
  "Glutamatergic",  "SUB-ProS",           "excitatory",
  "Glutamatergic",  "CT SUB",             "excitatory",
  
  # --- Inhibitory neurons (GABAergic) ---
  "GABAergic",      "Pvalb",              "inhibitory",
  "GABAergic",      "Sst",                "inhibitory",
  "GABAergic",      "Sst Chodl",          "inhibitory",
  "GABAergic",      "Vip",                "inhibitory",
  "GABAergic",      "Lamp5",              "inhibitory",
  "GABAergic",      "Sncg",               "inhibitory",
  "GABAergic",      "Meis2",              "inhibitory",
  
  # --- Non-neuronal ---
  "Non-Neuronal",   "Micro-PVM",          "microglia",
  "Non-Neuronal",   "Astro",              "astrocyte",
  "Non-Neuronal",   "Oligo",              "oligodendrocyte",
  "Non-Neuronal",   "Endo",               "endothelial_pericyte",
  "Non-Neuronal",   "SMC-Peri",           "endothelial_pericyte",
  "Non-Neuronal",   "VLMC",               "endothelial_pericyte"
)

# --- Apply mapping ---
ref_meta_mapped <- ref_meta_filtered %>%
  inner_join(subgroup_mapping,
             by = c("class_label", "subclass_label"))

cat("\nAfter subgroup mapping: ", nrow(ref_meta_mapped), " cells\n")
cat("\nCells per group:\n")
print(ref_meta_mapped %>% count(our_group, sort = TRUE))

cat("\nCells per group × subclass:\n")
group_subclass_counts <- ref_meta_mapped %>%
  count(our_group, subclass_label, sort = TRUE)
print(group_subclass_counts, n = 60)
write_csv(group_subclass_counts,
          file.path(SUBSET_DIR, "allen_cells_per_group_subclass.csv"))

# --- Check for unmapped subclasses ---
unmapped <- ref_meta_filtered %>%
  anti_join(subgroup_mapping,
            by = c("class_label", "subclass_label")) %>%
  count(class_label, subclass_label, sort = TRUE)

if (nrow(unmapped) > 0) {
  cat("\n*** Unmapped subclasses (not assigned to any group): ***\n")
  print(unmapped, n = 40)
  write_csv(unmapped,
            file.path(SUBSET_DIR, "allen_unmapped_subclasses.csv"))
}

# =============================================================
# 5. Downsample per subclass
# =============================================================
# Target: keep proportional representation of each subclass
# within each group, but cap total cells per group.

# Max cells per reference group (adjust based on your RAM)
MAX_CELLS_PER_GROUP <- 50000

# Minimum cells to keep per subclass (even tiny ones get at least this)
MIN_CELLS_PER_SUBCLASS <- 50

downsample_group <- function(meta_group, max_total, min_per_subclass) {
  n_total <- nrow(meta_group)
  if (n_total <= max_total) {
    message("  Group has ", n_total, " cells ≤ max (", max_total,
            ") — keeping all")
    return(meta_group)
  }
  
  subclass_counts <- meta_group %>%
    count(subclass_label, name = "n_orig")
  
  # Proportional allocation
  subclass_counts <- subclass_counts %>%
    mutate(
      prop       = n_orig / sum(n_orig),
      n_target   = pmax(round(prop * max_total), min_per_subclass),
      n_keep     = pmin(n_target, n_orig)
    )
  
  # Adjust if total exceeds max (can happen with min_per_subclass bumps)
  total_keep <- sum(subclass_counts$n_keep)
  if (total_keep > max_total * 1.1) {
    # Scale back proportionally, but keep minimums
    excess <- total_keep - max_total
    can_reduce <- subclass_counts %>%
      mutate(reducible = n_keep - min_per_subclass) %>%
      filter(reducible > 0)
    scale_factor <- 1 - (excess / sum(can_reduce$reducible))
    subclass_counts <- subclass_counts %>%
      mutate(
        reducible = n_keep - min_per_subclass,
        n_keep    = ifelse(reducible > 0,
                           min_per_subclass + round(reducible * scale_factor),
                           n_keep)
      ) %>%
      select(-reducible)
  }
  
  message("  Downsampling ", n_total, " → ~", sum(subclass_counts$n_keep),
          " cells across ", nrow(subclass_counts), " subclasses")
  
  # Sample within each subclass
  sample_plan <- subclass_counts %>% select(subclass_label, n_keep)
  
  sampled_parts <- list()
  for (i in seq_len(nrow(sample_plan))) {
    sc   <- sample_plan$subclass_label[i]
    n_k  <- sample_plan$n_keep[i]
    rows <- meta_group %>% filter(subclass_label == sc)
    sampled_parts[[i]] <- slice_sample(rows, n = n_k)
  }
  
  bind_rows(sampled_parts)
}

# Apply downsampling per group
groups_to_build <- unique(ref_meta_mapped$our_group)
downsampled_list <- list()

for (grp in groups_to_build) {
  message("\nDownsampling group: ", grp)
  grp_meta <- ref_meta_mapped %>% filter(our_group == grp)
  downsampled_list[[grp]] <- downsample_group(
    grp_meta, MAX_CELLS_PER_GROUP, MIN_CELLS_PER_SUBCLASS
  )
}

# Summary
downsample_summary <- bind_rows(lapply(names(downsampled_list), function(grp) {
  downsampled_list[[grp]] %>%
    count(our_group, subclass_label, name = "n_downsampled")
}))
write_csv(downsample_summary,
          file.path(SUBSET_DIR, "allen_downsample_summary.csv"))

cat("\n\nDownsample summary (cells per group):\n")
print(downsample_summary %>%
        group_by(our_group) %>%
        summarise(n_cells = sum(n_downsampled),
                  n_subclasses = n(), .groups = "drop"))

# =============================================================
# 6. Read gene names + cell IDs from HDF5
# =============================================================

message("\n\nReading gene names and cell ID index from HDF5...")

h5_genes   <- h5read(EXPRESSION_H5, "/data/gene")
h5_samples <- h5read(EXPRESSION_H5, "/data/samples")

cat("HDF5: ", length(h5_genes), " genes × ",
    length(h5_samples), " cells\n")

# Quick sanity check: do metadata cell IDs match HDF5 cell IDs?
all_meta_cells <- unique(unlist(lapply(downsampled_list, function(x) x$sample_name)))
n_found <- sum(all_meta_cells %in% h5_samples)
cat("Cell ID check: ", n_found, " / ", length(all_meta_cells),
    " metadata cells found in HDF5\n")
if (n_found == 0) {
  cat("  First 3 metadata IDs: ",
      paste(head(all_meta_cells, 3), collapse = ", "), "\n")
  cat("  First 3 HDF5 IDs:    ",
      paste(head(h5_samples, 3), collapse = ", "), "\n")
  stop("No cell IDs match — check if metadata sample_name matches HDF5 /data/samples")
}

# Build cell_id → column index lookup
# (HDF5 is 1-indexed in R via rhdf5)
sample_to_idx <- setNames(seq_along(h5_samples), h5_samples)

# =============================================================
# 7. Build Seurat reference for each group
# =============================================================

build_reference <- function(grp_name, grp_meta, h5_path,
                            gene_names, sample_to_idx) {
  message("\n========== Building reference: ", grp_name, " ==========")
  
  cell_ids <- grp_meta$sample_name
  n_cells  <- length(cell_ids)
  
  # Find column indices in HDF5
  col_idx <- sample_to_idx[cell_ids]
  missing  <- sum(is.na(col_idx))
  if (missing > 0) {
    warning("  ", missing, " cell IDs not found in HDF5 — dropping")
    keep     <- !is.na(col_idx)
    cell_ids <- cell_ids[keep]
    col_idx  <- col_idx[keep]
    grp_meta <- grp_meta[keep, ]
    n_cells  <- length(cell_ids)
  }
  
  cat("  Reading ", n_cells, " cells from HDF5...\n")
  cat("  Index range: ", min(col_idx), " to ", max(col_idx),
      " (HDF5 dim 2 = ", length(sample_to_idx), ")\n")
  
  # --- Read expression in chunks to avoid memory spikes ---
  # HDF5 layout: /data/counts is (31053 genes × 1169320 cells)
  # BUT rhdf5 reverses dimension order, so in R:
  #   dim[1] = 1169320 (cells)
  #   dim[2] = 31053   (genes)
  # So index = list(cell_indices, NULL) to select cells × all genes,
  # then transpose to get genes × cells for Seurat.
  
  # Sort column indices for efficient HDF5 access
  # unname() is critical — rhdf5 misinterprets named vectors
  col_idx <- unname(col_idx)
  sort_order <- order(col_idx)
  col_idx_sorted <- col_idx[sort_order]
  
  # Read in chunks of ~5000 cells to keep memory manageable
  CHUNK_SIZE <- 5000
  n_chunks <- ceiling(n_cells / CHUNK_SIZE)
  count_chunks <- list()
  
  for (i in seq_len(n_chunks)) {
    start_i <- (i - 1) * CHUNK_SIZE + 1
    end_i   <- min(i * CHUNK_SIZE, n_cells)
    chunk_idx <- col_idx_sorted[start_i:end_i]
    
    if (i %% 5 == 1 || i == n_chunks) {
      message("  chunk ", i, "/", n_chunks,
              " (cells ", start_i, "-", end_i, ")")
    }
    
    # cells in dim 1, genes in dim 2 (rhdf5 reverses HDF5 dims)
    chunk_data <- h5read(
      EXPRESSION_H5, "/data/counts",
      index = list(as.integer(chunk_idx), NULL)
    )
    
    # Transpose: (cells × genes) → (genes × cells), then sparse
    count_chunks[[i]] <- Matrix(t(chunk_data), sparse = TRUE)
    rm(chunk_data); gc(verbose = FALSE)
  }
  
  # Combine chunks
  counts_sparse <- do.call(cbind, count_chunks)
  rm(count_chunks); gc(verbose = FALSE)
  
  # Restore original cell order (undo sort)
  unsort_order <- order(sort_order)
  counts_sparse <- counts_sparse[, unsort_order]
  
  # Set dimnames
  # Allen reference can have duplicate gene symbols — make unique
  gene_names_unique <- make.unique(gene_names, sep = ".")
  n_dup <- sum(duplicated(gene_names))
  if (n_dup > 0) {
    message("  Note: ", n_dup, " duplicate gene names made unique")
  }
  rownames(counts_sparse) <- gene_names_unique
  colnames(counts_sparse) <- cell_ids
  
  cat("  Expression matrix: ", nrow(counts_sparse), " genes × ",
      ncol(counts_sparse), " cells\n")
  
  # --- Create Seurat object ---
  obj <- CreateSeuratObject(
    counts   = counts_sparse,
    project  = paste0("allen_", grp_name),
    min.cells    = 0,
    min.features = 0
  )
  
  # Attach Allen metadata
  meta_cols <- c(
    "class_label", "subclass_label", "cluster_label",
    "neighborhood_label", "region_label",
    "cell_type_alias_label", "cell_type_designation_label",
    "our_group"
  )
  meta_cols_present <- intersect(meta_cols, colnames(grp_meta))
  
  md_to_add <- grp_meta %>%
    select(sample_name, all_of(meta_cols_present)) %>%
    column_to_rownames("sample_name")
  
  # Match order to Seurat cells
  md_to_add <- md_to_add[colnames(obj), , drop = FALSE]
  for (col in colnames(md_to_add)) {
    obj[[col]] <- md_to_add[[col]]
  }
  
  rm(counts_sparse); gc(verbose = FALSE)
  
  # --- Standard Seurat processing ---
  message("  Normalizing + PCA + UMAP...")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  
  # --- Save ---
  out_path <- file.path(REF_DIR,
                        paste0("allen_", grp_name, "_reference.rds"))
  saveRDS(obj, out_path)
  message("  Saved: ", out_path)
  
  # --- QC plots ---
  p_sub <- DimPlot(obj, group.by = "subclass_label",
                   label = TRUE, repel = TRUE) +
    ggtitle(paste0("Allen ", grp_name, " reference — subclass")) +
    theme(legend.position = "right")
  ggsave(file.path(SUBSET_DIR,
                   paste0("allen_", grp_name, "_umap_subclass.png")),
         plot = p_sub, width = 10, height = 7, dpi = 200)
  
  if ("cluster_label" %in% colnames(obj@meta.data)) {
    p_cl <- DimPlot(obj, group.by = "cluster_label",
                    label = FALSE) +
      ggtitle(paste0("Allen ", grp_name, " reference — cluster")) +
      NoLegend()
    ggsave(file.path(SUBSET_DIR,
                     paste0("allen_", grp_name, "_umap_cluster.png")),
           plot = p_cl, width = 8, height = 7, dpi = 200)
  }
  
  p_reg <- DimPlot(obj, group.by = "region_label",
                   label = FALSE) +
    ggtitle(paste0("Allen ", grp_name, " reference — region"))
  ggsave(file.path(SUBSET_DIR,
                   paste0("allen_", grp_name, "_umap_region.png")),
         plot = p_reg, width = 10, height = 7, dpi = 200)
  
  # Composition table
  comp <- obj@meta.data %>%
    count(subclass_label, cluster_label, name = "n_cells") %>%
    arrange(subclass_label, desc(n_cells))
  write_csv(comp,
            file.path(SUBSET_DIR,
                      paste0("allen_", grp_name, "_subclass_cluster_counts.csv")))
  
  cat("  Done: ", ncol(obj), " cells, ",
      n_distinct(obj$subclass_label), " subclasses\n")
  
  obj
}

# =============================================================
# 8. Run all groups
# =============================================================

# Which groups to build? (comment out any you want to skip)
groups_to_build <- c(
  "excitatory",
  "inhibitory",
  "microglia",
  "astrocyte",
  "oligodendrocyte",
  "endothelial_pericyte"
)

allen_refs <- list()

for (grp in groups_to_build) {
  if (!grp %in% names(downsampled_list)) {
    warning("No cells for group '", grp, "' — skipping")
    next
  }
  
  allen_refs[[grp]] <- build_reference(
    grp_name      = grp,
    grp_meta      = downsampled_list[[grp]],
    h5_path       = EXPRESSION_H5,
    gene_names    = h5_genes,
    sample_to_idx = sample_to_idx
  )
  
  # Free memory between groups
  gc(verbose = FALSE)
}

# =============================================================
# 9. Summary
# =============================================================

cat("\n\n==================== SUMMARY ====================\n\n")

for (grp in names(allen_refs)) {
  obj <- allen_refs[[grp]]
  cat(grp, ": ",
      ncol(obj), " cells, ",
      nrow(obj), " genes, ",
      n_distinct(obj$subclass_label), " subclasses\n")
}

cat("\nReference objects saved in: ", REF_DIR, "\n")
cat("QC outputs saved in:       ", SUBSET_DIR, "\n")

cat("\n\nNext step: run 07_transfer_labels.R to map your query\n")
cat("cells onto these references using Seurat's TransferData.\n")

# =============================================================
# 10. Load dendrogram (for later use in transfer / plotting)
# =============================================================

# The dendrogram encodes the Allen taxonomy hierarchy.
# It's not needed for building references but will be useful
# for ordering labels and plotting in 07.

dendro_env <- new.env()
load(DENDRO_RDATA, envir = dendro_env)
cat("\nDendrogram objects loaded: ",
    paste(ls(dendro_env), collapse = ", "), "\n")

# Save the dendrogram alongside references for easy access
saveRDS(
  as.list(dendro_env),
  file.path(REF_DIR, "allen_dendrogram.rds")
)