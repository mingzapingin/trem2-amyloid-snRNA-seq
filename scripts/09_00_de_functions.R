# ============================================================
# 09_00_de_functions.R
#
# Shared setup for all pseudobulk DE scripts (09-01 to 09-06).
# Source this at the top of every DE script:
#   source("./scripts/09_00_de_functions.R")
#
# Contains:
#   - Library loading
#   - Global settings (thresholds, output dirs)
#   - Object loading functions
#   - Pseudobulk aggregation function
#   - DESeq2 wrapper
#   - Results formatting + saving
#   - DE-ready subtype lists for 7m and 15m
#
# Change thresholds HERE — all downstream scripts inherit them.
# ============================================================

# =============================================================
# Libraries
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
})

# Fix namespace conflicts: DESeq2 loads AnnotationDbi which masks dplyr
select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
count   <- dplyr::count
slice   <- dplyr::slice
desc    <- dplyr::desc
first   <- dplyr::first
collapse <- dplyr::collapse

set.seed(1234)

# =============================================================
# Global settings — EDIT THESE AS NEEDED
# =============================================================

# DE significance thresholds (GLOBAL DEFAULTS)
# These are fallbacks — can be overridden per-comparison or per-subtype.
PADJ_CUTOFF  <- 0.05
LFC_CUTOFF   <- 0.25    # log2 fold change

# =============================================================
# Color palettes — used by all plots
# =============================================================

# Volcano plot: UP / DOWN / NS
VOLCANO_COLORS <- c(
  "UP"   = "#E64B35",   # vermillion red
  "DOWN" = "#4DBBD5",   # cyan blue
  "NS"   = "#CCCCCC"    # light grey
)

# Cross-comparison categories (used in script 10)
CATEGORY_COLORS <- c(
  "Cat1"  = "#F39B7F",  # salmon — pure amyloid
  "Cat2"  = "#E64B35",  # red — Trem2-dependent
  "Cat3"  = "#4DBBD5",  # blue — Trem2-independent
  "Cat4"  = "#F28E2B",  # orange — exacerbated
  "Cat5"  = "#B07AA1",  # purple — Trem2-autonomous
  "Cat6"  = "#59A14F",  # green — redirected
  "NS"    = "#E5E5E5"   # very light grey
)

# Regional categories (15m Q4a-Q4d cross-comparison)
REGIONAL_COLORS <- c(
  "Reg1"  = "#636363",  # dark grey — constitutive regional
  "Reg2"  = "#E64B35",  # red — amyloid-emergent regional
  "Reg3"  = "#4DBBD5",  # blue — amyloid-lost regional
  "Reg4"  = "#F28E2B",  # orange — Trem2-dependent regional
  "Reg5"  = "#B07AA1",  # purple — Trem2-emergent regional
  "Reg6"  = "#59A14F",  # green — compound (amyloid + KO) regional
  "NS"    = "#E5E5E5"
)

# Broad type palette (for cell-type heatmaps / summaries)
BROAD_TYPE_COLORS <- c(
  "Excitatory_neuron"    = "#3182BD",  # blue
  "Inhibitory_neuron"    = "#E6550D",  # orange
  "Astrocyte"            = "#31A354",  # green
  "Oligodendrocyte"      = "#756BB1",  # purple
  "OPC"                  = "#9E9AC8",  # light purple
  "Microglia"            = "#E41A1C",  # red
  "Endothelial_Pericyte" = "#636363",  # grey
  "Unknown"              = "#BDBDBD"   # light grey
)

# Subtype palette — tab20-inspired, 25 distinct colors
# Ordered: excitatory shades → inhibitory shades → glia → vascular
SUBTYPE_COLORS <- c(
  # Excitatory (blues / teals)
  "L2/3 IT CTX"  = "#1F77B4",
  "L4/5 IT CTX"  = "#AEC7E8",
  "L5 IT CTX"    = "#3B8ABF",
  "L5 PT CTX"    = "#6BAED6",
  "L5/6 NP CTX"  = "#9ECAE1",
  "L6 IT CTX"    = "#2171B5",
  "L6 CT CTX"    = "#08519C",
  "L6b CTX"      = "#08306B",
  "L2/3 IT PPP"  = "#4292C6",
  "L2/3 IT ENTl" = "#C6DBEF",
  "L6b/CT ENT"   = "#2C6DAA",
  "L2 IT ENTl"   = "#5BA3CF",
  "L3 IT ENT"    = "#85BCDB",
  "Car3"         = "#17BECF",
  # Inhibitory (oranges / reds)
  "Pvalb"        = "#FF7F0E",
  "Sst"          = "#D62728",
  "Meis2"        = "#E377C2",
  "Lamp5"        = "#FFBB78",
  "Vip"          = "#FF9896",
  "Sncg"         = "#F7B6D2",
  # Glia
  "Micro-PVM"    = "#E41A1C",
  "Astro"        = "#2CA02C",
  "Oligo"        = "#9467BD",
  # Vascular
  "VLMC"         = "#8C564B",
  "Endo"         = "#7F7F7F"
)

# =============================================================
# Threshold resolution: subtype override → comparison override → global
# =============================================================
# Usage in child scripts:
#
#   # Per-comparison override (all subtypes in this comparison):
#   COMPARISONS$Q1$padj_cutoff <- 0.01
#   COMPARISONS$Q1$lfc_cutoff  <- 0.5
#
#   # Per-subtype override (finest control):
#   COMPARISONS$Q1$subtype_thresholds <- list(
#     "Micro-PVM" = list(padj = 0.01, lfc = 0.5),
#     "Pvalb"     = list(padj = 0.1,  lfc = 0.1)
#   )
#
# Resolution order: subtype_thresholds > comparison-level > global

resolve_thresholds <- function(subtype = NULL, comp = NULL) {
  padj <- PADJ_CUTOFF
  lfc  <- LFC_CUTOFF
  
  # Comparison-level override
  if (!is.null(comp$padj_cutoff)) padj <- comp$padj_cutoff
  if (!is.null(comp$lfc_cutoff))  lfc  <- comp$lfc_cutoff
  
  # Subtype-level override (most specific wins)
  if (!is.null(subtype) && !is.null(comp$subtype_thresholds)) {
    st <- comp$subtype_thresholds[[subtype]]
    if (!is.null(st$padj)) padj <- st$padj
    if (!is.null(st$lfc))  lfc  <- st$lfc
  }
  
  list(padj = padj, lfc = lfc)
}

# Pseudobulk filters
MIN_CELLS_PER_PSEUDOBULK <- 10   # min cells in a sample×subtype to include
MIN_GENES_DETECTED       <- 200  # min genes with nonzero counts in pseudobulk

# DESeq2 settings
DESEQ2_TEST     <- "Wald"
DESEQ2_FIT_TYPE <- "local"   # "parametric", "local", or "mean"

# =============================================================
# Output directories
# =============================================================

DE_BASE_DIR <- "./result/de_results"

dir.create(DE_BASE_DIR, recursive = TRUE, showWarnings = FALSE)

make_de_dir <- function(comparison_name) {
  d <- file.path(DE_BASE_DIR, comparison_name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# =============================================================
# Object paths
# =============================================================

ANNOTATED_PATHS <- list(
  "7m"  = "./processed/seurat_7m_final_annotated.rds",
  "15m" = "./processed/seurat_15m_final_annotated.rds"
)

# =============================================================
# 1. Load annotated object (cached to avoid reloading)
# =============================================================

.de_cache <- new.env(parent = emptyenv())

load_annotated <- function(cohort_name) {
  if (!exists(cohort_name, envir = .de_cache)) {
    path <- ANNOTATED_PATHS[[cohort_name]]
    if (!file.exists(path)) stop("File not found: ", path)
    message("Loading ", cohort_name, " from ", path)
    obj <- readRDS(path)
    assign(cohort_name, obj, envir = .de_cache)
  }
  get(cohort_name, envir = .de_cache)
}

# =============================================================
# 2. Subset cells for a comparison
# =============================================================

#' Subset to cells matching a cell type label + set of conditions
#'
#' @param obj Seurat object (final annotated)
#' @param cell_type Label to match (subtype name or broad type name)
#' @param conditions Character vector of condition values to keep
#' @param filter_col Metadata column to filter on
#'   "allen_subclass_filtered" for subtype, "broad_type" for broad
#' @param min_cells Minimum cells per sample to include that sample
#' @return Named list: cells, meta, sample_table, valid, reason
subset_for_de <- function(obj, cell_type, conditions,
                          filter_col = "allen_subclass_filtered",
                          min_cells = MIN_CELLS_PER_PSEUDOBULK) {
  
  md <- obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    filter(
      .data[[filter_col]] == cell_type,
      condition %in% conditions
    )
  
  if (nrow(md) == 0) {
    return(list(cells = character(0), meta = md,
                sample_table = NULL, valid = FALSE,
                reason = "no cells"))
  }
  
  # Count per sample
  sample_counts <- md %>%
    count(condition, sample_id, name = "n_cells")
  
  # Drop samples with too few cells
  adequate <- sample_counts %>% filter(n_cells >= min_cells)
  
  if (nrow(adequate) == 0) {
    return(list(cells = character(0), meta = md,
                sample_table = sample_counts, valid = FALSE,
                reason = "no samples with enough cells"))
  }
  
  # Check: at least 2 samples per condition
  cond_coverage <- adequate %>%
    count(condition, name = "n_samples")
  
  missing_conds <- setdiff(conditions, cond_coverage$condition)
  low_conds     <- cond_coverage %>% filter(n_samples < 2)
  
  if (length(missing_conds) > 0 || nrow(low_conds) > 0) {
    reason <- paste0(
      if (length(missing_conds) > 0)
        paste0("missing conditions: ", paste(missing_conds, collapse = ", "))
      else "",
      if (nrow(low_conds) > 0)
        paste0("; <2 samples in: ",
               paste(low_conds$condition, collapse = ", "))
      else ""
    )
    return(list(cells = character(0), meta = md,
                sample_table = sample_counts, valid = FALSE,
                reason = reason))
  }
  
  # Keep only cells from adequate samples
  keep_samples <- adequate$sample_id
  md_filtered  <- md %>% filter(sample_id %in% keep_samples)
  
  list(
    cells        = md_filtered$cell_id,
    meta         = md_filtered,
    sample_table = sample_counts,
    valid        = TRUE,
    reason       = "OK"
  )
}

# =============================================================
# 3. Pseudobulk aggregation
# =============================================================

#' Aggregate single-cell counts to pseudobulk (one column per sample)
#'
#' @param obj Seurat object
#' @param cell_ids Cell IDs to include
#' @param group_by Column to group by (default: "sample_id")
#' @return List: counts (genes × samples matrix), 
#'         coldata (sample-level metadata)
make_pseudobulk <- function(obj, cell_ids, group_by = "sample_id") {
  
  # Get raw counts
  DefaultAssay(obj) <- "RNA"
  counts_mat <- GetAssayData(obj, layer = "counts")
  
  # Subset to cells
  counts_mat <- counts_mat[, cell_ids, drop = FALSE]
  
  # Get grouping vector
  groups <- obj@meta.data[cell_ids, group_by]
  group_levels <- sort(unique(groups))
  
  # Aggregate: sum counts per group
  pb_list <- lapply(group_levels, function(g) {
    cell_idx <- which(groups == g)
    if (length(cell_idx) == 1) {
      counts_mat[, cell_idx, drop = TRUE]
    } else {
      Matrix::rowSums(counts_mat[, cell_idx, drop = FALSE])
    }
  })
  
  pb_mat <- do.call(cbind, pb_list)
  colnames(pb_mat) <- group_levels
  
  # Convert to regular matrix for DESeq2
  if (inherits(pb_mat, "sparseMatrix")) {
    pb_mat <- as.matrix(pb_mat)
  }
  
  # Filter genes: keep genes detected in at least half the samples
  genes_detected <- rowSums(pb_mat > 0)
  min_samples    <- max(2, floor(ncol(pb_mat) / 2))
  keep_genes     <- genes_detected >= min_samples
  pb_mat         <- pb_mat[keep_genes, , drop = FALSE]
  
  # Build coldata from one representative cell per sample
  md <- obj@meta.data[cell_ids, ]
  coldata <- md %>%
    rownames_to_column("cell_id") %>%
    group_by(!!sym(group_by)) %>%
    summarise(
      condition      = condition[1],
      sample_id      = sample_id[1],
      genotype       = genotype[1],
      trem2_status   = trem2_status[1],
      amyloid_status = amyloid_status[1],
      n_cells        = n(),
      .groups = "drop"
    ) %>%
    column_to_rownames(group_by)
  
  # Add region if present (for 15m)
  if ("region" %in% colnames(md)) {
    region_info <- md %>%
      rownames_to_column("cell_id") %>%
      group_by(!!sym(group_by)) %>%
      summarise(region = region[1], .groups = "drop") %>%
      column_to_rownames(group_by)
    coldata$region <- region_info[rownames(coldata), "region"]
  }
  
  # Ensure column order matches
  coldata <- coldata[colnames(pb_mat), , drop = FALSE]
  
  list(counts = pb_mat, coldata = coldata)
}

# =============================================================
# 4. DESeq2 wrapper
# =============================================================

#' Run DESeq2 on pseudobulk data
#'
#' @param pb_counts Pseudobulk count matrix (genes × samples)
#' @param pb_coldata Sample-level metadata (must have 'condition' column)
#' @param contrast Character vector of length 3: c("condition", "test", "ref")
#' @param design Formula (default: ~ condition)
#' @param test DESeq2 test type
#' @param fit_type DESeq2 fit type
#' @return DESeq2 results as a data.frame, or NULL if it fails
run_deseq2 <- function(pb_counts, pb_coldata, contrast,
                       design      = ~ condition,
                       test        = DESEQ2_TEST,
                       fit_type    = DESEQ2_FIT_TYPE,
                       padj_cutoff = PADJ_CUTOFF,
                       lfc_cutoff  = LFC_CUTOFF) {
  
  # Ensure condition is a factor with reference level
  ref_level <- contrast[3]
  pb_coldata$condition <- factor(pb_coldata$condition)
  pb_coldata$condition <- relevel(pb_coldata$condition, ref = ref_level)
  
  # Create DESeq2 dataset
  dds <- tryCatch({
    DESeqDataSetFromMatrix(
      countData = pb_counts,
      colData   = pb_coldata,
      design    = design
    )
  }, error = function(e) {
    message("    DESeqDataSet creation failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(dds)) return(NULL)
  
  # Filter lowly expressed genes
  keep <- rowSums(counts(dds) >= 5) >= 2
  dds  <- dds[keep, ]
  
  if (nrow(dds) < 50) {
    message("    Too few genes after filtering (", nrow(dds), ")")
    return(NULL)
  }
  
  # Run DESeq2
  dds <- tryCatch({
    DESeq(dds, test = test, fitType = fit_type, quiet = TRUE)
  }, error = function(e) {
    message("    DESeq2 failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(dds)) return(NULL)
  
  # Extract results
  res <- tryCatch({
    results(dds, contrast = contrast, alpha = padj_cutoff)
  }, error = function(e) {
    message("    Results extraction failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(res)) return(NULL)
  
  # Convert to data.frame
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(
      significant = !is.na(padj) & padj < padj_cutoff &
        abs(log2FoldChange) >= lfc_cutoff,
      direction   = case_when(
        !significant          ~ "NS",
        log2FoldChange > 0    ~ "UP",
        log2FoldChange < 0    ~ "DOWN",
        TRUE                  ~ "NS"
      )
    )
  
  res_df
}

# =============================================================
# 5. Run DE for one subtype
# =============================================================

#' Full pipeline for one subtype: subset → pseudobulk → DESeq2 → save
#'
#' @param obj Seurat object
#' @param subtype Allen subclass label
#' @param conditions 2-element vector: c("test_condition", "ref_condition")
#' @param comparison_name For file naming (e.g., "Q1_WT_vs_WT5XFAD")
#' @param design Formula (default: ~ condition)
#' @param out_dir Output directory
#' @return List with results df, summary stats, or NULL
run_de_one_subtype <- function(obj, subtype, conditions,
                               comparison_name, design = ~ condition,
                               filter_col = "allen_subclass_filtered",
                               comp = NULL, out_dir = NULL) {
  
  if (is.null(out_dir)) out_dir <- make_de_dir(comparison_name)
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
  
  # --- Resolve thresholds for this subtype ---
  thresh <- resolve_thresholds(subtype, comp)
  
  test_cond <- conditions[1]
  ref_cond  <- conditions[2]
  contrast  <- c("condition", test_cond, ref_cond)
  
  message("  ", subtype, ": ", test_cond, " vs ", ref_cond,
          "  [padj<", thresh$padj, ", |lfc|≥", thresh$lfc, "]")
  
  # --- Subset ---
  sub_info <- subset_for_de(obj, subtype, conditions, filter_col = filter_col)
  
  if (!sub_info$valid) {
    message("    SKIP — ", sub_info$reason)
    return(list(
      subtype    = subtype,
      comparison = comparison_name,
      status     = "SKIPPED",
      reason     = sub_info$reason,
      results    = NULL,
      summary    = NULL
    ))
  }
  
  n_cells <- length(sub_info$cells)
  message("    ", n_cells, " cells across ",
          n_distinct(sub_info$meta$sample_id), " samples")
  
  # --- Pseudobulk ---
  pb <- make_pseudobulk(obj, sub_info$cells, group_by = "sample_id")
  
  message("    Pseudobulk: ", nrow(pb$counts), " genes × ",
          ncol(pb$counts), " samples")
  message("    Samples per condition: ",
          paste(names(table(pb$coldata$condition)),
                table(pb$coldata$condition),
                sep = "=", collapse = ", "))
  
  # --- DESeq2 ---
  res_df <- run_deseq2(pb$counts, pb$coldata, contrast,
                       design = design,
                       padj_cutoff = thresh$padj,
                       lfc_cutoff  = thresh$lfc)
  
  if (is.null(res_df)) {
    return(list(
      subtype    = subtype,
      comparison = comparison_name,
      status     = "FAILED",
      reason     = "DESeq2 error",
      results    = NULL,
      summary    = NULL
    ))
  }
  
  # --- Summary stats ---
  n_up   <- sum(res_df$direction == "UP",   na.rm = TRUE)
  n_down <- sum(res_df$direction == "DOWN", na.rm = TRUE)
  n_sig  <- n_up + n_down
  n_tested <- sum(!is.na(res_df$padj))
  
  summary_row <- data.frame(
    subtype      = subtype,
    comparison   = comparison_name,
    test_cond    = test_cond,
    ref_cond     = ref_cond,
    n_cells      = n_cells,
    n_samples    = ncol(pb$counts),
    n_genes_tested = n_tested,
    n_sig        = n_sig,
    n_up         = n_up,
    n_down       = n_down,
    padj_cutoff  = thresh$padj,
    lfc_cutoff   = thresh$lfc,
    stringsAsFactors = FALSE
  )
  
  message("    Results: ", n_sig, " DE genes (",
          n_up, " up, ", n_down, " down) / ",
          n_tested, " tested")
  
  # --- Save per-subtype results ---
  write_csv(
    res_df,
    file.path(out_dir,
              paste0(comparison_name, "_", safe_name, "_full_results.csv"))
  )
  
  # Save significant only
  sig_df <- res_df %>% filter(significant)
  if (nrow(sig_df) > 0) {
    write_csv(
      sig_df,
      file.path(out_dir,
                paste0(comparison_name, "_", safe_name, "_sig_only.csv"))
    )
  }
  
  list(
    subtype    = subtype,
    comparison = comparison_name,
    status     = "OK",
    reason     = NA_character_,
    results    = res_df,
    summary    = summary_row
  )
}

# =============================================================
# 6. Run DE across all subtypes for one comparison
# =============================================================

#' Loop over subtypes, run DE, combine summaries
#'
#' @param obj Seurat object
#' @param subtypes Character vector of Allen subclass labels
#' @param conditions c("test", "reference")
#' @param comparison_name For file naming
#' @param design Formula
#' @return List of all results + combined summary table
run_de_comparison <- function(obj, subtypes, conditions,
                              comparison_name, design = ~ condition,
                              filter_col = "allen_subclass_filtered",
                              comp = NULL) {
  
  out_dir <- make_de_dir(comparison_name)
  all_results <- list()
  all_summaries <- list()
  all_skipped <- list()
  
  for (subtype in subtypes) {
    res <- run_de_one_subtype(
      obj, subtype, conditions, comparison_name,
      design = design, filter_col = filter_col,
      comp = comp, out_dir = out_dir
    )
    
    if (res$status == "OK") {
      all_results[[subtype]]  <- res$results
      all_summaries[[subtype]] <- res$summary
    } else {
      all_skipped[[subtype]] <- data.frame(
        subtype    = subtype,
        comparison = comparison_name,
        status     = res$status,
        reason     = res$reason,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # --- Combined summary ---
  summary_df <- bind_rows(all_summaries)
  if (nrow(summary_df) > 0) {
    summary_df <- summary_df %>% arrange(desc(n_sig))
    write_csv(summary_df,
              file.path(out_dir,
                        paste0(comparison_name, "_summary.csv")))
    
    cat("\n\n=== ", comparison_name, " Summary ===\n")
    print(as.data.frame(summary_df), row.names = FALSE)
  }
  
  # Skipped subtypes
  skipped_df <- bind_rows(all_skipped)
  if (nrow(skipped_df) > 0) {
    write_csv(skipped_df,
              file.path(out_dir,
                        paste0(comparison_name, "_skipped.csv")))
    cat("\nSkipped subtypes:\n")
    print(as.data.frame(skipped_df), row.names = FALSE)
  }
  
  # --- Volcano plots for each subtype ---
  for (subtype in names(all_results)) {
    res_df <- all_results[[subtype]]
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    thresh <- resolve_thresholds(subtype, comp)
    
    p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),
                            color = direction)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = VOLCANO_COLORS) +
      geom_hline(yintercept = -log10(thresh$padj),
                 linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-thresh$lfc, thresh$lfc),
                 linetype = "dashed", color = "grey40") +
      theme_bw() +
      ggtitle(paste0(comparison_name, " — ", subtype)) +
      labs(subtitle = paste0("padj<", thresh$padj, ", |lfc|≥", thresh$lfc)) +
      theme(legend.position = "none",
            panel.background = element_rect(fill = "white"),
            plot.background  = element_rect(fill = "white"))
    
    ggsave(
      file.path(out_dir,
                paste0(comparison_name, "_", safe_name, "_volcano.png")),
      plot = p, width = 7, height = 6, dpi = 200
    )
  }
  
  list(
    results  = all_results,
    summary  = summary_df,
    skipped  = skipped_df,
    out_dir  = out_dir
  )
}

# =============================================================
# 7. DE-ready subtypes
# =============================================================
# Filled from 08 output: 7m_de_readiness.csv, 15m_de_readiness.csv
# Only subtypes with de_ready=TRUE (7m) or enough total cells (15m)
# are included. Downstream functions will still skip per-comparison
# if a specific sample lacks cells.

# --- 7m: pseudobulk DESeq2 (3 replicates per condition) ---
# de_ready=TRUE means all 4 conditions have ≥2 samples with ≥10 cells

DE_READY_7m_EXCITATORY <- c(
  "L2/3 IT CTX",     # 6234 cells
  "L4/5 IT CTX",     # 6028
  "L6 CT CTX",       # 4628
  "L6 IT CTX",       # 1639
  "L5 IT CTX",       # 1616
  "L5 PT CTX",       # 1581
  "L2/3 IT PPP",     #  997
  "L5/6 NP CTX",     #  824
  "Car3",            #  700
  "L6b CTX",         #  570
  "L2/3 IT ENTl",    #  443
  "L6b/CT ENT",      #  424
  "L2 IT ENTl",      #  401  (SOME_SAMPLES_LOW but de_ready=TRUE)
  "L3 IT ENT"        #  252  (SOME_SAMPLES_LOW but de_ready=TRUE)
)

DE_READY_7m_INHIBITORY <- c(
  "Pvalb",           # 2473
  "Sst",             # 1593
  "Meis2",           # 1100
  "Lamp5",           # 1008
  "Vip",             # 1002
  "Sncg"             #  322
)

DE_READY_7m_NEURONAL <- c(DE_READY_7m_EXCITATORY, DE_READY_7m_INHIBITORY)

DE_READY_7m_GLIAL <- c(
  "Micro-PVM",       # 3991
  "Astro",           # 4021
  "Oligo"            # 16180
)

DE_READY_7m_VASCULAR <- c(
  "VLMC",            #  557
  "Endo"             #  317
  # SMC-Peri (196) — de_ready=FALSE, TOO_FEW_SAMPLES in 1 condition
)

DE_READY_7m_NONNEURONAL <- c(DE_READY_7m_GLIAL, DE_READY_7m_VASCULAR)

DE_READY_7m_ALL <- c(
  DE_READY_7m_NEURONAL,
  DE_READY_7m_GLIAL,
  DE_READY_7m_VASCULAR
)

# --- 15m: cell-level FindMarkers (no replicates per comparison) ---
# For cell-level DE, we include any subtype with enough total cells
# that individual samples likely have ≥30 cells. The run_cellevel
# functions will skip per-comparison if a sample is too sparse.
# Being inclusive here — the function handles filtering.

DE_READY_15m_ALL <- c(
  # Excitatory (only large subtypes — most are too sparse per sample)
  "L2/3 IT CTX",     # 1303
  "DG",              # 1148
  "L6 CT CTX",       # 1119
  "L4/5 IT CTX",     #  861
  "L2/3 IT PPP",     #  458
  "CR",              #  286
  
  # Inhibitory
  "Pvalb",           #  854
  "Sst",             #  659
  "Vip",             #  523
  "Lamp5",           #  466
  
  # Non-neuronal
  "Micro-PVM",       #  713
  "Astro",           # 2493
  "Oligo",           # 2373
  "VLMC"             #  591
)

# --- Broad cell type DE ---
# Uses broad_type column instead of allen_subclass_filtered.
# All cells of a type are pooled together — maximum power,
# gives the "headline" result per cell class.

BROAD_TYPES_7m <- c(
  "Excitatory_neuron",     # 35878
  "Inhibitory_neuron",     # 24186
  "Microglia",             #  3991
  "Astrocyte",             #  4021
  "Oligodendrocyte",       # 13628  (includes OPC cluster)
  "Endothelial_Pericyte"   #  1488
)

BROAD_TYPES_15m <- c(
  "Excitatory_neuron",
  "Inhibitory_neuron",
  "Microglia",
  "Astrocyte",
  "Oligodendrocyte",
  "Endothelial_Pericyte"
)

# =============================================================
# 8. Cell-level DE (FindMarkers) for unreplicated designs (15m)
# =============================================================
# For 15m, each comparison is n=1 vs n=1. DESeq2 pseudobulk
# can't estimate dispersion. Instead, use Seurat's FindMarkers
# (Wilcoxon) which treats cells as observations.
#
# IMPORTANT: p-values are anti-conservative due to pseudo-
# replication. These are EXPLORATORY results. Use 7m pseudobulk
# DE for primary inference.

#' Subset cells for a cell-level comparison (by sample_id)
#'
#' @param obj Seurat object
#' @param subtype Allen subclass label
#' @param samples_test Sample IDs for the test group
#' @param samples_ref  Sample IDs for the reference group
#' @param min_cells Minimum cells per group
subset_for_cellevel <- function(obj, subtype, samples_test, samples_ref,
                                filter_col = "allen_subclass_filtered",
                                min_cells = 30) {
  md <- obj@meta.data
  
  cells_test <- rownames(md)[
    md[[filter_col]] == subtype &
      md$sample_id %in% samples_test
  ]
  cells_ref <- rownames(md)[
    md[[filter_col]] == subtype &
      md$sample_id %in% samples_ref
  ]
  
  # Remove NAs
  cells_test <- cells_test[!is.na(cells_test)]
  cells_ref  <- cells_ref[!is.na(cells_ref)]
  
  valid <- length(cells_test) >= min_cells & length(cells_ref) >= min_cells
  
  reason <- if (valid) {
    "OK"
  } else {
    paste0("too few cells: test=", length(cells_test),
           ", ref=", length(cells_ref))
  }
  
  list(
    cells_test = cells_test,
    cells_ref  = cells_ref,
    n_test     = length(cells_test),
    n_ref      = length(cells_ref),
    valid      = valid,
    reason     = reason
  )
}

#' Run cell-level FindMarkers for one subtype
run_cellevel_one_subtype <- function(obj, subtype, samples_test, samples_ref,
                                     test_label, ref_label,
                                     comparison_name, out_dir = NULL,
                                     filter_col = "allen_subclass_filtered",
                                     comp = NULL,
                                     test.use = "wilcox",
                                     min_cells = 30) {
  
  if (is.null(out_dir)) out_dir <- make_de_dir(comparison_name)
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
  
  # --- Resolve thresholds ---
  thresh <- resolve_thresholds(subtype, comp)
  
  message("  ", subtype, ": ", test_label, " vs ", ref_label,
          "  [padj<", thresh$padj, ", |lfc|≥", thresh$lfc, "]")
  
  # --- Subset ---
  sub_info <- subset_for_cellevel(obj, subtype, samples_test, samples_ref,
                                  filter_col = filter_col,
                                  min_cells = min_cells)
  
  if (!sub_info$valid) {
    message("    SKIP — ", sub_info$reason)
    return(list(subtype = subtype, comparison = comparison_name,
                status = "SKIPPED", reason = sub_info$reason,
                results = NULL, summary = NULL))
  }
  
  message("    test: ", sub_info$n_test, " cells (",
          paste(samples_test, collapse = "+"), ")")
  message("    ref:  ", sub_info$n_ref,  " cells (",
          paste(samples_ref, collapse = "+"), ")")
  
  # --- Set up idents ---
  all_cells <- c(sub_info$cells_test, sub_info$cells_ref)
  obj_sub <- subset(obj, cells = all_cells)
  obj_sub <- JoinLayers(obj_sub, assay = "RNA")
  
  obj_sub$de_group <- ifelse(
    colnames(obj_sub) %in% sub_info$cells_test,
    test_label, ref_label
  )
  Idents(obj_sub) <- "de_group"
  
  # --- FindMarkers ---
  res_df <- tryCatch({
    markers <- FindMarkers(
      obj_sub,
      ident.1    = test_label,
      ident.2    = ref_label,
      test.use   = test.use,
      min.pct    = 0.1,
      logfc.threshold = 0,
      verbose    = FALSE
    )
    markers %>%
      rownames_to_column("gene") %>%
      rename(padj = p_val_adj, log2FoldChange = avg_log2FC, pvalue = p_val) %>%
      mutate(
        significant = !is.na(padj) & padj < thresh$padj &
          abs(log2FoldChange) >= thresh$lfc,
        direction = case_when(
          !significant          ~ "NS",
          log2FoldChange > 0    ~ "UP",
          log2FoldChange < 0    ~ "DOWN",
          TRUE                  ~ "NS"
        )
      ) %>%
      arrange(padj, desc(abs(log2FoldChange)))
  }, error = function(e) {
    message("    FindMarkers failed: ", e$message)
    NULL
  })
  
  if (is.null(res_df)) {
    return(list(subtype = subtype, comparison = comparison_name,
                status = "FAILED", reason = "FindMarkers error",
                results = NULL, summary = NULL))
  }
  
  # --- Summary ---
  n_up     <- sum(res_df$direction == "UP",   na.rm = TRUE)
  n_down   <- sum(res_df$direction == "DOWN", na.rm = TRUE)
  n_sig    <- n_up + n_down
  n_tested <- sum(!is.na(res_df$padj))
  
  summary_row <- data.frame(
    subtype      = subtype,
    comparison   = comparison_name,
    test_label   = test_label,
    ref_label    = ref_label,
    n_cells_test = sub_info$n_test,
    n_cells_ref  = sub_info$n_ref,
    method       = test.use,
    n_genes_tested = n_tested,
    n_sig        = n_sig,
    n_up         = n_up,
    n_down       = n_down,
    padj_cutoff  = thresh$padj,
    lfc_cutoff   = thresh$lfc,
    stringsAsFactors = FALSE
  )
  
  message("    Results: ", n_sig, " DE genes (",
          n_up, " up, ", n_down, " down) / ", n_tested, " tested")
  message("    NOTE: p-values are exploratory (no biological replicates)")
  
  # --- Save ---
  write_csv(res_df,
            file.path(out_dir,
                      paste0(comparison_name, "_", safe_name, "_full_results.csv")))
  
  sig_df <- res_df %>% filter(significant)
  if (nrow(sig_df) > 0) {
    write_csv(sig_df,
              file.path(out_dir,
                        paste0(comparison_name, "_", safe_name, "_sig_only.csv")))
  }
  
  list(subtype = subtype, comparison = comparison_name,
       status = "OK", reason = NA_character_,
       results = res_df, summary = summary_row)
}

#' Loop over subtypes for a cell-level comparison
run_cellevel_comparison <- function(obj, subtypes, samples_test, samples_ref,
                                    test_label, ref_label,
                                    comparison_name,
                                    filter_col = "allen_subclass_filtered",
                                    comp = NULL) {
  
  out_dir <- make_de_dir(comparison_name)
  all_results  <- list()
  all_summaries <- list()
  all_skipped  <- list()
  
  for (subtype in subtypes) {
    res <- run_cellevel_one_subtype(
      obj, subtype, samples_test, samples_ref,
      test_label, ref_label, comparison_name,
      filter_col = filter_col, comp = comp, out_dir = out_dir
    )
    
    if (res$status == "OK") {
      all_results[[subtype]]  <- res$results
      all_summaries[[subtype]] <- res$summary
    } else {
      all_skipped[[subtype]] <- data.frame(
        subtype = subtype, comparison = comparison_name,
        status = res$status, reason = res$reason,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # --- Combined summary ---
  summary_df <- bind_rows(all_summaries)
  if (nrow(summary_df) > 0) {
    summary_df <- summary_df %>% arrange(desc(n_sig))
    write_csv(summary_df,
              file.path(out_dir, paste0(comparison_name, "_summary.csv")))
    cat("\n\n=== ", comparison_name, " Summary (cell-level, exploratory) ===\n")
    print(as.data.frame(summary_df), row.names = FALSE)
  }
  
  skipped_df <- bind_rows(all_skipped)
  if (nrow(skipped_df) > 0) {
    write_csv(skipped_df,
              file.path(out_dir, paste0(comparison_name, "_skipped.csv")))
    cat("\nSkipped subtypes:\n")
    print(as.data.frame(skipped_df), row.names = FALSE)
  }
  
  # --- Volcano plots ---
  for (subtype in names(all_results)) {
    res_df <- all_results[[subtype]]
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    thresh <- resolve_thresholds(subtype, comp)
    
    p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),
                            color = direction)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = VOLCANO_COLORS) +
      geom_hline(yintercept = -log10(thresh$padj),
                 linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-thresh$lfc, thresh$lfc),
                 linetype = "dashed", color = "grey40") +
      theme_bw() +
      ggtitle(paste0(comparison_name, " — ", subtype, " (cell-level)")) +
      labs(subtitle = paste0("EXPLORATORY | padj<", thresh$padj,
                             ", |lfc|≥", thresh$lfc)) +
      theme(legend.position = "none",
            panel.background = element_rect(fill = "white"),
            plot.background  = element_rect(fill = "white"))
    
    ggsave(file.path(out_dir,
                     paste0(comparison_name, "_", safe_name, "_volcano.png")),
           plot = p, width = 7, height = 6, dpi = 200)
  }
  
  list(results = all_results, summary = summary_df,
       skipped = skipped_df, out_dir = out_dir)
}

# =============================================================
# 9. Comparison definitions
# =============================================================
# method = "pseudobulk" → uses DESeq2 (7m, replicated)
# method = "cellevel"   → uses FindMarkers (15m, unreplicated)

COMPARISONS <- list(
  
  # --- 7m: pseudobulk DESeq2 (3 biological replicates per condition) ---
  
  Q1 = list(
    name        = "Q1_WT_vs_WT5XFAD",
    cohort      = "7m",
    method      = "pseudobulk",
    conditions  = c("WT_5XFAD", "WT"),
    design      = ~ condition,
    subtypes    = DE_READY_7m_ALL,
    description = "7m: Effect of 5XFAD amyloid pathology (WT background)"
  ),
  Q2 = list(
    name        = "Q2_WT5XFAD_vs_Trem2KO5XFAD",
    cohort      = "7m",
    method      = "pseudobulk",
    conditions  = c("Trem2_KO_5XFAD", "WT_5XFAD"),
    design      = ~ condition,
    subtypes    = DE_READY_7m_ALL,
    description = "7m: Effect of Trem2 KO in 5XFAD background"
  ),
  Q3 = list(
    name        = "Q3_WT_vs_Trem2KO",
    cohort      = "7m",
    method      = "pseudobulk",
    conditions  = c("Trem2_KO", "WT"),
    design      = ~ condition,
    subtypes    = DE_READY_7m_ALL,
    description = "7m: Effect of Trem2 KO alone (no amyloid)"
  ),
  
  # --- 15m: cell-level FindMarkers (no replicates, exploratory) ---
  # Region comparisons: are neuronal effects region-specific?
  
  Q4a = list(
    name          = "Q4a_15m_WT_Cor_vs_Hip",
    cohort        = "15m",
    method        = "cellevel",
    samples_test  = c("WT_Cor"),
    samples_ref   = c("WT_Hip"),
    test_label    = "WT_Cor",
    ref_label     = "WT_Hip",
    subtypes      = DE_READY_15m_ALL,
    description   = "15m: Cortex vs Hippocampus in WT (baseline regional differences)"
  ),
  Q4b = list(
    name          = "Q4b_15m_WT5XFAD_Cor_vs_Hip",
    cohort        = "15m",
    method        = "cellevel",
    samples_test  = c("WT_5XFAD_Cor"),
    samples_ref   = c("WT_5XFAD_Hip"),
    test_label    = "WT_5XFAD_Cor",
    ref_label     = "WT_5XFAD_Hip",
    subtypes      = DE_READY_15m_ALL,
    description   = "15m: Cortex vs Hippocampus in 5XFAD (regional response to amyloid)"
  ),
  
  # Within-region condition comparisons
  
  Q4c = list(
    name          = "Q4c_15m_Trem2KO_Cor_vs_Hip",
    cohort        = "15m",
    method        = "cellevel",
    samples_test  = c("Trem2_KO_Cor"),
    samples_ref   = c("Trem2_KO_Hip"),
    test_label    = "Trem2_KO_Cor",
    ref_label     = "Trem2_KO_Hip",
    subtypes      = DE_READY_15m_ALL,
    description   = "15m: Cortex vs Hippocampus in Trem2 KO (regional diff without amyloid)"
  ),
  Q4d = list(
    name          = "Q4d_15m_Trem2KO5XFAD_Cor_vs_Hip",
    cohort        = "15m",
    method        = "cellevel",
    samples_test  = c("Trem2_KO_5XFAD_Cor"),
    samples_ref   = c("Trem2_KO_5XFAD_Hip"),
    test_label    = "Trem2_KO_5XFAD_Cor",
    ref_label     = "Trem2_KO_5XFAD_Hip",
    subtypes      = DE_READY_15m_ALL,
    description   = "15m: Cortex vs Hippocampus in Trem2 KO + 5XFAD (regional diff with amyloid + KO)"
  ),
  
  # --- 7m BROAD cell type DE (pseudobulk) ---
  # Same questions as Q1-Q3 but pooling all cells per broad type.
  # filter_col = "broad_type" tells the pipeline to use broad_type
  # instead of allen_subclass_filtered.
  
  Q1_broad = list(
    name        = "Q1_broad_WT_vs_WT5XFAD",
    cohort      = "7m",
    method      = "pseudobulk",
    conditions  = c("WT_5XFAD", "WT"),
    design      = ~ condition,
    subtypes    = BROAD_TYPES_7m,
    filter_col  = "broad_type",
    description = "7m BROAD: 5XFAD effect per cell class"
  ),
  Q2_broad = list(
    name        = "Q2_broad_WT5XFAD_vs_Trem2KO5XFAD",
    cohort      = "7m",
    method      = "pseudobulk",
    conditions  = c("Trem2_KO_5XFAD", "WT_5XFAD"),
    design      = ~ condition,
    subtypes    = BROAD_TYPES_7m,
    filter_col  = "broad_type",
    description = "7m BROAD: Trem2 KO in 5XFAD per cell class"
  ),
  Q3_broad = list(
    name        = "Q3_broad_WT_vs_Trem2KO",
    cohort      = "7m",
    method      = "pseudobulk",
    conditions  = c("Trem2_KO", "WT"),
    design      = ~ condition,
    subtypes    = BROAD_TYPES_7m,
    filter_col  = "broad_type",
    description = "7m BROAD: Trem2 KO alone per cell class"
  ),
  
  # --- 15m BROAD cell type DE (cell-level) ---
  
  Q4a_broad = list(
    name          = "Q4a_broad_15m_WT_Cor_vs_Hip",
    cohort        = "15m",
    method        = "cellevel",
    samples_test  = c("WT_Cor"),
    samples_ref   = c("WT_Hip"),
    test_label    = "WT_Cor",
    ref_label     = "WT_Hip",
    subtypes      = BROAD_TYPES_15m,
    filter_col    = "broad_type",
    description   = "15m BROAD: Cor vs Hip in WT"
  ),
  Q4b_broad = list(
    name          = "Q4b_broad_15m_WT5XFAD_Cor_vs_Hip",
    cohort        = "15m",
    method        = "cellevel",
    samples_test  = c("WT_5XFAD_Cor"),
    samples_ref   = c("WT_5XFAD_Hip"),
    test_label    = "WT_5XFAD_Cor",
    ref_label     = "WT_5XFAD_Hip",
    subtypes      = BROAD_TYPES_15m,
    filter_col    = "broad_type",
    description   = "15m BROAD: Cor vs Hip in 5XFAD"
  ),
  Q4c_broad = list(
    name          = "Q4c_broad_15m_Trem2KO_Cor_vs_Hip",
    cohort        = "15m",
    method        = "cellevel",
    samples_test  = c("Trem2_KO_Cor"),
    samples_ref   = c("Trem2_KO_Hip"),
    test_label    = "Trem2_KO_Cor",
    ref_label     = "Trem2_KO_Hip",
    subtypes      = BROAD_TYPES_15m,
    filter_col    = "broad_type",
    description   = "15m BROAD: Cor vs Hip in Trem2 KO"
  ),
  Q4d_broad = list(
    name          = "Q4d_broad_15m_Trem2KO5XFAD_Cor_vs_Hip",
    cohort        = "15m",
    method        = "cellevel",
    samples_test  = c("Trem2_KO_5XFAD_Cor"),
    samples_ref   = c("Trem2_KO_5XFAD_Hip"),
    test_label    = "Trem2_KO_5XFAD_Cor",
    ref_label     = "Trem2_KO_5XFAD_Hip",
    subtypes      = BROAD_TYPES_15m,
    filter_col    = "broad_type",
    description   = "15m BROAD: Cor vs Hip in Trem2 KO + 5XFAD"
  ),
  
  # --- 7m GLOBAL DE (pseudobulk across ALL cells per sample) ---
  # "Global" = treat whole population as one pseudo cell type.
  # filter_col = "global_type" needs to exist in the Seurat metadata
  # and be set to "Global" for every cell (done in 09_07 script).
  
  Q1_global = list(
    name          = "Q1_global_WT_vs_WT5XFAD",
    cohort        = "7m",
    method        = "pseudobulk",
    conditions    = c("WT_5XFAD", "WT"),
    subtypes      = c("Global"),
    design        = ~ condition,
    filter_col    = "global_type",
    description   = "7m GLOBAL: whole population, WT_5XFAD vs WT"
  ),
  Q2_global = list(
    name          = "Q2_global_WT5XFAD_vs_Trem2KO5XFAD",
    cohort        = "7m",
    method        = "pseudobulk",
    conditions    = c("Trem2_KO_5XFAD", "WT_5XFAD"),
    subtypes      = c("Global"),
    design        = ~ condition,
    filter_col    = "global_type",
    description   = "7m GLOBAL: whole population, Trem2_KO_5XFAD vs WT_5XFAD"
  ),
  Q3_global = list(
    name          = "Q3_global_WT_vs_Trem2KO",
    cohort        = "7m",
    method        = "pseudobulk",
    conditions    = c("Trem2_KO", "WT"),
    subtypes      = c("Global"),
    design        = ~ condition,
    filter_col    = "global_type",
    description   = "7m GLOBAL: whole population, Trem2_KO vs WT"
  )
)

# =============================================================
# 10. Dispatcher: run any comparison from its definition
# =============================================================
# Reads method, filter_col, etc. from the COMPARISONS entry
# so downstream scripts are just:
#   run_comparison(obj, COMPARISONS$Q1)

run_comparison <- function(obj, comp) {
  filter_col <- comp$filter_col %||% "allen_subclass_filtered"
  default_thresh <- resolve_thresholds(NULL, comp)
  
  cat("\n", comp$description, "\n")
  cat("  method: ", comp$method, " | filter_col: ", filter_col, "\n")
  cat("  default thresholds: padj<", default_thresh$padj,
      ", |lfc|≥", default_thresh$lfc)
  if (!is.null(comp$subtype_thresholds)) {
    cat("  (", length(comp$subtype_thresholds),
        " subtype-specific overrides)")
  }
  cat("\n\n")
  
  if (comp$method == "pseudobulk") {
    run_de_comparison(
      obj             = obj,
      subtypes        = comp$subtypes,
      conditions      = comp$conditions,
      comparison_name = comp$name,
      design          = comp$design,
      filter_col      = filter_col,
      comp            = comp
    )
  } else if (comp$method == "cellevel") {
    run_cellevel_comparison(
      obj             = obj,
      subtypes        = comp$subtypes,
      samples_test    = comp$samples_test,
      samples_ref     = comp$samples_ref,
      test_label      = comp$test_label,
      ref_label       = comp$ref_label,
      comparison_name = comp$name,
      filter_col      = filter_col,
      comp            = comp
    )
  } else {
    stop("Unknown method: ", comp$method)
  }
}

# =============================================================
# 11. Utility: quick summary print
# =============================================================

print_de_settings <- function() {
  cat("============ DE Settings ============\n")
  cat("  padj cutoff:           ", PADJ_CUTOFF, "\n")
  cat("  log2FC cutoff:         ", LFC_CUTOFF, "\n")
  cat("  Min cells/pseudobulk:  ", MIN_CELLS_PER_PSEUDOBULK, "\n")
  cat("  DESeq2 test:           ", DESEQ2_TEST, "\n")
  cat("  DESeq2 fit type:       ", DESEQ2_FIT_TYPE, "\n")
  cat("  Output base dir:       ", DE_BASE_DIR, "\n")
  cat("=====================================\n\n")
}

message("09_00_de_functions.R loaded successfully")
print_de_settings()