# ============================================================
# 13_00_propeller_functions.R
#
# Shared infrastructure for compositional abundance analysis
# using propeller (from speckle package).
#
# Tests whether cell type proportions differ between conditions,
# using arcsine-square-root-transformed per-sample proportions
# with empirical Bayes moderation.
#
# Source from any propeller script:
#   source("./scripts/13_00_propeller_functions.R")
#
# Requires:
#   BiocManager::install("speckle")
#   BiocManager::install("limma")   # dependency
# ============================================================

suppressPackageStartupMessages({
  library(speckle)
  library(limma)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(Seurat)
})

# Fix namespace conflicts
select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

# =============================================================
# 1. Settings
# =============================================================

# Transformation for propeller
PROP_TRANSFORM <- "asin"   # arcsine square root, default for proportions

# Filters for subtype analysis
MIN_CELLS_PER_SUBTYPE  <- 50    # min total cells for a subtype to be tested
MIN_SAMPLES_WITH_CELLS <- 2     # min samples with cells present for a subtype

# Color palette (matches 09_00 where possible)
COND_COLORS <- c(
  "WT"                 = "#4DBBD5",
  "WT_5XFAD"           = "#E64B35",
  "Trem2_KO"           = "#00A087",
  "Trem2_KO_5XFAD"     = "#B07AA1",
  "WT_Cor"             = "#4DBBD5",
  "WT_Hip"             = "#3182BD",
  "WT_5XFAD_Cor"       = "#E64B35",
  "WT_5XFAD_Hip"       = "#CB3A27",
  "Trem2_KO_Cor"       = "#00A087",
  "Trem2_KO_Hip"       = "#007A65",
  "Trem2_KO_5XFAD_Cor" = "#B07AA1",
  "Trem2_KO_5XFAD_Hip" = "#8B5A80"
)

ABUNDANCE_BASE_DIR <- "./result/abundance"

# =============================================================
# 2. Calculate per-sample proportions
# =============================================================

#' Build a sample × cluster proportion table
#'
#' @param meta data.frame with columns: sample, condition, cluster
#' @return list with:
#'   - props:    numeric matrix, samples × clusters (proportions)
#'   - counts:   integer matrix, samples × clusters (raw cell counts)
#'   - sample_info: data.frame with one row per sample (sample, condition, total_cells)
calculate_proportions_per_sample <- function(meta) {
  
  required <- c("sample", "condition", "cluster")
  missing  <- setdiff(required, colnames(meta))
  if (length(missing) > 0)
    stop("Missing columns: ", paste(missing, collapse = ", "))
  
  # Counts per sample × cluster
  counts <- meta %>%
    count(sample, cluster, name = "n") %>%
    pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
  
  # Keep sample column separate
  sample_ids <- counts$sample
  count_mat  <- as.matrix(counts %>% select(-sample))
  rownames(count_mat) <- sample_ids
  
  # Proportions (row-normalized)
  total_per_sample <- rowSums(count_mat)
  prop_mat <- sweep(count_mat, 1, total_per_sample, FUN = "/")
  
  # Sample info table (sample → condition)
  sample_info <- meta %>%
    distinct(sample, condition) %>%
    mutate(total_cells = total_per_sample[sample])
  
  list(
    props       = prop_mat,
    counts      = count_mat,
    sample_info = sample_info
  )
}

# =============================================================
# 3. Filter subtypes for propeller
# =============================================================

#' Filter subtypes that have enough cells and sample representation
filter_subtypes_for_abundance <- function(counts_mat,
                                          min_cells   = MIN_CELLS_PER_SUBTYPE,
                                          min_samples = MIN_SAMPLES_WITH_CELLS) {
  total_cells   <- colSums(counts_mat)
  samples_with  <- colSums(counts_mat > 0)
  
  keep <- total_cells >= min_cells & samples_with >= min_samples
  
  msg <- paste0(
    "  Filter: ", sum(keep), "/", length(keep),
    " subtypes pass (min ", min_cells, " cells, min ",
    min_samples, " samples with cells)"
  )
  message(msg)
  
  if (sum(!keep) > 0) {
    dropped <- colnames(counts_mat)[!keep]
    message("  Dropped: ", paste(dropped, collapse = ", "))
  }
  
  names(keep)[keep]
}

# =============================================================
# 4. Run propeller
# =============================================================

#' Run propeller compositional test
#'
#' @param counts_mat samples × clusters matrix of counts
#' @param sample_info data.frame with sample, condition columns
#' @param formula formula for the design matrix (e.g., ~ condition or ~ region + condition)
#' @param reference_cond reference level for condition factor
#' @param transform "asin" or "logit"
#' @return data.frame with propeller results
run_propeller <- function(counts_mat, sample_info,
                          formula_str = "~ condition",
                          reference_cond = NULL,
                          transform = PROP_TRANSFORM) {
  
  # Ensure sample order matches
  sample_info <- sample_info %>% arrange(match(sample, rownames(counts_mat)))
  stopifnot(all(sample_info$sample == rownames(counts_mat)))
  
  # Condition as factor
  if (!is.null(reference_cond))
    sample_info$condition <- relevel(factor(sample_info$condition),
                                     ref = reference_cond)
  else
    sample_info$condition <- factor(sample_info$condition)
  
  # Design matrix
  design <- model.matrix(as.formula(formula_str), data = sample_info)
  
  # Propeller needs cluster × sample matrix of proportions
  # Re-derive from counts for consistency
  props <- t(counts_mat / rowSums(counts_mat))
  
  # Run propeller — use the generic propeller() with transformed props
  result <- tryCatch({
    propeller(
      clusters = rep(colnames(counts_mat), ncol(props)),
      sample   = rep(rownames(counts_mat), each = ncol(counts_mat)),
      group    = rep(sample_info$condition[match(rownames(counts_mat),
                                                 sample_info$sample)],
                     each = ncol(counts_mat)),
      transform = transform
    )
  }, error = function(e) {
    # Fallback: use propeller.ttest() directly with transformed proportions
    message("  propeller() error: ", e$message, " — trying manual path")
    NULL
  })
  
  # If propeller's convenience function failed, build it manually
  if (is.null(result)) {
    # Transform proportions
    if (transform == "asin") {
      trans_props <- asin(sqrt(props))
    } else {
      trans_props <- log(props / (1 - props))
    }
    
    # Run propeller.ttest for two-group, propeller.anova for >2
    n_groups <- length(levels(sample_info$condition))
    if (n_groups == 2) {
      result <- propeller.ttest(
        prop.list = list(Proportions = props,
                         TransformedProps = trans_props),
        design    = design,
        contrasts = c(0, 1),
        robust    = TRUE, trend = FALSE, sort = TRUE
      )
    } else {
      result <- propeller.anova(
        prop.list = list(Proportions = props,
                         TransformedProps = trans_props),
        design    = design,
        coef      = 2:n_groups,
        robust    = TRUE, trend = FALSE, sort = TRUE
      )
    }
  }
  
  as.data.frame(result) %>% rownames_to_column("cluster")
}

# =============================================================
# 5. Summary table with outlier detection
# =============================================================

#' Format propeller output into a publication-ready summary
#'
#' Adds per-condition mean proportion, difference, max per-sample
#' deviation, and full statistics.
make_summary_table <- function(prop_res,
                               props_mat,
                               sample_info,
                               comparison_name,
                               test_cond, ref_cond) {
  
  # Per-condition mean proportions
  cond_means <- as.data.frame(props_mat) %>%
    rownames_to_column("sample") %>%
    left_join(sample_info %>% select(sample, condition), by = "sample") %>%
    pivot_longer(-c(sample, condition),
                 names_to = "cluster", values_to = "prop") %>%
    group_by(cluster, condition) %>%
    summarise(mean_prop = mean(prop),
              sd_prop   = sd(prop),
              .groups   = "drop") %>%
    pivot_wider(names_from = condition,
                values_from = c(mean_prop, sd_prop),
                names_sep = "_")
  
  # Max per-sample deviation from its group mean (outlier flag)
  outlier_check <- as.data.frame(props_mat) %>%
    rownames_to_column("sample") %>%
    left_join(sample_info %>% select(sample, condition), by = "sample") %>%
    pivot_longer(-c(sample, condition),
                 names_to = "cluster", values_to = "prop") %>%
    group_by(cluster, condition) %>%
    mutate(group_mean = mean(prop),
           abs_dev = abs(prop - group_mean)) %>%
    group_by(cluster) %>%
    summarise(max_sample_deviation = round(max(abs_dev), 4),
              .groups = "drop")
  
  # Number of samples per condition
  n_per_cond <- sample_info %>%
    count(condition, name = "n_samples") %>%
    pivot_wider(names_from = condition, values_from = n_samples,
                names_prefix = "n_")
  
  # Merge with propeller results
  summary <- prop_res %>%
    rename_with(~ gsub("^P\\.", "p_", .x)) %>%
    left_join(cond_means, by = "cluster") %>%
    left_join(outlier_check, by = "cluster") %>%
    mutate(
      comparison = comparison_name,
      test_cond  = test_cond,
      ref_cond   = ref_cond,
      .before    = 1
    )
  
  # Add n_samples columns if we can
  if (ncol(n_per_cond) > 0) {
    for (col in colnames(n_per_cond))
      summary[[col]] <- n_per_cond[[col]][1]
  }
  
  summary
}

# =============================================================
# 6. Plotting
# =============================================================

#' Per-sample proportion plot — points colored by condition, one
#' facet per cluster. Shows outliers immediately.
plot_proportions_per_sample <- function(props_mat, sample_info,
                                        title = "Per-sample proportions") {
  
  plot_df <- as.data.frame(props_mat) %>%
    rownames_to_column("sample") %>%
    left_join(sample_info %>% select(sample, condition), by = "sample") %>%
    pivot_longer(-c(sample, condition),
                 names_to = "cluster", values_to = "proportion")
  
  cond_levels <- unique(plot_df$condition)
  cond_palette <- COND_COLORS[intersect(names(COND_COLORS), cond_levels)]
  missing_cols <- setdiff(cond_levels, names(cond_palette))
  if (length(missing_cols) > 0) {
    extra <- rep(c("#FF7F0E", "#2CA02C", "#9467BD"), length.out = length(missing_cols))
    cond_palette <- c(cond_palette, setNames(extra, missing_cols))
  }
  
  ggplot(plot_df, aes(x = condition, y = proportion, color = condition)) +
    geom_jitter(width = 0.15, size = 2.2, alpha = 0.85) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.5, linewidth = 0.4, color = "grey30") +
    facet_wrap(~ cluster, scales = "free_y") +
    scale_color_manual(values = cond_palette) +
    theme_bw(base_size = 11) +
    ggtitle(title) +
    xlab(NULL) + ylab("Proportion of cells") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none",
          strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
}

#' Boxplot of proportions per condition per cluster
plot_proportions_boxplot <- function(props_mat, sample_info,
                                     title = "Proportions") {
  
  plot_df <- as.data.frame(props_mat) %>%
    rownames_to_column("sample") %>%
    left_join(sample_info %>% select(sample, condition), by = "sample") %>%
    pivot_longer(-c(sample, condition),
                 names_to = "cluster", values_to = "proportion")
  
  cond_levels <- unique(plot_df$condition)
  cond_palette <- COND_COLORS[intersect(names(COND_COLORS), cond_levels)]
  missing_cols <- setdiff(cond_levels, names(cond_palette))
  if (length(missing_cols) > 0) {
    extra <- rep(c("#FF7F0E", "#2CA02C", "#9467BD"), length.out = length(missing_cols))
    cond_palette <- c(cond_palette, setNames(extra, missing_cols))
  }
  
  ggplot(plot_df, aes(x = condition, y = proportion, fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
    facet_wrap(~ cluster, scales = "free_y") +
    scale_fill_manual(values = cond_palette) +
    theme_bw(base_size = 11) +
    ggtitle(title) +
    xlab(NULL) + ylab("Proportion of cells") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none",
          strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
}

# =============================================================
# 7. Full pipeline: metadata → propeller → outputs
# =============================================================

#' Run one propeller comparison end-to-end
#'
#' @param meta     metadata with sample, condition, cluster columns
#' @param test_cond label for test condition
#' @param ref_cond  label for reference condition
#' @param comparison_name short name (e.g., "Q1")
#' @param out_dir   output directory
#' @param filter_min_cells whether to apply subtype cell filter
#' @param formula_str design formula
run_propeller_comparison <- function(meta,
                                     test_cond,
                                     ref_cond,
                                     comparison_name,
                                     out_dir,
                                     filter_min_cells = FALSE,
                                     formula_str      = "~ condition") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n=== ", comparison_name, ": ", test_cond, " vs ", ref_cond, " ===")
  
  # Subset to the two conditions
  meta_sub <- meta %>% filter(condition %in% c(test_cond, ref_cond))
  n_samples <- n_distinct(meta_sub$sample)
  message("  Samples: ", n_samples, ", cells: ", nrow(meta_sub))
  
  if (n_samples < 2) {
    message("  Too few samples — skipping")
    return(NULL)
  }
  
  # Build proportion matrices
  prop_data <- calculate_proportions_per_sample(meta_sub)
  counts_mat <- prop_data$counts
  props_mat  <- prop_data$props
  sample_info <- prop_data$sample_info
  
  # Optional filter for subtypes
  if (filter_min_cells) {
    keep_clusters <- filter_subtypes_for_abundance(counts_mat)
    if (length(keep_clusters) < 2) {
      message("  Too few clusters after filtering — skipping")
      return(NULL)
    }
    counts_mat <- counts_mat[, keep_clusters, drop = FALSE]
    props_mat  <- props_mat[, keep_clusters, drop = FALSE]
  }
  
  # Run propeller
  prop_res <- tryCatch({
    run_propeller(counts_mat, sample_info,
                  formula_str = formula_str,
                  reference_cond = ref_cond)
  }, error = function(e) {
    message("  propeller failed: ", e$message)
    NULL
  })
  
  if (is.null(prop_res)) return(NULL)
  
  # Build summary
  summary <- make_summary_table(prop_res, props_mat, sample_info,
                                comparison_name, test_cond, ref_cond)
  
  # Save
  write_csv(summary,
            file.path(out_dir,
                      paste0(comparison_name, "_propeller_results.csv")))
  
  # Per-sample plot
  p1 <- plot_proportions_per_sample(
    props_mat, sample_info,
    title = paste0(comparison_name, ": per-sample proportions")
  )
  n_clust <- ncol(props_mat)
  ncol_facet <- min(4, n_clust)
  nrow_facet <- ceiling(n_clust / ncol_facet)
  ggsave(file.path(out_dir,
                   paste0(comparison_name, "_per_sample_proportions.png")),
         plot = p1, width = 3 + ncol_facet * 2.2,
         height = 2.5 + nrow_facet * 2.2, dpi = 200)
  
  # Boxplot
  p2 <- plot_proportions_boxplot(
    props_mat, sample_info,
    title = paste0(comparison_name, ": proportions boxplot")
  )
  ggsave(file.path(out_dir,
                   paste0(comparison_name, "_boxplot.png")),
         plot = p2, width = 3 + ncol_facet * 2.2,
         height = 2.5 + nrow_facet * 2.2, dpi = 200)
  
  n_sig <- sum(get_fdr_values(summary) < 0.05, na.rm = TRUE)
  message("  Saved. ", nrow(summary), " clusters tested, ",
          n_sig, " significant (FDR<0.05)")
  
  summary
}

# =============================================================
# 7b. Helper: extract FDR/padj column (handles varying propeller output)
# =============================================================

#' propeller returns different column names depending on test type:
#'   propeller.ttest → FDR
#'   propeller.anova → FDR (sometimes "fdr")
#'   Some speckle versions → "adj.P.Val"
#' This helper picks whichever exists and returns it as a numeric vector.
get_fdr_values <- function(df) {
  candidates <- c("FDR", "fdr", "adj.P.Val", "adj_P_Val", "padj")
  found <- intersect(candidates, colnames(df))
  if (length(found) == 0) return(rep(NA_real_, nrow(df)))
  df[[found[1]]]
}

# =============================================================
# 8. Metadata loader helper
# =============================================================

#' Extract a clean propeller-ready metadata table from a Seurat object
#'
#' @param obj Seurat object with sample_id, condition, and a cluster column
#' @param cluster_col Name of the column holding cluster/cell-type labels
build_meta_from_seurat <- function(obj, cluster_col = "broad_type") {
  
  md <- obj@meta.data
  if (!cluster_col %in% colnames(md))
    stop("Column not found: ", cluster_col)
  
  data.frame(
    sample    = md$sample_id,
    condition = md$condition,
    cluster   = md[[cluster_col]],
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(cluster), cluster != "Unknown", cluster != "")
}

message("13_00_propeller_functions.R loaded successfully")
message("  Transform: ", PROP_TRANSFORM,
        " | Min cells/subtype: ", MIN_CELLS_PER_SUBTYPE)