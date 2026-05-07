# ============================================================
# 14_01_partition_ora_broad.R
#
# Partition-based ORA on broad cell types.
# For each broad type, builds gene set partitions from Q1/Q2/Q3
# overlaps and runs ORA on each non-empty partition.
#
# Also builds cross-cell-type pathway heatmaps showing which
# pathways are shared vs specific across broad types for each
# partition type.
#
# Outputs: ./result/pathway/partition_ora/broad/
# ============================================================

source("./scripts/14_00_partition_functions.R")

OUT_BASE <- "./result/pathway/partition_ora/broad"
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. Run partition analysis per broad type
# =============================================================

message("\n========== Broad Type Partition ORA ==========\n")

all_results <- list()

for (bt in BROAD_TYPES_7m) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
  out_dir   <- file.path(OUT_BASE, safe_name)
  
  # Load Q1/Q2/Q3 DE tables for this broad type
  de_tables <- list(
    Q1 = load_de_for_partition("Q1_broad_WT_vs_WT5XFAD", bt),
    Q2 = load_de_for_partition("Q2_broad_WT5XFAD_vs_Trem2KO5XFAD", bt),
    Q3 = load_de_for_partition("Q3_broad_WT_vs_Trem2KO", bt)
  )
  
  res <- run_partition_analysis(
    de_tables = de_tables,
    label     = bt,
    out_dir   = out_dir
  )
  
  if (!is.null(res))
    all_results[[bt]] <- res
}

# =============================================================
# 2. Cross-cell-type summaries
# =============================================================

if (length(all_results) > 0) {
  
  CROSS_DIR <- file.path(OUT_BASE, "cross_celltype")
  dir.create(CROSS_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Partition size summary
  p_sizes <- plot_partition_sizes(all_results,
                                  title = "Broad types: partition sizes (Q1/Q2/Q3)")
  if (!is.null(p_sizes))
    ggsave(file.path(CROSS_DIR, "partition_size_summary.png"),
           plot = p_sizes, width = 10, height = 8, dpi = 200)
  
  # Pathway heatmaps for key partitions
  partition_patterns <- c(
    "Q1_only"       = "Q1_only",
    "Q2_only"       = "Q2_only",
    "Q3_only"       = "Q3_only",
    "Q1_AND_Q2"     = "Q1_AND_Q2",
    "Q1_AND_Q3"     = "Q1_AND_Q3",
    "Q2_AND_Q3"     = "Q2_AND_Q3",
    "all_Q1_Q2_Q3"  = "all_Q1_Q2_Q3"
  )
  
  for (pname in names(partition_patterns)) {
    p_hm <- plot_partition_pathway_heatmap(
      all_results,
      partition_pattern = partition_patterns[[pname]],
      title = paste0("Broad types: ", pname, " pathways")
    )
    if (!is.null(p_hm)) {
      safe_pname <- gsub("[^A-Za-z0-9_-]", "_", pname)
      ggsave(file.path(CROSS_DIR,
                       paste0(safe_pname, "_pathway_heatmap.png")),
             plot = p_hm, width = 10, height = 8, dpi = 200)
    }
  }
  
  # Grand summary: all partitions × all cell types
  grand <- bind_rows(lapply(names(all_results), function(ct) {
    res <- all_results[[ct]]
    if (!is.null(res$summary) && nrow(res$summary) > 0)
      res$summary %>% mutate(cell_type = ct)
    else NULL
  }))
  
  if (nrow(grand) > 0)
    write_csv(grand,
              file.path(CROSS_DIR, "broad_partition_ORA_grand_summary.csv"))
}

cat("\n\n==================== BROAD PARTITION ORA COMPLETE ====================\n")
cat("Outputs in: ", OUT_BASE, "\n")
cat("  <broad_type>/  — per-partition ORA results + gene lists\n")
cat("  cross_celltype/ — pathway heatmaps + grand summary\n")