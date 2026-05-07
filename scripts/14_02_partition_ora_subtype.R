# ============================================================
# 14_02_partition_ora_subtype.R
#
# Partition-based ORA on neuronal subtypes.
# For each subtype, builds gene set partitions from Q1/Q2/Q3
# overlaps and runs ORA on each non-empty partition.
#
# Many subtypes will have small or empty partitions —
# handled gracefully with minimum-size filtering.
#
# Outputs: ./result/pathway/partition_ora/subtype/
# ============================================================

source("./scripts/14_00_partition_functions.R")

OUT_BASE <- "./result/pathway/partition_ora/subtype"
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. Run partition analysis per subtype
# =============================================================

message("\n========== Subtype Partition ORA ==========\n")

all_results <- list()

for (subtype in DE_READY_7m_ALL) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
  out_dir   <- file.path(OUT_BASE, safe_name)
  
  de_tables <- list(
    Q1 = load_de_for_partition("Q1_WT_vs_WT5XFAD", subtype),
    Q2 = load_de_for_partition("Q2_WT5XFAD_vs_Trem2KO5XFAD", subtype),
    Q3 = load_de_for_partition("Q3_WT_vs_Trem2KO", subtype)
  )
  
  res <- run_partition_analysis(
    de_tables = de_tables,
    label     = subtype,
    out_dir   = out_dir
  )
  
  if (!is.null(res))
    all_results[[subtype]] <- res
}

# =============================================================
# 2. Cross-subtype summaries
# =============================================================

if (length(all_results) > 0) {
  
  CROSS_DIR <- file.path(OUT_BASE, "cross_subtype")
  dir.create(CROSS_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Partition size heatmap — which subtypes have which partitions
  p_sizes <- plot_partition_sizes(
    all_results,
    title = "Subtypes: partition sizes (Q1/Q2/Q3)"
  )
  if (!is.null(p_sizes))
    ggsave(file.path(CROSS_DIR, "partition_size_summary.png"),
           plot = p_sizes, width = 14, height = 10, dpi = 200)
  
  # Pathway heatmaps per partition type
  partition_patterns <- c(
    "Q1_only"       = "Q1_only",
    "Q2_only"       = "Q2_only",
    "Q1_AND_Q2"     = "Q1_AND_Q2",
    "Q1_AND_Q3"     = "Q1_AND_Q3",
    "all_Q1_Q2_Q3"  = "all_Q1_Q2_Q3"
  )
  
  for (pname in names(partition_patterns)) {
    p_hm <- plot_partition_pathway_heatmap(
      all_results,
      partition_pattern = partition_patterns[[pname]],
      title = paste0("Subtypes: ", pname, " pathways")
    )
    if (!is.null(p_hm)) {
      safe_pname <- gsub("[^A-Za-z0-9_-]", "_", pname)
      ggsave(file.path(CROSS_DIR,
                       paste0(safe_pname, "_pathway_heatmap.png")),
             plot = p_hm, width = 14, height = 10, dpi = 200)
    }
  }
  
  # Grand summary
  grand <- bind_rows(lapply(names(all_results), function(ct) {
    res <- all_results[[ct]]
    if (!is.null(res$summary) && nrow(res$summary) > 0)
      res$summary %>% mutate(cell_type = ct)
    else NULL
  }))
  
  if (nrow(grand) > 0)
    write_csv(grand,
              file.path(CROSS_DIR,
                        "subtype_partition_ORA_grand_summary.csv"))
  
  # Which subtypes have triple-overlap genes?
  triple_subtypes <- bind_rows(lapply(names(all_results), function(ct) {
    info <- all_results[[ct]]$partition_info
    if (is.null(info)) return(NULL)
    triple <- info %>% filter(grepl("all_", partition))
    if (nrow(triple) > 0)
      triple %>% mutate(cell_type = ct)
    else NULL
  }))
  
  if (nrow(triple_subtypes) > 0) {
    cat("\n\n=== Subtypes with Q1 ∩ Q2 ∩ Q3 overlap ===\n")
    print(as.data.frame(triple_subtypes %>%
                          select(cell_type, partition, n_genes)),
          row.names = FALSE)
  }
}

cat("\n\n==================== SUBTYPE PARTITION ORA COMPLETE ====================\n")
cat("Outputs in: ", OUT_BASE, "\n")
cat("  <subtype>/  — per-partition ORA results + gene lists\n")
cat("  cross_subtype/ — heatmaps + grand summary\n")