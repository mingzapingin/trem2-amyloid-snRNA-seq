# ============================================================
# 10_03_cross_comparison_broad_vs_global.R
#
# For each broad type × comparison (Q1/Q2/Q3), classify every
# DE gene as:
#   - Cell_type_specific (DE broad, NOT global)
#   - Shared_with_global (DE both, same direction)
#   - Global_only (DE global, NOT broad)
#   - Direction_mismatch (DE both, opposite direction)
#
# Outputs: ./result/cross_comparison/broad_vs_global/Q{1,2,3}/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/10_00_cross_comparison_functions.R")

OUT_BASE <- "./result/cross_comparison/broad_vs_global"
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. Comparisons map: broad comparison → matching global
# =============================================================

comps <- list(
  Q1 = list(broad  = "Q1_broad_WT_vs_WT5XFAD",
            global = "Q1_global_WT_vs_WT5XFAD"),
  Q2 = list(broad  = "Q2_broad_WT5XFAD_vs_Trem2KO5XFAD",
            global = "Q2_global_WT5XFAD_vs_Trem2KO5XFAD"),
  Q3 = list(broad  = "Q3_broad_WT_vs_Trem2KO",
            global = "Q3_global_WT_vs_Trem2KO")
)

# =============================================================
# 2. Run per comparison
# =============================================================

for (q_name in names(comps)) {
  q <- comps[[q_name]]
  out_dir <- file.path(OUT_BASE, q_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\n\n========== ", q_name, " ==========")
  
  # Load global DE
  global_de <- load_de_table(q$global, "Global")
  if (is.null(global_de)) {
    message("  Missing global DE — run 09_07 first")
    next
  }
  
  all_counts <- list()
  all_classified <- list()
  
  for (bt in BROAD_TYPES_7m) {
    ct_de <- load_de_table(q$broad, bt)
    if (is.null(ct_de)) {
      message("  ", bt, ": no DE table — skipping")
      next
    }
    
    message("  ", bt)
    
    classified <- classify_vs_global(
      ct_de_table     = ct_de,
      global_de_table = global_de,
      cell_label      = bt
    )
    
    all_classified[[bt]] <- classified
    all_counts[[bt]]     <- count_vs_global_categories(classified)
    
    # Per-cell-type output
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
    write_csv(classified,
              file.path(out_dir,
                        paste0(safe_name, "_vs_global_categories.csv")))
    
    # LFC scatter
    p_sc <- plot_vs_global_scatter(
      classified,
      title = paste0(bt, " vs Global — ", q_name)
    )
    if (!is.null(p_sc))
      ggsave(file.path(out_dir,
                       paste0(safe_name, "_vs_global_scatter.png")),
             plot = p_sc, width = 8, height = 7, dpi = 200)
  }
  
  if (length(all_counts) == 0) next
  
  # Combined category counts
  counts_df <- bind_rows(all_counts)
  write_csv(counts_df, file.path(out_dir, "category_counts.csv"))
  
  # Heatmap
  p_hm <- plot_vs_global_heatmap(
    counts_df,
    title = paste0(q_name, ": broad types vs global")
  )
  if (!is.null(p_hm))
    ggsave(file.path(out_dir, "category_counts_heatmap.png"),
           plot = p_hm, width = 8, height = 5, dpi = 200)
  
  # Shared genes summary — one row per gene per cell type
  shared_genes <- bind_rows(all_classified) %>%
    filter(vs_global_category != "Not_DE") %>%
    select(cell_type, gene, vs_global_category,
           ct_lfc, ct_padj, gl_lfc, gl_padj) %>%
    arrange(vs_global_category, cell_type, gene)
  
  write_csv(shared_genes, file.path(out_dir, "shared_genes_summary.csv"))
  
  cat("\n=== ", q_name, " category counts ===\n")
  print(as.data.frame(counts_df %>%
                        pivot_wider(names_from = vs_global_category, values_from = n_genes,
                                    values_fill = 0)), row.names = FALSE)
}

cat("\n\n==================== BROAD vs GLOBAL COMPLETE ====================\n")
cat("Outputs in: ", OUT_BASE, "\n")