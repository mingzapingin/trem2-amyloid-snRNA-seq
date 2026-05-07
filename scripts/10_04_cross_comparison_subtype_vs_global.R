# ============================================================
# 10_04_cross_comparison_subtype_vs_global.R
#
# Subtype × global cross-comparison.
# Same logic as 10_03 but with all 25 subtypes instead of 6
# broad types. Subtype-specific genes here are the strongest
# cell-type-specific findings — they survived two filters:
# DE in the subtype AND not explained by global changes.
#
# Outputs: ./result/cross_comparison/subtype_vs_global/Q{1,2,3}/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/10_00_cross_comparison_functions.R")

OUT_BASE <- "./result/cross_comparison/subtype_vs_global"
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

comps <- list(
  Q1 = list(subtype = "Q1_WT_vs_WT5XFAD",
            global  = "Q1_global_WT_vs_WT5XFAD"),
  Q2 = list(subtype = "Q2_WT5XFAD_vs_Trem2KO5XFAD",
            global  = "Q2_global_WT5XFAD_vs_Trem2KO5XFAD"),
  Q3 = list(subtype = "Q3_WT_vs_Trem2KO",
            global  = "Q3_global_WT_vs_Trem2KO")
)

for (q_name in names(comps)) {
  q <- comps[[q_name]]
  out_dir <- file.path(OUT_BASE, q_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\n\n========== ", q_name, " ==========")
  
  global_de <- load_de_table(q$global, "Global")
  if (is.null(global_de)) {
    message("  Missing global DE — run 09_07 first")
    next
  }
  
  all_counts <- list()
  all_classified <- list()
  
  for (subtype in DE_READY_7m_ALL) {
    ct_de <- load_de_table(q$subtype, subtype)
    if (is.null(ct_de)) next
    
    message("  ", subtype)
    
    classified <- classify_vs_global(
      ct_de_table     = ct_de,
      global_de_table = global_de,
      cell_label      = subtype
    )
    
    all_classified[[subtype]] <- classified
    all_counts[[subtype]]     <- count_vs_global_categories(classified)
    
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    write_csv(classified,
              file.path(out_dir,
                        paste0(safe_name, "_vs_global_categories.csv")))
    
    p_sc <- plot_vs_global_scatter(
      classified,
      title = paste0(subtype, " vs Global — ", q_name)
    )
    if (!is.null(p_sc))
      ggsave(file.path(out_dir,
                       paste0(safe_name, "_vs_global_scatter.png")),
             plot = p_sc, width = 8, height = 7, dpi = 200)
    
    # Venn diagram
    plot_vs_global_venn(
      classified,
      title    = paste0(subtype, " vs Global — ", q_name),
      filepath = file.path(out_dir,
                           paste0(safe_name, "_vs_global_venn.png"))
    )
  }
  
  if (length(all_counts) == 0) next
  
  counts_df <- bind_rows(all_counts)
  write_csv(counts_df, file.path(out_dir, "category_counts.csv"))
  
  p_hm <- plot_vs_global_heatmap(
    counts_df,
    title = paste0(q_name, ": subtypes vs global")
  )
  if (!is.null(p_hm))
    ggsave(file.path(out_dir, "category_counts_heatmap.png"),
           plot = p_hm, width = 9, height = 10, dpi = 200)
  
  # Subtype-specificity summary — for each gene, which subtype
  # called it specific / shared / mismatch
  specificity <- bind_rows(all_classified) %>%
    filter(vs_global_category != "Not_DE") %>%
    select(gene, cell_type, vs_global_category, ct_lfc, gl_lfc) %>%
    arrange(gene, vs_global_category)
  
  # Also: wide format — one row per gene, columns = subtype × category
  specificity_wide <- specificity %>%
    select(gene, cell_type, vs_global_category) %>%
    distinct() %>%
    mutate(value = vs_global_category) %>%
    pivot_wider(names_from = cell_type, values_from = value,
                values_fill = "")
  
  write_csv(specificity,
            file.path(out_dir, "subtype_specificity_long.csv"))
  write_csv(specificity_wide,
            file.path(out_dir, "subtype_specificity_wide.csv"))
  
  cat("\n=== ", q_name, " subtype category counts ===\n")
  print(as.data.frame(counts_df %>%
                        pivot_wider(names_from = vs_global_category, values_from = n_genes,
                                    values_fill = 0) %>%
                        arrange(desc(rowSums(across(where(is.numeric)))))),
        row.names = FALSE)
}

cat("\n\n==================== SUBTYPE vs GLOBAL COMPLETE ====================\n")
cat("Outputs in: ", OUT_BASE, "\n")