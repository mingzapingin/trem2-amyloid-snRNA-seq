# ============================================================
# 11_01_ora_subtype.R
#
# Direction-aware ORA on subtype-level DE results.
# Runs up and down separately, combines into summary table.
#
# Outputs: ./result/ora_results/subtype/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/11_00_pathway_functions.R")

ORA_OUT <- "./result/ora_results/subtype"
dir.create(ORA_OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. Q1/Q2/Q3 per subtype
# =============================================================

message("\n========== Subtype ORA: Q1/Q2/Q3 ==========\n")

comps_to_ora <- list(
  list(name = "Q1_WT_vs_WT5XFAD",
       test = "WT_5XFAD", ref = "WT",
       subtypes = DE_READY_7m_ALL),
  list(name = "Q2_WT5XFAD_vs_Trem2KO5XFAD",
       test = "Trem2_KO_5XFAD", ref = "WT_5XFAD",
       subtypes = DE_READY_7m_ALL),
  list(name = "Q3_WT_vs_Trem2KO",
       test = "Trem2_KO", ref = "WT",
       subtypes = DE_READY_7m_ALL)
)

all_comp_summaries <- list()

for (comp in comps_to_ora) {
  comp_dir <- file.path(DE_BASE_DIR, comp$name)
  out_dir  <- file.path(ORA_OUT, comp$name)
  
  message("\n--- ", comp$name, " ---")
  
  for (subtype in comp$subtypes) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    filepath  <- file.path(comp_dir,
                           paste0(comp$name, "_", safe_name,
                                  "_full_results.csv"))
    
    res <- run_ora_on_de_file(
      filepath    = filepath,
      label       = subtype,
      comparison  = comp$name,
      out_dir     = out_dir,
      test_cond   = comp$test,
      ref_cond    = comp$ref
    )
    
    if (!is.null(res$summary) && nrow(res$summary) > 0)
      all_comp_summaries[[paste0(comp$name, "_", safe_name)]] <- res$summary
  }
}

# Save grand summary across all comparisons
grand_summary <- bind_rows(all_comp_summaries)
if (nrow(grand_summary) > 0) {
  write_csv(grand_summary,
            file.path(ORA_OUT, "subtype_Q1Q2Q3_ORA_grand_summary.csv"))
  cat("\nGrand summary: ", nrow(grand_summary), " enriched terms across all\n")
}

# =============================================================
# 2. Cross-comparison categories per subtype
# =============================================================

message("\n\n========== Subtype ORA: Categories ==========\n")

CROSS_DIR  <- file.path(DE_BASE_DIR, "cross_comparison")
CAT_ORA_OUT <- file.path(ORA_OUT, "categories")

for (subtype in DE_READY_7m_ALL) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
  cat_file  <- file.path(CROSS_DIR,
                         paste0(safe_name, "_gene_categories.csv"))
  
  if (!file.exists(cat_file)) next
  
  run_ora_for_categories(
    categories_file = cat_file,
    label           = subtype,
    out_dir         = CAT_ORA_OUT,
    comparison      = "cross_comparison"
  )
}

# =============================================================
# 3. Regional categories per subtype (15m)
# =============================================================

message("\n\n========== Subtype ORA: Regional (15m) ==========\n")

REG_DIR     <- file.path(CROSS_DIR, "regional")
REG_ORA_OUT <- file.path(ORA_OUT, "regional_categories")

for (subtype in DE_READY_15m_ALL) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
  reg_file  <- file.path(REG_DIR,
                         paste0(safe_name, "_regional_categories.csv"))
  
  if (!file.exists(reg_file)) next
  
  run_ora_for_regional(
    categories_file = reg_file,
    label           = subtype,
    out_dir         = REG_ORA_OUT,
    comparison      = "regional"
  )
}

# =============================================================
# 4. Summary
# =============================================================

cat("\n\n==================== SUBTYPE ORA COMPLETE ====================\n")
cat("Outputs in: ", ORA_OUT, "\n")
cat("  <comp>/<subtype>_<dir>_<DB>.csv          — per direction per database\n")
cat("  <comp>/<subtype>_ORA_ALL_databases.csv   — combined summary with direction\n")
cat("  <comp>/<subtype>_ORA_direction_plot.png   — two-panel up/down plot\n")
cat("  <comp>/<subtype>_syngo_mapping_report.csv — SynGO QC (neuronal only)\n")
cat("  subtype_Q1Q2Q3_ORA_grand_summary.csv     — all comparisons combined\n")