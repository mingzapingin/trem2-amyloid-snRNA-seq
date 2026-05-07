# ============================================================
# 12_01_gsea_subtype.R
#
# Direction-aware GSEA on subtype-level DE results.
# Q1/Q2/Q3 (7m) + Q4a-Q4d (15m).
#
# Outputs: ./result/gsea_results/subtype/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/12_00_gsea_functions.R")

GSEA_OUT <- "./result/gsea_results/subtype"
dir.create(GSEA_OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. 7m subtype GSEA: Q1/Q2/Q3
# =============================================================

message("\n========== Subtype GSEA: 7m Q1/Q2/Q3 ==========\n")

comps_7m <- list(
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

all_summaries <- list()

for (comp in comps_7m) {
  comp_dir <- file.path(DE_BASE_DIR, comp$name)
  out_dir  <- file.path(GSEA_OUT, comp$name)
  
  message("\n--- ", comp$name, " ---")
  
  for (subtype in comp$subtypes) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    filepath  <- file.path(comp_dir,
                           paste0(comp$name, "_", safe_name,
                                  "_full_results.csv"))
    
    res <- run_gsea_on_de_file(
      filepath   = filepath,
      label      = subtype,
      comparison = comp$name,
      out_dir    = out_dir,
      test_cond  = comp$test,
      ref_cond   = comp$ref
    )
    
    if (!is.null(res$combined) && nrow(res$combined) > 0)
      all_summaries[[paste0(comp$name, "_", safe_name)]] <- res$combined
  }
}

# =============================================================
# 2. 15m subtype GSEA: Q4a-Q4d
# =============================================================

message("\n\n========== Subtype GSEA: 15m Q4a-Q4d ==========\n")

comps_15m <- list(
  list(name = "Q4a_15m_WT_Cor_vs_Hip",
       test = "WT_Cor", ref = "WT_Hip",
       subtypes = DE_READY_15m_ALL),
  list(name = "Q4b_15m_WT5XFAD_Cor_vs_Hip",
       test = "WT_5XFAD_Cor", ref = "WT_5XFAD_Hip",
       subtypes = DE_READY_15m_ALL),
  list(name = "Q4c_15m_Trem2KO_Cor_vs_Hip",
       test = "Trem2_KO_Cor", ref = "Trem2_KO_Hip",
       subtypes = DE_READY_15m_ALL),
  list(name = "Q4d_15m_Trem2KO5XFAD_Cor_vs_Hip",
       test = "Trem2_KO_5XFAD_Cor", ref = "Trem2_KO_5XFAD_Hip",
       subtypes = DE_READY_15m_ALL)
)

for (comp in comps_15m) {
  comp_dir <- file.path(DE_BASE_DIR, comp$name)
  out_dir  <- file.path(GSEA_OUT, comp$name)
  
  message("\n--- ", comp$name, " ---")
  
  for (subtype in comp$subtypes) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    filepath  <- file.path(comp_dir,
                           paste0(comp$name, "_", safe_name,
                                  "_full_results.csv"))
    
    res <- run_gsea_on_de_file(
      filepath   = filepath,
      label      = subtype,
      comparison = comp$name,
      out_dir    = out_dir,
      test_cond  = comp$test,
      ref_cond   = comp$ref
    )
    
    if (!is.null(res$combined) && nrow(res$combined) > 0)
      all_summaries[[paste0(comp$name, "_", safe_name)]] <- res$combined
  }
}

# =============================================================
# 3. Grand summary
# =============================================================

grand_summary <- bind_rows(all_summaries)
if (nrow(grand_summary) > 0) {
  write_csv(grand_summary,
            file.path(GSEA_OUT, "subtype_GSEA_grand_summary.csv"))
  
  n_act <- sum(grand_summary$direction == "Activated")
  n_sup <- sum(grand_summary$direction == "Suppressed")
  cat("\nGrand summary: ", nrow(grand_summary), " terms (",
      n_act, " activated, ", n_sup, " suppressed)\n")
}

cat("\n\n==================== SUBTYPE GSEA COMPLETE ====================\n")
cat("Outputs in: ", GSEA_OUT, "\n")
cat("  <comp>/<subtype>_GSEA_<DB>.csv              — per database\n")
cat("  <comp>/<subtype>_GSEA_ALL_databases.csv      — combined with direction\n")
cat("  <comp>/<subtype>_GSEA_activated/suppressed.csv — split by NES sign\n")
cat("  <comp>/<subtype>_GSEA_nes_dot/bar.png        — balanced top-N plots\n")
cat("  <comp>/<subtype>_syngo_mapping_report.csv    — SynGO QC\n")
cat("  subtype_GSEA_grand_summary.csv               — everything combined\n")