# ============================================================
# 12_02_gsea_broad.R
#
# Direction-aware GSEA on broad cell-type level DE results.
# Q1-Q3 broad (7m) + Q4a-Q4d broad (15m).
#
# Outputs: ./result/gsea_results/broad/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/12_00_gsea_functions.R")

GSEA_OUT <- "./result/gsea_results/broad"
dir.create(GSEA_OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. 7m broad GSEA: Q1/Q2/Q3
# =============================================================

message("\n========== Broad GSEA: 7m Q1/Q2/Q3 ==========\n")

broad_7m <- list(
  list(name = "Q1_broad_WT_vs_WT5XFAD",
       test = "WT_5XFAD", ref = "WT",
       types = BROAD_TYPES_7m),
  list(name = "Q2_broad_WT5XFAD_vs_Trem2KO5XFAD",
       test = "Trem2_KO_5XFAD", ref = "WT_5XFAD",
       types = BROAD_TYPES_7m),
  list(name = "Q3_broad_WT_vs_Trem2KO",
       test = "Trem2_KO", ref = "WT",
       types = BROAD_TYPES_7m)
)

all_summaries <- list()

for (comp in broad_7m) {
  comp_dir <- file.path(DE_BASE_DIR, comp$name)
  out_dir  <- file.path(GSEA_OUT, comp$name)
  
  message("\n--- ", comp$name, " ---")
  
  for (bt in comp$types) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
    filepath  <- file.path(comp_dir,
                           paste0(comp$name, "_", safe_name,
                                  "_full_results.csv"))
    
    res <- run_gsea_on_de_file(
      filepath   = filepath,
      label      = bt,
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
# 2. 15m broad GSEA: Q4a-Q4d
# =============================================================

message("\n\n========== Broad GSEA: 15m Q4a-Q4d ==========\n")

broad_15m <- list(
  list(name = "Q4a_broad_15m_WT_Cor_vs_Hip",
       test = "WT_Cor", ref = "WT_Hip",
       types = BROAD_TYPES_15m),
  list(name = "Q4b_broad_15m_WT5XFAD_Cor_vs_Hip",
       test = "WT_5XFAD_Cor", ref = "WT_5XFAD_Hip",
       types = BROAD_TYPES_15m),
  list(name = "Q4c_broad_15m_Trem2KO_Cor_vs_Hip",
       test = "Trem2_KO_Cor", ref = "Trem2_KO_Hip",
       types = BROAD_TYPES_15m),
  list(name = "Q4d_broad_15m_Trem2KO5XFAD_Cor_vs_Hip",
       test = "Trem2_KO_5XFAD_Cor", ref = "Trem2_KO_5XFAD_Hip",
       types = BROAD_TYPES_15m)
)

for (comp in broad_15m) {
  comp_dir <- file.path(DE_BASE_DIR, comp$name)
  out_dir  <- file.path(GSEA_OUT, comp$name)
  
  message("\n--- ", comp$name, " ---")
  
  for (bt in comp$types) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
    filepath  <- file.path(comp_dir,
                           paste0(comp$name, "_", safe_name,
                                  "_full_results.csv"))
    
    res <- run_gsea_on_de_file(
      filepath   = filepath,
      label      = bt,
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
            file.path(GSEA_OUT, "broad_GSEA_grand_summary.csv"))
  
  n_act <- sum(grand_summary$direction == "Activated")
  n_sup <- sum(grand_summary$direction == "Suppressed")
  cat("\nGrand summary: ", nrow(grand_summary), " terms (",
      n_act, " activated, ", n_sup, " suppressed)\n")
}

cat("\n\n==================== BROAD GSEA COMPLETE ====================\n")
cat("Outputs in: ", GSEA_OUT, "\n")
cat("  <comp>/<broad_type>_GSEA_<DB>.csv              — per database\n")
cat("  <comp>/<broad_type>_GSEA_ALL_databases.csv      — combined with direction\n")
cat("  <comp>/<broad_type>_GSEA_activated/suppressed.csv\n")
cat("  <comp>/<broad_type>_GSEA_nes_dot/bar.png        — balanced plots\n")
cat("  broad_GSEA_grand_summary.csv                    — everything combined\n")