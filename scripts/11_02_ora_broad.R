# ============================================================
# 11_02_ora_broad.R
#
# Direction-aware ORA on broad cell-type level DE results.
#
# Outputs: ./result/ora_results/broad/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/11_00_pathway_functions.R")

ORA_OUT <- "./result/ora_results/broad"
dir.create(ORA_OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. 7m broad Q1/Q2/Q3
# =============================================================

message("\n========== Broad ORA: 7m Q1/Q2/Q3 ==========\n")

broad_comps <- list(
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

for (comp in broad_comps) {
  comp_dir <- file.path(DE_BASE_DIR, comp$name)
  out_dir  <- file.path(ORA_OUT, comp$name)
  
  message("\n--- ", comp$name, " ---")
  
  for (bt in comp$types) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
    filepath  <- file.path(comp_dir,
                           paste0(comp$name, "_", safe_name,
                                  "_full_results.csv"))
    
    res <- run_ora_on_de_file(
      filepath    = filepath,
      label       = bt,
      comparison  = comp$name,
      out_dir     = out_dir,
      test_cond   = comp$test,
      ref_cond    = comp$ref
    )
    
    if (!is.null(res$summary) && nrow(res$summary) > 0)
      all_summaries[[paste0(comp$name, "_", safe_name)]] <- res$summary
  }
}

# =============================================================
# 2. 15m broad Q4a-Q4d
# =============================================================

message("\n\n========== Broad ORA: 15m Q4a-Q4d ==========\n")

broad_15m_comps <- list(
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

for (comp in broad_15m_comps) {
  comp_dir <- file.path(DE_BASE_DIR, comp$name)
  out_dir  <- file.path(ORA_OUT, comp$name)
  
  message("\n--- ", comp$name, " ---")
  
  for (bt in comp$types) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
    filepath  <- file.path(comp_dir,
                           paste0(comp$name, "_", safe_name,
                                  "_full_results.csv"))
    
    res <- run_ora_on_de_file(
      filepath    = filepath,
      label       = bt,
      comparison  = comp$name,
      out_dir     = out_dir,
      test_cond   = comp$test,
      ref_cond    = comp$ref
    )
    
    if (!is.null(res$summary) && nrow(res$summary) > 0)
      all_summaries[[paste0(comp$name, "_", safe_name)]] <- res$summary
  }
}

# =============================================================
# 3. Broad cross-comparison categories
# =============================================================

message("\n\n========== Broad ORA: Categories ==========\n")

CROSS_BROAD_DIR <- file.path(DE_BASE_DIR, "cross_comparison_broad")
CAT_ORA_OUT     <- file.path(ORA_OUT, "categories")

for (bt in BROAD_TYPES_7m) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
  cat_file  <- file.path(CROSS_BROAD_DIR,
                         paste0(safe_name, "_gene_categories.csv"))
  
  if (!file.exists(cat_file)) next
  
  run_ora_for_categories(
    categories_file = cat_file,
    label           = bt,
    out_dir         = CAT_ORA_OUT,
    comparison      = "cross_comparison_broad"
  )
}

# =============================================================
# 4. Broad regional categories
# =============================================================

message("\n\n========== Broad ORA: Regional ==========\n")

REG_BROAD_DIR <- file.path(CROSS_BROAD_DIR, "regional")
REG_ORA_OUT   <- file.path(ORA_OUT, "regional_categories")

for (bt in BROAD_TYPES_15m) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
  reg_file  <- file.path(REG_BROAD_DIR,
                         paste0(safe_name, "_regional_categories.csv"))
  
  if (!file.exists(reg_file)) next
  
  run_ora_for_regional(
    categories_file = reg_file,
    label           = bt,
    out_dir         = REG_ORA_OUT,
    comparison      = "regional_broad"
  )
}

# =============================================================
# 5. Grand summary
# =============================================================

grand_summary <- bind_rows(all_summaries)
if (nrow(grand_summary) > 0) {
  write_csv(grand_summary,
            file.path(ORA_OUT, "broad_ORA_grand_summary.csv"))
}

cat("\n\n==================== BROAD ORA COMPLETE ====================\n")
cat("Outputs in: ", ORA_OUT, "\n")
cat("  <comp>/<broad_type>_<dir>_<DB>.csv        — per direction per database\n")
cat("  <comp>/<broad_type>_ORA_ALL_databases.csv  — combined summary\n")
cat("  <comp>/<broad_type>_ORA_direction_plot.png  — two-panel plot\n")
cat("  broad_ORA_grand_summary.csv                — all comparisons combined\n")