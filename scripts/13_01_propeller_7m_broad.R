# ============================================================
# 13_01_propeller_7m_broad.R
#
# Compositional analysis on 7m broad cell types.
# Runs Q1, Q2, Q3 — same comparisons as DE, but asking whether
# cell type ABUNDANCE changes.
#
# Outputs: ./result/abundance/7m/broad/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/13_00_propeller_functions.R")

OUT_DIR <- file.path(ABUNDANCE_BASE_DIR, "7m", "broad")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. Load 7m object and build metadata
# =============================================================

obj <- load_annotated("7m")
meta <- build_meta_from_seurat(obj, cluster_col = "broad_type")

cat("7m metadata: ", nrow(meta), " cells across ",
    n_distinct(meta$sample), " samples × ",
    n_distinct(meta$cluster), " broad types\n\n")

# =============================================================
# 2. Run three comparisons
# =============================================================

all_summaries <- list()

# Q1: Amyloid effect
res_Q1 <- run_propeller_comparison(
  meta            = meta,
  test_cond       = "WT_5XFAD",
  ref_cond        = "WT",
  comparison_name = "Q1_WT_vs_WT5XFAD",
  out_dir         = OUT_DIR,
  filter_min_cells = FALSE
)
if (!is.null(res_Q1)) all_summaries[["Q1"]] <- res_Q1

# Q2: Trem2 KO in amyloid
res_Q2 <- run_propeller_comparison(
  meta            = meta,
  test_cond       = "Trem2_KO_5XFAD",
  ref_cond        = "WT_5XFAD",
  comparison_name = "Q2_WT5XFAD_vs_Trem2KO5XFAD",
  out_dir         = OUT_DIR,
  filter_min_cells = FALSE
)
if (!is.null(res_Q2)) all_summaries[["Q2"]] <- res_Q2

# Q3: Trem2 KO alone
res_Q3 <- run_propeller_comparison(
  meta            = meta,
  test_cond       = "Trem2_KO",
  ref_cond        = "WT",
  comparison_name = "Q3_WT_vs_Trem2KO",
  out_dir         = OUT_DIR,
  filter_min_cells = FALSE
)
if (!is.null(res_Q3)) all_summaries[["Q3"]] <- res_Q3

# =============================================================
# 3. Combined summary across all comparisons
# =============================================================

grand <- bind_rows(all_summaries)
if (nrow(grand) > 0) {
  write_csv(grand, file.path(OUT_DIR, "summary_all_comparisons.csv"))
  cat("\n\n=== 7m Broad Compositional Summary ===\n")
  sig <- grand %>% filter(get_fdr_values(.) < 0.05)
  cat("Total tests: ", nrow(grand),
      " | Significant (FDR<0.05): ", nrow(sig), "\n")
  if (nrow(sig) > 0) print(as.data.frame(sig), row.names = FALSE)
}

cat("\nOutputs in: ", OUT_DIR, "\n")