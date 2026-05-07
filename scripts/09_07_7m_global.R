# ============================================================
# 09_07_7m_global.R
#
# Whole-population (global) pseudobulk DE for 7m.
# Pools ALL cells per sample — ignores cell type entirely.
# Useful as a reference for cross-comparison with broad/subtype
# DE in scripts 10_03 and 10_04.
#
# Design: aggregate all cells per sample, run DESeq2 same as
# broad-type pseudobulk but with one pseudo-cluster "Global".
#
# Outputs: ./result/de_results/Q{1,2,3}_global_*/
# ============================================================

source("./scripts/09_00_de_functions.R")

# --- Optional threshold overrides ---
# COMPARISONS$Q1_global$padj_cutoff <- 0.01
# COMPARISONS$Q2_global$lfc_cutoff  <- 0.5

# --- Load + add global_type column ---
obj <- load_annotated("7m")
obj$global_type <- "Global"

# --- Run ---
res_Q1g <- run_comparison(obj, COMPARISONS$Q1_global)
cat("\n\nQ1_global complete. Results in: ", res_Q1g$out_dir, "\n")

res_Q2g <- run_comparison(obj, COMPARISONS$Q2_global)
cat("\n\nQ2_global complete. Results in: ", res_Q2g$out_dir, "\n")

res_Q3g <- run_comparison(obj, COMPARISONS$Q3_global)
cat("\n\nQ3_global complete. Results in: ", res_Q3g$out_dir, "\n")