# 09_03_Q3_WT_vs_Trem2KO.R
# 7m pseudobulk: Trem2 KO alone, no amyloid (subtype-level)
# See 09_01 for threshold override examples.

source("./scripts/09_00_de_functions.R")

# --- Optional overrides ---
# COMPARISONS$Q3$padj_cutoff <- 0.01
# COMPARISONS$Q3$subtype_thresholds <- list(...)

obj <- load_annotated("7m")
res_Q3 <- run_comparison(obj, COMPARISONS$Q3)
cat("\n\nQ3 complete. Results in: ", res_Q3$out_dir, "\n")