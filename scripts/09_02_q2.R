# 09_02_Q2_WT5XFAD_vs_Trem2KO5XFAD.R
# 7m pseudobulk: Trem2 KO effect in 5XFAD background (subtype-level)
# See 09_01 for threshold override examples.

source("./scripts/09_00_de_functions.R")

# --- Optional overrides ---
# COMPARISONS$Q2$padj_cutoff <- 0.01
# COMPARISONS$Q2$lfc_cutoff  <- 0.5
# COMPARISONS$Q2$subtype_thresholds <- list(
#   "Micro-PVM" = list(padj = 0.01, lfc = 0.5)
# )

obj <- load_annotated("7m")
res_Q2 <- run_comparison(obj, COMPARISONS$Q2)
cat("\n\nQ2 complete. Results in: ", res_Q2$out_dir, "\n")