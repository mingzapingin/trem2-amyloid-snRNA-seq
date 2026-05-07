# 09_06_Q4_broad_15m_regional.R
# 15m cell-level DE at broad cell-type level (6 cell classes).
# Same 4 regional comparisons (Cor vs Hip) as Q4a-Q4d but pooling
# all cells per broad type. EXPLORATORY — no biological replicates.
# Q4a: WT | Q4b: WT_5XFAD | Q4c: Trem2_KO | Q4d: Trem2_KO_5XFAD
# See 09_01 for threshold override examples.

source("./scripts/09_00_de_functions.R")

# --- Optional overrides ---
# COMPARISONS$Q4a_broad$padj_cutoff <- 0.01
# COMPARISONS$Q4c_broad$subtype_thresholds <- list(
#   "Microglia" = list(padj = 0.01, lfc = 0.5)
# )

obj <- load_annotated("15m")

res_Q4ab <- run_comparison(obj, COMPARISONS$Q4a_broad)
cat("\n\nQ4a_broad complete. Results in: ", res_Q4ab$out_dir, "\n")

res_Q4bb <- run_comparison(obj, COMPARISONS$Q4b_broad)
cat("\n\nQ4b_broad complete. Results in: ", res_Q4bb$out_dir, "\n")

res_Q4cb <- run_comparison(obj, COMPARISONS$Q4c_broad)
cat("\n\nQ4c_broad complete. Results in: ", res_Q4cb$out_dir, "\n")

res_Q4db <- run_comparison(obj, COMPARISONS$Q4d_broad)
cat("\n\nQ4d_broad complete. Results in: ", res_Q4db$out_dir, "\n")