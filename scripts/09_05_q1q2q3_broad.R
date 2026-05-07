# 09_05_Q1Q2Q3_broad.R
# 7m pseudobulk DE at broad cell-type level (6 cell classes).
# Same 3 comparisons as Q1-Q3 but pooling all cells per broad type.
# See 09_01 for threshold override examples.

source("./scripts/09_00_de_functions.R")

# --- Optional: broad-type specific overrides ---
# COMPARISONS$Q1_broad$subtype_thresholds <- list(
#   "Excitatory_neuron"    = list(padj = 0.01, lfc = 0.5),
#   "Microglia"            = list(padj = 0.01, lfc = 0.25),
#   "Endothelial_Pericyte" = list(padj = 0.1,  lfc = 0.15)
# )

obj <- load_annotated("7m")

res_Q1b <- run_comparison(obj, COMPARISONS$Q1_broad)
cat("\n\nQ1_broad complete. Results in: ", res_Q1b$out_dir, "\n")

res_Q2b <- run_comparison(obj, COMPARISONS$Q2_broad)
cat("\n\nQ2_broad complete. Results in: ", res_Q2b$out_dir, "\n")

res_Q3b <- run_comparison(obj, COMPARISONS$Q3_broad)
cat("\n\nQ3_broad complete. Results in: ", res_Q3b$out_dir, "\n")