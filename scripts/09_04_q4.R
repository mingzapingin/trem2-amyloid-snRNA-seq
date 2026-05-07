# 09_04_Q4_15m_regional.R
# 15m cell-level DE: regional comparisons (subtype-level)
# Q4a: WT Cor vs Hip              — baseline regional differences
# Q4b: WT_5XFAD Cor vs Hip       — regional response to amyloid
# Q4c: Trem2_KO Cor vs Hip       — regional diff without amyloid
# Q4d: Trem2_KO_5XFAD Cor vs Hip — regional diff with amyloid + KO
#
# All EXPLORATORY — no biological replicates.
# All four are Cor vs Hip, one per genotype×amyloid condition.
# See 09_01 for threshold override examples.

source("./scripts/09_00_de_functions.R")

# --- Optional: per-comparison overrides ---
# COMPARISONS$Q4a$padj_cutoff <- 0.01
# COMPARISONS$Q4c$subtype_thresholds <- list(
#   "Micro-PVM" = list(padj = 0.01, lfc = 0.5)
# )

obj <- load_annotated("15m")

res_Q4a <- run_comparison(obj, COMPARISONS$Q4a)
cat("\n\nQ4a complete. Results in: ", res_Q4a$out_dir, "\n")

res_Q4b <- run_comparison(obj, COMPARISONS$Q4b)
cat("\n\nQ4b complete. Results in: ", res_Q4b$out_dir, "\n")

res_Q4c <- run_comparison(obj, COMPARISONS$Q4c)
cat("\n\nQ4c complete. Results in: ", res_Q4c$out_dir, "\n")

res_Q4d <- run_comparison(obj, COMPARISONS$Q4d)
cat("\n\nQ4d complete. Results in: ", res_Q4d$out_dir, "\n")