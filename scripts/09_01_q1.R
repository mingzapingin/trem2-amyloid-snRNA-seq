# ============================================================
# 09_01_Q1_WT_vs_WT5XFAD.R
#
# 7m pseudobulk DE: WT_5XFAD vs WT (subtype-level)
# Question: What is the effect of 5XFAD amyloid pathology?
#
# Threshold hierarchy (edit below as needed):
#   Global defaults:  padj < 0.05, |lfc| >= 0.25  (in 09_00)
#   Comparison-level: override for ALL subtypes in Q1
#   Subtype-level:    override for specific subtypes only
# ============================================================

source("./scripts/09_00_de_functions.R")

# --- Optional: comparison-level override ---
# Uncomment to tighten/relax ALL subtypes in Q1:
# COMPARISONS$Q1$padj_cutoff <- 0.05
# COMPARISONS$Q1$lfc_cutoff  <- 0.25

# --- Optional: per-subtype overrides ---
# More stringent for large excitatory populations,
# more relaxed for small/rare subtypes:
# COMPARISONS$Q1$subtype_thresholds <- list(
#   "L2/3 IT CTX" = list(padj = 0.01, lfc = 0.5),
#   "L4/5 IT CTX" = list(padj = 0.01, lfc = 0.5),
#   "L6 CT CTX"   = list(padj = 0.01, lfc = 0.5),
#   "Micro-PVM"   = list(padj = 0.01, lfc = 0.25),
#   "Sncg"        = list(padj = 0.1,  lfc = 0.15),
#   "L3 IT ENT"   = list(padj = 0.1,  lfc = 0.15)
# )

# --- Load + Run ---
obj <- load_annotated("7m")
res_Q1 <- run_comparison(obj, COMPARISONS$Q1)
cat("\n\nQ1 complete. Results in: ", res_Q1$out_dir, "\n")