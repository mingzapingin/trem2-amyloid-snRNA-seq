# ============================================================
# 13_02_propeller_7m_subtype.R
#
# Compositional analysis on 7m neuronal subtypes.
# Runs excitatory and inhibitory subtypes SEPARATELY to avoid
# mixing abundance ranges (propeller's empirical Bayes moderation
# works better on groups with similar variance structure).
#
# Only DE-ready subtypes are tested — consistency with DE analysis.
#
# Outputs: ./result/abundance/7m/subtype/{excitatory,inhibitory}/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/13_00_propeller_functions.R")

# =============================================================
# 1. Load + classify subtypes into excitatory vs inhibitory
# =============================================================

obj <- load_annotated("7m")
meta_full <- build_meta_from_seurat(obj, cluster_col = "allen_subclass_filtered")

# Classify subtypes using Allen class labels
# (The Seurat object should have `allen_class_label` from label transfer)
class_lookup <- obj@meta.data %>%
  distinct(allen_subclass_filtered, allen_class_label) %>%
  filter(!is.na(allen_subclass_filtered)) %>%
  rename(cluster = allen_subclass_filtered, class = allen_class_label)

meta_full <- meta_full %>%
  left_join(class_lookup, by = "cluster")

excitatory_subtypes <- class_lookup %>%
  filter(class == "Glutamatergic") %>% pull(cluster) %>%
  intersect(DE_READY_7m_ALL)

inhibitory_subtypes <- class_lookup %>%
  filter(class == "GABAergic") %>% pull(cluster) %>%
  intersect(DE_READY_7m_ALL)

cat("7m subtype compositional analysis:\n")
cat("  Excitatory (DE-ready): ", length(excitatory_subtypes), "\n")
cat("  Inhibitory (DE-ready): ", length(inhibitory_subtypes), "\n\n")

# =============================================================
# 2. Helper: run all three comparisons on a subtype set
# =============================================================

run_all_comparisons <- function(meta, out_dir, label) {
  
  cat("\n\n########## ", label, " ##########\n")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  summaries <- list()
  
  summaries[["Q1"]] <- run_propeller_comparison(
    meta, "WT_5XFAD", "WT",
    "Q1_WT_vs_WT5XFAD", out_dir, filter_min_cells = TRUE
  )
  summaries[["Q2"]] <- run_propeller_comparison(
    meta, "Trem2_KO_5XFAD", "WT_5XFAD",
    "Q2_WT5XFAD_vs_Trem2KO5XFAD", out_dir, filter_min_cells = TRUE
  )
  summaries[["Q3"]] <- run_propeller_comparison(
    meta, "Trem2_KO", "WT",
    "Q3_WT_vs_Trem2KO", out_dir, filter_min_cells = TRUE
  )
  
  grand <- bind_rows(summaries)
  if (nrow(grand) > 0) {
    write_csv(grand, file.path(out_dir, "summary_all_comparisons.csv"))
    cat("\n=== ", label, " summary ===\n")
    sig <- grand %>% filter(get_fdr_values(.) < 0.05)
    cat("Total: ", nrow(grand), " | Sig (FDR<0.05): ", nrow(sig), "\n")
    if (nrow(sig) > 0) print(as.data.frame(sig), row.names = FALSE)
  }
  
  grand
}

# =============================================================
# 3. Excitatory subtypes
# =============================================================

meta_exc <- meta_full %>% filter(cluster %in% excitatory_subtypes)
out_exc  <- file.path(ABUNDANCE_BASE_DIR, "7m", "subtype", "excitatory")
res_exc  <- run_all_comparisons(meta_exc, out_exc, "EXCITATORY subtypes")

# =============================================================
# 4. Inhibitory subtypes
# =============================================================

meta_inh <- meta_full %>% filter(cluster %in% inhibitory_subtypes)
out_inh  <- file.path(ABUNDANCE_BASE_DIR, "7m", "subtype", "inhibitory")
res_inh  <- run_all_comparisons(meta_inh, out_inh, "INHIBITORY subtypes")

cat("\n\n==================== 7m SUBTYPE COMPLETE ====================\n")
cat("Outputs in: ", file.path(ABUNDANCE_BASE_DIR, "7m", "subtype"), "\n")