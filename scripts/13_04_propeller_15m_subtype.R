# ============================================================
# 13_04_propeller_15m_subtype.R
#
# Compositional analysis on 15m neuronal subtypes.
# Uses the pooled approach (Cor+Hip per condition, region as
# covariate) from 13_03. Excitatory and inhibitory tested
# separately.
#
# Outputs: ./result/abundance/15m/subtype/{excitatory,inhibitory}/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/13_00_propeller_functions.R")

# =============================================================
# 1. Load + classify subtypes
# =============================================================

obj <- load_annotated("15m")

meta_full <- data.frame(
  sample    = obj$sample_id,
  condition = obj$condition,
  cluster   = obj$allen_subclass_filtered,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(cluster), cluster != "", cluster != "Unknown")

meta_full$region <- ifelse(grepl("_Cor$", meta_full$sample), "Cor", "Hip")

# Pool condition (strip region suffix)
meta_full$condition <- gsub("_(Cor|Hip)$", "", meta_full$condition)

# Class lookup from Seurat metadata
class_lookup <- obj@meta.data %>%
  distinct(allen_subclass_filtered, allen_class_label) %>%
  filter(!is.na(allen_subclass_filtered)) %>%
  rename(cluster = allen_subclass_filtered, class = allen_class_label)

meta_full <- meta_full %>% left_join(class_lookup, by = "cluster")

excitatory_subtypes <- class_lookup %>%
  filter(class == "Glutamatergic") %>% pull(cluster) %>%
  intersect(DE_READY_15m_ALL)

inhibitory_subtypes <- class_lookup %>%
  filter(class == "GABAergic") %>% pull(cluster) %>%
  intersect(DE_READY_15m_ALL)

cat("15m subtype compositional analysis (pooled):\n")
cat("  Excitatory (DE-ready): ", length(excitatory_subtypes), "\n")
cat("  Inhibitory (DE-ready): ", length(inhibitory_subtypes), "\n\n")

# =============================================================
# 2. Helper: pooled propeller with region covariate
# =============================================================

run_pooled_subtype_comparison <- function(meta, test_cond, ref_cond,
                                          comparison_name, out_dir) {
  
  sub_meta <- meta %>% filter(condition %in% c(test_cond, ref_cond))
  if (n_distinct(sub_meta$sample) < 3) {
    message("\n=== ", comparison_name, " skipped (too few samples) ===")
    return(NULL)
  }
  
  message("\n=== ", comparison_name, ": ", test_cond, " vs ", ref_cond, " ===")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  prop_data <- calculate_proportions_per_sample(sub_meta)
  
  # Filter low-count subtypes
  keep <- filter_subtypes_for_abundance(prop_data$counts)
  if (length(keep) < 2) {
    message("  Too few subtypes after filtering — skipping")
    return(NULL)
  }
  prop_data$counts <- prop_data$counts[, keep, drop = FALSE]
  prop_data$props  <- prop_data$props[, keep, drop = FALSE]
  
  # Add region to sample_info
  region_lookup <- sub_meta %>% distinct(sample, region)
  prop_data$sample_info <- prop_data$sample_info %>%
    left_join(region_lookup, by = "sample")
  
  prop_res <- tryCatch({
    run_propeller(prop_data$counts, prop_data$sample_info,
                  formula_str = "~ region + condition",
                  reference_cond = ref_cond)
  }, error = function(e) {
    message("  propeller failed: ", e$message); NULL
  })
  
  if (is.null(prop_res)) return(NULL)
  
  summary <- make_summary_table(prop_res, prop_data$props,
                                prop_data$sample_info,
                                comparison_name, test_cond, ref_cond)
  write_csv(summary,
            file.path(out_dir,
                      paste0(comparison_name, "_propeller_results.csv")))
  
  p1 <- plot_proportions_per_sample(
    prop_data$props, prop_data$sample_info,
    title = paste0(comparison_name, ": per-sample")
  )
  n_c <- ncol(prop_data$props); nc <- min(4, n_c); nr <- ceiling(n_c / nc)
  ggsave(file.path(out_dir,
                   paste0(comparison_name, "_per_sample_proportions.png")),
         plot = p1, width = 3 + nc * 2.2, height = 2.5 + nr * 2.2, dpi = 200)
  
  p2 <- plot_proportions_boxplot(prop_data$props, prop_data$sample_info,
                                 title = paste0(comparison_name, ": boxplot"))
  ggsave(file.path(out_dir, paste0(comparison_name, "_boxplot.png")),
         plot = p2, width = 3 + nc * 2.2, height = 2.5 + nr * 2.2, dpi = 200)
  
  summary
}

# =============================================================
# 3. Run for both neuronal classes
# =============================================================

run_class <- function(meta, out_dir, label) {
  cat("\n\n########## ", label, " ##########\n")
  
  summaries <- list()
  summaries[["P1"]] <- run_pooled_subtype_comparison(
    meta, "WT_5XFAD", "WT", "P1_WT_vs_WT5XFAD", out_dir
  )
  summaries[["P2"]] <- run_pooled_subtype_comparison(
    meta, "Trem2_KO_5XFAD", "WT_5XFAD",
    "P2_WT5XFAD_vs_Trem2KO5XFAD", out_dir
  )
  summaries[["P3"]] <- run_pooled_subtype_comparison(
    meta, "Trem2_KO", "WT", "P3_WT_vs_Trem2KO", out_dir
  )
  
  grand <- bind_rows(summaries)
  if (nrow(grand) > 0) {
    write_csv(grand, file.path(out_dir, "summary_all_comparisons.csv"))
    cat("\n=== ", label, " summary ===\n")
    sig <- grand %>% filter(FDR < 0.05 | fdr < 0.05)
    cat("Total: ", nrow(grand), " | Sig (FDR<0.05): ", nrow(sig), "\n")
    if (nrow(sig) > 0) print(as.data.frame(sig), row.names = FALSE)
  }
  grand
}

# Excitatory
meta_exc <- meta_full %>% filter(cluster %in% excitatory_subtypes)
out_exc  <- file.path(ABUNDANCE_BASE_DIR, "15m", "subtype", "excitatory")
res_exc  <- run_class(meta_exc, out_exc, "EXCITATORY subtypes (pooled)")

# Inhibitory
meta_inh <- meta_full %>% filter(cluster %in% inhibitory_subtypes)
out_inh  <- file.path(ABUNDANCE_BASE_DIR, "15m", "subtype", "inhibitory")
res_inh  <- run_class(meta_inh, out_inh, "INHIBITORY subtypes (pooled)")

cat("\n\n==================== 15m SUBTYPE COMPLETE ====================\n")
cat("Outputs in: ", file.path(ABUNDANCE_BASE_DIR, "15m", "subtype"), "\n")