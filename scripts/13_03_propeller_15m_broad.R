# ============================================================
# 13_03_propeller_15m_broad.R
#
# Compositional analysis on 15m broad cell types.
#
# DESIGN NOTE: 15m has n=1 per region per condition. Splitting
# by region (Cor vs Hip) gives zero degrees of freedom for
# propeller. Two approaches here:
#
#   Approach A (main): Pool Cor + Hip within each condition so
#                       n=2 per condition, treat region as a
#                       covariate in the formula.
#
#   Approach B (exploratory per-region): Q4a-Q4d comparisons use
#                       a single sample per region per condition
#                       — clearly underpowered. Reported as
#                       "EXPLORATORY — no replicates" in output.
#
# Outputs: ./result/abundance/15m/broad/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/13_00_propeller_functions.R")

OUT_DIR <- file.path(ABUNDANCE_BASE_DIR, "15m", "broad")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. Load 15m object
# =============================================================

obj <- load_annotated("15m")

# Base metadata
meta_base <- data.frame(
  sample    = obj$sample_id,
  condition = obj$condition,
  cluster   = obj$broad_type,
  stringsAsFactors = FALSE
) %>% filter(!is.na(cluster), cluster != "Unknown")

# Extract region from sample name
meta_base$region <- ifelse(grepl("_Cor$", meta_base$sample), "Cor", "Hip")

# "Base" condition without region suffix (WT_5XFAD_Cor -> WT_5XFAD)
meta_base$base_condition <- gsub("_(Cor|Hip)$", "", meta_base$sample)
meta_base$base_condition <- gsub("_\\d+$", "", meta_base$base_condition)

cat("15m broad compositional analysis:\n")
cat("  Samples: ", n_distinct(meta_base$sample),
    " | Cells: ", nrow(meta_base),
    " | Broad types: ", n_distinct(meta_base$cluster), "\n\n")

# =============================================================
# 2. Approach A: Pooled across regions (recommended)
# =============================================================
# Use `condition` as the core grouping (pools Cor+Hip), with
# region as a covariate in the design.

message("\n########## APPROACH A: POOLED (region as covariate) ##########\n")

pooled_meta <- meta_base %>%
  mutate(condition = gsub("_(Cor|Hip)$", "", condition))
# Now conditions are: WT, WT_5XFAD, Trem2_KO, Trem2_KO_5XFAD (each n=2)

POOLED_DIR <- file.path(OUT_DIR, "pooled")
dir.create(POOLED_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper that passes `region` as covariate
run_pooled_comparison <- function(test_cond, ref_cond, comparison_name) {
  sub_meta <- pooled_meta %>% filter(condition %in% c(test_cond, ref_cond))
  
  if (n_distinct(sub_meta$sample) < 3) {
    message("\n=== ", comparison_name, " skipped (too few samples) ===")
    return(NULL)
  }
  
  message("\n=== ", comparison_name, ": ", test_cond, " vs ", ref_cond,
          " (pooled, region covariate) ===")
  
  prop_data <- calculate_proportions_per_sample(sub_meta)
  
  # Add region info to sample_info for the design
  region_lookup <- sub_meta %>% distinct(sample, region)
  prop_data$sample_info <- prop_data$sample_info %>%
    left_join(region_lookup, by = "sample")
  
  # Use formula with region covariate
  prop_res <- tryCatch({
    run_propeller(prop_data$counts, prop_data$sample_info,
                  formula_str = "~ region + condition",
                  reference_cond = ref_cond)
  }, error = function(e) {
    message("  propeller failed: ", e$message)
    NULL
  })
  
  if (is.null(prop_res)) return(NULL)
  
  summary <- make_summary_table(prop_res, prop_data$props,
                                prop_data$sample_info,
                                comparison_name, test_cond, ref_cond)
  write_csv(summary,
            file.path(POOLED_DIR,
                      paste0(comparison_name, "_propeller_results.csv")))
  
  p1 <- plot_proportions_per_sample(
    prop_data$props, prop_data$sample_info,
    title = paste0(comparison_name, ": pooled per-sample proportions")
  )
  n_c <- ncol(prop_data$props); nc <- min(4, n_c); nr <- ceiling(n_c / nc)
  ggsave(file.path(POOLED_DIR,
                   paste0(comparison_name, "_per_sample_proportions.png")),
         plot = p1, width = 3 + nc * 2.2, height = 2.5 + nr * 2.2, dpi = 200)
  
  p2 <- plot_proportions_boxplot(prop_data$props, prop_data$sample_info,
                                 title = paste0(comparison_name, ": pooled boxplot"))
  ggsave(file.path(POOLED_DIR,
                   paste0(comparison_name, "_boxplot.png")),
         plot = p2, width = 3 + nc * 2.2, height = 2.5 + nr * 2.2, dpi = 200)
  
  summary
}

pooled_summaries <- list()

pooled_summaries[["P1"]] <- run_pooled_comparison(
  "WT_5XFAD", "WT", "P1_WT_vs_WT5XFAD"
)
pooled_summaries[["P2"]] <- run_pooled_comparison(
  "Trem2_KO_5XFAD", "WT_5XFAD", "P2_WT5XFAD_vs_Trem2KO5XFAD"
)
pooled_summaries[["P3"]] <- run_pooled_comparison(
  "Trem2_KO", "WT", "P3_WT_vs_Trem2KO"
)

pooled_grand <- bind_rows(pooled_summaries)
if (nrow(pooled_grand) > 0)
  write_csv(pooled_grand,
            file.path(POOLED_DIR, "summary_all_pooled.csv"))

# =============================================================
# 3. Approach B: Regional Q4 (EXPLORATORY, n=1 per region)
# =============================================================

message("\n########## APPROACH B: REGIONAL Q4a-d (EXPLORATORY) ##########\n")
message("WARNING: n=1 per region per condition — p-values unreliable\n")

REG_DIR <- file.path(OUT_DIR, "regional_q4")
dir.create(REG_DIR, recursive = TRUE, showWarnings = FALSE)

q4_comps <- list(
  list(test = "WT_Cor",             ref = "WT_Hip",             name = "Q4a_WT_Cor_vs_Hip"),
  list(test = "WT_5XFAD_Cor",       ref = "WT_5XFAD_Hip",       name = "Q4b_WT5XFAD_Cor_vs_Hip"),
  list(test = "Trem2_KO_Cor",       ref = "Trem2_KO_Hip",       name = "Q4c_Trem2KO_Cor_vs_Hip"),
  list(test = "Trem2_KO_5XFAD_Cor", ref = "Trem2_KO_5XFAD_Hip", name = "Q4d_Trem2KO5XFAD_Cor_vs_Hip")
)

# For n=1 per group, propeller won't produce meaningful stats,
# but we still save the proportion tables for visual inspection.

q4_summaries <- list()

for (comp in q4_comps) {
  sub_meta <- meta_base %>% filter(condition %in% c(comp$test, comp$ref))
  
  if (n_distinct(sub_meta$sample) < 2) next
  
  message("\n--- ", comp$name, " (EXPLORATORY) ---")
  
  prop_data <- calculate_proportions_per_sample(sub_meta)
  
  # Save the proportions directly — no meaningful p-value with n=1
  summary <- data.frame(
    cluster = colnames(prop_data$props),
    stringsAsFactors = FALSE
  )
  for (s in rownames(prop_data$props)) {
    summary[[paste0("prop_", s)]] <- round(prop_data$props[s, ], 5)
  }
  summary$comparison <- comp$name
  summary$test_cond  <- comp$test
  summary$ref_cond   <- comp$ref
  summary$note       <- "EXPLORATORY: n=1 per region, no p-value"
  
  write_csv(summary,
            file.path(REG_DIR, paste0(comp$name, "_proportions.csv")))
  
  p <- plot_proportions_per_sample(prop_data$props, prop_data$sample_info,
                                   title = paste0(comp$name,
                                                  " (EXPLORATORY n=1)"))
  n_c <- ncol(prop_data$props); nc <- min(4, n_c); nr <- ceiling(n_c / nc)
  ggsave(file.path(REG_DIR, paste0(comp$name, "_proportions.png")),
         plot = p, width = 3 + nc * 2.2, height = 2.5 + nr * 2.2, dpi = 200)
  
  q4_summaries[[comp$name]] <- summary
}

if (length(q4_summaries) > 0)
  write_csv(bind_rows(q4_summaries),
            file.path(REG_DIR, "summary_all_regional.csv"))

# =============================================================
# 4. Final summary
# =============================================================

cat("\n\n==================== 15m BROAD COMPLETE ====================\n")
cat("Primary (pooled): ", POOLED_DIR, "\n")
cat("  P1/P2/P3_propeller_results.csv — with region as covariate\n")
cat("Exploratory (regional): ", REG_DIR, "\n")
cat("  Q4a-Q4d_proportions.csv — descriptive only, no p-values\n")