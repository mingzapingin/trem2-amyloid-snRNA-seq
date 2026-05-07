# ============================================================
# 08_filter_merge_labels.R
#
# Post-transfer QC:
#   1. Load transferred objects from 07
#   2. Filter out low-confidence predictions
#   3. Summarize before/after cell counts
#   4. Merge transferred labels back into the full broad-typed
#      objects (so every cell has its Allen label or NA)
#   5. Report cell counts per subtype × condition × sample
#   6. Flag subtypes too small for pseudobulk DE
#
# Inputs:
#   ./processed/transfer_objects/seurat_<cohort>_<group>_transferred.rds
#   ./processed/seurat_7m_broad_typed.rds
#   ./processed/seurat_15m_broad_typed.rds
#
# Outputs:
#   ./processed/seurat_7m_final_annotated.rds
#   ./processed/seurat_15m_final_annotated.rds
#   ./result/transfer_qc/*.csv
#   ./result/transfer_qc/*.png
# ============================================================

library(Seurat)
library(SeuratObject)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(ggplot2)

set.seed(1234)

dir.create("./result/transfer_qc", recursive = TRUE, showWarnings = FALSE)

TRANSFER_DIR <- "./processed/transfer_objects"
OUT_DIR      <- "./result/transfer_qc"

# =============================================================
# 1. Settings
# =============================================================

# Minimum prediction score to keep a label
# Cells below this get allen_subclass_label set to NA
SCORE_THRESHOLD <- 0.9

# Minimum cells per subtype × sample to be usable for pseudobulk DE
# (DESeq2/edgeR need enough cells to aggregate into a reliable pseudobulk)
MIN_CELLS_PER_SAMPLE <- 10

# Minimum number of samples (replicates) per condition with enough cells
# to run pseudobulk DE for that subtype × condition
MIN_SAMPLES_PER_CONDITION <- 2

# =============================================================
# 2. Cohort + group registry
# =============================================================

cohort_inputs <- list(
  "7m"  = "./processed/seurat_7m_broad_typed.rds",
  "15m" = "./processed/seurat_15m_broad_typed.rds"
)

group_names <- c(
  "excitatory", "inhibitory", "microglia",
  "astrocyte", "oligodendrocyte", "endothelial_pericyte"
)

cohorts_to_run <- c("7m", "15m")

# =============================================================
# 3. Process one cohort
# =============================================================

process_cohort <- function(cohort_name) {
  
  message("\n\n========== ", cohort_name, " ==========")
  
  # --- Load full broad-typed object ---
  obj_full <- readRDS(cohort_inputs[[cohort_name]])
  n_total  <- ncol(obj_full)
  cat("Full object: ", n_total, " cells\n")
  
  # Initialize columns for merged labels
  obj_full$allen_subclass_label       <- NA_character_
  obj_full$allen_subclass_label_score <- NA_real_
  obj_full$allen_cluster_label        <- NA_character_
  obj_full$allen_cluster_label_score  <- NA_real_
  obj_full$allen_subclass_filtered    <- NA_character_  # post-filter label
  
  # --- Collect per-group summaries ---
  all_before_after <- list()
  all_subtype_detail <- list()
  
  for (group_name in group_names) {
    rds_path <- file.path(
      TRANSFER_DIR,
      paste0("seurat_", cohort_name, "_", group_name, "_transferred.rds")
    )
    
    if (!file.exists(rds_path)) {
      message("  ", group_name, ": transferred object not found — skipping")
      next
    }
    
    message("\n  --- ", group_name, " ---")
    obj_sub <- readRDS(rds_path)
    n_sub   <- ncol(obj_sub)
    
    # Check required columns
    if (!"allen_subclass_label" %in% colnames(obj_sub@meta.data)) {
      message("  No allen_subclass_label — skipping")
      next
    }
    
    md <- obj_sub@meta.data
    
    # -------------------------------------------------------
    # A. Before/after filtering summary
    # -------------------------------------------------------
    
    # Score distribution
    score_stats <- md %>%
      summarise(
        n_cells     = n(),
        mean_score  = round(mean(allen_subclass_label_score, na.rm = TRUE), 3),
        median_score = round(median(allen_subclass_label_score, na.rm = TRUE), 3),
        pct_above_threshold = round(
          100 * mean(allen_subclass_label_score >= SCORE_THRESHOLD, na.rm = TRUE), 1
        ),
        n_pass = sum(allen_subclass_label_score >= SCORE_THRESHOLD, na.rm = TRUE),
        n_fail = sum(allen_subclass_label_score < SCORE_THRESHOLD, na.rm = TRUE)
      ) %>%
      mutate(group = group_name)
    
    cat("    Total: ", score_stats$n_cells, " cells\n")
    cat("    Pass (score ≥ ", SCORE_THRESHOLD, "): ", score_stats$n_pass,
        " (", score_stats$pct_above_threshold, "%)\n")
    cat("    Fail: ", score_stats$n_fail, "\n")
    cat("    Median score: ", score_stats$median_score, "\n")
    
    # Per-subclass before/after
    subclass_ba <- md %>%
      group_by(allen_subclass_label) %>%
      summarise(
        n_before      = n(),
        mean_score    = round(mean(allen_subclass_label_score), 3),
        median_score  = round(median(allen_subclass_label_score), 3),
        n_after       = sum(allen_subclass_label_score >= SCORE_THRESHOLD),
        n_removed     = n_before - n_after,
        pct_kept      = round(100 * n_after / n_before, 1),
        .groups = "drop"
      ) %>%
      mutate(group = group_name) %>%
      arrange(desc(n_before))
    
    all_before_after[[group_name]] <- subclass_ba
    
    # -------------------------------------------------------
    # B. Apply filter: set filtered label
    # -------------------------------------------------------
    
    filtered_label <- ifelse(
      md$allen_subclass_label_score >= SCORE_THRESHOLD,
      md$allen_subclass_label,
      NA_character_
    )
    
    # -------------------------------------------------------
    # C. Merge into full object
    # -------------------------------------------------------
    
    cells_sub <- rownames(md)
    obj_full$allen_subclass_label[cells_sub]       <- md$allen_subclass_label
    obj_full$allen_subclass_label_score[cells_sub]  <- md$allen_subclass_label_score
    obj_full$allen_subclass_filtered[cells_sub]     <- filtered_label
    
    if ("allen_cluster_label" %in% colnames(md)) {
      obj_full$allen_cluster_label[cells_sub] <- md$allen_cluster_label
    }
    if ("allen_cluster_label_score" %in% colnames(md)) {
      obj_full$allen_cluster_label_score[cells_sub] <- md$allen_cluster_label_score
    }
    
    # -------------------------------------------------------
    # D. Detailed counts: subtype × condition × sample
    # -------------------------------------------------------
    
    # Use filtered labels
    md$allen_subclass_filtered <- filtered_label
    
    detail <- md %>%
      filter(!is.na(allen_subclass_filtered)) %>%
      group_by(allen_subclass_filtered) %>%
      mutate(subtype_total = n()) %>%
      ungroup()
    
    if ("condition" %in% colnames(md) && "orig.ident" %in% colnames(md)) {
      subtype_detail <- detail %>%
        count(allen_subclass_filtered, condition, orig.ident,
              name = "n_cells") %>%
        arrange(allen_subclass_filtered, condition, orig.ident) %>%
        mutate(group = group_name)
    } else {
      subtype_detail <- detail %>%
        count(allen_subclass_filtered, orig.ident,
              name = "n_cells") %>%
        arrange(allen_subclass_filtered, orig.ident) %>%
        mutate(group = group_name)
    }
    
    all_subtype_detail[[group_name]] <- subtype_detail
    
    rm(obj_sub, md); gc(verbose = FALSE)
  }
  
  # -------------------------------------------------------
  # 4. Save before/after summary
  # -------------------------------------------------------
  
  ba_combined <- bind_rows(all_before_after)
  if (nrow(ba_combined) > 0) {
    write_csv(ba_combined,
              file.path(OUT_DIR,
                        paste0(cohort_name, "_subclass_before_after_filter.csv")))
    
    cat("\n\n  === Before/After Filter Summary ===\n")
    ba_group <- ba_combined %>%
      group_by(group) %>%
      summarise(
        n_subtypes    = n(),
        total_before  = sum(n_before),
        total_after   = sum(n_after),
        total_removed = sum(n_removed),
        pct_kept      = round(100 * total_after / total_before, 1),
        .groups = "drop"
      )
    print(as.data.frame(ba_group), row.names = FALSE)
    write_csv(ba_group,
              file.path(OUT_DIR,
                        paste0(cohort_name, "_group_filter_summary.csv")))
    
    # --- Plot: score distribution per group ---
    # Rebuild from full object
    score_df <- obj_full@meta.data %>%
      filter(!is.na(allen_subclass_label_score)) %>%
      select(broad_type, allen_subclass_label, allen_subclass_label_score)
    
    if (nrow(score_df) > 0) {
      p_hist <- ggplot(score_df,
                       aes(x = allen_subclass_label_score)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white") +
        geom_vline(xintercept = SCORE_THRESHOLD,
                   linetype = "dashed", color = "red", linewidth = 0.8) +
        facet_wrap(~ broad_type, scales = "free_y") +
        theme_bw() +
        ggtitle(paste0(cohort_name,
                       " — prediction score distribution (threshold = ",
                       SCORE_THRESHOLD, ")")) +
        xlab("Prediction score") + ylab("Cells")
      ggsave(file.path(OUT_DIR,
                       paste0(cohort_name, "_score_distribution.png")),
             plot = p_hist, width = 12, height = 8, dpi = 200)
    }
  }
  
  # -------------------------------------------------------
  # 5. Detailed subtype × condition × sample table
  # -------------------------------------------------------
  
  detail_combined <- bind_rows(all_subtype_detail)
  if (nrow(detail_combined) > 0) {
    write_csv(detail_combined,
              file.path(OUT_DIR,
                        paste0(cohort_name,
                               "_subtype_by_condition_sample_detail.csv")))
    
    # --- Overall subtype counts (after filter) ---
    subtype_totals <- detail_combined %>%
      group_by(group, allen_subclass_filtered) %>%
      summarise(total_cells = sum(n_cells), .groups = "drop") %>%
      arrange(group, desc(total_cells))
    write_csv(subtype_totals,
              file.path(OUT_DIR,
                        paste0(cohort_name, "_subtype_total_counts.csv")))
    
    # --- Detailed wide table: subtype × sample ---
    subtype_by_sample_wide <- detail_combined %>%
      select(group, allen_subclass_filtered, orig.ident, n_cells) %>%
      pivot_wider(
        names_from  = orig.ident,
        values_from = n_cells,
        values_fill = 0
      ) %>%
      mutate(total = rowSums(across(where(is.numeric)))) %>%
      arrange(group, desc(total))
    
    write_csv(subtype_by_sample_wide,
              file.path(OUT_DIR,
                        paste0(cohort_name,
                               "_subtype_by_sample_wide.csv")))
    
    cat("\n\n  === Subtype × Sample (", cohort_name, ") ===\n")
    print(as.data.frame(subtype_by_sample_wide), row.names = FALSE)
    
    # --- Flag subtypes too small for pseudobulk ---
    if ("condition" %in% colnames(detail_combined)) {
      # For each subtype × condition, count how many samples have ≥ MIN_CELLS
      sample_adequacy <- detail_combined %>%
        mutate(adequate = n_cells >= MIN_CELLS_PER_SAMPLE) %>%
        group_by(group, allen_subclass_filtered, condition) %>%
        summarise(
          n_samples       = n(),
          n_adequate      = sum(adequate),
          total_cells     = sum(n_cells),
          min_per_sample  = min(n_cells),
          max_per_sample  = max(n_cells),
          .groups = "drop"
        ) %>%
        mutate(
          usable_for_de = n_adequate >= MIN_SAMPLES_PER_CONDITION,
          flag = case_when(
            total_cells == 0                         ~ "NO_CELLS",
            n_adequate < MIN_SAMPLES_PER_CONDITION   ~ "TOO_FEW_SAMPLES",
            min_per_sample < MIN_CELLS_PER_SAMPLE    ~ "SOME_SAMPLES_LOW",
            TRUE                                     ~ "OK"
          )
        ) %>%
        arrange(group, allen_subclass_filtered, condition)
      
      write_csv(sample_adequacy,
                file.path(OUT_DIR,
                          paste0(cohort_name,
                                 "_subtype_sample_adequacy.csv")))
      
      # Concise DE-readiness table
      de_ready <- sample_adequacy %>%
        group_by(group, allen_subclass_filtered) %>%
        summarise(
          total_cells       = sum(total_cells),
          n_conditions_ok   = sum(usable_for_de),
          n_conditions_total = n(),
          worst_flag        = case_when(
            any(flag == "NO_CELLS")        ~ "NO_CELLS",
            any(flag == "TOO_FEW_SAMPLES") ~ "TOO_FEW_SAMPLES",
            any(flag == "SOME_SAMPLES_LOW") ~ "SOME_SAMPLES_LOW",
            TRUE                            ~ "OK"
          ),
          .groups = "drop"
        ) %>%
        mutate(
          de_ready = n_conditions_ok == n_conditions_total
        ) %>%
        arrange(group, desc(total_cells))
      
      write_csv(de_ready,
                file.path(OUT_DIR,
                          paste0(cohort_name, "_de_readiness.csv")))
      
      # Print summary
      cat("\n\n  === DE Readiness (", cohort_name, ") ===\n")
      cat("  Score threshold: ", SCORE_THRESHOLD, "\n")
      cat("  Min cells/sample: ", MIN_CELLS_PER_SAMPLE, "\n")
      cat("  Min adequate samples/condition: ", MIN_SAMPLES_PER_CONDITION, "\n\n")
      
      n_ready   <- sum(de_ready$de_ready)
      n_total_s <- nrow(de_ready)
      cat("  ", n_ready, " / ", n_total_s,
          " subtypes are DE-ready across all conditions\n\n")
      
      # Show problem subtypes
      problems <- de_ready %>% filter(!de_ready)
      if (nrow(problems) > 0) {
        cat("  Subtypes NOT ready for DE:\n")
        print(as.data.frame(problems), row.names = FALSE)
      }
      
      cat("\n  Subtypes ready for DE:\n")
      ready <- de_ready %>% filter(de_ready)
      if (nrow(ready) > 0) {
        print(as.data.frame(ready), row.names = FALSE)
      }
    }
  }
  
  # -------------------------------------------------------
  # 6. Summary of full annotated object
  # -------------------------------------------------------
  
  cat("\n\n  === Final Annotated Object ===\n")
  cat("  Total cells: ", ncol(obj_full), "\n")
  
  labeled <- sum(!is.na(obj_full$allen_subclass_filtered))
  unlabeled <- sum(is.na(obj_full$allen_subclass_filtered))
  cat("  With confident Allen label:    ", labeled, "\n")
  cat("  Without (low score / Unknown): ", unlabeled, "\n")
  
  anno_summary <- obj_full@meta.data %>%
    count(broad_type, allen_subclass_filtered, name = "n_cells") %>%
    arrange(broad_type, desc(n_cells))
  write_csv(anno_summary,
            file.path(OUT_DIR,
                      paste0(cohort_name, "_final_annotation_summary.csv")))
  
  # --- Save final object ---
  out_path <- paste0("./processed/seurat_", cohort_name,
                     "_final_annotated.rds")
  saveRDS(obj_full, out_path)
  message("  Saved: ", out_path)
  
  obj_full
}

# =============================================================
# 4. Run all cohorts
# =============================================================

final_objects <- list()

for (cohort_name in cohorts_to_run) {
  final_objects[[cohort_name]] <- process_cohort(cohort_name)
}

# =============================================================
# 5. Overall summary
# =============================================================

cat("\n\n==================== OVERALL SUMMARY ====================\n")

for (cohort_name in names(final_objects)) {
  obj <- final_objects[[cohort_name]]
  cat("\n", cohort_name, ":\n")
  cat("  Total cells: ", ncol(obj), "\n")
  
  # Broad type × filtered subclass counts
  bt_sub <- obj@meta.data %>%
    filter(!is.na(allen_subclass_filtered)) %>%
    count(broad_type, name = "n_labeled") %>%
    arrange(broad_type)
  print(as.data.frame(bt_sub), row.names = FALSE)
}

cat("\n\nFinal objects saved:\n")
for (cohort_name in names(final_objects)) {
  cat("  ./processed/seurat_", cohort_name, "_final_annotated.rds\n", sep = "")
}
