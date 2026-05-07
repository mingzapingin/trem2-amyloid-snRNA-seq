# ============================================================
# 10_cross_comparison.R
#
# Cross-comparison analysis across Q1, Q2, Q3 for 7m subtypes.
#
# For each subtype, classifies every DE gene into biological
# categories based on the pattern across comparisons:
#
#   Cat 1: Pure amyloid        — DE in Q1, not Q3
#   Cat 2: Trem2-dependent     — DE in Q1, lost/reversed in Q2
#   Cat 3: Trem2-independent   — DE in Q1, not Q2
#   Cat 4: Trem2-exacerbated   — DE in Q1, stronger in Q2
#   Cat 5: Trem2-autonomous    — DE in Q3, not Q1
#   Cat 6: Redirected          — DE in Q2, not Q1
#
# Outputs:
#   ./result/de_results/cross_comparison/
#     *_gene_categories.csv         per subtype
#     *_category_summary.csv        counts per category per subtype
#     *_lfc_scatter.png             Q1 vs Q2 LFC per subtype
#     *_volcano_category.png        volcano colored by category
#     *_upset.png                   UpSet plot per subtype
#     *_venn.png                    Venn diagram per subtype
#     summary_matrix.csv            subtypes × comparisons
#     summary_matrix_heatmap.png
#     category_counts_all.csv       all subtypes × categories
#     category_counts_heatmap.png
#
# Requires: UpSetR, VennDiagram (install if needed)
# ============================================================

source("./scripts/09_00_de_functions.R")

suppressPackageStartupMessages({
  if (!requireNamespace("UpSetR", quietly = TRUE))
    install.packages("UpSetR", quiet = TRUE)
  if (!requireNamespace("VennDiagram", quietly = TRUE))
    install.packages("VennDiagram", quiet = TRUE)
  library(UpSetR)
  library(VennDiagram)
})

set.seed(1234)

CROSS_DIR <- make_de_dir("cross_comparison")

# =============================================================
# 1. Settings
# =============================================================

# Directories where Q1/Q2/Q3 results live
Q1_DIR <- file.path(DE_BASE_DIR, "Q1_WT_vs_WT5XFAD")
Q2_DIR <- file.path(DE_BASE_DIR, "Q2_WT5XFAD_vs_Trem2KO5XFAD")
Q3_DIR <- file.path(DE_BASE_DIR, "Q3_WT_vs_Trem2KO")

# LFC ratio threshold for "exacerbated" vs "dependent"
# If Q2 LFC / Q1 LFC > this AND same direction, call it exacerbated
EXACERBATION_RATIO <- 1.5

# Subtypes to analyze (intersection of what ran successfully in Q1/Q2/Q3)
SUBTYPES_TO_COMPARE <- DE_READY_7m_ALL

# =============================================================
# 2. Load DE results for one subtype across Q1/Q2/Q3
# =============================================================

load_de_results <- function(subtype, q_dir, q_name) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
  path <- file.path(q_dir, paste0(q_name, "_", safe_name, "_full_results.csv"))
  
  if (!file.exists(path)) {
    message("  No results for ", subtype, " in ", q_name)
    return(NULL)
  }
  
  read_csv(path, show_col_types = FALSE)
}

# =============================================================
# 3. Merge Q1/Q2/Q3 and classify genes
# =============================================================

classify_genes <- function(subtype) {
  message("\n--- Classifying: ", subtype, " ---")
  
  q1 <- load_de_results(subtype, Q1_DIR, "Q1_WT_vs_WT5XFAD")
  q2 <- load_de_results(subtype, Q2_DIR, "Q2_WT5XFAD_vs_Trem2KO5XFAD")
  q3 <- load_de_results(subtype, Q3_DIR, "Q3_WT_vs_Trem2KO")
  
  if (is.null(q1)) return(NULL)
  
  # Build merged table — all genes from any comparison
  all_genes <- unique(c(q1$gene,
                        if (!is.null(q2)) q2$gene else character(0),
                        if (!is.null(q3)) q3$gene else character(0)))
  
  merged <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  
  # Q1 columns
  q1_slim <- q1 %>%
    select(gene,
           q1_lfc = log2FoldChange, q1_padj = padj,
           q1_sig = significant, q1_dir = direction)
  merged <- merged %>% left_join(q1_slim, by = "gene")
  
  # Q2 columns
  if (!is.null(q2)) {
    q2_slim <- q2 %>%
      select(gene,
             q2_lfc = log2FoldChange, q2_padj = padj,
             q2_sig = significant, q2_dir = direction)
    merged <- merged %>% left_join(q2_slim, by = "gene")
  } else {
    merged$q2_lfc <- NA_real_
    merged$q2_padj <- NA_real_
    merged$q2_sig <- FALSE
    merged$q2_dir <- "NS"
  }
  
  # Q3 columns
  if (!is.null(q3)) {
    q3_slim <- q3 %>%
      select(gene,
             q3_lfc = log2FoldChange, q3_padj = padj,
             q3_sig = significant, q3_dir = direction)
    merged <- merged %>% left_join(q3_slim, by = "gene")
  } else {
    merged$q3_lfc <- NA_real_
    merged$q3_padj <- NA_real_
    merged$q3_sig <- FALSE
    merged$q3_dir <- "NS"
  }
  
  # Replace NAs in sig columns with FALSE
  merged <- merged %>%
    mutate(
      q1_sig = replace_na(q1_sig, FALSE),
      q2_sig = replace_na(q2_sig, FALSE),
      q3_sig = replace_na(q3_sig, FALSE),
      q1_dir = replace_na(q1_dir, "NS"),
      q2_dir = replace_na(q2_dir, "NS"),
      q3_dir = replace_na(q3_dir, "NS")
    )
  
  # --- Classify ---
  merged <- merged %>%
    mutate(
      # Direction comparison between Q1 and Q2
      same_direction_q1q2 = (q1_dir == q2_dir) & q1_dir != "NS",
      reversed_q1q2 = (q1_dir == "UP" & q2_dir == "DOWN") |
        (q1_dir == "DOWN" & q2_dir == "UP"),
      
      # LFC ratio for exacerbation detection
      lfc_ratio_q2q1 = ifelse(
        !is.na(q1_lfc) & !is.na(q2_lfc) & abs(q1_lfc) > 0,
        abs(q2_lfc) / abs(q1_lfc),
        NA_real_
      ),
      
      category = case_when(
        # Cat 6: DE in Q2 but NOT Q1 — redirected
        !q1_sig & q2_sig                                      ~ "Cat6_Redirected",
        
        # Cat 5: DE in Q3 but NOT Q1 — Trem2-autonomous
        !q1_sig & q3_sig                                      ~ "Cat5_Trem2_autonomous",
        
        # Cat 4: DE in Q1 AND Q2, same direction, Q2 stronger
        q1_sig & q2_sig & same_direction_q1q2 &
          !is.na(lfc_ratio_q2q1) &
          lfc_ratio_q2q1 > EXACERBATION_RATIO                 ~ "Cat4_Trem2_exacerbated",
        
        # Cat 2: DE in Q1, AND (DE in Q2 reversed OR DE in Q2 lost effect)
        # "Trem2-dependent" = amyloid response requires Trem2
        # Lost: DE in Q1, Q2 sig with reversed direction
        q1_sig & q2_sig & reversed_q1q2                       ~ "Cat2_Trem2_dependent",
        
        # Cat 2 alternate: DE in Q1, Q2 sig same direction but weaker
        q1_sig & q2_sig & same_direction_q1q2 &
          !is.na(lfc_ratio_q2q1) &
          lfc_ratio_q2q1 <= EXACERBATION_RATIO                ~ "Cat2_Trem2_dependent",
        
        # Cat 1: DE in Q1, not Q2, not Q3 — pure amyloid, Trem2 irrelevant
        q1_sig & !q2_sig & !q3_sig                            ~ "Cat1_Pure_amyloid",
        
        # Cat 3: DE in Q1 and Q3, but not Q2 — both affect gene independently
        q1_sig & !q2_sig & q3_sig                             ~ "Cat3_Trem2_independent",
        
        # Not DE anywhere interesting
        TRUE                                                   ~ "Not_DE"
      )
    )
  
  # Refine Cat1 vs Cat3: Cat1 requires not Q3 sig,
  # but if already classified as Cat2/3/4 from Q2 logic, keep that
  # (the case_when priority handles this — earlier matches win)
  
  merged$subtype <- subtype
  merged
}

# =============================================================
# 4. Run classification for all subtypes
# =============================================================

message("\n========== Gene Classification ==========\n")

all_classified <- list()

for (subtype in SUBTYPES_TO_COMPARE) {
  result <- tryCatch(
    classify_genes(subtype),
    error = function(e) {
      message("  ERROR: ", e$message)
      NULL
    }
  )
  
  if (!is.null(result)) {
    all_classified[[subtype]] <- result
    
    # Save per-subtype
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    write_csv(result,
              file.path(CROSS_DIR,
                        paste0(safe_name, "_gene_categories.csv")))
  }
}

# =============================================================
# 5. Category summary table
# =============================================================

cat_counts_list <- lapply(names(all_classified), function(subtype) {
  df <- all_classified[[subtype]]
  df %>%
    filter(category != "Not_DE") %>%
    count(subtype, category, name = "n_genes") %>%
    mutate(subtype = subtype)
})

cat_counts <- bind_rows(cat_counts_list)

if (nrow(cat_counts) > 0) {
  # Wide format
  cat_wide <- cat_counts %>%
    pivot_wider(names_from = category, values_from = n_genes,
                values_fill = 0) %>%
    mutate(total_de = rowSums(across(where(is.numeric)))) %>%
    arrange(desc(total_de))
  
  write_csv(cat_wide,
            file.path(CROSS_DIR, "category_counts_all.csv"))
  
  cat("\n\n=== Category Counts by Subtype ===\n")
  print(as.data.frame(cat_wide), row.names = FALSE)
  
  # --- Heatmap of categories × subtypes ---
  cat_long <- cat_counts %>%
    mutate(subtype = factor(subtype,
                            levels = rev(cat_wide$subtype)))
  
  p_cat <- ggplot(cat_long, aes(x = category, y = subtype,
                                fill = n_genes)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_genes), size = 3) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) +
    ggtitle("Gene categories across subtypes") +
    xlab(NULL) + ylab(NULL)
  
  ggsave(file.path(CROSS_DIR, "category_counts_heatmap.png"),
         plot = p_cat, width = 12, height = 10, dpi = 200)
}

# =============================================================
# 6. Summary matrix: subtypes × comparisons (n DE genes)
# =============================================================

summary_rows <- list()

for (subtype in names(all_classified)) {
  df <- all_classified[[subtype]]
  summary_rows[[subtype]] <- data.frame(
    subtype    = subtype,
    Q1_up      = sum(df$q1_dir == "UP" & df$q1_sig, na.rm = TRUE),
    Q1_down    = sum(df$q1_dir == "DOWN" & df$q1_sig, na.rm = TRUE),
    Q2_up      = sum(df$q2_dir == "UP" & df$q2_sig, na.rm = TRUE),
    Q2_down    = sum(df$q2_dir == "DOWN" & df$q2_sig, na.rm = TRUE),
    Q3_up      = sum(df$q3_dir == "UP" & df$q3_sig, na.rm = TRUE),
    Q3_down    = sum(df$q3_dir == "DOWN" & df$q3_sig, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summary_matrix <- bind_rows(summary_rows) %>%
  mutate(Q1_total = Q1_up + Q1_down,
         Q2_total = Q2_up + Q2_down,
         Q3_total = Q3_up + Q3_down) %>%
  arrange(desc(Q1_total))

write_csv(summary_matrix,
          file.path(CROSS_DIR, "summary_matrix.csv"))

cat("\n\n=== Summary Matrix ===\n")
print(as.data.frame(summary_matrix), row.names = FALSE)

# --- Summary matrix heatmap ---
sm_long <- summary_matrix %>%
  select(subtype, Q1_total, Q2_total, Q3_total) %>%
  pivot_longer(-subtype, names_to = "comparison", values_to = "n_de") %>%
  mutate(subtype = factor(subtype, levels = rev(summary_matrix$subtype)),
         comparison = gsub("_total", "", comparison))

p_sm <- ggplot(sm_long, aes(x = comparison, y = subtype, fill = n_de)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_de), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) +
  ggtitle("Total DE genes per subtype × comparison") +
  xlab(NULL) + ylab(NULL)

ggsave(file.path(CROSS_DIR, "summary_matrix_heatmap.png"),
       plot = p_sm, width = 7, height = 10, dpi = 200)

# =============================================================
# 7. Per-subtype plots: LFC scatter, UpSet, Venn, volcano
# =============================================================

for (subtype in names(all_classified)) {
  df <- all_classified[[subtype]]
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
  
  # ---------------------------------------------------------
  # A. LFC scatter: Q1 vs Q2
  # ---------------------------------------------------------
  plot_df <- df %>%
    filter(!is.na(q1_lfc) & !is.na(q2_lfc)) %>%
    mutate(
      highlight = case_when(
        category %in% c("Cat2_Trem2_dependent")    ~ "Trem2-dependent",
        category %in% c("Cat4_Trem2_exacerbated")  ~ "Trem2-exacerbated",
        category %in% c("Cat3_Trem2_independent")  ~ "Trem2-independent",
        category %in% c("Cat6_Redirected")          ~ "Redirected",
        TRUE                                         ~ "Other"
      )
    )
  
  if (nrow(plot_df) > 100) {
    p_lfc <- ggplot(plot_df, aes(x = q1_lfc, y = q2_lfc, color = highlight)) +
      geom_point(alpha = 0.4, size = 0.8) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, color = "grey80") +
      geom_vline(xintercept = 0, color = "grey80") +
      scale_color_manual(values = c(
        "Trem2-dependent"    = CATEGORY_COLORS["Cat2"],
        "Trem2-exacerbated"  = CATEGORY_COLORS["Cat4"],
        "Trem2-independent"  = CATEGORY_COLORS["Cat3"],
        "Redirected"         = CATEGORY_COLORS["Cat6"],
        "Other"              = CATEGORY_COLORS["NS"]
      )) +
      theme_bw() +
      ggtitle(paste0(subtype, " — Q1 LFC vs Q2 LFC")) +
      xlab("log2FC: WT_5XFAD vs WT (Q1)") +
      ylab("log2FC: Trem2KO_5XFAD vs WT_5XFAD (Q2)") +
      theme(panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"))
    
    ggsave(file.path(CROSS_DIR,
                     paste0(safe_name, "_lfc_scatter_Q1vQ2.png")),
           plot = p_lfc, width = 8, height = 7, dpi = 200)
  }
  
  # ---------------------------------------------------------
  # B. UpSet plot: Q1 ∩ Q2 ∩ Q3 sig gene overlaps
  # ---------------------------------------------------------
  upset_input <- data.frame(
    gene = df$gene,
    Q1   = as.integer(df$q1_sig),
    Q2   = as.integer(df$q2_sig),
    Q3   = as.integer(df$q3_sig)
  ) %>%
    filter(Q1 == 1 | Q2 == 1 | Q3 == 1)
  
  if (nrow(upset_input) >= 5) {
    tryCatch({
      png(file.path(CROSS_DIR,
                    paste0(safe_name, "_upset.png")),
          width = 8, height = 5, units = "in", res = 200)
      print(upset(upset_input,
                  sets = c("Q1", "Q2", "Q3"),
                  order.by = "freq",
                  sets.bar.color = "steelblue",
                  main.bar.color = "grey30",
                  text.scale = 1.3,
                  mainbar.y.label = paste0(subtype, " — DE gene intersections")))
      dev.off()
    }, error = function(e) {
      message("  UpSet failed for ", subtype, ": ", e$message)
      try(dev.off(), silent = TRUE)
    })
  }
  
  # ---------------------------------------------------------
  # C. Venn diagram: Q1 ∩ Q2 ∩ Q3
  # ---------------------------------------------------------
  q1_genes <- df$gene[df$q1_sig]
  q2_genes <- df$gene[df$q2_sig]
  q3_genes <- df$gene[df$q3_sig]
  
  if (length(q1_genes) + length(q2_genes) + length(q3_genes) >= 3) {
    tryCatch({
      venn_list <- list(
        Q1_Amyloid  = q1_genes,
        Q2_Trem2_5X = q2_genes,
        Q3_Trem2    = q3_genes
      )
      # Remove empty sets
      venn_list <- venn_list[lengths(venn_list) > 0]
      
      if (length(venn_list) >= 2) {
        futile.logger::flog.threshold(futile.logger::ERROR)
        venn.diagram(
          x = venn_list,
          filename = file.path(CROSS_DIR,
                               paste0(safe_name, "_venn.png")),
          output = TRUE,
          imagetype = "png",
          height = 2400, width = 2800, resolution = 300,
          fill = c("firebrick", "steelblue", "forestgreen")[seq_along(venn_list)],
          alpha = 0.3,
          main = paste0(subtype, " — DE gene overlap"),
          main.cex = 1.2,
          cat.cex = 0.9,
          cex = 1.0
        )
      }
    }, error = function(e) {
      message("  Venn failed for ", subtype, ": ", e$message)
    })
  }
  
  # ---------------------------------------------------------
  # D. Volcano colored by category (using Q1 as base)
  # ---------------------------------------------------------
  vol_df <- df %>%
    filter(!is.na(q1_lfc) & !is.na(q1_padj)) %>%
    mutate(
      cat_color = case_when(
        category == "Cat1_Pure_amyloid"               ~ "Cat1",
        category == "Cat2_Trem2_dependent"            ~ "Cat2",
        category == "Cat3_Trem2_independent"          ~ "Cat3",
        category == "Cat4_Trem2_exacerbated"          ~ "Cat4",
        category == "Cat5_Trem2_autonomous"           ~ "Cat5",
        category == "Cat6_Redirected"                 ~ "Cat6",
        TRUE                                           ~ "NS"
      )
    )
  
  cat_colors <- CATEGORY_COLORS
  
  p_vol <- ggplot(vol_df, aes(x = q1_lfc, y = -log10(q1_padj),
                              color = cat_color)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = cat_colors, name = "Category") +
    theme_bw() +
    ggtitle(paste0(subtype, " — Q1 volcano colored by cross-comparison category")) +
    xlab("log2FC (WT_5XFAD vs WT)") +
    ylab("-log10(padj)") +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
  
  ggsave(file.path(CROSS_DIR,
                   paste0(safe_name, "_volcano_category.png")),
         plot = p_vol, width = 9, height = 7, dpi = 200)
}

# =============================================================
# 8. Top candidates table
# =============================================================
# Pull the most interesting genes: Cat2 (Trem2-dependent) and
# Cat4 (exacerbated) with strongest effects

top_candidates <- bind_rows(all_classified) %>%
  filter(category %in% c("Cat2_Trem2_dependent",
                         "Cat4_Trem2_exacerbated",
                         "Cat6_Redirected")) %>%
  mutate(q1_abs_lfc = abs(q1_lfc)) %>%
  arrange(category, desc(q1_abs_lfc)) %>%
  select(subtype, gene, category,
         q1_lfc, q1_padj, q2_lfc, q2_padj, q3_lfc, q3_padj,
         lfc_ratio_q2q1)

write_csv(top_candidates,
          file.path(CROSS_DIR, "top_candidates_cat2_cat4_cat6.csv"))

cat("\n\n=== Top Trem2-Dependent Genes (Cat2, top 20) ===\n")
top_cat2 <- top_candidates %>%
  filter(category == "Cat2_Trem2_dependent") %>%
  head(20)
print(as.data.frame(top_cat2), row.names = FALSE)

cat("\n\n=== Top Trem2-Exacerbated Genes (Cat4, top 20) ===\n")
top_cat4 <- top_candidates %>%
  filter(category == "Cat4_Trem2_exacerbated") %>%
  head(20)
print(as.data.frame(top_cat4), row.names = FALSE)

cat("\n\n=== Top Redirected Genes (Cat6, top 20) ===\n")
top_cat6 <- top_candidates %>%
  filter(category == "Cat6_Redirected") %>%
  head(20)
print(as.data.frame(top_cat6), row.names = FALSE)

# =============================================================
# 9. REGIONAL ANALYSIS (15m Q4a–Q4d)
# =============================================================
# All four Q4 comparisons are Cortex vs Hippocampus under
# different genotype × amyloid conditions:
#   Q4a: WT           — baseline regional differences
#   Q4b: WT_5XFAD     — regional response to amyloid
#   Q4c: Trem2_KO     — regional diff without amyloid
#   Q4d: Trem2_KO_5XFAD — regional diff with amyloid + KO
#
# Regional categories:
#   Reg1: Constitutive      — DE in Q4a, stable across conditions
#   Reg2: Amyloid-emergent  — DE in Q4b, NOT Q4a (amyloid creates regional diff)
#   Reg3: Amyloid-lost      — DE in Q4a, NOT Q4b (amyloid erases regional diff)
#   Reg4: Trem2-dependent   — DE in Q4b, NOT Q4d (regional amyloid response needs Trem2)
#   Reg5: Trem2-emergent    — DE in Q4c or Q4d, NOT Q4a (Trem2 KO alters regional pattern)
#   Reg6: Compound          — DE in Q4d only (needs both amyloid + KO)

message("\n\n========== Regional Cross-Comparison (15m) ==========\n")

REGIONAL_DIR <- file.path(CROSS_DIR, "regional")
dir.create(REGIONAL_DIR, recursive = TRUE, showWarnings = FALSE)

Q4A_DIR <- file.path(DE_BASE_DIR, "Q4a_15m_WT_Cor_vs_Hip")
Q4B_DIR <- file.path(DE_BASE_DIR, "Q4b_15m_WT5XFAD_Cor_vs_Hip")
Q4C_DIR <- file.path(DE_BASE_DIR, "Q4c_15m_Trem2KO_Cor_vs_Hip")
Q4D_DIR <- file.path(DE_BASE_DIR, "Q4d_15m_Trem2KO5XFAD_Cor_vs_Hip")

SUBTYPES_REGIONAL <- DE_READY_15m_ALL

# --- Load + classify for one subtype ---
classify_regional <- function(subtype) {
  message("\n--- Regional: ", subtype, " ---")
  
  q4a <- load_de_results(subtype, Q4A_DIR, "Q4a_15m_WT_Cor_vs_Hip")
  q4b <- load_de_results(subtype, Q4B_DIR, "Q4b_15m_WT5XFAD_Cor_vs_Hip")
  q4c <- load_de_results(subtype, Q4C_DIR, "Q4c_15m_Trem2KO_Cor_vs_Hip")
  q4d <- load_de_results(subtype, Q4D_DIR, "Q4d_15m_Trem2KO5XFAD_Cor_vs_Hip")
  
  # Need at least Q4a
  if (is.null(q4a) && is.null(q4b) && is.null(q4c) && is.null(q4d)) {
    message("  No regional results — skipping")
    return(NULL)
  }
  
  # Collect all genes
  all_genes <- unique(c(
    if (!is.null(q4a)) q4a$gene else character(0),
    if (!is.null(q4b)) q4b$gene else character(0),
    if (!is.null(q4c)) q4c$gene else character(0),
    if (!is.null(q4d)) q4d$gene else character(0)
  ))
  
  merged <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  
  # Helper to join one comparison
  join_q <- function(merged, q_df, prefix) {
    if (!is.null(q_df)) {
      slim <- q_df %>%
        select(gene,
               !!paste0(prefix, "_lfc") := log2FoldChange,
               !!paste0(prefix, "_padj") := padj,
               !!paste0(prefix, "_sig") := significant,
               !!paste0(prefix, "_dir") := direction)
      merged <- merged %>% left_join(slim, by = "gene")
    } else {
      merged[[paste0(prefix, "_lfc")]]  <- NA_real_
      merged[[paste0(prefix, "_padj")]] <- NA_real_
      merged[[paste0(prefix, "_sig")]]  <- FALSE
      merged[[paste0(prefix, "_dir")]]  <- "NS"
    }
    merged
  }
  
  merged <- join_q(merged, q4a, "q4a")
  merged <- join_q(merged, q4b, "q4b")
  merged <- join_q(merged, q4c, "q4c")
  merged <- join_q(merged, q4d, "q4d")
  
  # Fill NAs
  merged <- merged %>%
    mutate(across(ends_with("_sig"), ~ replace_na(.x, FALSE)),
           across(ends_with("_dir"), ~ replace_na(.x, "NS")))
  
  # --- Classify ---
  merged <- merged %>%
    mutate(
      reg_category = case_when(
        # Reg6: ONLY in Q4d (compound: needs amyloid + Trem2 KO)
        !q4a_sig & !q4b_sig & !q4c_sig & q4d_sig  ~ "Reg6_Compound",
        
        # Reg5: DE in Q4c or Q4d, NOT Q4a (Trem2 KO alters regional pattern)
        !q4a_sig & (q4c_sig | q4d_sig)             ~ "Reg5_Trem2_emergent",
        
        # Reg2: DE in Q4b, NOT Q4a (amyloid creates new regional diff)
        !q4a_sig & q4b_sig                         ~ "Reg2_Amyloid_emergent",
        
        # Reg4: DE in Q4b (amyloid regional), NOT Q4d (lost when KO)
        q4b_sig & !q4d_sig                         ~ "Reg4_Trem2_dependent_regional",
        
        # Reg3: DE in Q4a (baseline), NOT Q4b (amyloid erases it)
        q4a_sig & !q4b_sig                         ~ "Reg3_Amyloid_lost",
        
        # Reg1: DE in Q4a baseline — present across conditions
        q4a_sig                                     ~ "Reg1_Constitutive",
        
        TRUE                                        ~ "Not_regional_DE"
      ),
      subtype = subtype
    )
  
  merged
}

# --- Run for all subtypes ---
all_regional <- list()

for (subtype in SUBTYPES_REGIONAL) {
  result <- tryCatch(
    classify_regional(subtype),
    error = function(e) {
      message("  ERROR: ", e$message)
      NULL
    }
  )
  if (!is.null(result)) {
    all_regional[[subtype]] <- result
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    write_csv(result,
              file.path(REGIONAL_DIR,
                        paste0(safe_name, "_regional_categories.csv")))
  }
}

# --- Regional category summary ---
if (length(all_regional) > 0) {
  
  reg_counts_list <- lapply(names(all_regional), function(subtype) {
    all_regional[[subtype]] %>%
      filter(reg_category != "Not_regional_DE") %>%
      count(subtype, reg_category, name = "n_genes")
  })
  
  reg_counts <- bind_rows(reg_counts_list)
  
  if (nrow(reg_counts) > 0) {
    reg_wide <- reg_counts %>%
      pivot_wider(names_from = reg_category, values_from = n_genes,
                  values_fill = 0) %>%
      mutate(total_reg_de = rowSums(across(where(is.numeric)))) %>%
      arrange(desc(total_reg_de))
    
    write_csv(reg_wide,
              file.path(REGIONAL_DIR, "regional_category_counts.csv"))
    
    cat("\n\n=== Regional Category Counts ===\n")
    print(as.data.frame(reg_wide), row.names = FALSE)
    
    # --- Heatmap ---
    reg_long <- reg_counts %>%
      mutate(subtype = factor(subtype, levels = rev(reg_wide$subtype)))
    
    p_reg_cat <- ggplot(reg_long,
                        aes(x = reg_category, y = subtype, fill = n_genes)) +
      geom_tile(color = "white") +
      geom_text(aes(label = n_genes), size = 3) +
      scale_fill_gradient(low = "white", high = "#E64B35") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white")) +
      ggtitle("15m regional gene categories across subtypes") +
      xlab(NULL) + ylab(NULL)
    
    ggsave(file.path(REGIONAL_DIR, "regional_category_heatmap.png"),
           plot = p_reg_cat, width = 12, height = 8, dpi = 200)
  }
  
  # --- Regional summary matrix: subtypes × Q4a/Q4b/Q4c/Q4d DE counts ---
  reg_summary_rows <- list()
  for (subtype in names(all_regional)) {
    df <- all_regional[[subtype]]
    reg_summary_rows[[subtype]] <- data.frame(
      subtype   = subtype,
      Q4a_up    = sum(df$q4a_dir == "UP" & df$q4a_sig, na.rm = TRUE),
      Q4a_down  = sum(df$q4a_dir == "DOWN" & df$q4a_sig, na.rm = TRUE),
      Q4b_up    = sum(df$q4b_dir == "UP" & df$q4b_sig, na.rm = TRUE),
      Q4b_down  = sum(df$q4b_dir == "DOWN" & df$q4b_sig, na.rm = TRUE),
      Q4c_up    = sum(df$q4c_dir == "UP" & df$q4c_sig, na.rm = TRUE),
      Q4c_down  = sum(df$q4c_dir == "DOWN" & df$q4c_sig, na.rm = TRUE),
      Q4d_up    = sum(df$q4d_dir == "UP" & df$q4d_sig, na.rm = TRUE),
      Q4d_down  = sum(df$q4d_dir == "DOWN" & df$q4d_sig, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  reg_summary <- bind_rows(reg_summary_rows) %>%
    mutate(Q4a_total = Q4a_up + Q4a_down,
           Q4b_total = Q4b_up + Q4b_down,
           Q4c_total = Q4c_up + Q4c_down,
           Q4d_total = Q4d_up + Q4d_down) %>%
    arrange(desc(Q4a_total + Q4b_total + Q4c_total + Q4d_total))
  
  write_csv(reg_summary,
            file.path(REGIONAL_DIR, "regional_summary_matrix.csv"))
  
  cat("\n\n=== Regional Summary Matrix ===\n")
  print(as.data.frame(
    reg_summary %>% select(subtype, Q4a_total, Q4b_total, Q4c_total, Q4d_total)
  ), row.names = FALSE)
  
  # --- Regional summary heatmap ---
  rsm_long <- reg_summary %>%
    select(subtype, Q4a_total, Q4b_total, Q4c_total, Q4d_total) %>%
    pivot_longer(-subtype, names_to = "comparison", values_to = "n_de") %>%
    mutate(subtype = factor(subtype, levels = rev(reg_summary$subtype)),
           comparison = gsub("_total", "", comparison))
  
  p_rsm <- ggplot(rsm_long, aes(x = comparison, y = subtype, fill = n_de)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_de), size = 3) +
    scale_fill_gradient(low = "white", high = "#4DBBD5") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) +
    ggtitle("15m: regional DE genes per subtype × condition") +
    labs(subtitle = "Q4a=WT | Q4b=5XFAD | Q4c=Trem2KO | Q4d=Trem2KO_5XFAD") +
    xlab(NULL) + ylab(NULL)
  
  ggsave(file.path(REGIONAL_DIR, "regional_summary_heatmap.png"),
         plot = p_rsm, width = 8, height = 8, dpi = 200)
  
  # --- Per-subtype plots ---
  for (subtype in names(all_regional)) {
    df <- all_regional[[subtype]]
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", subtype)
    
    # LFC scatter: Q4a (baseline) vs Q4b (amyloid)
    scat_df <- df %>% filter(!is.na(q4a_lfc) & !is.na(q4b_lfc))
    if (nrow(scat_df) > 50) {
      scat_df <- scat_df %>%
        mutate(highlight = case_when(
          reg_category == "Reg2_Amyloid_emergent"           ~ "Amyloid-emergent",
          reg_category == "Reg3_Amyloid_lost"               ~ "Amyloid-lost",
          reg_category == "Reg4_Trem2_dependent_regional"   ~ "Trem2-dependent",
          reg_category == "Reg1_Constitutive"               ~ "Constitutive",
          TRUE                                               ~ "Other"
        ))
      
      p_scat <- ggplot(scat_df, aes(x = q4a_lfc, y = q4b_lfc,
                                    color = highlight)) +
        geom_point(alpha = 0.4, size = 0.8) +
        geom_abline(slope = 1, intercept = 0,
                    linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = 0, color = "grey80") +
        geom_vline(xintercept = 0, color = "grey80") +
        scale_color_manual(values = c(
          "Constitutive"     = REGIONAL_COLORS["Reg1"],
          "Amyloid-emergent" = REGIONAL_COLORS["Reg2"],
          "Amyloid-lost"     = REGIONAL_COLORS["Reg3"],
          "Trem2-dependent"  = REGIONAL_COLORS["Reg4"],
          "Other"            = REGIONAL_COLORS["NS"]
        )) +
        theme_bw() +
        ggtitle(paste0(subtype, " — Q4a (WT) vs Q4b (5XFAD) regional LFC")) +
        xlab("Cor vs Hip LFC in WT (Q4a)") +
        ylab("Cor vs Hip LFC in 5XFAD (Q4b)") +
        theme(panel.background = element_rect(fill = "white"),
              plot.background = element_rect(fill = "white"))
      
      ggsave(file.path(REGIONAL_DIR,
                       paste0(safe_name, "_lfc_scatter_Q4avQ4b.png")),
             plot = p_scat, width = 8, height = 7, dpi = 200)
    }
    
    # UpSet: Q4a ∩ Q4b ∩ Q4c ∩ Q4d
    upset_df <- data.frame(
      gene = df$gene,
      Q4a_WT        = as.integer(df$q4a_sig),
      Q4b_5XFAD     = as.integer(df$q4b_sig),
      Q4c_Trem2KO   = as.integer(df$q4c_sig),
      Q4d_KO_5XFAD  = as.integer(df$q4d_sig)
    ) %>%
      filter(Q4a_WT == 1 | Q4b_5XFAD == 1 | Q4c_Trem2KO == 1 | Q4d_KO_5XFAD == 1)
    
    if (nrow(upset_df) >= 5) {
      tryCatch({
        png(file.path(REGIONAL_DIR,
                      paste0(safe_name, "_regional_upset.png")),
            width = 9, height = 5, units = "in", res = 200)
        print(upset(upset_df,
                    sets = c("Q4a_WT", "Q4b_5XFAD", "Q4c_Trem2KO", "Q4d_KO_5XFAD"),
                    order.by = "freq",
                    sets.bar.color = "#4DBBD5",
                    main.bar.color = "grey30",
                    text.scale = 1.2,
                    mainbar.y.label = paste0(subtype, " — regional DE intersections")))
        dev.off()
      }, error = function(e) {
        message("  UpSet failed for ", subtype, ": ", e$message)
        try(dev.off(), silent = TRUE)
      })
    }
    
    # Venn: Q4a ∩ Q4b ∩ Q4c ∩ Q4d
    venn_list <- list(
      Q4a_WT       = df$gene[df$q4a_sig],
      Q4b_5XFAD    = df$gene[df$q4b_sig],
      Q4c_Trem2KO  = df$gene[df$q4c_sig],
      Q4d_KO_5XFAD = df$gene[df$q4d_sig]
    )
    venn_list <- venn_list[lengths(venn_list) > 0]
    
    if (length(venn_list) >= 2) {
      tryCatch({
        futile.logger::flog.threshold(futile.logger::ERROR)
        venn.diagram(
          x = venn_list,
          filename = file.path(REGIONAL_DIR,
                               paste0(safe_name, "_regional_venn.png")),
          output = TRUE,
          imagetype = "png",
          height = 2600, width = 3000, resolution = 300,
          fill = c("#636363", "#E64B35", "#B07AA1", "#59A14F")[seq_along(venn_list)],
          alpha = 0.3,
          main = paste0(subtype, " — regional DE overlap"),
          main.cex = 1.1,
          cat.cex = 0.8,
          cex = 0.9
        )
      }, error = function(e) {
        message("  Venn failed for ", subtype, ": ", e$message)
      })
    }
  }
  
  # --- Top regional candidates ---
  reg_top <- bind_rows(all_regional) %>%
    filter(reg_category %in% c("Reg2_Amyloid_emergent",
                               "Reg4_Trem2_dependent_regional",
                               "Reg5_Trem2_emergent",
                               "Reg6_Compound")) %>%
    mutate(max_abs_lfc = pmax(abs(q4a_lfc), abs(q4b_lfc),
                              abs(q4c_lfc), abs(q4d_lfc),
                              na.rm = TRUE)) %>%
    arrange(reg_category, desc(max_abs_lfc)) %>%
    select(subtype, gene, reg_category,
           q4a_lfc, q4a_padj, q4b_lfc, q4b_padj,
           q4c_lfc, q4c_padj, q4d_lfc, q4d_padj)
  
  write_csv(reg_top,
            file.path(REGIONAL_DIR, "regional_top_candidates.csv"))
  
  cat("\n\n=== Top Amyloid-Emergent Regional Genes (Reg2, top 15) ===\n")
  print(as.data.frame(
    reg_top %>% filter(reg_category == "Reg2_Amyloid_emergent") %>% head(15)
  ), row.names = FALSE)
  
  cat("\n\n=== Top Trem2-Dependent Regional Genes (Reg4, top 15) ===\n")
  print(as.data.frame(
    reg_top %>% filter(reg_category == "Reg4_Trem2_dependent_regional") %>% head(15)
  ), row.names = FALSE)
}

# =============================================================
# 10. Final Summary
# =============================================================

cat("\n\n==================== CROSS-COMPARISON COMPLETE ====================\n")
cat("\n--- 7m Trem2 × Amyloid (Q1/Q2/Q3) ---\n")
cat("Subtypes analyzed: ", length(all_classified), "\n")

cat("\n--- 15m Regional (Q4a/Q4b/Q4c/Q4d) ---\n")
cat("Subtypes analyzed: ", length(all_regional), "\n")

cat("\nOutput directory:  ", CROSS_DIR, "\n\n")
cat("Key files (7m cross-comparison):\n")
cat("  category_counts_all.csv          — genes per category per subtype\n")
cat("  category_counts_heatmap.png      — which subtypes are Trem2-sensitive\n")
cat("  summary_matrix.csv               — DE counts Q1/Q2/Q3 per subtype\n")
cat("  top_candidates_cat2_cat4_cat6.csv — follow-up gene list\n")
cat("  <subtype>_lfc_scatter_Q1vQ2.png  — direction-aware LFC comparison\n")
cat("  <subtype>_upset.png / _venn.png  — set intersections\n")
cat("  <subtype>_volcano_category.png   — volcano colored by biology\n")
cat("\nKey files (15m regional):\n")
cat("  regional/regional_category_counts.csv    — regional categories per subtype\n")
cat("  regional/regional_category_heatmap.png   — which subtypes show regional effects\n")
cat("  regional/regional_summary_matrix.csv     — DE counts Q4a/Q4b/Q4c/Q4d\n")
cat("  regional/regional_top_candidates.csv     — interesting region-specific genes\n")
cat("  regional/<subtype>_lfc_scatter_Q4avQ4b.png — baseline vs amyloid regional LFC\n")
cat("  regional/<subtype>_regional_upset.png    — 4-way intersection\n")