# ============================================================
# 10_01_cross_comparison_broad.R
#
# Same cross-comparison framework as script 10 but at the
# broad cell-type level (6 cell classes) using results from
# 09_05 (Q1-Q3 broad) and 09_06 (Q4a-Q4d broad).
#
# 7m: Cat1–Cat6 classification across Q1_broad/Q2_broad/Q3_broad
# 15m: Reg1–Reg6 classification across Q4a_broad–Q4d_broad
#
# Outputs:
#   ./result/de_results/cross_comparison_broad/
#   ./result/de_results/cross_comparison_broad/regional/
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

CROSS_DIR <- make_de_dir("cross_comparison_broad")

# =============================================================
# 1. Settings
# =============================================================

Q1B_DIR <- file.path(DE_BASE_DIR, "Q1_broad_WT_vs_WT5XFAD")
Q2B_DIR <- file.path(DE_BASE_DIR, "Q2_broad_WT5XFAD_vs_Trem2KO5XFAD")
Q3B_DIR <- file.path(DE_BASE_DIR, "Q3_broad_WT_vs_Trem2KO")

Q4AB_DIR <- file.path(DE_BASE_DIR, "Q4a_broad_15m_WT_Cor_vs_Hip")
Q4BB_DIR <- file.path(DE_BASE_DIR, "Q4b_broad_15m_WT5XFAD_Cor_vs_Hip")
Q4CB_DIR <- file.path(DE_BASE_DIR, "Q4c_broad_15m_Trem2KO_Cor_vs_Hip")
Q4DB_DIR <- file.path(DE_BASE_DIR, "Q4d_broad_15m_Trem2KO5XFAD_Cor_vs_Hip")

EXACERBATION_RATIO <- 1.5

BROAD_TYPES_7m_COMPARE  <- BROAD_TYPES_7m
BROAD_TYPES_15m_COMPARE <- BROAD_TYPES_15m

# =============================================================
# 2. Helper: load DE results
# =============================================================

load_de_results <- function(cell_type, q_dir, q_name) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", cell_type)
  path <- file.path(q_dir, paste0(q_name, "_", safe_name, "_full_results.csv"))
  if (!file.exists(path)) {
    message("  No results for ", cell_type, " in ", q_name)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE)
}

# =============================================================
# 3. Helper: merge + join one comparison
# =============================================================

join_comparison <- function(merged, q_df, prefix) {
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

fill_na_sig <- function(merged) {
  merged %>%
    mutate(across(ends_with("_sig"), ~ replace_na(.x, FALSE)),
           across(ends_with("_dir"), ~ replace_na(.x, "NS")))
}

# =============================================================
# 4. 7m BROAD: Cat1–Cat6 classification
# =============================================================

message("\n\n========== 7m Broad Cross-Comparison (Q1/Q2/Q3) ==========\n")

all_classified <- list()

for (bt in BROAD_TYPES_7m_COMPARE) {
  message("\n--- ", bt, " ---")
  
  q1 <- load_de_results(bt, Q1B_DIR, "Q1_broad_WT_vs_WT5XFAD")
  q2 <- load_de_results(bt, Q2B_DIR, "Q2_broad_WT5XFAD_vs_Trem2KO5XFAD")
  q3 <- load_de_results(bt, Q3B_DIR, "Q3_broad_WT_vs_Trem2KO")
  
  if (is.null(q1)) next
  
  all_genes <- unique(c(q1$gene,
                        if (!is.null(q2)) q2$gene else character(0),
                        if (!is.null(q3)) q3$gene else character(0)))
  
  merged <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  merged <- join_comparison(merged, q1, "q1")
  merged <- join_comparison(merged, q2, "q2")
  merged <- join_comparison(merged, q3, "q3")
  merged <- fill_na_sig(merged)
  
  merged <- merged %>%
    mutate(
      same_direction_q1q2 = (q1_dir == q2_dir) & q1_dir != "NS",
      reversed_q1q2 = (q1_dir == "UP" & q2_dir == "DOWN") |
        (q1_dir == "DOWN" & q2_dir == "UP"),
      lfc_ratio_q2q1 = ifelse(
        !is.na(q1_lfc) & !is.na(q2_lfc) & abs(q1_lfc) > 0,
        abs(q2_lfc) / abs(q1_lfc), NA_real_),
      
      category = case_when(
        !q1_sig & q2_sig                                        ~ "Cat6_Redirected",
        !q1_sig & q3_sig                                        ~ "Cat5_Trem2_autonomous",
        q1_sig & q2_sig & same_direction_q1q2 &
          !is.na(lfc_ratio_q2q1) &
          lfc_ratio_q2q1 > EXACERBATION_RATIO                   ~ "Cat4_Trem2_exacerbated",
        q1_sig & q2_sig & reversed_q1q2                         ~ "Cat2_Trem2_dependent",
        q1_sig & q2_sig & same_direction_q1q2 &
          !is.na(lfc_ratio_q2q1) &
          lfc_ratio_q2q1 <= EXACERBATION_RATIO                  ~ "Cat2_Trem2_dependent",
        q1_sig & !q2_sig & !q3_sig                              ~ "Cat1_Pure_amyloid",
        q1_sig & !q2_sig & q3_sig                               ~ "Cat3_Trem2_independent",
        TRUE                                                     ~ "Not_DE"
      ),
      broad_type = bt
    )
  
  all_classified[[bt]] <- merged
  
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
  write_csv(merged, file.path(CROSS_DIR,
                              paste0(safe_name, "_gene_categories.csv")))
}

# --- Category summary ---
cat_counts <- bind_rows(lapply(names(all_classified), function(bt) {
  all_classified[[bt]] %>%
    filter(category != "Not_DE") %>%
    count(broad_type, category, name = "n_genes")
}))

if (nrow(cat_counts) > 0) {
  cat_wide <- cat_counts %>%
    pivot_wider(names_from = category, values_from = n_genes,
                values_fill = 0) %>%
    mutate(total_de = rowSums(across(where(is.numeric)))) %>%
    arrange(desc(total_de))
  
  write_csv(cat_wide, file.path(CROSS_DIR, "category_counts_broad.csv"))
  
  cat("\n\n=== Broad Type Category Counts ===\n")
  print(as.data.frame(cat_wide), row.names = FALSE)
  
  # Heatmap
  cat_long <- cat_counts %>%
    mutate(broad_type = factor(broad_type, levels = rev(cat_wide$broad_type)))
  
  p_cat <- ggplot(cat_long, aes(x = category, y = broad_type, fill = n_genes)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_genes), size = 4) +
    scale_fill_gradient(low = "white", high = "#E64B35") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) +
    ggtitle("7m BROAD: gene categories across cell classes") +
    xlab(NULL) + ylab(NULL)
  
  ggsave(file.path(CROSS_DIR, "category_counts_broad_heatmap.png"),
         plot = p_cat, width = 11, height = 6, dpi = 200)
}

# --- Summary matrix ---
summary_rows <- list()
for (bt in names(all_classified)) {
  df <- all_classified[[bt]]
  summary_rows[[bt]] <- data.frame(
    broad_type = bt,
    Q1_up   = sum(df$q1_dir == "UP" & df$q1_sig, na.rm = TRUE),
    Q1_down = sum(df$q1_dir == "DOWN" & df$q1_sig, na.rm = TRUE),
    Q2_up   = sum(df$q2_dir == "UP" & df$q2_sig, na.rm = TRUE),
    Q2_down = sum(df$q2_dir == "DOWN" & df$q2_sig, na.rm = TRUE),
    Q3_up   = sum(df$q3_dir == "UP" & df$q3_sig, na.rm = TRUE),
    Q3_down = sum(df$q3_dir == "DOWN" & df$q3_sig, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

summary_matrix <- bind_rows(summary_rows) %>%
  mutate(Q1_total = Q1_up + Q1_down,
         Q2_total = Q2_up + Q2_down,
         Q3_total = Q3_up + Q3_down) %>%
  arrange(desc(Q1_total))

write_csv(summary_matrix, file.path(CROSS_DIR, "summary_matrix_broad.csv"))

cat("\n\n=== Broad Summary Matrix ===\n")
print(as.data.frame(summary_matrix), row.names = FALSE)

# Summary heatmap
sm_long <- summary_matrix %>%
  select(broad_type, Q1_total, Q2_total, Q3_total) %>%
  pivot_longer(-broad_type, names_to = "comparison", values_to = "n_de") %>%
  mutate(broad_type = factor(broad_type, levels = rev(summary_matrix$broad_type)),
         comparison = gsub("_total", "", comparison))

p_sm <- ggplot(sm_long, aes(x = comparison, y = broad_type, fill = n_de)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_de), size = 4) +
  scale_fill_gradient(low = "white", high = "#4DBBD5") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) +
  ggtitle("7m BROAD: total DE genes per cell class × comparison") +
  xlab(NULL) + ylab(NULL)

ggsave(file.path(CROSS_DIR, "summary_matrix_broad_heatmap.png"),
       plot = p_sm, width = 7, height = 5, dpi = 200)

# --- Per broad-type plots ---
for (bt in names(all_classified)) {
  df <- all_classified[[bt]]
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
  
  # LFC scatter Q1 vs Q2
  scat_df <- df %>%
    filter(!is.na(q1_lfc) & !is.na(q2_lfc)) %>%
    mutate(highlight = case_when(
      category %in% c("Cat2_Trem2_dependent")   ~ "Trem2-dependent",
      category %in% c("Cat4_Trem2_exacerbated") ~ "Trem2-exacerbated",
      category %in% c("Cat3_Trem2_independent") ~ "Trem2-independent",
      category %in% c("Cat6_Redirected")         ~ "Redirected",
      TRUE                                        ~ "Other"))
  
  if (nrow(scat_df) > 50) {
    p_lfc <- ggplot(scat_df, aes(x = q1_lfc, y = q2_lfc, color = highlight)) +
      geom_point(alpha = 0.5, size = 1.0) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, color = "grey80") +
      geom_vline(xintercept = 0, color = "grey80") +
      scale_color_manual(values = c(
        "Trem2-dependent"   = CATEGORY_COLORS["Cat2"],
        "Trem2-exacerbated" = CATEGORY_COLORS["Cat4"],
        "Trem2-independent" = CATEGORY_COLORS["Cat3"],
        "Redirected"        = CATEGORY_COLORS["Cat6"],
        "Other"             = CATEGORY_COLORS["NS"])) +
      theme_bw() +
      ggtitle(paste0(bt, " — Q1 vs Q2 LFC (broad)")) +
      xlab("log2FC: WT_5XFAD vs WT") +
      ylab("log2FC: Trem2KO_5XFAD vs WT_5XFAD") +
      theme(panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"))
    
    ggsave(file.path(CROSS_DIR, paste0(safe_name, "_lfc_scatter_Q1vQ2.png")),
           plot = p_lfc, width = 8, height = 7, dpi = 200)
  }
  
  # UpSet
  upset_input <- data.frame(
    gene = df$gene,
    Q1 = as.integer(df$q1_sig),
    Q2 = as.integer(df$q2_sig),
    Q3 = as.integer(df$q3_sig)) %>%
    filter(Q1 == 1 | Q2 == 1 | Q3 == 1)
  
  if (nrow(upset_input) >= 5) {
    tryCatch({
      png(file.path(CROSS_DIR, paste0(safe_name, "_upset.png")),
          width = 8, height = 5, units = "in", res = 200)
      print(upset(upset_input, sets = c("Q1", "Q2", "Q3"),
                  order.by = "freq", sets.bar.color = "#4DBBD5",
                  main.bar.color = "grey30", text.scale = 1.3,
                  mainbar.y.label = paste0(bt, " — DE gene intersections")))
      dev.off()
    }, error = function(e) {
      message("  UpSet failed: ", e$message); try(dev.off(), silent = TRUE)
    })
  }
  
  # Venn
  venn_list <- list(
    Q1_Amyloid  = df$gene[df$q1_sig],
    Q2_Trem2_5X = df$gene[df$q2_sig],
    Q3_Trem2    = df$gene[df$q3_sig])
  venn_list <- venn_list[lengths(venn_list) > 0]
  
  if (length(venn_list) >= 2) {
    tryCatch({
      futile.logger::flog.threshold(futile.logger::ERROR)
      venn.diagram(
        x = venn_list,
        filename = file.path(CROSS_DIR, paste0(safe_name, "_venn.png")),
        output = TRUE, imagetype = "png",
        height = 2400, width = 2800, resolution = 300,
        fill = c("#E64B35", "#4DBBD5", "#59A14F")[seq_along(venn_list)],
        alpha = 0.3,
        main = paste0(bt, " — DE gene overlap (broad)"),
        main.cex = 1.2, cat.cex = 0.9, cex = 1.0)
    }, error = function(e) message("  Venn failed: ", e$message))
  }
  
  # Volcano colored by category
  vol_df <- df %>%
    filter(!is.na(q1_lfc) & !is.na(q1_padj)) %>%
    mutate(cat_color = case_when(
      grepl("Cat1", category) ~ "Cat1",
      category == "Cat2_Trem2_dependent"   ~ "Cat2",
      category == "Cat3_Trem2_independent" ~ "Cat3",
      category == "Cat4_Trem2_exacerbated" ~ "Cat4",
      category == "Cat5_Trem2_autonomous"  ~ "Cat5",
      category == "Cat6_Redirected"        ~ "Cat6",
      TRUE ~ "NS"))
  
  p_vol <- ggplot(vol_df, aes(x = q1_lfc, y = -log10(q1_padj), color = cat_color)) +
    geom_point(alpha = 0.5, size = 1.0) +
    scale_color_manual(values = CATEGORY_COLORS, name = "Category") +
    theme_bw() +
    ggtitle(paste0(bt, " — Q1 volcano (broad, colored by category)")) +
    xlab("log2FC (WT_5XFAD vs WT)") + ylab("-log10(padj)") +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
  
  ggsave(file.path(CROSS_DIR, paste0(safe_name, "_volcano_category.png")),
         plot = p_vol, width = 9, height = 7, dpi = 200)
}

# --- Top candidates ---
top_broad <- bind_rows(all_classified) %>%
  filter(category %in% c("Cat2_Trem2_dependent",
                         "Cat4_Trem2_exacerbated",
                         "Cat6_Redirected")) %>%
  mutate(q1_abs_lfc = abs(q1_lfc)) %>%
  arrange(category, desc(q1_abs_lfc)) %>%
  select(broad_type, gene, category,
         q1_lfc, q1_padj, q2_lfc, q2_padj, q3_lfc, q3_padj,
         lfc_ratio_q2q1)

write_csv(top_broad, file.path(CROSS_DIR, "top_candidates_broad.csv"))

cat("\n\n=== Top Broad Cat2 (Trem2-dependent, top 20) ===\n")
print(as.data.frame(top_broad %>% filter(category == "Cat2_Trem2_dependent") %>% head(20)),
      row.names = FALSE)

# =============================================================
# 5. 15m BROAD REGIONAL: Reg1–Reg6 classification
# =============================================================

message("\n\n========== 15m Broad Regional Cross-Comparison ==========\n")

REGIONAL_DIR <- file.path(CROSS_DIR, "regional")
dir.create(REGIONAL_DIR, recursive = TRUE, showWarnings = FALSE)

all_regional <- list()

for (bt in BROAD_TYPES_15m_COMPARE) {
  message("\n--- Regional: ", bt, " ---")
  
  q4a <- load_de_results(bt, Q4AB_DIR, "Q4a_broad_15m_WT_Cor_vs_Hip")
  q4b <- load_de_results(bt, Q4BB_DIR, "Q4b_broad_15m_WT5XFAD_Cor_vs_Hip")
  q4c <- load_de_results(bt, Q4CB_DIR, "Q4c_broad_15m_Trem2KO_Cor_vs_Hip")
  q4d <- load_de_results(bt, Q4DB_DIR, "Q4d_broad_15m_Trem2KO5XFAD_Cor_vs_Hip")
  
  if (is.null(q4a) && is.null(q4b) && is.null(q4c) && is.null(q4d)) next
  
  all_genes <- unique(c(
    if (!is.null(q4a)) q4a$gene else character(0),
    if (!is.null(q4b)) q4b$gene else character(0),
    if (!is.null(q4c)) q4c$gene else character(0),
    if (!is.null(q4d)) q4d$gene else character(0)))
  
  merged <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  merged <- join_comparison(merged, q4a, "q4a")
  merged <- join_comparison(merged, q4b, "q4b")
  merged <- join_comparison(merged, q4c, "q4c")
  merged <- join_comparison(merged, q4d, "q4d")
  merged <- fill_na_sig(merged)
  
  merged <- merged %>%
    mutate(
      reg_category = case_when(
        !q4a_sig & !q4b_sig & !q4c_sig & q4d_sig ~ "Reg6_Compound",
        !q4a_sig & (q4c_sig | q4d_sig)            ~ "Reg5_Trem2_emergent",
        !q4a_sig & q4b_sig                        ~ "Reg2_Amyloid_emergent",
        q4b_sig & !q4d_sig                        ~ "Reg4_Trem2_dependent_regional",
        q4a_sig & !q4b_sig                        ~ "Reg3_Amyloid_lost",
        q4a_sig                                    ~ "Reg1_Constitutive",
        TRUE                                       ~ "Not_regional_DE"),
      broad_type = bt)
  
  all_regional[[bt]] <- merged
  
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
  write_csv(merged, file.path(REGIONAL_DIR,
                              paste0(safe_name, "_regional_categories.csv")))
}

if (length(all_regional) > 0) {
  
  # --- Regional category summary ---
  reg_counts <- bind_rows(lapply(names(all_regional), function(bt) {
    all_regional[[bt]] %>%
      filter(reg_category != "Not_regional_DE") %>%
      count(broad_type, reg_category, name = "n_genes")
  }))
  
  if (nrow(reg_counts) > 0) {
    reg_wide <- reg_counts %>%
      pivot_wider(names_from = reg_category, values_from = n_genes,
                  values_fill = 0) %>%
      mutate(total = rowSums(across(where(is.numeric)))) %>%
      arrange(desc(total))
    
    write_csv(reg_wide, file.path(REGIONAL_DIR, "regional_category_counts_broad.csv"))
    
    cat("\n\n=== Broad Regional Category Counts ===\n")
    print(as.data.frame(reg_wide), row.names = FALSE)
    
    # Heatmap
    reg_long <- reg_counts %>%
      mutate(broad_type = factor(broad_type, levels = rev(reg_wide$broad_type)))
    
    p_reg <- ggplot(reg_long, aes(x = reg_category, y = broad_type, fill = n_genes)) +
      geom_tile(color = "white") +
      geom_text(aes(label = n_genes), size = 4) +
      scale_fill_gradient(low = "white", high = "#E64B35") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white")) +
      ggtitle("15m BROAD: regional gene categories") +
      xlab(NULL) + ylab(NULL)
    
    ggsave(file.path(REGIONAL_DIR, "regional_category_broad_heatmap.png"),
           plot = p_reg, width = 11, height = 5, dpi = 200)
  }
  
  # --- Regional summary matrix ---
  reg_sm_rows <- list()
  for (bt in names(all_regional)) {
    df <- all_regional[[bt]]
    reg_sm_rows[[bt]] <- data.frame(
      broad_type = bt,
      Q4a_total = sum(df$q4a_sig, na.rm = TRUE),
      Q4b_total = sum(df$q4b_sig, na.rm = TRUE),
      Q4c_total = sum(df$q4c_sig, na.rm = TRUE),
      Q4d_total = sum(df$q4d_sig, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
  
  reg_sm <- bind_rows(reg_sm_rows) %>%
    arrange(desc(Q4a_total + Q4b_total + Q4c_total + Q4d_total))
  
  write_csv(reg_sm, file.path(REGIONAL_DIR, "regional_summary_matrix_broad.csv"))
  
  cat("\n\n=== Broad Regional Summary Matrix ===\n")
  print(as.data.frame(reg_sm), row.names = FALSE)
  
  rsm_long <- reg_sm %>%
    pivot_longer(-broad_type, names_to = "comparison", values_to = "n_de") %>%
    mutate(broad_type = factor(broad_type, levels = rev(reg_sm$broad_type)),
           comparison = gsub("_total", "", comparison))
  
  p_rsm <- ggplot(rsm_long, aes(x = comparison, y = broad_type, fill = n_de)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_de), size = 4) +
    scale_fill_gradient(low = "white", high = "#4DBBD5") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) +
    ggtitle("15m BROAD: regional DE genes per cell class × condition") +
    labs(subtitle = "Q4a=WT | Q4b=5XFAD | Q4c=Trem2KO | Q4d=Trem2KO_5XFAD") +
    xlab(NULL) + ylab(NULL)
  
  ggsave(file.path(REGIONAL_DIR, "regional_summary_broad_heatmap.png"),
         plot = p_rsm, width = 8, height = 5, dpi = 200)
  
  # --- Per broad-type regional plots ---
  for (bt in names(all_regional)) {
    df <- all_regional[[bt]]
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", bt)
    
    # LFC scatter Q4a vs Q4b
    scat_df <- df %>% filter(!is.na(q4a_lfc) & !is.na(q4b_lfc))
    if (nrow(scat_df) > 50) {
      scat_df <- scat_df %>%
        mutate(highlight = case_when(
          reg_category == "Reg2_Amyloid_emergent"         ~ "Amyloid-emergent",
          reg_category == "Reg3_Amyloid_lost"             ~ "Amyloid-lost",
          reg_category == "Reg4_Trem2_dependent_regional" ~ "Trem2-dependent",
          reg_category == "Reg1_Constitutive"             ~ "Constitutive",
          TRUE ~ "Other"))
      
      p_scat <- ggplot(scat_df, aes(x = q4a_lfc, y = q4b_lfc, color = highlight)) +
        geom_point(alpha = 0.5, size = 1.0) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = 0, color = "grey80") +
        geom_vline(xintercept = 0, color = "grey80") +
        scale_color_manual(values = c(
          "Constitutive"     = REGIONAL_COLORS["Reg1"],
          "Amyloid-emergent" = REGIONAL_COLORS["Reg2"],
          "Amyloid-lost"     = REGIONAL_COLORS["Reg3"],
          "Trem2-dependent"  = REGIONAL_COLORS["Reg4"],
          "Other"            = REGIONAL_COLORS["NS"])) +
        theme_bw() +
        ggtitle(paste0(bt, " — Q4a vs Q4b regional LFC (broad)")) +
        xlab("Cor vs Hip LFC in WT") + ylab("Cor vs Hip LFC in 5XFAD") +
        theme(panel.background = element_rect(fill = "white"),
              plot.background = element_rect(fill = "white"))
      
      ggsave(file.path(REGIONAL_DIR, paste0(safe_name, "_lfc_scatter_Q4avQ4b.png")),
             plot = p_scat, width = 8, height = 7, dpi = 200)
    }
    
    # UpSet Q4a–Q4d
    upset_df <- data.frame(
      gene = df$gene,
      Q4a_WT = as.integer(df$q4a_sig),
      Q4b_5XFAD = as.integer(df$q4b_sig),
      Q4c_Trem2KO = as.integer(df$q4c_sig),
      Q4d_KO_5XFAD = as.integer(df$q4d_sig)) %>%
      filter(Q4a_WT == 1 | Q4b_5XFAD == 1 | Q4c_Trem2KO == 1 | Q4d_KO_5XFAD == 1)
    
    if (nrow(upset_df) >= 5) {
      tryCatch({
        png(file.path(REGIONAL_DIR, paste0(safe_name, "_regional_upset.png")),
            width = 9, height = 5, units = "in", res = 200)
        print(upset(upset_df,
                    sets = c("Q4a_WT", "Q4b_5XFAD", "Q4c_Trem2KO", "Q4d_KO_5XFAD"),
                    order.by = "freq", sets.bar.color = "#4DBBD5",
                    main.bar.color = "grey30", text.scale = 1.2,
                    mainbar.y.label = paste0(bt, " — regional DE intersections")))
        dev.off()
      }, error = function(e) {
        message("  UpSet failed: ", e$message); try(dev.off(), silent = TRUE)
      })
    }
    
    # Venn Q4a–Q4d
    venn_list <- list(
      Q4a_WT = df$gene[df$q4a_sig],
      Q4b_5XFAD = df$gene[df$q4b_sig],
      Q4c_Trem2KO = df$gene[df$q4c_sig],
      Q4d_KO_5XFAD = df$gene[df$q4d_sig])
    venn_list <- venn_list[lengths(venn_list) > 0]
    
    if (length(venn_list) >= 2) {
      tryCatch({
        futile.logger::flog.threshold(futile.logger::ERROR)
        venn.diagram(
          x = venn_list,
          filename = file.path(REGIONAL_DIR, paste0(safe_name, "_regional_venn.png")),
          output = TRUE, imagetype = "png",
          height = 2600, width = 3000, resolution = 300,
          fill = c("#636363", "#E64B35", "#B07AA1", "#59A14F")[seq_along(venn_list)],
          alpha = 0.3,
          main = paste0(bt, " — regional DE overlap (broad)"),
          main.cex = 1.1, cat.cex = 0.8, cex = 0.9)
      }, error = function(e) message("  Venn failed: ", e$message))
    }
  }
  
  # --- Regional top candidates ---
  reg_top <- bind_rows(all_regional) %>%
    filter(reg_category %in% c("Reg2_Amyloid_emergent",
                               "Reg4_Trem2_dependent_regional",
                               "Reg5_Trem2_emergent",
                               "Reg6_Compound")) %>%
    mutate(max_abs_lfc = pmax(abs(q4a_lfc), abs(q4b_lfc),
                              abs(q4c_lfc), abs(q4d_lfc), na.rm = TRUE)) %>%
    arrange(reg_category, desc(max_abs_lfc)) %>%
    select(broad_type, gene, reg_category,
           q4a_lfc, q4a_padj, q4b_lfc, q4b_padj,
           q4c_lfc, q4c_padj, q4d_lfc, q4d_padj)
  
  write_csv(reg_top, file.path(REGIONAL_DIR, "regional_top_candidates_broad.csv"))
  
  cat("\n\n=== Broad Regional Top Candidates (top 15) ===\n")
  print(as.data.frame(head(reg_top, 15)), row.names = FALSE)
}

# =============================================================
# 6. Summary
# =============================================================

cat("\n\n==================== BROAD CROSS-COMPARISON COMPLETE ====================\n")
cat("\n7m broad types analyzed:  ", length(all_classified), "\n")
cat("15m broad types analyzed: ", length(all_regional), "\n")
cat("\nOutputs in: ", CROSS_DIR, "\n")
cat("  category_counts_broad.csv / _heatmap.png  — 7m Cat1–6\n")
cat("  summary_matrix_broad.csv / _heatmap.png   — 7m DE counts\n")
cat("  top_candidates_broad.csv                  — 7m follow-up genes\n")
cat("  regional/regional_category_counts_broad.csv / _heatmap.png — 15m Reg1–6\n")
cat("  regional/regional_summary_matrix_broad.csv / _heatmap.png  — 15m DE counts\n")
cat("  regional/regional_top_candidates_broad.csv                 — 15m follow-up genes\n")
cat("  <broad_type>_lfc_scatter_*.png, _upset.png, _venn.png, _volcano_category.png\n")