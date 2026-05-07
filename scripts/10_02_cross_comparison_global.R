# ============================================================
# 10_02_cross_comparison_global.R
#
# Cat1–6 cross-comparison on GLOBAL DE (Q1/Q2/Q3 global).
# Treats "Global" as a single pseudo-cell-type.
#
# This is a simpler sibling of script 10 — same categorization
# logic but with one row instead of per-subtype.
#
# Outputs: ./result/de_results/cross_comparison_global/
# ============================================================

source("./scripts/09_00_de_functions.R")

suppressPackageStartupMessages({
  if (!requireNamespace("UpSetR", quietly = TRUE))
    install.packages("UpSetR", quiet = TRUE)
  library(UpSetR)
})

set.seed(1234)

CROSS_DIR <- make_de_dir("cross_comparison_global")
EXACERBATION_RATIO <- 1.5

# =============================================================
# 1. Load global Q1/Q2/Q3
# =============================================================

load_global <- function(q_name) {
  path <- file.path(DE_BASE_DIR, q_name,
                    paste0(q_name, "_Global_full_results.csv"))
  if (!file.exists(path)) {
    message("  Missing: ", path)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE)
}

q1 <- load_global("Q1_global_WT_vs_WT5XFAD")
q2 <- load_global("Q2_global_WT5XFAD_vs_Trem2KO5XFAD")
q3 <- load_global("Q3_global_WT_vs_Trem2KO")

if (is.null(q1) || is.null(q2) || is.null(q3))
  stop("One or more global DE tables missing — run 09_07 first")

# =============================================================
# 2. Merge + classify
# =============================================================

all_genes <- unique(c(q1$gene, q2$gene, q3$gene))
merged <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

merged <- merged %>%
  left_join(q1 %>% select(gene, q1_lfc = log2FoldChange, q1_padj = padj,
                          q1_sig = significant, q1_dir = direction),
            by = "gene") %>%
  left_join(q2 %>% select(gene, q2_lfc = log2FoldChange, q2_padj = padj,
                          q2_sig = significant, q2_dir = direction),
            by = "gene") %>%
  left_join(q3 %>% select(gene, q3_lfc = log2FoldChange, q3_padj = padj,
                          q3_sig = significant, q3_dir = direction),
            by = "gene") %>%
  mutate(across(ends_with("_sig"), ~ replace_na(.x, FALSE)),
         across(ends_with("_dir"), ~ replace_na(.x, "NS")))

merged <- merged %>%
  mutate(
    same_direction_q1q2 = (q1_dir == q2_dir) & q1_dir != "NS",
    reversed_q1q2 = (q1_dir == "UP" & q2_dir == "DOWN") |
      (q1_dir == "DOWN" & q2_dir == "UP"),
    lfc_ratio_q2q1 = ifelse(
      !is.na(q1_lfc) & !is.na(q2_lfc) & abs(q1_lfc) > 0,
      abs(q2_lfc) / abs(q1_lfc), NA_real_
    ),
    category = case_when(
      !q1_sig & q2_sig                                      ~ "Cat6_Redirected",
      !q1_sig & q3_sig                                      ~ "Cat5_Trem2_autonomous",
      q1_sig & q2_sig & same_direction_q1q2 &
        !is.na(lfc_ratio_q2q1) &
        lfc_ratio_q2q1 > EXACERBATION_RATIO                 ~ "Cat4_Trem2_exacerbated",
      q1_sig & q2_sig & reversed_q1q2                       ~ "Cat2_Trem2_dependent",
      q1_sig & q2_sig & same_direction_q1q2 &
        !is.na(lfc_ratio_q2q1) &
        lfc_ratio_q2q1 <= EXACERBATION_RATIO                ~ "Cat2_Trem2_dependent",
      q1_sig & !q2_sig & !q3_sig                              ~ "Cat1_Pure_amyloid",
      q1_sig & !q2_sig & q3_sig                               ~ "Cat3_Trem2_independent",
      TRUE                                                     ~ "Not_DE"
    ),
    cell_type = "Global"
  )

write_csv(merged, file.path(CROSS_DIR, "Global_gene_categories.csv"))

# =============================================================
# 3. Category summary
# =============================================================

cat_counts <- merged %>%
  filter(category != "Not_DE") %>%
  count(category, name = "n_genes") %>%
  arrange(desc(n_genes))

write_csv(cat_counts,
          file.path(CROSS_DIR, "Global_category_counts.csv"))

cat("\n=== Global Cross-Comparison Categories ===\n")
print(as.data.frame(cat_counts), row.names = FALSE)

# Bar plot
p <- ggplot(cat_counts,
            aes(x = reorder(category, n_genes), y = n_genes,
                fill = category)) +
  geom_col(alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(values = CATEGORY_COLORS) +
  theme_bw(base_size = 11) +
  ggtitle("Global DE: gene categories") +
  xlab(NULL) + ylab("Number of genes") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

ggsave(file.path(CROSS_DIR, "Global_category_bar.png"),
       plot = p, width = 9, height = 5, dpi = 200)

# =============================================================
# 4. LFC scatter Q1 vs Q2
# =============================================================

scat_df <- merged %>%
  filter(!is.na(q1_lfc) & !is.na(q2_lfc)) %>%
  mutate(cat_color = case_when(
    grepl("Cat1", category) ~ "Cat1",
    category == "Cat2_Trem2_dependent"   ~ "Cat2",
    category == "Cat3_Trem2_independent" ~ "Cat3",
    category == "Cat4_Trem2_exacerbated" ~ "Cat4",
    category == "Cat5_Trem2_autonomous"  ~ "Cat5",
    category == "Cat6_Redirected"        ~ "Cat6",
    TRUE ~ "NS"
  ))

p_sc <- ggplot(scat_df,
               aes(x = q1_lfc, y = q2_lfc, color = cat_color)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  scale_color_manual(values = CATEGORY_COLORS, name = "Category") +
  theme_bw(base_size = 11) +
  ggtitle("Global Q1 vs Q2 LFC") +
  xlab("log2FC: WT_5XFAD vs WT (Q1)") +
  ylab("log2FC: Trem2KO_5XFAD vs WT_5XFAD (Q2)") +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

ggsave(file.path(CROSS_DIR, "Global_lfc_scatter_Q1vQ2.png"),
       plot = p_sc, width = 8, height = 7, dpi = 200)

cat("\n\nOutputs in: ", CROSS_DIR, "\n")