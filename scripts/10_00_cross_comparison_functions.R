# ============================================================
# 10_00_cross_comparison_functions.R
#
# Shared functions for cross-comparison between cell-type-level
# DE and whole-population ("global") DE.
#
# Used by:
#   10_03_cross_comparison_broad_vs_global.R
#   10_04_cross_comparison_subtype_vs_global.R
#
# Core function: classify_vs_global()
#   For a given cell type's DE table + the matching global DE
#   table (same comparison), classifies every gene as:
#     - Cell-type-specific (DE here, NOT global)
#     - Shared with global (DE both, same direction)
#     - Global-only (DE global, NOT here)
#     - Direction mismatch (DE both, opposite directions)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  if (!requireNamespace("VennDiagram", quietly = TRUE))
    install.packages("VennDiagram", quiet = TRUE)
  library(VennDiagram)
})

select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
count   <- dplyr::count
desc    <- dplyr::desc

# =============================================================
# 1. Category colors for vs-global analysis
# =============================================================

VS_GLOBAL_COLORS <- c(
  "Cell_type_specific"   = "#E64B35",  # red — the most interesting
  "Shared_with_global"   = "#4DBBD5",  # blue — robust across scales
  "Global_only"          = "#8C8C8C",  # grey — driven by other cells
  "Direction_mismatch"   = "#B07AA1",  # purple — opposite direction
  "Not_DE"               = "#E5E5E5"
)

# =============================================================
# 2. Classify vs global
# =============================================================

#' Classify every gene as specific/shared/global-only/mismatch
#'
#' @param ct_de_table    Full DE table for a cell type (from 09_xx output)
#' @param global_de_table Full DE table for global (from 09_07 output)
#' @param cell_label     String label for this cell type (e.g., "Pvalb")
#' @param padj_cutoff    Significance threshold
#' @param lfc_cutoff     Min absolute log2FC
#' @return data.frame with gene-level classification
classify_vs_global <- function(ct_de_table, global_de_table,
                               cell_label,
                               padj_cutoff = 0.05,
                               lfc_cutoff  = 0.25) {
  
  # Extract just gene, lfc, padj, and derive significance
  ct <- ct_de_table %>%
    select(gene, ct_lfc = log2FoldChange, ct_padj = padj) %>%
    mutate(
      ct_sig = !is.na(ct_padj) & ct_padj < padj_cutoff &
        abs(ct_lfc) >= lfc_cutoff,
      ct_dir = case_when(
        !ct_sig            ~ "NS",
        ct_lfc > 0         ~ "UP",
        ct_lfc < 0         ~ "DOWN",
        TRUE               ~ "NS"
      )
    )
  
  gl <- global_de_table %>%
    select(gene, gl_lfc = log2FoldChange, gl_padj = padj) %>%
    mutate(
      gl_sig = !is.na(gl_padj) & gl_padj < padj_cutoff &
        abs(gl_lfc) >= lfc_cutoff,
      gl_dir = case_when(
        !gl_sig            ~ "NS",
        gl_lfc > 0         ~ "UP",
        gl_lfc < 0         ~ "DOWN",
        TRUE               ~ "NS"
      )
    )
  
  merged <- full_join(ct, gl, by = "gene") %>%
    mutate(
      ct_sig = replace_na(ct_sig, FALSE),
      gl_sig = replace_na(gl_sig, FALSE),
      ct_dir = replace_na(ct_dir, "NS"),
      gl_dir = replace_na(gl_dir, "NS")
    )
  
  merged <- merged %>%
    mutate(
      vs_global_category = case_when(
        ct_sig & gl_sig & ct_dir == gl_dir             ~ "Shared_with_global",
        ct_sig & gl_sig & ct_dir != gl_dir             ~ "Direction_mismatch",
        ct_sig & !gl_sig                                ~ "Cell_type_specific",
        !ct_sig & gl_sig                                ~ "Global_only",
        TRUE                                             ~ "Not_DE"
      ),
      cell_type = cell_label
    )
  
  merged
}

# =============================================================
# 3. Category counts per cell type
# =============================================================

count_vs_global_categories <- function(classified_df) {
  classified_df %>%
    filter(vs_global_category != "Not_DE") %>%
    count(cell_type, vs_global_category, name = "n_genes")
}

# =============================================================
# 4. Plots
# =============================================================

#' Heatmap: cell types (y) × categories (x), cells show gene counts
plot_vs_global_heatmap <- function(count_df, title = "vs Global") {
  if (nrow(count_df) == 0) return(NULL)
  
  cat_levels <- c("Cell_type_specific", "Shared_with_global",
                  "Global_only", "Direction_mismatch")
  
  count_wide <- count_df %>%
    pivot_wider(names_from = vs_global_category, values_from = n_genes,
                values_fill = 0) %>%
    mutate(total = rowSums(across(where(is.numeric)))) %>%
    arrange(desc(total))
  
  plot_df <- count_df %>%
    mutate(cell_type = factor(cell_type,
                              levels = rev(count_wide$cell_type)),
           vs_global_category = factor(
             vs_global_category,
             levels = intersect(cat_levels, unique(vs_global_category))
           ))
  
  ggplot(plot_df, aes(x = vs_global_category, y = cell_type,
                      fill = n_genes)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_genes), size = 3.5) +
    scale_fill_gradient(low = "white", high = "#E64B35") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) +
    ggtitle(title) +
    xlab(NULL) + ylab(NULL)
}

#' LFC scatter: cell type LFC (x) vs global LFC (y), colored by category
plot_vs_global_scatter <- function(classified_df, title = "Cell-type vs Global LFC") {
  
  plot_df <- classified_df %>%
    filter(!is.na(ct_lfc) & !is.na(gl_lfc)) %>%
    filter(vs_global_category != "Not_DE") %>%
    mutate(vs_global_category = factor(
      vs_global_category,
      levels = names(VS_GLOBAL_COLORS)
    ))
  
  if (nrow(plot_df) == 0) return(NULL)
  
  ggplot(plot_df, aes(x = ct_lfc, y = gl_lfc, color = vs_global_category)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    scale_color_manual(values = VS_GLOBAL_COLORS, name = "Category") +
    theme_bw(base_size = 11) +
    ggtitle(title) +
    xlab("log2FC (cell-type)") + ylab("log2FC (global)") +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
}

# =============================================================
# 5. Venn diagram: cell-type vs global gene overlap
# =============================================================

#' Draw a 2-way Venn for cell-type DE genes vs global DE genes
#'
#' @param classified_df Output from classify_vs_global()
#' @param title  Plot title
#' @param filepath Output PNG path
plot_vs_global_venn <- function(classified_df, title, filepath) {
  ct_label <- unique(classified_df$cell_type)[1]
  
  ct_genes <- classified_df %>%
    filter(ct_sig) %>% pull(gene)
  gl_genes <- classified_df %>%
    filter(gl_sig) %>% pull(gene)
  
  if (length(ct_genes) + length(gl_genes) < 2) return(invisible(NULL))
  
  venn_list <- list(ct_genes, gl_genes)
  names(venn_list) <- c(ct_label, "Global")
  venn_list <- venn_list[lengths(venn_list) > 0]
  
  if (length(venn_list) < 2) return(invisible(NULL))
  
  tryCatch({
    futile.logger::flog.threshold(futile.logger::ERROR)
    venn.diagram(
      x         = venn_list,
      filename  = filepath,
      output    = TRUE,
      imagetype = "png",
      height    = 2000,
      width     = 2400,
      resolution = 300,
      fill      = c("#E64B35", "#4DBBD5"),
      alpha     = 0.35,
      main      = title,
      main.cex  = 1.1,
      cat.cex   = 0.9,
      cex       = 1.0
    )
  }, error = function(e) {
    message("  Venn failed: ", e$message)
  })
}

# =============================================================
# 6. Load DE table helper
# =============================================================

load_de_table <- function(comparison_name, label) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", label)
  path <- file.path(DE_BASE_DIR, comparison_name,
                    paste0(comparison_name, "_", safe_name,
                           "_full_results.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}

message("10_00_cross_comparison_functions.R loaded")