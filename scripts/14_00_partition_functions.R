# ============================================================
# 14_00_partition_functions.R
#
# Build gene set partitions from multi-comparison DE overlaps
# and run ORA on each non-empty partition.
#
# For any cell type with DE across Q1/Q2/Q3 (or Q4a-d), creates
# up to 2^n - 1 partitions from the n-way intersection, then
# runs multi-database ORA on each partition that meets a minimum
# gene count threshold.
#
# Reuses 11_00_pathway_functions.R for all ORA calls.
#
# Source:
#   source("./scripts/14_00_partition_functions.R")
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/11_00_pathway_functions.R")

suppressPackageStartupMessages({
  if (!requireNamespace("VennDiagram", quietly = TRUE))
    install.packages("VennDiagram", quiet = TRUE)
  library(VennDiagram)
})

# =============================================================
# 1. Settings
# =============================================================

PARTITION_MIN_GENES <- 5     # minimum genes per partition to run ORA
PARTITION_MIN_ORA   <- 10    # minimum for "interpretable" ORA

# =============================================================
# 2. Build partitions from DE tables
# =============================================================

#' Build gene set partitions from multiple DE comparisons
#'
#' Takes a named list of DE tables (e.g., list(Q1=df1, Q2=df2, Q3=df3)),
#' extracts significant genes from each, and creates all non-empty
#' intersection partitions.
#'
#' Works with any number of inputs (2, 3, 4, ...).
#'
#' @param de_tables Named list of DE data.frames (full results)
#' @param padj_cutoff Significance threshold
#' @param lfc_cutoff  Minimum absolute log2FC
#' @return list with:
#'   - partitions: named list of character vectors (gene names)
#'   - summary: data.frame with partition name, member comparisons, gene count
#'   - sig_genes: named list of sig gene sets per comparison (for reference)
build_partitions <- function(de_tables,
                             padj_cutoff = ORA_PADJ_CUTOFF,
                             lfc_cutoff  = ORA_LFC_CUTOFF) {
  
  comp_names <- names(de_tables)
  n_comps    <- length(comp_names)
  
  # Extract significant genes from each comparison
  sig_genes <- lapply(de_tables, function(df) {
    if (is.null(df)) return(character(0))
    df %>%
      filter(!is.na(padj), padj < padj_cutoff,
             abs(log2FoldChange) >= lfc_cutoff) %>%
      pull(gene) %>% unique()
  })
  
  # All unique genes across comparisons
  all_sig <- unique(unlist(sig_genes))
  if (length(all_sig) == 0) {
    message("    No significant genes in any comparison")
    return(list(partitions = list(), summary = data.frame(), sig_genes = sig_genes))
  }
  
  # Build membership matrix: gene × comparison (TRUE/FALSE)
  membership <- sapply(sig_genes, function(g) all_sig %in% g)
  
  # Force matrix when only 1 gene (sapply returns a vector)
  if (!is.matrix(membership)) {
    membership <- matrix(membership, nrow = length(all_sig),
                         ncol = length(sig_genes),
                         dimnames = list(all_sig, comp_names))
  } else {
    rownames(membership) <- all_sig
  }
  
  # Generate partition labels from membership patterns
  # Each unique combination of TRUE/FALSE across comparisons = one partition
  patterns <- apply(membership, 1, function(row) {
    paste(comp_names[row], collapse = " ∩ ")
  })
  
  # Remove genes not in any comparison (shouldn't happen but safety)
  patterns <- patterns[patterns != ""]
  
  # Split genes by pattern
  partitions <- split(names(patterns), patterns)
  
  # Also create "X_only" labels for single-comparison partitions
  partition_labels <- sapply(names(partitions), function(p) {
    members <- strsplit(p, " ∩ ")[[1]]
    if (length(members) == 1) {
      paste0(members, "_only")
    } else if (length(members) == n_comps) {
      paste0("all_", paste(comp_names, collapse = "_"))
    } else {
      gsub(" ∩ ", "_AND_", p)
    }
  })
  
  names(partitions) <- partition_labels
  
  # Summary table
  summary_df <- data.frame(
    partition    = partition_labels,
    members      = names(partition_labels),
    n_genes      = sapply(partitions, length),
    can_run_ora  = sapply(partitions, length) >= PARTITION_MIN_GENES,
    confidence   = sapply(partitions, function(g) {
      n <- length(g)
      if (n < PARTITION_MIN_GENES) "SKIP"
      else if (n < PARTITION_MIN_ORA) "exploratory"
      else "interpretable"
    }),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  summary_df <- summary_df %>% arrange(desc(n_genes))
  
  message("    Partitions: ", nrow(summary_df), " total, ",
          sum(summary_df$can_run_ora), " with enough genes for ORA")
  
  list(partitions = partitions, summary = summary_df, sig_genes = sig_genes)
}

# =============================================================
# 3. Run ORA on all partitions for one cell type
# =============================================================

#' @param partitions Output from build_partitions()
#' @param universe   Background gene universe (all tested genes)
#' @param label      Cell type label
#' @param out_dir    Output directory
#' @param is_neuronal Whether to run SynGO
#' @return list of ORA results per partition
run_ora_on_partitions <- function(partitions,
                                  universe,
                                  label,
                                  out_dir,
                                  is_neuronal = NULL) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  safe_label <- gsub("[^A-Za-z0-9_-]", "_", label)
  
  if (is.null(is_neuronal))
    is_neuronal <- label %in% SYNGO_NEURONAL_TYPES
  
  part_list <- partitions$partitions
  part_summary <- partitions$summary
  
  # Save partition summary
  write_csv(part_summary,
            file.path(out_dir,
                      paste0(safe_label, "_partition_summary.csv")))
  
  # Save gene lists for all partitions (even tiny ones)
  for (pname in names(part_list)) {
    safe_pname <- gsub("[^A-Za-z0-9_-]", "_", pname)
    writeLines(part_list[[pname]],
               file.path(out_dir,
                         paste0(safe_label, "_", safe_pname, "_genes.txt")))
  }
  
  # Run ORA on partitions meeting threshold
  all_ora_results <- list()
  all_summaries   <- list()
  
  for (pname in names(part_list)) {
    genes <- part_list[[pname]]
    n_g   <- length(genes)
    
    if (n_g < PARTITION_MIN_GENES) {
      message("    ", pname, ": ", n_g, " genes — too few, skipping ORA")
      next
    }
    
    safe_pname <- gsub("[^A-Za-z0-9_-]", "_", pname)
    message("\n    ", pname, ": ", n_g, " genes — running ORA")
    
    ora_out <- run_all_ora(genes, universe, run_syngo_flag = is_neuronal)
    ora_res <- ora_out$results
    
    if (length(ora_res) > 0) {
      all_ora_results[[pname]] <- ora_res
      
      # Save per-database
      for (db in names(ora_res))
        write_csv(as.data.frame(ora_res[[db]]),
                  file.path(out_dir,
                            paste0(safe_label, "_", safe_pname,
                                   "_", db, ".csv")))
      
      # Standardized summary
      std <- standardize_ora_result(
        ora_res,
        label       = label,
        comparison  = pname,
        direction   = "partition",
        n_input_genes = n_g
      )
      if (!is.null(std) && nrow(std) > 0)
        all_summaries[[pname]] <- std
    }
  }
  
  # Combined ORA summary across all partitions
  combined <- bind_rows(all_summaries)
  if (nrow(combined) > 0) {
    write_csv(combined,
              file.path(out_dir,
                        paste0(safe_label, "_partition_ORA_ALL.csv")))
  }
  
  list(
    ora_results = all_ora_results,
    summary     = combined,
    partition_info = part_summary
  )
}

# =============================================================
# 3b. Venn diagram of comparison overlaps
# =============================================================

#' Draw Venn diagram showing gene count overlaps across comparisons
#'
#' @param sig_genes Named list of significant gene vectors per comparison
#' @param title  Plot title
#' @param filepath Output PNG path
plot_partition_venn <- function(sig_genes, title, filepath) {
  # Remove empty sets
  sig_genes <- sig_genes[lengths(sig_genes) > 0]
  if (length(sig_genes) < 2) return(invisible(NULL))
  
  venn_colors <- c("#E64B35", "#4DBBD5", "#59A14F", "#B07AA1")[seq_along(sig_genes)]
  
  tryCatch({
    futile.logger::flog.threshold(futile.logger::ERROR)
    venn.diagram(
      x         = sig_genes,
      filename  = filepath,
      output    = TRUE,
      imagetype = "png",
      height    = 2400,
      width     = 2800,
      resolution = 300,
      fill      = venn_colors,
      alpha     = 0.3,
      main      = title,
      main.cex  = 1.1,
      cat.cex   = 0.9,
      cex       = 1.0
    )
    message("    Venn saved: ", filepath)
  }, error = function(e) {
    message("    Venn failed: ", e$message)
  })
}

# =============================================================
# 3c. ORA summary with confidence coloring
# =============================================================

CONFIDENCE_COLORS <- c(
  "interpretable" = "#2CA02C",   # green
  "exploratory"   = "#FFB347",   # amber/orange
  "SKIP"          = "#CCCCCC"    # grey
)

#' Combined ORA summary plot: one panel per partition, confidence
#' as background color, top ORA terms listed.
#'
#' @param partition_info partition summary data.frame
#' @param ora_summary combined ORA results (from standardize_ora_result)
#' @param title Plot title
#' @param top_n Top terms per partition
plot_partition_ora_with_confidence <- function(ora_summary,
                                               title = "Partition ORA",
                                               top_n = 8) {
  
  if (is.null(ora_summary) || nrow(ora_summary) == 0) return(NULL)
  
  # Get top terms per partition
  top_terms <- ora_summary %>%
    group_by(comparison) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  if (nrow(top_terms) == 0) return(NULL)
  
  # confidence and n_input_genes are already in ora_summary
  # from standardize_ora_result() — no join needed
  top_terms <- top_terms %>%
    mutate(
      facet_label = paste0(comparison, "\n(",
                           n_input_genes, " genes, ", confidence, ")"),
      Description = ifelse(nchar(Description) > 50,
                           paste0(substr(Description, 1, 47), "..."),
                           Description)
    )
  
  # Order facets by gene count (largest first)
  facet_order <- top_terms %>%
    distinct(comparison, n_input_genes) %>%
    arrange(desc(n_input_genes)) %>%
    pull(comparison)
  
  top_terms <- top_terms %>%
    mutate(
      comparison = factor(comparison, levels = facet_order),
      facet_label = factor(facet_label,
                           levels = unique(facet_label[order(match(comparison, facet_order))]))
    )
  
  # Database shapes
  db_shapes <- c(GO_BP = 16, KEGG = 17, Reactome = 15, SynGO = 18)
  
  p <- ggplot(top_terms, aes(x = Count, y = reorder(Description, Count))) +
    geom_col(aes(fill = confidence), alpha = 0.75) +
    geom_point(aes(shape = database), size = 2.5, color = "grey20") +
    facet_wrap(~ facet_label, scales = "free", ncol = 1) +
    scale_fill_manual(values = CONFIDENCE_COLORS, name = "Confidence") +
    scale_shape_manual(values = db_shapes, name = "Database") +
    theme_bw(base_size = 10) +
    ggtitle(title) +
    xlab("Gene count") + ylab(NULL) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold", size = 10),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.position = "right"
    )
  
  p
}

# =============================================================
# 4. Full pipeline: cell type + DE tables → partitions → ORA
# =============================================================

#' @param de_tables Named list (e.g., list(Q1=df, Q2=df, Q3=df))
#' @param label     Cell type name
#' @param out_dir   Output directory
#' @param universe  Background genes (if NULL, union of all tested genes)
run_partition_analysis <- function(de_tables,
                                   label,
                                   out_dir,
                                   universe     = NULL,
                                   is_neuronal  = NULL,
                                   padj_cutoff  = ORA_PADJ_CUTOFF,
                                   lfc_cutoff   = ORA_LFC_CUTOFF) {
  
  message("\n=== Partition ORA: ", label, " ===")
  safe_label <- gsub("[^A-Za-z0-9_-]", "_", label)
  
  # Remove NULL entries
  de_tables <- de_tables[!sapply(de_tables, is.null)]
  if (length(de_tables) < 2) {
    message("  Need at least 2 comparisons — skipping")
    return(NULL)
  }
  
  # Build universe from union of all tested genes
  if (is.null(universe)) {
    universe <- unique(unlist(lapply(de_tables, function(df) df$gene)))
  }
  
  # Build partitions
  parts <- build_partitions(de_tables,
                            padj_cutoff = padj_cutoff,
                            lfc_cutoff  = lfc_cutoff)
  
  if (length(parts$partitions) == 0) {
    message("  No partitions to analyze")
    return(NULL)
  }
  
  # Only run ORA if there are partitions with >= 2 distinct sets
  # (single-comparison partitions are already covered by script 11)
  multi_part_count <- sum(grepl("_AND_|^all_", names(parts$partitions)))
  single_part_count <- sum(grepl("_only$", names(parts$partitions)))
  
  message("  ", single_part_count, " single-comparison partitions, ",
          multi_part_count, " intersection partitions")
  
  # Run ORA on all qualifying partitions
  results <- run_ora_on_partitions(
    partitions  = parts,
    universe    = universe,
    label       = label,
    out_dir     = out_dir,
    is_neuronal = is_neuronal
  )
  
  # --- Venn diagram of comparison overlaps ---
  plot_partition_venn(
    sig_genes = parts$sig_genes,
    title     = paste0(label, " — DE gene overlaps"),
    filepath  = file.path(out_dir,
                          paste0(safe_label, "_partition_venn.png"))
  )
  
  # --- ORA summary with confidence coloring ---
  p_ora <- plot_partition_ora_with_confidence(
    ora_summary    = results$summary,
    title          = paste0(label, " — Partition ORA")
  )
  if (!is.null(p_ora)) {
    n_parts <- sum(parts$summary$can_run_ora)
    ggsave(file.path(out_dir,
                     paste0(safe_label, "_partition_ora_confidence.png")),
           plot = p_ora,
           width = 12, height = 2.5 + n_parts * 3, dpi = 200)
  }
  
  results
}

# =============================================================
# 5. Cross-cell-type pathway heatmap
# =============================================================

#' For a given partition type (e.g., "Q1_only"), collect ORA results
#' across all cell types that had that partition and make a heatmap
#' showing top pathways × cell types.
plot_partition_pathway_heatmap <- function(all_results,
                                           partition_pattern,
                                           title = NULL,
                                           top_n = 10) {
  
  # Collect ORA summaries matching the partition pattern
  matching <- list()
  
  for (ct in names(all_results)) {
    res <- all_results[[ct]]
    if (is.null(res$summary) || nrow(res$summary) == 0) next
    
    rows <- res$summary %>%
      filter(grepl(partition_pattern, comparison, fixed = TRUE))
    
    if (nrow(rows) > 0) {
      matching[[ct]] <- rows %>%
        mutate(cell_type = ct)
    }
  }
  
  combined <- bind_rows(matching)
  if (nrow(combined) == 0) return(NULL)
  
  # Get top N pathways by frequency across cell types
  top_terms <- combined %>%
    count(Description, sort = TRUE) %>%
    head(top_n) %>%
    pull(Description)
  
  plot_df <- combined %>%
    filter(Description %in% top_terms) %>%
    mutate(
      Description = ifelse(nchar(Description) > 50,
                           paste0(substr(Description, 1, 47), "..."),
                           Description),
      neg_log_padj = -log10(p.adjust)
    ) %>%
    select(cell_type, Description, neg_log_padj) %>%
    distinct(cell_type, Description, .keep_all = TRUE)
  
  if (is.null(title))
    title <- paste0("Partition: ", partition_pattern)
  
  ggplot(plot_df, aes(x = cell_type, y = Description, fill = neg_log_padj)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "#E64B35",
                        name = "-log10(padj)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) +
    ggtitle(title) +
    xlab(NULL) + ylab(NULL)
}

#' Partition size summary across cell types
plot_partition_sizes <- function(all_results, title = "Partition sizes") {
  
  sizes <- bind_rows(lapply(names(all_results), function(ct) {
    info <- all_results[[ct]]$partition_info
    if (is.null(info) || nrow(info) == 0) return(NULL)
    info %>% mutate(cell_type = ct)
  }))
  
  if (nrow(sizes) == 0) return(NULL)
  
  ggplot(sizes, aes(x = cell_type, y = partition, fill = n_genes)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_genes), size = 2.8) +
    scale_fill_gradient(low = "white", high = "#4DBBD5", name = "Genes") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")) +
    ggtitle(title) +
    xlab(NULL) + ylab(NULL)
}

# =============================================================
# 6. Helper: load DE table
# =============================================================

load_de_for_partition <- function(comparison_name, label) {
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", label)
  path <- file.path(DE_BASE_DIR, comparison_name,
                    paste0(comparison_name, "_", safe_name,
                           "_full_results.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}

message("14_00_partition_functions.R loaded")
message("  Min genes for ORA: ", PARTITION_MIN_GENES,
        " | Interpretable: ", PARTITION_MIN_ORA)