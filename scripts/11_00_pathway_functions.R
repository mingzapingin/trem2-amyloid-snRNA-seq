# ============================================================
# 11_00_pathway_functions.R
#
# Direction-aware multi-database ORA (GO BP, KEGG, Reactome, SynGO).
#
# Design rules:
#   1. Never run ORA on mixed "all" genes — always up and down separately
#   2. Skip ORA when a direction has < 10 genes
#   3. Mark confidence: <10 skip, 10–19 exploratory, 20+ interpretable
#   4. Prioritize enriched terms with Count ≥ 3
#   5. SynGO: mouse→human orthologs, report mapping loss, neuronal only
#   6. Combined summary table with direction column
#
# Requires:
#   BiocManager::install(c("clusterProfiler", "org.Mm.eg.db",
#                          "org.Hs.eg.db", "ReactomePA"))
#
# SynGO files:
#   ./metadata/syngo_TERM2GENE.csv
#   ./metadata/syngo_TERM2NAME.csv
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

HAS_REACTOME <- requireNamespace("ReactomePA", quietly = TRUE)
if (HAS_REACTOME) suppressPackageStartupMessages(library(ReactomePA))

HAS_ORG_HS <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
if (HAS_ORG_HS) suppressPackageStartupMessages(library(org.Hs.eg.db))

# Fix namespace conflicts
select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
count   <- dplyr::count
desc    <- dplyr::desc

# =============================================================
# 1. Settings
# =============================================================

ORA_PADJ_CUTOFF   <- 0.05    # for reading DE genes
ORA_LFC_CUTOFF    <- 0.25

ORA_PVALUE_CUTOFF <- 0.05    # enrichment significance
ORA_QVALUE_CUTOFF <- 0.1
ORA_MIN_GENE_SIZE <- 10
ORA_MAX_GENE_SIZE <- 500
ORA_TOP_N         <- 15

KEGG_ORGANISM <- "mmu"

# Confidence bins (gene list size)
MIN_GENES_TO_RUN  <- 10      # below this: skip ORA entirely
MIN_GENES_STRONG  <- 20      # above this: "interpretable"
MIN_TERM_COUNT    <- 3       # minimum gene overlap for a strong term

SYNGO_TERM2GENE <- "./metadata/syngo_TERM2GENE.csv"
SYNGO_TERM2NAME <- "./metadata/syngo_TERM2NAME.csv"

# Neuronal broad types — SynGO only runs for these
SYNGO_NEURONAL_TYPES <- c(
  "Excitatory_neuron", "Inhibitory_neuron",
  # Allen subclass names that are neuronal
  "L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L5 PT CTX",
  "L5/6 NP CTX", "L6 IT CTX", "L6 CT CTX", "L6b CTX",
  "L2/3 IT PPP", "L2/3 IT ENTl", "L6b/CT ENT", "L2 IT ENTl",
  "L3 IT ENT", "Car3", "DG", "CA1-ProS",
  "Pvalb", "Sst", "Meis2", "Lamp5", "Vip", "Sncg"
)

# =============================================================
# 2. Confidence classification
# =============================================================

classify_confidence <- function(n_genes) {
  case_when(
    n_genes < MIN_GENES_TO_RUN ~ "SKIP",
    n_genes < MIN_GENES_STRONG ~ "exploratory",
    TRUE                       ~ "interpretable"
  )
}

flag_term_strength <- function(count) {
  case_when(
    count >= 4 ~ "strong",
    count >= 3 ~ "usable",
    count >= 2 ~ "weak",
    TRUE       ~ "very_weak"
  )
}

# =============================================================
# 3. Mouse → Human ortholog conversion
# =============================================================

.ortho_cache <- new.env(parent = emptyenv())

build_ortholog_table <- function() {
  if (exists("table", envir = .ortho_cache))
    return(get("table", envir = .ortho_cache))
  
  message("  Building mouse→human ortholog table...")
  
  if (HAS_ORG_HS) {
    mm_map <- tryCatch({
      bitr(keys(org.Mm.eg.db, keytype = "SYMBOL"),
           fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Mm.eg.db) %>%
        distinct(SYMBOL, .keep_all = TRUE)
    }, error = function(e) NULL)
    
    hs_map <- tryCatch({
      bitr(keys(org.Hs.eg.db, keytype = "SYMBOL"),
           fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Hs.eg.db) %>%
        distinct(SYMBOL, .keep_all = TRUE)
    }, error = function(e) NULL)
    
    if (!is.null(mm_map) && !is.null(hs_map)) {
      mm_map$upper <- toupper(mm_map$SYMBOL)
      hs_map$upper <- toupper(hs_map$SYMBOL)
      
      ortho <- mm_map %>%
        select(mouse_symbol = SYMBOL, upper) %>%
        inner_join(hs_map %>% select(human_symbol = SYMBOL, upper),
                   by = "upper") %>%
        select(mouse_symbol, human_symbol) %>%
        distinct(mouse_symbol, .keep_all = TRUE)
      
      message("    Ortholog table: ", nrow(ortho), " mappings")
      assign("table", ortho, envir = .ortho_cache)
      return(ortho)
    }
  }
  
  message("    Using toupper() fallback")
  assign("table", NULL, envir = .ortho_cache)
  NULL
}

mouse_to_human <- function(mouse_genes) {
  ortho_table <- build_ortholog_table()
  
  if (!is.null(ortho_table)) {
    mapped <- ortho_table %>% filter(mouse_symbol %in% mouse_genes)
    unmapped <- setdiff(mouse_genes, mapped$mouse_symbol)
    fallback <- data.frame(mouse_symbol = unmapped,
                           human_symbol = toupper(unmapped),
                           stringsAsFactors = FALSE)
    all_mapped <- bind_rows(mapped, fallback)
  } else {
    all_mapped <- data.frame(mouse_symbol = mouse_genes,
                             human_symbol = toupper(mouse_genes),
                             stringsAsFactors = FALSE)
  }
  setNames(all_mapped$human_symbol, all_mapped$mouse_symbol)
}

# =============================================================
# 4. Mouse Entrez ID conversion
# =============================================================

symbols_to_entrez <- function(symbols) {
  mapping <- tryCatch({
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Mm.eg.db)
  }, error = function(e) {
    data.frame(SYMBOL = character(0), ENTREZID = character(0))
  })
  if (nrow(mapping) == 0) return(character(0))
  mapping <- mapping %>% distinct(SYMBOL, .keep_all = TRUE)
  setNames(mapping$ENTREZID, mapping$SYMBOL)
}

# =============================================================
# 5. Read DE genes from CSV (up and down only)
# =============================================================

read_de_genes <- function(filepath,
                          gene_col    = "gene",
                          padj_col    = "padj",
                          lfc_col     = "log2FoldChange",
                          padj_cutoff = ORA_PADJ_CUTOFF,
                          lfc_cutoff  = ORA_LFC_CUTOFF) {
  
  if (!file.exists(filepath)) {
    warning("File not found: ", filepath)
    return(list(up = character(0), down = character(0), full_table = NULL))
  }
  
  df <- read_csv(filepath, show_col_types = FALSE)
  
  needed <- c(gene_col, padj_col, lfc_col)
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    warning("Missing columns: ", paste(missing, collapse = ", "))
    return(list(up = character(0), down = character(0), full_table = df))
  }
  
  sig <- df %>%
    filter(!is.na(.data[[padj_col]]),
           .data[[padj_col]] < padj_cutoff,
           abs(.data[[lfc_col]]) >= lfc_cutoff)
  
  up   <- sig %>% filter(.data[[lfc_col]] > 0) %>% pull(.data[[gene_col]])
  down <- sig %>% filter(.data[[lfc_col]] < 0) %>% pull(.data[[gene_col]])
  
  message("  ", nrow(df), " total -> ", length(up), " up, ", length(down), " down")
  
  list(up = up, down = down, full_table = df)
}

# =============================================================
# 6. Individual ORA functions
# =============================================================

run_go_bp <- function(genes, universe = NULL) {
  if (length(genes) < 3) return(NULL)
  tryCatch({
    enrichGO(gene = genes, universe = universe,
             OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = ORA_PVALUE_CUTOFF, qvalueCutoff = ORA_QVALUE_CUTOFF,
             minGSSize = ORA_MIN_GENE_SIZE, maxGSSize = ORA_MAX_GENE_SIZE,
             readable = TRUE)
  }, error = function(e) { message("    GO BP failed: ", e$message); NULL })
}

run_kegg <- function(entrez_ids, universe_entrez = NULL) {
  if (length(entrez_ids) < 3) return(NULL)
  tryCatch({
    enrichKEGG(gene = entrez_ids, universe = universe_entrez,
               organism = KEGG_ORGANISM, keyType = "kegg",
               pAdjustMethod = "BH",
               pvalueCutoff = ORA_PVALUE_CUTOFF, qvalueCutoff = ORA_QVALUE_CUTOFF,
               minGSSize = ORA_MIN_GENE_SIZE, maxGSSize = ORA_MAX_GENE_SIZE)
  }, error = function(e) { message("    KEGG failed: ", e$message); NULL })
}

run_reactome <- function(entrez_ids, universe_entrez = NULL) {
  if (!HAS_REACTOME || length(entrez_ids) < 3) return(NULL)
  tryCatch({
    enrichPathway(gene = entrez_ids, universe = universe_entrez,
                  organism = "mouse", pAdjustMethod = "BH",
                  pvalueCutoff = ORA_PVALUE_CUTOFF, qvalueCutoff = ORA_QVALUE_CUTOFF,
                  minGSSize = ORA_MIN_GENE_SIZE, maxGSSize = ORA_MAX_GENE_SIZE,
                  readable = TRUE)
  }, error = function(e) { message("    Reactome failed: ", e$message); NULL })
}

# --- SynGO ---

.syngo_cache <- new.env(parent = emptyenv())

load_syngo <- function() {
  if (!file.exists(SYNGO_TERM2GENE)) return(NULL)
  
  t2g <- read_csv(SYNGO_TERM2GENE, show_col_types = FALSE)
  t2n <- if (file.exists(SYNGO_TERM2NAME))
    read_csv(SYNGO_TERM2NAME, show_col_types = FALSE) else NULL
  
  colnames(t2g) <- c("term", "gene")
  if (!is.null(t2n)) colnames(t2n) <- c("term", "name")
  
  t2g <- t2g %>% filter(!is.na(gene), gene != "") %>% distinct()
  universe <- unique(t2g$gene)
  
  message("    SynGO: ", n_distinct(t2g$term), " terms, ",
          length(universe), " genes in universe")
  
  list(term2gene = t2g, term2name = t2n, universe = universe)
}

get_syngo <- function() {
  if (!exists("data", envir = .syngo_cache))
    assign("data", load_syngo(), envir = .syngo_cache)
  get("data", envir = .syngo_cache)
}

#' Run SynGO with mapping report
#' @param mouse_genes Mouse gene symbols
#' @param mouse_universe Mouse background genes
#' @param return_report If TRUE, return mapping stats alongside results
run_syngo <- function(mouse_genes, mouse_universe = NULL,
                      return_report = TRUE) {
  sg <- get_syngo()
  if (is.null(sg)) return(if (return_report) list(result = NULL, report = NULL) else NULL)
  if (length(mouse_genes) < 3) return(if (return_report) list(result = NULL, report = NULL) else NULL)
  
  # Mouse → Human
  ortho_map     <- mouse_to_human(mouse_genes)
  human_genes   <- unname(ortho_map)
  human_in_syngo <- human_genes[human_genes %in% sg$universe]
  
  # Background
  if (!is.null(mouse_universe)) {
    bg_ortho        <- mouse_to_human(mouse_universe)
    human_bg        <- unname(bg_ortho)
    human_bg_syngo  <- intersect(human_bg, sg$universe)
  } else {
    human_bg        <- NULL
    human_bg_syngo  <- sg$universe
  }
  
  # Mapping report
  report <- data.frame(
    mouse_input        = length(mouse_genes),
    human_orthologs    = length(human_genes),
    human_in_syngo     = length(human_in_syngo),
    genes_lost_mapping = length(mouse_genes) - length(unique(human_genes)),
    genes_lost_syngo   = length(unique(human_genes)) - length(human_in_syngo),
    pct_retained       = round(100 * length(human_in_syngo) / max(length(mouse_genes), 1), 1),
    bg_mouse           = if (!is.null(mouse_universe)) length(mouse_universe) else NA_integer_,
    bg_human_mapped    = if (!is.null(human_bg)) length(unique(human_bg)) else NA_integer_,
    bg_in_syngo        = length(human_bg_syngo),
    stringsAsFactors   = FALSE
  )
  
  message("    SynGO mapping: ", report$mouse_input, " mouse → ",
          report$human_orthologs, " human → ",
          report$human_in_syngo, " in SynGO (",
          report$pct_retained, "% retained)")
  
  if (length(human_in_syngo) < 3) {
    message("    Too few genes in SynGO — skipping")
    return(if (return_report) list(result = NULL, report = report) else NULL)
  }
  
  result <- tryCatch({
    enricher(
      gene = human_in_syngo, universe = human_bg_syngo,
      TERM2GENE = sg$term2gene, TERM2NAME = sg$term2name,
      pAdjustMethod = "BH",
      pvalueCutoff = ORA_PVALUE_CUTOFF, qvalueCutoff = ORA_QVALUE_CUTOFF,
      minGSSize = ORA_MIN_GENE_SIZE, maxGSSize = ORA_MAX_GENE_SIZE)
  }, error = function(e) {
    message("    SynGO enrichment failed: ", e$message)
    NULL
  })
  
  if (!is.null(result))
    message("    SynGO: ", nrow(as.data.frame(result)), " terms")
  
  if (return_report) list(result = result, report = report) else result
}

# =============================================================
# 7. Run all databases for one direction
# =============================================================

#' @param genes Mouse gene symbols (one direction)
#' @param universe Mouse background
#' @param run_syngo_flag Whether to run SynGO (FALSE for non-neuronal)
run_all_ora <- function(genes, universe = NULL, run_syngo_flag = TRUE) {
  results <- list()
  syngo_report <- NULL
  
  # GO BP
  message("    GO_BP...")
  res_go <- run_go_bp(genes, universe)
  if (!is.null(res_go) && nrow(as.data.frame(res_go)) > 0)
    results[["GO_BP"]] <- res_go
  
  # KEGG + Reactome (need entrez)
  entrez_map <- symbols_to_entrez(genes)
  entrez_univ <- if (!is.null(universe)) symbols_to_entrez(universe) else NULL
  
  message("    KEGG...")
  res_kegg <- run_kegg(unname(entrez_map),
                       if (!is.null(entrez_univ)) unname(entrez_univ) else NULL)
  if (!is.null(res_kegg) && nrow(as.data.frame(res_kegg)) > 0)
    results[["KEGG"]] <- res_kegg
  
  if (HAS_REACTOME) {
    message("    Reactome...")
    res_react <- run_reactome(unname(entrez_map),
                              if (!is.null(entrez_univ)) unname(entrez_univ) else NULL)
    if (!is.null(res_react) && nrow(as.data.frame(res_react)) > 0)
      results[["Reactome"]] <- res_react
  }
  
  # SynGO (only if flagged)
  if (run_syngo_flag) {
    message("    SynGO...")
    syngo_out <- run_syngo(genes, universe, return_report = TRUE)
    syngo_report <- syngo_out$report
    if (!is.null(syngo_out$result) && nrow(as.data.frame(syngo_out$result)) > 0)
      results[["SynGO"]] <- syngo_out$result
  }
  
  list(results = results, syngo_report = syngo_report)
}

# =============================================================
# 8. Standardize and annotate enrichment results
# =============================================================

#' Convert enrichResult to annotated data.frame with metadata
standardize_ora_result <- function(ora_results, label, comparison,
                                   direction, n_input_genes,
                                   test_cond = NA, ref_cond = NA) {
  if (length(ora_results) == 0) return(NULL)
  
  confidence <- classify_confidence(n_input_genes)
  
  bind_rows(lapply(names(ora_results), function(db) {
    df <- as.data.frame(ora_results[[db]])
    if (nrow(df) == 0) return(NULL)
    
    df %>%
      mutate(
        subtype       = label,
        comparison    = comparison,
        test_cond     = test_cond,
        ref_cond      = ref_cond,
        direction     = direction,
        database      = db,
        n_input_genes = n_input_genes,
        confidence    = confidence,
        term_strength = flag_term_strength(Count),
        .before = 1
      )
  }))
}

# =============================================================
# 9. Plotting: direction-aware two-panel
# =============================================================

plot_ora_direction <- function(up_results, down_results,
                               title = "ORA", top_n = ORA_TOP_N) {
  
  up_df <- if (length(up_results) > 0) {
    bind_rows(lapply(names(up_results), function(db) {
      as.data.frame(up_results[[db]]) %>%
        arrange(p.adjust) %>% head(top_n) %>%
        mutate(database = db, direction = "UP")
    }))
  } else { NULL }
  
  down_df <- if (length(down_results) > 0) {
    bind_rows(lapply(names(down_results), function(db) {
      as.data.frame(down_results[[db]]) %>%
        arrange(p.adjust) %>% head(top_n) %>%
        mutate(database = db, direction = "DOWN")
    }))
  } else { NULL }
  
  combined <- bind_rows(up_df, down_df)
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  
  db_shapes <- c(GO_BP = 16, KEGG = 17, Reactome = 15, SynGO = 18)
  
  combined <- combined %>%
    mutate(Description = ifelse(nchar(Description) > 55,
                                paste0(substr(Description, 1, 52), "..."),
                                Description)) %>%
    arrange(direction, database, p.adjust) %>%
    mutate(Description = factor(Description,
                                levels = rev(unique(Description))))
  
  p <- ggplot(combined, aes(x = Count, y = Description,
                            color = -log10(p.adjust), shape = database)) +
    geom_point(size = 3) +
    facet_wrap(~ direction, scales = "free_y", ncol = 1) +
    scale_color_gradient(low = "grey70", high = "#E64B35",
                         name = "-log10(padj)") +
    scale_shape_manual(values = db_shapes, name = "Database") +
    theme_bw(base_size = 10) +
    ggtitle(title) +
    xlab("Gene count") + ylab(NULL) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold", size = 11),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
  
  p
}

# =============================================================
# 10. Full pipeline: DE CSV → direction-aware ORA
# =============================================================

#' @param filepath DE results CSV
#' @param label Subtype or broad type name
#' @param comparison Comparison name (e.g., "Q1_WT_vs_WT5XFAD")
#' @param out_dir Output directory
#' @param test_cond Test condition label (for summary table)
#' @param ref_cond Reference condition label
#' @param is_neuronal Whether to run SynGO
run_ora_on_de_file <- function(filepath, label, comparison, out_dir,
                               test_cond  = NA, ref_cond = NA,
                               is_neuronal = NULL,
                               padj_cutoff = ORA_PADJ_CUTOFF,
                               lfc_cutoff  = ORA_LFC_CUTOFF) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  safe_label <- gsub("[^A-Za-z0-9_-]", "_", label)
  
  message("\n=== ORA: ", label, " (", comparison, ") ===")
  
  # Auto-detect neuronal if not specified
  if (is.null(is_neuronal)) {
    is_neuronal <- label %in% SYNGO_NEURONAL_TYPES
  }
  
  gene_lists <- read_de_genes(filepath,
                              padj_cutoff = padj_cutoff,
                              lfc_cutoff  = lfc_cutoff)
  
  if (is.null(gene_lists$full_table)) return(NULL)
  
  universe <- gene_lists$full_table$gene
  
  # --- Process each direction ---
  all_summaries     <- list()
  all_syngo_reports <- list()
  direction_results <- list()
  
  for (dir_name in c("up", "down")) {
    genes <- gene_lists[[dir_name]]
    n_genes <- length(genes)
    conf <- classify_confidence(n_genes)
    
    message("\n  ", toupper(dir_name), ": ", n_genes, " genes [", conf, "]")
    
    if (conf == "SKIP") {
      message("    Skipping ORA (< ", MIN_GENES_TO_RUN, " genes)")
      next
    }
    
    # Run all databases
    ora_out <- run_all_ora(genes, universe, run_syngo_flag = is_neuronal)
    ora_res <- ora_out$results
    
    if (!is.null(ora_out$syngo_report)) {
      all_syngo_reports[[dir_name]] <- ora_out$syngo_report %>%
        mutate(direction = dir_name, subtype = label, comparison = comparison)
    }
    
    if (length(ora_res) > 0) {
      direction_results[[dir_name]] <- ora_res
      
      # Save per-database CSVs
      for (db in names(ora_res))
        write_csv(as.data.frame(ora_res[[db]]),
                  file.path(out_dir,
                            paste0(safe_label, "_", dir_name, "_", db, ".csv")))
      
      # Standardized summary
      std <- standardize_ora_result(
        ora_res, label = label, comparison = comparison,
        direction = dir_name, n_input_genes = n_genes,
        test_cond = test_cond, ref_cond = ref_cond
      )
      if (!is.null(std)) all_summaries[[dir_name]] <- std
    }
  }
  
  # --- Combined direction summary CSV ---
  summary_df <- bind_rows(all_summaries)
  if (nrow(summary_df) > 0) {
    write_csv(summary_df,
              file.path(out_dir,
                        paste0(safe_label, "_ORA_ALL_databases.csv")))
  }
  
  # --- SynGO mapping report ---
  syngo_report_df <- bind_rows(all_syngo_reports)
  if (nrow(syngo_report_df) > 0) {
    write_csv(syngo_report_df,
              file.path(out_dir,
                        paste0(safe_label, "_syngo_mapping_report.csv")))
  }
  
  # --- Two-panel direction plot ---
  p <- plot_ora_direction(
    up_results   = direction_results[["up"]],
    down_results = direction_results[["down"]],
    title = paste0(label, " — ", comparison)
  )
  if (!is.null(p)) {
    ggsave(file.path(out_dir,
                     paste0(safe_label, "_ORA_direction_plot.png")),
           plot = p, width = 12, height = 10, dpi = 200)
  }
  
  list(
    direction_results = direction_results,
    summary           = summary_df,
    syngo_report      = syngo_report_df
  )
}

# =============================================================
# 11. ORA on cross-comparison categories
# =============================================================

run_ora_for_categories <- function(categories_file, label, out_dir,
                                   comparison = "cross_comparison",
                                   categories_to_run = c(
                                     "Cat2_Trem2_dependent",
                                     "Cat3_Trem2_independent",
                                     "Cat4_Trem2_exacerbated",
                                     "Cat6_Redirected"),
                                   cat_col = "category",
                                   is_neuronal = NULL) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  safe_label <- gsub("[^A-Za-z0-9_-]", "_", label)
  
  if (!file.exists(categories_file)) {
    warning("Not found: ", categories_file); return(NULL)
  }
  
  if (is.null(is_neuronal))
    is_neuronal <- label %in% SYNGO_NEURONAL_TYPES
  
  df <- read_csv(categories_file, show_col_types = FALSE)
  universe <- unique(df$gene)
  
  # For categories we need LFC to determine direction
  has_lfc <- "q1_lfc" %in% colnames(df)
  
  all_summaries     <- list()
  all_syngo_reports <- list()
  results           <- list()
  
  for (cat_name in categories_to_run) {
    cat_genes <- df %>%
      filter(.data[[cat_col]] == cat_name) %>%
      pull(gene) %>% unique()
    
    if (length(cat_genes) < 3) {
      message("  ", label, " ", cat_name, ": ", length(cat_genes), " genes — skip")
      next
    }
    
    cat_safe <- gsub("[^A-Za-z0-9_-]", "_", cat_name)
    
    # Split by direction if LFC available
    if (has_lfc) {
      cat_df <- df %>% filter(.data[[cat_col]] == cat_name)
      up_genes   <- cat_df %>% filter(q1_lfc > 0) %>% pull(gene) %>% unique()
      down_genes <- cat_df %>% filter(q1_lfc < 0) %>% pull(gene) %>% unique()
    } else {
      # Can't split — run as single set labeled "all"
      up_genes   <- cat_genes
      down_genes <- character(0)
    }
    
    cat_results <- list()
    
    for (dir_info in list(
      list(name = "up",   genes = up_genes),
      list(name = "down", genes = down_genes)
    )) {
      n_g <- length(dir_info$genes)
      conf <- classify_confidence(n_g)
      
      message("\n  ", label, " ", cat_name, " ", toupper(dir_info$name),
              ": ", n_g, " genes [", conf, "]")
      
      if (conf == "SKIP") next
      
      ora_out <- run_all_ora(dir_info$genes, universe,
                             run_syngo_flag = is_neuronal)
      
      if (!is.null(ora_out$syngo_report)) {
        all_syngo_reports[[paste0(cat_name, "_", dir_info$name)]] <-
          ora_out$syngo_report %>%
          mutate(direction = dir_info$name, subtype = label,
                 comparison = comparison, category = cat_name)
      }
      
      if (length(ora_out$results) > 0) {
        cat_results[[dir_info$name]] <- ora_out$results
        
        for (db in names(ora_out$results))
          write_csv(as.data.frame(ora_out$results[[db]]),
                    file.path(out_dir,
                              paste0(safe_label, "_", cat_safe, "_",
                                     dir_info$name, "_", db, ".csv")))
        
        std <- standardize_ora_result(
          ora_out$results, label = label,
          comparison = paste0(comparison, "_", cat_name),
          direction = dir_info$name, n_input_genes = n_g
        )
        if (!is.null(std)) all_summaries[[paste0(cat_name, "_", dir_info$name)]] <- std
      }
    }
    
    results[[cat_name]] <- cat_results
    
    # Direction plot per category
    p <- plot_ora_direction(
      cat_results[["up"]], cat_results[["down"]],
      title = paste0(label, " — ", cat_name)
    )
    if (!is.null(p))
      ggsave(file.path(out_dir,
                       paste0(safe_label, "_", cat_safe, "_direction_plot.png")),
             plot = p, width = 12, height = 10, dpi = 200)
  }
  
  # Combined summary
  summary_df <- bind_rows(all_summaries)
  if (nrow(summary_df) > 0)
    write_csv(summary_df,
              file.path(out_dir,
                        paste0(safe_label, "_categories_ALL_databases.csv")))
  
  # SynGO reports
  syngo_df <- bind_rows(all_syngo_reports)
  if (nrow(syngo_df) > 0)
    write_csv(syngo_df,
              file.path(out_dir,
                        paste0(safe_label, "_categories_syngo_reports.csv")))
  
  results
}

# =============================================================
# 12. ORA on regional categories
# =============================================================

run_ora_for_regional <- function(categories_file, label, out_dir,
                                 comparison = "regional",
                                 categories_to_run = c(
                                   "Reg1_Constitutive",
                                   "Reg2_Amyloid_emergent",
                                   "Reg4_Trem2_dependent_regional",
                                   "Reg5_Trem2_emergent"),
                                 is_neuronal = NULL) {
  
  run_ora_for_categories(
    categories_file   = categories_file,
    label             = label,
    out_dir           = out_dir,
    comparison        = comparison,
    categories_to_run = categories_to_run,
    cat_col           = "reg_category",
    is_neuronal       = is_neuronal
  )
}

# =============================================================
# 13. Startup
# =============================================================

message("11_00_pathway_functions.R loaded successfully")
message("  Databases: GO_BP=YES, KEGG=YES",
        ", Reactome=", ifelse(HAS_REACTOME, "YES", "NO"),
        ", SynGO=", ifelse(file.exists(SYNGO_TERM2GENE), "YES", "NO"))
message("  SynGO neuronal-only: ", length(SYNGO_NEURONAL_TYPES), " types")
message("  Confidence: skip<", MIN_GENES_TO_RUN,
        ", exploratory=", MIN_GENES_TO_RUN, "-", MIN_GENES_STRONG - 1,
        ", interpretable>=", MIN_GENES_STRONG)