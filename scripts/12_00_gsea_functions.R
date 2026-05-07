# ============================================================
# 12_00_gsea_functions.R
#
# Multi-database GSEA: Hallmark, GO BP, Reactome, SynGO.
#
# Design rules (matching 11_00 ORA improvements):
#   1. Run GSEA once on full ranked list, then split results
#      by NES direction (positive = higher in test, negative = lower)
#   2. Summary table always has a direction column
#   3. Select top pathways separately from each side
#   4. SynGO: mouse→human, neuronal only, mapping report saved
#   5. Combined CSV with: subtype, comparison, test/ref cond,
#      database, NES, direction, padj, set size
#
# Ranking stat: sign(log2FC) × -log10(pvalue)
#
# Requires:
#   BiocManager::install(c("clusterProfiler", "org.Mm.eg.db",
#                          "org.Hs.eg.db", "ReactomePA", "fgsea"))
#   install.packages("msigdbr")
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

HAS_MSIGDBR <- requireNamespace("msigdbr", quietly = TRUE)
if (HAS_MSIGDBR) suppressPackageStartupMessages(library(msigdbr))

HAS_FGSEA <- requireNamespace("fgsea", quietly = TRUE)
if (HAS_FGSEA) suppressPackageStartupMessages(library(fgsea))

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

GSEA_PVALUE_CUTOFF <- 0.05
GSEA_MIN_SIZE      <- 10
GSEA_MAX_SIZE      <- 500
GSEA_NPERMS        <- 10000

GSEA_TOP_N_PER_SIDE <- 10   # top N from positive AND negative separately

KEGG_ORGANISM <- "mmu"

SYNGO_TERM2GENE <- "./metadata/syngo_TERM2GENE.csv"
SYNGO_TERM2NAME <- "./metadata/syngo_TERM2NAME.csv"

# Neuronal types (SynGO only runs for these — same list as 11_00)
SYNGO_NEURONAL_TYPES <- c(
  "Excitatory_neuron", "Inhibitory_neuron",
  "L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L5 PT CTX",
  "L5/6 NP CTX", "L6 IT CTX", "L6 CT CTX", "L6b CTX",
  "L2/3 IT PPP", "L2/3 IT ENTl", "L6b/CT ENT", "L2 IT ENTl",
  "L3 IT ENT", "Car3", "DG", "CA1-ProS",
  "Pvalb", "Sst", "Meis2", "Lamp5", "Vip", "Sncg"
)

# =============================================================
# 2. Build ranked gene list
# =============================================================

build_rank_vector <- function(filepath,
                              gene_col = "gene",
                              lfc_col  = "log2FoldChange",
                              pval_col = "pvalue",
                              padj_col = "padj") {
  
  if (!file.exists(filepath)) {
    warning("File not found: ", filepath)
    return(list(rank_vector = NULL, full_table = NULL))
  }
  
  df <- read_csv(filepath, show_col_types = FALSE)
  
  # Fall back to padj if pvalue not present
  if (!pval_col %in% colnames(df) && padj_col %in% colnames(df))
    pval_col <- padj_col
  
  needed <- c(gene_col, lfc_col, pval_col)
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    warning("Missing columns: ", paste(missing, collapse = ", "))
    return(list(rank_vector = NULL, full_table = df))
  }
  
  ranked <- df %>%
    filter(!is.na(.data[[lfc_col]]), !is.na(.data[[pval_col]])) %>%
    mutate(pval_safe = ifelse(.data[[pval_col]] == 0,
                              min(.data[[pval_col]][.data[[pval_col]] > 0],
                                  na.rm = TRUE),
                              .data[[pval_col]])) %>%
    mutate(rank_stat = sign(.data[[lfc_col]]) * -log10(pval_safe)) %>%
    distinct(.data[[gene_col]], .keep_all = TRUE) %>%
    arrange(desc(rank_stat))
  
  rank_vec <- setNames(ranked$rank_stat, ranked[[gene_col]])
  
  message("  Ranked ", length(rank_vec), " genes",
          " (range: ", round(min(rank_vec), 2), " to ",
          round(max(rank_vec), 2), ")")
  
  list(rank_vector = rank_vec, full_table = df)
}

# =============================================================
# 3. Gene ID conversions
# =============================================================

symbols_to_entrez <- function(symbols) {
  mapping <- tryCatch({
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Mm.eg.db)
  }, error = function(e) {
    data.frame(SYMBOL = character(0), ENTREZID = character(0))
  })
  if (nrow(mapping) == 0) return(NULL)
  mapping <- mapping %>% distinct(SYMBOL, .keep_all = TRUE)
  setNames(mapping$ENTREZID, mapping$SYMBOL)
}

rank_to_entrez <- function(rank_vec) {
  sym2ent <- symbols_to_entrez(names(rank_vec))
  if (is.null(sym2ent)) return(NULL)
  keep <- names(rank_vec) %in% names(sym2ent)
  entrez_rank <- rank_vec[keep]
  names(entrez_rank) <- unname(sym2ent[names(entrez_rank)])
  entrez_rank <- entrez_rank[!duplicated(names(entrez_rank))]
  sort(entrez_rank, decreasing = TRUE)
}

# --- Mouse → Human (same as 11_00) ---

.ortho_cache <- new.env(parent = emptyenv())

build_ortholog_table <- function() {
  if (exists("table", envir = .ortho_cache))
    return(get("table", envir = .ortho_cache))
  
  if (HAS_ORG_HS) {
    mm_map <- tryCatch({
      bitr(keys(org.Mm.eg.db, keytype = "SYMBOL"),
           fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Mm.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)
    }, error = function(e) NULL)
    
    hs_map <- tryCatch({
      bitr(keys(org.Hs.eg.db, keytype = "SYMBOL"),
           fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)
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
      assign("table", ortho, envir = .ortho_cache)
      return(ortho)
    }
  }
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

rank_to_human <- function(rank_vec) {
  ortho <- mouse_to_human(names(rank_vec))
  human_rank <- rank_vec
  names(human_rank) <- unname(ortho[names(human_rank)])
  human_rank <- human_rank[!duplicated(names(human_rank))]
  sort(human_rank, decreasing = TRUE)
}

# =============================================================
# 4. Individual GSEA functions
# =============================================================

# --- Hallmark ---
run_gsea_hallmark <- function(rank_vec) {
  if (!HAS_MSIGDBR) return(NULL)
  h_sets <- msigdbr(db_species = "MM", category = "H") %>%
    select(gs_name, gene_symbol)
  if (nrow(h_sets) == 0) return(NULL)
  tryCatch({
    GSEA(geneList = rank_vec, TERM2GENE = h_sets,
         pvalueCutoff = GSEA_PVALUE_CUTOFF,
         minGSSize = GSEA_MIN_SIZE, maxGSSize = GSEA_MAX_SIZE,
         pAdjustMethod = "BH", verbose = FALSE)
  }, error = function(e) { message("    Hallmark failed: ", e$message); NULL })
}

# --- GO BP ---
run_gsea_go <- function(rank_vec) {
  tryCatch({
    gseGO(geneList = rank_vec, OrgDb = org.Mm.eg.db,
          keyType = "SYMBOL", ont = "BP",
          pvalueCutoff = GSEA_PVALUE_CUTOFF,
          minGSSize = GSEA_MIN_SIZE, maxGSSize = GSEA_MAX_SIZE,
          pAdjustMethod = "BH", verbose = FALSE)
  }, error = function(e) { message("    GO BP failed: ", e$message); NULL })
}

# --- Reactome ---
run_gsea_reactome <- function(entrez_rank) {
  if (!HAS_REACTOME || is.null(entrez_rank)) return(NULL)
  tryCatch({
    gsePathway(geneList = entrez_rank, organism = "mouse",
               pvalueCutoff = GSEA_PVALUE_CUTOFF,
               minGSSize = GSEA_MIN_SIZE, maxGSSize = GSEA_MAX_SIZE,
               pAdjustMethod = "BH", verbose = FALSE)
  }, error = function(e) { message("    Reactome failed: ", e$message); NULL })
}

# --- SynGO (fgsea, human) ---

.syngo_cache <- new.env(parent = emptyenv())

load_syngo_genesets <- function() {
  if (exists("data", envir = .syngo_cache))
    return(get("data", envir = .syngo_cache))
  
  if (!file.exists(SYNGO_TERM2GENE)) {
    assign("data", NULL, envir = .syngo_cache)
    return(NULL)
  }
  
  t2g <- read_csv(SYNGO_TERM2GENE, show_col_types = FALSE)
  colnames(t2g) <- c("term", "gene")
  t2g <- t2g %>% filter(!is.na(gene), gene != "") %>% distinct()
  
  t2n <- if (file.exists(SYNGO_TERM2NAME)) {
    tmp <- read_csv(SYNGO_TERM2NAME, show_col_types = FALSE)
    colnames(tmp) <- c("term", "name")
    tmp %>% distinct()
  } else { NULL }
  
  pathways <- split(t2g$gene, t2g$term)
  if (!is.null(t2n)) {
    term_names <- setNames(t2n$name, t2n$term)
    names(pathways) <- ifelse(names(pathways) %in% names(term_names),
                              term_names[names(pathways)], names(pathways))
  }
  
  message("    SynGO: ", length(pathways), " gene sets loaded")
  assign("data", pathways, envir = .syngo_cache)
  pathways
}

#' Run SynGO GSEA with mapping report
run_gsea_syngo <- function(rank_vec) {
  if (!HAS_FGSEA) return(list(result = NULL, report = NULL))
  
  pathways <- load_syngo_genesets()
  if (is.null(pathways)) return(list(result = NULL, report = NULL))
  
  message("    Converting mouse → human for SynGO...")
  human_rank <- rank_to_human(rank_vec)
  
  # Filter pathways
  pathways_filt <- lapply(pathways, function(gs) intersect(gs, names(human_rank)))
  pathways_filt <- pathways_filt[lengths(pathways_filt) >= GSEA_MIN_SIZE]
  
  # Mapping report
  syngo_universe <- unique(unlist(pathways))
  human_in_syngo <- intersect(names(human_rank), syngo_universe)
  
  report <- data.frame(
    mouse_genes_ranked  = length(rank_vec),
    human_orthologs     = length(human_rank),
    human_in_syngo      = length(human_in_syngo),
    pct_in_syngo        = round(100 * length(human_in_syngo) / max(length(human_rank), 1), 1),
    syngo_sets_testable = length(pathways_filt),
    syngo_sets_total    = length(pathways),
    stringsAsFactors    = FALSE
  )
  
  message("    SynGO mapping: ", report$mouse_genes_ranked, " mouse → ",
          report$human_orthologs, " human → ",
          report$human_in_syngo, " in SynGO universe (",
          report$pct_in_syngo, "%), ",
          report$syngo_sets_testable, " testable sets")
  
  if (length(pathways_filt) == 0) {
    message("    No testable SynGO sets — skipping")
    return(list(result = NULL, report = report))
  }
  
  result <- tryCatch({
    fgsea(pathways = pathways_filt, stats = human_rank,
          minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE,
          nPermSimple = GSEA_NPERMS)
  }, error = function(e) {
    message("    fgsea failed: ", e$message); NULL
  })
  
  if (!is.null(result) && nrow(result) > 0) {
    result <- result %>% filter(padj < GSEA_PVALUE_CUTOFF) %>% arrange(padj)
    message("    SynGO: ", nrow(result), " significant sets")
  }
  
  list(result = result, report = report)
}

# =============================================================
# 5. Run all GSEA databases
# =============================================================

run_all_gsea <- function(rank_vec, run_syngo_flag = TRUE) {
  results <- list()
  syngo_report <- NULL
  
  message("    Hallmark...")
  res_h <- run_gsea_hallmark(rank_vec)
  if (!is.null(res_h) && nrow(as.data.frame(res_h)) > 0) {
    results[["Hallmark"]] <- res_h
    message("      ", nrow(as.data.frame(res_h)), " terms")
  }
  
  message("    GO_BP...")
  res_go <- run_gsea_go(rank_vec)
  if (!is.null(res_go) && nrow(as.data.frame(res_go)) > 0) {
    results[["GO_BP"]] <- res_go
    message("      ", nrow(as.data.frame(res_go)), " terms")
  }
  
  message("    Reactome...")
  entrez_rank <- rank_to_entrez(rank_vec)
  res_react <- run_gsea_reactome(entrez_rank)
  if (!is.null(res_react) && nrow(as.data.frame(res_react)) > 0) {
    results[["Reactome"]] <- res_react
    message("      ", nrow(as.data.frame(res_react)), " terms")
  }
  
  if (run_syngo_flag) {
    message("    SynGO...")
    syngo_out <- run_gsea_syngo(rank_vec)
    syngo_report <- syngo_out$report
    if (!is.null(syngo_out$result) && nrow(syngo_out$result) > 0)
      results[["SynGO"]] <- syngo_out$result
  }
  
  list(results = results, syngo_report = syngo_report)
}

# =============================================================
# 6. Standardize GSEA results into one table
# =============================================================

#' Normalize clusterProfiler + fgsea outputs into a common format,
#' split by NES direction, add metadata columns.
standardize_gsea_results <- function(gsea_results, label, comparison,
                                     test_cond = NA, ref_cond = NA) {
  if (length(gsea_results) == 0) return(NULL)
  
  bind_rows(lapply(names(gsea_results), function(db) {
    res <- gsea_results[[db]]
    
    if (db == "SynGO") {
      # fgsea data.table format
      df <- as.data.frame(res) %>%
        transmute(
          database    = db,
          ID          = pathway,
          Description = pathway,
          NES         = NES,
          pvalue      = pval,
          padj        = padj,
          setSize     = size,
          leading_edge = sapply(leadingEdge, function(x) paste(x, collapse = "/"))
        )
    } else {
      # clusterProfiler enrichResult format
      df <- as.data.frame(res) %>%
        transmute(
          database    = db,
          ID          = ID,
          Description = Description,
          NES         = NES,
          pvalue      = pvalue,
          padj        = p.adjust,
          setSize     = setSize,
          leading_edge = core_enrichment
        )
    }
    
    df %>%
      mutate(
        subtype    = label,
        comparison = comparison,
        test_cond  = test_cond,
        ref_cond   = ref_cond,
        direction  = ifelse(NES > 0, "Activated", "Suppressed"),
        .before    = 1
      )
  }))
}

# =============================================================
# 7. Plotting: direction-aware NES
# =============================================================

#' NES dot plot: top N from each side, faceted by database
plot_gsea_nes_dot <- function(combined_df, title = "GSEA",
                              top_n = GSEA_TOP_N_PER_SIDE) {
  if (is.null(combined_df) || nrow(combined_df) == 0) return(NULL)
  
  # Select top N per side per database
  top_pos <- combined_df %>% filter(NES > 0) %>%
    group_by(database) %>% arrange(padj) %>% slice_head(n = top_n) %>% ungroup()
  top_neg <- combined_df %>% filter(NES < 0) %>%
    group_by(database) %>% arrange(padj) %>% slice_head(n = top_n) %>% ungroup()
  
  plot_df <- bind_rows(top_pos, top_neg)
  if (nrow(plot_df) == 0) return(NULL)
  
  plot_df <- plot_df %>%
    mutate(Description = ifelse(nchar(Description) > 55,
                                paste0(substr(Description, 1, 52), "..."),
                                Description)) %>%
    arrange(database, NES) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  ggplot(plot_df, aes(x = NES, y = Description,
                      size = setSize, color = -log10(padj))) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ database, scales = "free_y", ncol = 1) +
    scale_color_gradient(low = "grey70", high = "#E64B35",
                         name = "-log10(padj)") +
    scale_size_continuous(range = c(2, 7), name = "Set size") +
    theme_bw(base_size = 10) +
    ggtitle(title) +
    xlab("Normalized Enrichment Score (NES)") + ylab(NULL) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold", size = 11),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
}

#' NES bar plot colored by direction
plot_gsea_nes_bar <- function(combined_df, title = "GSEA",
                              top_n = GSEA_TOP_N_PER_SIDE) {
  if (is.null(combined_df) || nrow(combined_df) == 0) return(NULL)
  
  top_pos <- combined_df %>% filter(NES > 0) %>%
    group_by(database) %>% arrange(padj) %>% slice_head(n = top_n) %>% ungroup()
  top_neg <- combined_df %>% filter(NES < 0) %>%
    group_by(database) %>% arrange(padj) %>% slice_head(n = top_n) %>% ungroup()
  
  plot_df <- bind_rows(top_pos, top_neg)
  if (nrow(plot_df) == 0) return(NULL)
  
  plot_df <- plot_df %>%
    mutate(Description = ifelse(nchar(Description) > 55,
                                paste0(substr(Description, 1, 52), "..."),
                                Description)) %>%
    arrange(database, NES) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  ggplot(plot_df, aes(x = NES, y = Description, fill = direction)) +
    geom_col(alpha = 0.85) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    facet_wrap(~ database, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Activated" = "#E64B35",
                                 "Suppressed" = "#4DBBD5")) +
    theme_bw(base_size = 10) +
    ggtitle(title) +
    xlab("Normalized Enrichment Score (NES)") + ylab(NULL) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold", size = 11),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
}

# =============================================================
# 8. Full pipeline: DE CSV → direction-aware GSEA
# =============================================================

run_gsea_on_de_file <- function(filepath, label, comparison, out_dir,
                                test_cond   = NA, ref_cond = NA,
                                is_neuronal = NULL) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  safe_label <- gsub("[^A-Za-z0-9_-]", "_", label)
  
  message("\n=== GSEA: ", label, " (", comparison, ") ===")
  
  if (is.null(is_neuronal))
    is_neuronal <- label %in% SYNGO_NEURONAL_TYPES
  
  # Build ranked list
  ranked <- build_rank_vector(filepath)
  if (is.null(ranked$rank_vector) || length(ranked$rank_vector) < 100) {
    message("  Too few ranked genes — skipping")
    return(NULL)
  }
  
  # Run all databases
  gsea_out <- run_all_gsea(ranked$rank_vector, run_syngo_flag = is_neuronal)
  gsea_res <- gsea_out$results
  
  if (length(gsea_res) == 0) {
    message("  No enriched gene sets")
    return(NULL)
  }
  
  # --- Save per-database CSVs ---
  for (db in names(gsea_res))
    write_csv(as.data.frame(gsea_res[[db]]),
              file.path(out_dir,
                        paste0(safe_label, "_GSEA_", db, ".csv")))
  
  # --- Standardize into combined table ---
  combined_df <- standardize_gsea_results(
    gsea_res, label = label, comparison = comparison,
    test_cond = test_cond, ref_cond = ref_cond
  )
  
  if (!is.null(combined_df) && nrow(combined_df) > 0) {
    # Save combined (all databases, both directions)
    write_csv(combined_df,
              file.path(out_dir,
                        paste0(safe_label, "_GSEA_ALL_databases.csv")))
    
    # Save split by direction
    activated  <- combined_df %>% filter(direction == "Activated")
    suppressed <- combined_df %>% filter(direction == "Suppressed")
    
    if (nrow(activated) > 0)
      write_csv(activated,
                file.path(out_dir,
                          paste0(safe_label, "_GSEA_activated.csv")))
    if (nrow(suppressed) > 0)
      write_csv(suppressed,
                file.path(out_dir,
                          paste0(safe_label, "_GSEA_suppressed.csv")))
    
    message("  Results: ", nrow(activated), " activated, ",
            nrow(suppressed), " suppressed pathways")
  }
  
  # --- SynGO mapping report ---
  if (!is.null(gsea_out$syngo_report)) {
    syngo_report <- gsea_out$syngo_report %>%
      mutate(subtype = label, comparison = comparison)
    write_csv(syngo_report,
              file.path(out_dir,
                        paste0(safe_label, "_syngo_mapping_report.csv")))
  }
  
  # --- Plots (balanced top N from each side) ---
  p_dot <- plot_gsea_nes_dot(combined_df,
                             title = paste0(label, " — ", comparison))
  if (!is.null(p_dot)) {
    n_db <- length(gsea_res)
    ggsave(file.path(out_dir,
                     paste0(safe_label, "_GSEA_nes_dot.png")),
           plot = p_dot, width = 12, height = 3 + n_db * 4, dpi = 200)
  }
  
  p_bar <- plot_gsea_nes_bar(combined_df,
                             title = paste0(label, " — ", comparison))
  if (!is.null(p_bar)) {
    n_db <- length(gsea_res)
    ggsave(file.path(out_dir,
                     paste0(safe_label, "_GSEA_nes_bar.png")),
           plot = p_bar, width = 12, height = 3 + n_db * 4, dpi = 200)
  }
  
  list(
    gsea_results = gsea_res,
    combined     = combined_df,
    syngo_report = gsea_out$syngo_report
  )
}

# =============================================================
# 9. Startup
# =============================================================

message("12_00_gsea_functions.R loaded successfully")
message("  Databases: Hallmark=", ifelse(HAS_MSIGDBR, "YES", "NO"),
        ", GO_BP=YES",
        ", Reactome=", ifelse(HAS_REACTOME, "YES", "NO"),
        ", SynGO=", ifelse(file.exists(SYNGO_TERM2GENE) & HAS_FGSEA, "YES", "NO"))
message("  Top N per side: ", GSEA_TOP_N_PER_SIDE,
        " (activated + suppressed separately)")