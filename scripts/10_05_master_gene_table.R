# ============================================================
# 10_05_master_gene_table.R
#
# Master integration table.
# For every gene that's DE anywhere, builds a wide table showing:
#   - which subtypes call it DE (and direction)
#   - which broad types call it DE (and direction)
#   - whether it's DE globally (and direction)
#   - a "call" column labeling it as global / broad-specific /
#     subtype-specific / mixed
#
# One row per gene per comparison (Q1/Q2/Q3).
#
# Outputs: ./result/cross_comparison/master_gene_table/
# ============================================================

source("./scripts/09_00_de_functions.R")
source("./scripts/10_00_cross_comparison_functions.R")

OUT_DIR <- "./result/cross_comparison/master_gene_table"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

comps <- list(
  Q1 = list(subtype = "Q1_WT_vs_WT5XFAD",
            broad   = "Q1_broad_WT_vs_WT5XFAD",
            global  = "Q1_global_WT_vs_WT5XFAD"),
  Q2 = list(subtype = "Q2_WT5XFAD_vs_Trem2KO5XFAD",
            broad   = "Q2_broad_WT5XFAD_vs_Trem2KO5XFAD",
            global  = "Q2_global_WT5XFAD_vs_Trem2KO5XFAD"),
  Q3 = list(subtype = "Q3_WT_vs_Trem2KO",
            broad   = "Q3_broad_WT_vs_Trem2KO",
            global  = "Q3_global_WT_vs_Trem2KO")
)

# =============================================================
# Helper: read a DE table, return gene → direction (UP/DOWN/blank)
# =============================================================

read_direction <- function(comparison, label) {
  df <- load_de_table(comparison, label)
  if (is.null(df)) return(NULL)
  
  df %>%
    filter(!is.na(significant) & significant) %>%
    mutate(dir = ifelse(log2FoldChange > 0, "up", "down")) %>%
    select(gene, dir)
}

# =============================================================
# Build one master table per comparison
# =============================================================

for (q_name in names(comps)) {
  q <- comps[[q_name]]
  message("\n========== Master table: ", q_name, " ==========")
  
  # Global
  gl <- read_direction(q$global, "Global")
  
  # Broad types
  broad_dirs <- list()
  for (bt in BROAD_TYPES_7m) {
    d <- read_direction(q$broad, bt)
    if (!is.null(d) && nrow(d) > 0)
      broad_dirs[[bt]] <- d %>% rename(!!paste0("broad_", bt) := dir)
  }
  
  # Subtypes
  subtype_dirs <- list()
  for (st in DE_READY_7m_ALL) {
    d <- read_direction(q$subtype, st)
    if (!is.null(d) && nrow(d) > 0) {
      safe <- gsub("[^A-Za-z0-9_-]", "_", st)
      subtype_dirs[[st]] <- d %>% rename(!!paste0("subtype_", safe) := dir)
    }
  }
  
  # Collect all genes that appear anywhere
  all_genes <- unique(c(
    if (!is.null(gl)) gl$gene else character(0),
    unlist(lapply(broad_dirs, function(x) x$gene)),
    unlist(lapply(subtype_dirs, function(x) x$gene))
  ))
  
  if (length(all_genes) == 0) {
    message("  No DE genes anywhere — skipping")
    next
  }
  
  master <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  
  # Join global
  if (!is.null(gl)) {
    master <- master %>%
      left_join(gl %>% rename(global = dir), by = "gene")
  } else {
    master$global <- NA_character_
  }
  
  # Join broad types
  for (bt in names(broad_dirs))
    master <- master %>% left_join(broad_dirs[[bt]], by = "gene")
  
  # Join subtypes
  for (st in names(subtype_dirs))
    master <- master %>% left_join(subtype_dirs[[st]], by = "gene")
  
  # Count hits per scale
  broad_cols   <- grep("^broad_",   colnames(master), value = TRUE)
  subtype_cols <- grep("^subtype_", colnames(master), value = TRUE)
  
  master$n_broad_hits   <- rowSums(!is.na(master[, broad_cols,   drop = FALSE]))
  master$n_subtype_hits <- rowSums(!is.na(master[, subtype_cols, drop = FALSE]))
  
  # Final call
  master <- master %>%
    mutate(call = case_when(
      !is.na(global) & n_broad_hits == 0 & n_subtype_hits == 0  ~ "global_only",
      !is.na(global) & n_broad_hits >= 1                        ~ "global_and_broad",
      is.na(global) & n_broad_hits >= 1 & n_subtype_hits == 0   ~ "broad_specific",
      is.na(global) & n_subtype_hits == 1                       ~ "subtype_specific",
      is.na(global) & n_subtype_hits >= 2                       ~ "multi_subtype",
      is.na(global) & n_broad_hits >= 1 & n_subtype_hits >= 1   ~ "broad_and_subtype",
      TRUE                                                       ~ "other"
    ))
  
  # Reorder columns: gene, global, broad_*, subtype_*, counts, call
  master <- master %>%
    select(gene, global, all_of(broad_cols), all_of(subtype_cols),
           n_broad_hits, n_subtype_hits, call)
  
  # Save
  write_csv(master, file.path(OUT_DIR, paste0(q_name, "_master_table.csv")))
  
  # Summary counts per call type
  call_summary <- master %>%
    count(call, name = "n_genes") %>%
    arrange(desc(n_genes))
  
  write_csv(call_summary,
            file.path(OUT_DIR, paste0(q_name, "_call_counts.csv")))
  
  cat("\n=== ", q_name, " call counts ===\n")
  print(as.data.frame(call_summary), row.names = FALSE)
  cat("Total rows: ", nrow(master), "\n")
}

cat("\n\n==================== MASTER TABLE COMPLETE ====================\n")
cat("Outputs in: ", OUT_DIR, "\n")
cat("  Q{1,2,3}_master_table.csv — one row per gene, columns per scale\n")
cat("  Q{1,2,3}_call_counts.csv  — summary of call types\n")