# =============================================================================
# Build SynGO Human Reference CSV
# =============================================================================
# Merges annotations.xlsx, genes.xlsx, and ontologies.xlsx into a single flat
# CSV file (syngo_human_reference.csv) ready for ORA with clusterProfiler.
#
# Input files (place in the same directory as this script, or adjust paths):
#   - annotations.xlsx
#   - genes.xlsx
#   - ontologies.xlsx
#
# Output:
#   - syngo_human_reference.csv   — one row per gene × GO-term annotation
#
# Install dependencies if needed:
#   install.packages(c("readxl", "dplyr", "tidyr", "writexl"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)

# =============================================================================
# 1. Load the three SynGO files
# =============================================================================

message("Loading SynGO files ...")

annotations <- read_excel("./metadata/synGO/annotations.xlsx")
genes       <- read_excel("./metadata/synGO/genes.xlsx")
ontologies  <- read_excel("./metadata/synGO/ontologies.xlsx")

message(sprintf("  annotations : %d rows x %d cols", nrow(annotations), ncol(annotations)))
message(sprintf("  genes       : %d rows x %d cols", nrow(genes),       ncol(genes)))
message(sprintf("  ontologies  : %d rows x %d cols", nrow(ontologies),  ncol(ontologies)))

# =============================================================================
# 2. Clean and prepare each table
# =============================================================================

# --- annotations: keep one row per unique gene × GO-term pair ---------------
# (multiple rows exist for the same pair when backed by different publications
#  or evidence codes; we collapse those into a single annotation row)

ann_clean <- annotations %>%
  select(
    hgnc_id, hgnc_symbol, uniprot_id,
    go_id, go_name, go_domain
  ) %>%
  distinct()

# --- genes: gene-level metadata (Ensembl, Entrez, full name, synonyms) ------

genes_clean <- genes %>%
  select(
    hgnc_id, hgnc_symbol,
    hgnc_name, hgnc_synonyms,
    ensembl_id, entrez_id
  )

# --- ontologies: term-level metadata (hierarchy, short name, parent) --------

onto_clean <- ontologies %>%
  select(
    go_id        = id,
    go_domain    = domain,
    go_name_full = name,
    go_name_hier = name_hierarchical,
    go_shortname = shortname,
    go_parent_id = parent_id
  ) %>%
  distinct(go_id, .keep_all = TRUE)   # ontologies has one row per term

# =============================================================================
# 3. Join all three tables
# =============================================================================

message("Joining tables ...")

syngo_full <- ann_clean %>%
  # Add gene metadata (Ensembl, Entrez, full name, synonyms)
  left_join(
    genes_clean %>% select(-hgnc_symbol),  # hgnc_symbol already in ann_clean
    by = "hgnc_id"
  ) %>%
  # Add ontology metadata (hierarchy, short name, parent term)
  left_join(
    onto_clean %>% select(-go_domain, -go_name_full),  # already present
    by = "go_id"
  ) %>%
  # Reorder columns logically
  select(
    # Gene identifiers
    hgnc_id, hgnc_symbol, hgnc_name, hgnc_synonyms,
    ensembl_id, entrez_id, uniprot_id,
    # GO term identifiers
    go_id, go_domain, go_name,
    # Ontology hierarchy extras
    go_shortname, go_parent_id, go_name_hier
  ) %>%
  arrange(go_domain, go_id, hgnc_symbol)

message(sprintf(
  "Final table: %d rows | %d unique genes | %d unique GO terms",
  nrow(syngo_full),
  n_distinct(syngo_full$hgnc_symbol),
  n_distinct(syngo_full$go_id)
))

# =============================================================================
# 4. Sanity checks
# =============================================================================

# Genes in annotations not found in genes table (should be 0 or very few)
missing_genes <- setdiff(ann_clean$hgnc_id, genes_clean$hgnc_id)
if (length(missing_genes) > 0) {
  message(sprintf(
    "Warning: %d hgnc_ids in annotations have no match in genes table.",
    length(missing_genes)
  ))
}

# GO terms in annotations not found in ontologies table
missing_terms <- setdiff(ann_clean$go_id, onto_clean$go_id)
if (length(missing_terms) > 0) {
  message(sprintf(
    "Warning: %d go_ids in annotations have no match in ontologies table: %s",
    length(missing_terms), paste(missing_terms, collapse = ", ")
  ))
}

# =============================================================================
# 5. Save the merged reference CSV
# =============================================================================

output_file <- "./metadata/synGO/syngo_human_reference.csv"
write.csv(syngo_full, output_file, row.names = FALSE)
message(sprintf("Saved: %s", output_file))

# =============================================================================
# 6. Also save a minimal TERM2GENE table (for clusterProfiler::enricher)
# =============================================================================
# enricher() requires exactly two columns: term ID and gene symbol.
# Two flavours are saved — one keyed on GO ID, one on human gene symbol.

term2gene <- syngo_full %>%
  select(go_id, hgnc_symbol) %>%
  distinct()

term2name <- syngo_full %>%
  mutate(go_label = paste0(go_name, " [", go_domain, "]")) %>%
  select(go_id, go_label) %>%
  distinct()

write.csv(term2gene, "./metadata/synGO/syngo_TERM2GENE.csv", row.names = FALSE)
write.csv(term2name, "./metadata/synGO/syngo_TERM2NAME.csv", row.names = FALSE)
message("Saved: syngo_TERM2GENE.csv")
message("Saved: syngo_TERM2NAME.csv")

# =============================================================================
# 7. Quick summary by domain
# =============================================================================

summary_tbl <- syngo_full %>%
  group_by(go_domain) %>%
  summarise(
    n_terms  = n_distinct(go_id),
    n_genes  = n_distinct(hgnc_symbol),
    n_rows   = n(),
    .groups  = "drop"
  )

message("\n--- Summary by domain ---")
print(as.data.frame(summary_tbl))

message("\nDone. Files ready for ORA:")
message("  syngo_human_reference.csv  — full merged reference")
message("  syngo_TERM2GENE.csv        — TERM2GENE for clusterProfiler::enricher()")
message("  syngo_TERM2NAME.csv        — TERM2NAME for clusterProfiler::enricher()")