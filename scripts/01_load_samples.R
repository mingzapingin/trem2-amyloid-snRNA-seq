# ============================================================
# 01_load_samples.R
# Read raw 10x matrices, merge into per-cohort Seurat objects.
#
# Inputs:  ./metadata/GSE140511_sample_metadata.csv
# Outputs: ./processed/seurat_7m_raw.rds
#          ./processed/seurat_15m_raw.rds
# ============================================================

library(Seurat)
library(readr)
library(dplyr)
library(purrr)

dir.create("./processed", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Read metadata
# -----------------------------
meta <- read_csv("./metadata/GSE140511_sample_metadata.csv", show_col_types = FALSE)

required_cols <- c("sample_id", "cohort",
                   "matrix_path", "barcodes_path", "features_path")
missing_cols  <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ",
       paste(missing_cols, collapse = ", "))
}

cat("Cohort breakdown:\n")
print(table(meta$cohort, useNA = "ifany"))

# -----------------------------
# 2. Read one sample into Seurat
# -----------------------------
read_one_sample <- function(row_df) {
  stopifnot(nrow(row_df) == 1)
  sid <- row_df$sample_id[[1]]
  message("Reading sample: ", sid)

  counts <- ReadMtx(
    mtx      = row_df$matrix_path[[1]],
    cells    = row_df$barcodes_path[[1]],
    features = row_df$features_path[[1]],
    feature.column = 2
  )

  obj <- CreateSeuratObject(counts   = counts,
                             project  = sid,
                             min.cells    = 3,
                             min.features = 200)

  # Attach all metadata except file paths
  md_cols <- setdiff(colnames(row_df),
                     c("matrix_path", "barcodes_path", "features_path"))
  for (nm in md_cols) obj[[nm]] <- row_df[[nm]][[1]]
  obj$sample_id <- sid
  obj
}

# -----------------------------
# 3. Read all samples
# -----------------------------
sample_rows <- split(meta, meta$sample_id)
obj_list    <- map(sample_rows, read_one_sample)

# -----------------------------
# 4. Merge by cohort
# -----------------------------
obj_list_7m  <- obj_list[meta$sample_id[meta$cohort == "7m"]]
obj_list_15m <- obj_list[meta$sample_id[meta$cohort == "15m"]]

seurat_7m <- merge(x = obj_list_7m[[1]], y = obj_list_7m[-1],
                   add.cell.ids = names(obj_list_7m),
                   project = "GSE140511_7m")

seurat_15m <- merge(x = obj_list_15m[[1]], y = obj_list_15m[-1],
                    add.cell.ids = names(obj_list_15m),
                    project = "GSE140511_15m")

# -----------------------------
# 5. Quick checks
# -----------------------------
cat("\n7m object:\n"); print(seurat_7m)
print(table(seurat_7m$sample_id, useNA = "ifany"))

cat("\n15m object:\n"); print(seurat_15m)
print(table(seurat_15m$sample_id, useNA = "ifany"))

# -----------------------------
# 6. Save
# -----------------------------
saveRDS(seurat_7m,  "./processed/seurat_7m_raw.rds")
saveRDS(seurat_15m, "./processed/seurat_15m_raw.rds")
