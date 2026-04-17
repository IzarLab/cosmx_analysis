# =============================================================================
# Step 03 — Batched PCA
# =============================================================================
#
# Purpose:
#   Compute study-wide PCA embeddings across all raw count matrices without
#   loading all cells into memory simultaneously.
#
#   Two methods available (config: pca_method):
#
#   "pearson" (default) — quasi-Poisson Pearson residual PCA (Lause 2021).
#     Works on raw counts. Uses per-batch gene frequencies for implicit
#     batch correction. Mathematically equivalent to scPearsonPCA.
#
#   "covariance" — standard covariance PCA on library-size-normalised data.
#     Applies normalisation on-the-fly from raw count matrices.
#     Mathematically equivalent to Seurat RunPCA on log1p-normalised data.
#
# Inputs:
#   hvg/hvg.rds              — character vector of HVG gene names
#   stats/global_stats.rds   — contains count_files and mean_totalcounts
#   stats/totalcounts.rds    — per-cell library sizes (Pearson only)
#
# Outputs:
#   pca/embeddings.rds   — cells × npcs matrix (rownames = cell IDs)
#   pca/loadings.rds     — genes × npcs matrix
#   pca/sdev.rds         — numeric vector of length npcs
#   pca/pca_params.rds   — list of parameters used (for reproducibility)
#
# Usage:
#   Rscript 03_run_pca.R       (standalone)
#   source("03_run_pca.R")     (from run_pipeline.R)
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(RSpectra)
})

if (!exists("out_dir")) source("config.R")

stopifnot(pca_method %in% c("pearson", "covariance"))

# ── Load inputs ───────────────────────────────────────────────────────────────

global_stats <- readRDS(file.path(out_dir, "stats", "global_stats.rds"))
count_files  <- global_stats$count_files
mean_totalcounts <- global_stats$mean_totalcounts

hvg          <- readRDS(file.path(out_dir, "hvg", "hvg.rds"))

# ── Output directory ──────────────────────────────────────────────────────────

dir_pca <- file.path(out_dir, "pca")
dir.create(dir_pca, showWarnings = FALSE, recursive = TRUE)

message("═══════════════════════════════════════════════════════════════")
message("Step 03 — Batched PCA")
message("═══════════════════════════════════════════════════════════════")
message("Method      : ", pca_method)
message("Count files : ", length(count_files))
message("HVGs        : ", length(hvg))
message("npcs        : ", npcs)
message("scale_max   : ", scale_max)
message("chunk_size  : ", chunk_size)
if (pca_method == "pearson") message("phi         : ", phi)
if (pca_method == "covariance") message("apply_log1p : ", apply_log1p)
message("═══════════════════════════════════════════════════════════════")

# ── Run PCA ───────────────────────────────────────────────────────────────────

if (pca_method == "pearson") {

  source("utils/pearson_pca_utils.R")

  totalcounts <- readRDS(file.path(out_dir, "stats", "totalcounts.rds"))

  pca_result <- pearson_pca_from_batches(
    batch_files = count_files,
    genes       = hvg,
    totalcounts = totalcounts,
    phi         = phi,
    npcs        = npcs,
    scale_max   = scale_max,
    do.center   = TRUE,
    do.scale    = TRUE,
    chunk_size  = chunk_size
  )

} else {   # "covariance"

  source("utils/cov_pca_utils.R")

  pca_result <- cov_pca_from_raw_batches(
    batch_files      = count_files,
    genes            = hvg,
    mean_totalcounts = mean_totalcounts,
    apply_log1p      = apply_log1p,
    npcs             = npcs,
    scale_max        = scale_max,
    chunk_size       = chunk_size
  )

}

# ── Attach cell IDs as rownames ───────────────────────────────────────────────
# pearson_pca_from_batches / cov_pca_from_raw_batches both preserve
# colnames from each batch file, so rownames(embeddings) == cell IDs.

n_cells <- nrow(pca_result$embeddings)
n_pcs   <- ncol(pca_result$embeddings)

message("\n── Summary ──")
message("Cells      : ", n_cells)
message("PCs        : ", n_pcs)
message("Var. exp. (PC1-5): ",
        paste(round(pca_result$sdev[1:min(5, n_pcs)]^2 /
                    sum(pca_result$sdev^2) * 100, 1),
              collapse = "%, "), "%")

# ── Save outputs ──────────────────────────────────────────────────────────────

saveRDS(pca_result$embeddings, file.path(dir_pca, "embeddings.rds"))
message("Saved: pca/embeddings.rds")

saveRDS(pca_result$loadings,   file.path(dir_pca, "loadings.rds"))
message("Saved: pca/loadings.rds")

saveRDS(pca_result$sdev,       file.path(dir_pca, "sdev.rds"))
message("Saved: pca/sdev.rds")

saveRDS(
  list(
    method     = pca_method,
    npcs       = npcs,
    scale_max  = scale_max,
    chunk_size = chunk_size,
    phi        = if (pca_method == "pearson") phi else NULL,
    apply_log1p = if (pca_method == "covariance") apply_log1p else NULL,
    hvg_n      = length(hvg),
    n_cells    = n_cells
  ),
  file.path(dir_pca, "pca_params.rds")
)
message("Saved: pca/pca_params.rds")

message("\n═══════════════════════════════════════════════════════════════")
message("Step 03 complete")
message("═══════════════════════════════════════════════════════════════")
