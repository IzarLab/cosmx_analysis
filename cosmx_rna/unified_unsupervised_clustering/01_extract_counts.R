# =============================================================================
# Step 01 — Extract Count Matrices from Seurat Objects
# =============================================================================
#
# Purpose:
#   Single pass over all Seurat objects in seurat_dir. For each object:
#   - Validate gene set consistency across slides
#   - Save the raw RNA count matrix as a sparse .rds file
#   - Accumulate global statistics (mean_totalcounts, per-cell library sizes)
#
# Inputs:
#   seurat_dir/  — directory of Seurat .rds objects (config: seurat_dir)
#
# Outputs:
#   counts/{slide_id}.rds          — sparse genes × cells raw count matrices
#   stats/global_stats.rds         — list(mean_totalcounts, total_cells,
#                                         gene_universe, count_files)
#   stats/totalcounts.rds          — named numeric vector of per-cell colSums
#                                     (from ALL genes; used by Pearson PCA)
#   stats/preprocessing_log.csv    — per-slide summary table
#   stats/cell_slide_lookup.csv    — cell_id, slide_id mapping (used by step 04)
#
# Usage:
#   Rscript 01_extract_counts.R       (standalone)
#   source("01_extract_counts.R")     (from run_pipeline.R)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

if (!exists("out_dir")) source("config.R")

# ── Directories ───────────────────────────────────────────────────────────────

dir_counts <- file.path(out_dir, "counts")
dir_stats  <- file.path(out_dir, "stats")
dir.create(dir_counts, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_stats,  showWarnings = FALSE, recursive = TRUE)

# ── Discover Seurat files ─────────────────────────────────────────────────────

seurat_files <- sort(list.files(seurat_dir, pattern = file_pattern,
                                full.names = TRUE, ignore.case = TRUE))
if (length(seurat_files) == 0L)
  stop("No Seurat files found in: ", seurat_dir)

slide_id_from_path <- function(f) tools::file_path_sans_ext(basename(f))

message("═══════════════════════════════════════════════════════════════")
message("Step 01 — Extract Count Matrices")
message("═══════════════════════════════════════════════════════════════")
message("Input  : ", seurat_dir)
message("Output : ", out_dir)
message("Slides : ", length(seurat_files))
message("═══════════════════════════════════════════════════════════════")

# ── Single pass ───────────────────────────────────────────────────────────────

total_counts_sum <- 0
total_cells      <- 0
gene_universe    <- NULL
totalcounts_list  <- vector("list", length(seurat_files))
count_files       <- character(length(seurat_files))
log_rows          <- vector("list", length(seurat_files))
cell_slide_rows   <- vector("list", length(seurat_files))

for (i in seq_along(seurat_files)) {
  f        <- seurat_files[[i]]
  slide_id <- slide_id_from_path(f)
  message("[", i, "/", length(seurat_files), "] ", basename(f))

  seu    <- readRDS(f)
  counts <- seu[[assay_name]]@counts   # sparse genes × cells, raw counts
  genes  <- rownames(counts)
  n_cells <- ncol(counts)
  tc      <- Matrix::colSums(counts)

  # ── Gene set consistency check ────────────────────────────────────────────
  if (is.null(gene_universe)) {
    gene_universe <- genes
    message("  Gene universe set: ", length(gene_universe), " genes")
  } else {
    if (!identical(genes, gene_universe)) {
      missing <- setdiff(gene_universe, genes)
      extra   <- setdiff(genes, gene_universe)
      if (length(missing) == 0L && length(extra) == 0L)
        stop("Gene mismatch in ", basename(f),
             ": same genes as reference but in a different order.")
      stop(
        "Gene mismatch in ", basename(f), "!\n",
        "  Missing (in reference but not here): ", length(missing), "\n",
        "  Extra   (here but not in reference): ", length(extra), "\n",
        "  First few missing: ", paste(head(missing, 5L), collapse = ", "), "\n",
        "  First few extra:   ", paste(head(extra,   5L), collapse = ", ")
      )
    }
  }

  # ── Accumulate global stats ───────────────────────────────────────────────
  total_counts_sum          <- total_counts_sum + sum(tc)
  total_cells               <- total_cells + n_cells
  totalcounts_list[[i]]     <- tc            # per-cell colSums for this slide

  # ── Save raw count matrix ─────────────────────────────────────────────────
  count_file <- file.path(dir_counts, paste0("counts_", slide_id, ".rds"))
  saveRDS(counts, count_file)
  count_files[[i]] <- count_file

  cell_slide_rows[[i]] <- data.table::data.table(
    cell_id  = colnames(counts),
    slide_id = slide_id
  )

  size_mb <- round(file.size(count_file) / 1e6, 1)
  message("  Cells: ", n_cells,
          " | Genes: ", length(genes),
          " | Mean counts/cell: ", round(mean(tc), 1),
          " | File: ", size_mb, " MB")

  log_rows[[i]] <- data.table(
    slide_id             = slide_id,
    source_file          = basename(f),
    n_cells              = n_cells,
    n_genes              = length(genes),
    mean_counts_per_cell = round(mean(tc), 2),
    count_file           = basename(count_file)
  )

  rm(seu, counts, tc); gc(verbose = FALSE)
}

# ── Finalise ──────────────────────────────────────────────────────────────────

mean_totalcounts <- floor(total_counts_sum / total_cells)
totalcounts      <- unlist(totalcounts_list)   # named by cell ID

message("\n── Summary ──")
message("Total slides  : ", length(seurat_files))
message("Total cells   : ", total_cells)
message("Total genes   : ", length(gene_universe))
message("mean_totalcounts : ", mean_totalcounts)

# ── Save outputs ──────────────────────────────────────────────────────────────

saveRDS(
  list(
    mean_totalcounts = mean_totalcounts,
    total_cells      = total_cells,
    gene_universe    = gene_universe,
    count_files      = count_files
  ),
  file = file.path(dir_stats, "global_stats.rds")
)
message("Saved: stats/global_stats.rds")

saveRDS(totalcounts, file = file.path(dir_stats, "totalcounts.rds"))
message("Saved: stats/totalcounts.rds (", length(totalcounts), " cells)")

log_dt <- data.table::rbindlist(log_rows)
data.table::fwrite(log_dt, file.path(dir_stats, "preprocessing_log.csv"))
message("Saved: stats/preprocessing_log.csv")

cell_slide_dt <- data.table::rbindlist(cell_slide_rows)
data.table::fwrite(
  cell_slide_dt,
  file.path(dir_stats, "cell_slide_lookup.csv")
)
message(
  "Saved: stats/cell_slide_lookup.csv (", nrow(cell_slide_dt), " cells)"
)

message("\n═══════════════════════════════════════════════════════════════")
message("Step 01 complete")
message("═══════════════════════════════════════════════════════════════")
