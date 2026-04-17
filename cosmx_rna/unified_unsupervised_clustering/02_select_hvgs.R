# =============================================================================
# Step 02 — Select Highly Variable Genes (Batched VST)
# =============================================================================
#
# Purpose:
#   One pass over the raw count matrices to accumulate per-gene sum and
#   sum-of-squares, then select study-wide HVGs via the Seurat VST approach:
#   fit LOESS on log10(variance) ~ log10(mean) and rank genes by
#   observed_variance / expected_variance.
#
# This step is modular and independent from Step 01 (count extraction).
# It re-reads the count matrices — a deliberate design choice that keeps each
# step self-contained and independently re-runnable.
#
# Inputs:
#   stats/global_stats.rds  — produced by 01_extract_counts.R
#   counts/*.rds            — raw count matrices
#
# Outputs:
#   hvg/hvg.rds       — character vector of selected HVG gene names
#   hvg/hvg_stats.csv — per-gene VST statistics table
#
# Usage:
#   Rscript 02_select_hvgs.R       (standalone)
#   source("02_select_hvgs.R")     (from run_pipeline.R)
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

if (!exists("out_dir")) source("config.R")
source("utils/hvg_utils.R")

# ── Load Step 01 outputs ──────────────────────────────────────────────────────

global_stats <- readRDS(file.path(out_dir, "stats", "global_stats.rds"))
count_files  <- global_stats$count_files
gene_universe <- global_stats$gene_universe

if (length(count_files) == 0L)
  stop("No count files found in global_stats.rds. Run 01_extract_counts.R first.")

# ── Output directory ──────────────────────────────────────────────────────────

dir_hvg <- file.path(out_dir, "hvg")
dir.create(dir_hvg, showWarnings = FALSE, recursive = TRUE)

message("═══════════════════════════════════════════════════════════════")
message("Step 02 — HVG Selection (Batched VST)")
message("═══════════════════════════════════════════════════════════════")
message("Count files : ", length(count_files))
message("Genes       : ", length(gene_universe))
message("n_hvg       : ", hvg_n)
message("loess_span  : ", hvg_loess_span)
message("═══════════════════════════════════════════════════════════════")

# ── Run batched HVG selection ─────────────────────────────────────────────────

hvg_stats <- batched_find_hvg(
  batch_files = count_files,
  genes       = gene_universe,
  n_hvg       = hvg_n,
  loess_span  = hvg_loess_span
)

hvg <- hvg_stats$gene[hvg_stats$highly.variable]

message("\n── Summary ──")
message("HVGs selected    : ", length(hvg))
message("Mean expression  : ", round(mean(hvg_stats$mean[hvg_stats$highly.variable]), 3))
message("Median std. var  : ", round(median(hvg_stats$variance.standardized[hvg_stats$highly.variable],
                                              na.rm = TRUE), 3))

# ── Save outputs ──────────────────────────────────────────────────────────────

saveRDS(hvg, file.path(dir_hvg, "hvg.rds"))
message("Saved: hvg/hvg.rds (", length(hvg), " genes)")

data.table::fwrite(
  data.table::as.data.table(hvg_stats),
  file.path(dir_hvg, "hvg_stats.csv")
)
message("Saved: hvg/hvg_stats.csv")

message("\n═══════════════════════════════════════════════════════════════")
message("Step 02 complete")
message("═══════════════════════════════════════════════════════════════")
