# =============================================================================
# Step 05 — Find Marker Genes (Batched Welch t-test)
# =============================================================================
#
# Purpose:
#   Run marker finding for every annotation column in cluster_labels.csv.
#   Each column defines a grouping: cells are assigned the label from that
#   column, and Welch t-tests are run across groups.
#
#   The special column "cluster_id" (always present) uses an identity lookup —
#   each cluster is its own group (numeric cluster-ID pass).
#   User-added columns map cluster_id → annotation string; blank entries are
#   excluded from the DE grouping.
#
#   Per-column caching: a column is skipped if its markers/ directory already
#   exists, unless the column name is listed in force_rerun.
#
# Inputs:
#   res{X.X}/cluster_labels.csv            — cluster_id [| ann_col1 | ...]
#   clustering/clustering.rds              — cluster assignments per cell
#   hvg/hvg.rds                            — genes to test
#   stats/totalcounts.rds                  — per-cell library sizes
#   stats/global_stats.rds                 — count_files list
#
# Outputs (per annotation column, in res{X.X}/{ann_col}/markers/):
#   markers_all.rds      — full marker data.table
#   markers_top.csv      — top n_top_markers per group
#   markers_params.rds   — parameters used
#
# Usage:
#   Rscript 05_find_markers.R          (standalone — processes all columns)
#   source("05_find_markers.R")        (from run_pipeline.R)
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(parallel)
})

if (!exists("out_dir")) source("config.R")
source("utils/marker_utils.R")

# ── Resolve shared inputs ──────────────────────────────────────────────────────

res_folder  <- file.path(out_dir, "annotation", paste0("res", active_resolution))
cluster_col <- paste0("cluster_res", active_resolution)
labels_file <- file.path(res_folder, "cluster_labels.csv")

if (!file.exists(labels_file))
  stop("cluster_labels.csv not found: ", labels_file,
       "\nRun step 04 first.")

labels_dt   <- data.table::fread(labels_file)
ann_cols    <- colnames(labels_dt)   # all columns, including "cluster_id"

message("═══════════════════════════════════════════════════════════════")
message("Step 05 — Find Marker Genes")
message("═══════════════════════════════════════════════════════════════")
message("Annotation columns: ", paste(ann_cols, collapse = ", "))
message("force_rerun       : ", if (length(force_rerun)) paste(force_rerun, collapse = ", ") else "(none)")
message("═══════════════════════════════════════════════════════════════")

# ── Load shared inputs (once) ─────────────────────────────────────────────────

global_stats  <- readRDS(file.path(out_dir, "stats", "global_stats.rds"))
count_files   <- global_stats$count_files
clustering_dt <- readRDS(file.path(out_dir, "clustering", "clustering.rds"))
hvg           <- readRDS(file.path(out_dir, "hvg", "hvg.rds"))
totalcounts   <- readRDS(file.path(out_dir, "stats", "totalcounts.rds"))

if (!cluster_col %in% colnames(clustering_dt))
  stop("Column '", cluster_col, "' not found in clustering.rds")

# ── Loop over annotation columns ──────────────────────────────────────────────

for (ann_col in ann_cols) {

  dir_markers <- file.path(res_folder, ann_col, "markers")

  if (dir.exists(dir_markers) && !ann_col %in% force_rerun) {
    message("\n[SKIP] '", ann_col, "'  (", ann_col, "/markers/ already exists)")
    next
  }

  message("\n── Annotation column: '", ann_col, "' ──────────────────────────")

  # Build per-cell group label vector
  if (ann_col == "cluster_id") {
    # Identity mapping: each cluster is its own group (numeric-cluster pass)
    lookup <- setNames(
      as.character(labels_dt$cluster_id),
      as.character(labels_dt$cluster_id)
    )
  } else {
    lookup_raw <- setNames(
      iconv(as.character(labels_dt[[ann_col]]), from = "UTF-8", to = "UTF-8", sub = ""),
      as.character(labels_dt$cluster_id)
    )
    lookup <- lookup_raw[nchar(trimws(lookup_raw)) > 0L]   # drop blank labels
  }

  n_groups <- length(unique(lookup))
  if (n_groups < 2L) {
    message("  [SKIP] Only ", n_groups, " non-blank group(s) — need at least 2")
    next
  }

  # Map cluster IDs → labels; NA for unlabeled clusters
  clust_dt <- data.table::copy(clustering_dt)
  clust_dt[, .group := lookup[as.character(get(cluster_col))]]

  n_labeled <- sum(!is.na(clust_dt$.group))
  message("  Groups   : ", n_groups)
  message("  Labeled  : ", format(n_labeled, big.mark = ","),
          " / ", format(nrow(clust_dt), big.mark = ","),
          " cells (", round(n_labeled / nrow(clust_dt) * 100, 1), "%)")
  if (ann_col != "cluster_id" && n_groups < length(lookup))
    message("  Merges   : ", length(lookup) - n_groups,
            " cluster(s) share the same label")

  dir.create(dir_markers, showWarnings = FALSE, recursive = TRUE)

  # ── Run batched Welch t-test ───────────────────────────────────────────────

  markers_all <- batched_find_markers(
    batch_files    = count_files,
    genes          = hvg,
    metadata       = clust_dt,
    cluster_column = ".group",
    cell_id_column = "cell_id",
    totalcounts    = totalcounts,
    ncores         = ncores_markers
  )

  markers_all[, priority_score := (cluster_expr + 0.025) /
                (clusterprime_expr + 0.025)]

  markers_all <- markers_all[order(cluster, -avg_log2FC, p_val)]

  if (marker_ranking == "priority") {
    markers_top <- markers_all[
      is.finite(priority_score) & !is.na(priority_score),
      .SD[order(-priority_score)][seq_len(min(.N, n_top_markers))],
      by = cluster
    ]
  } else {
    markers_top <- markers_all[
      is.finite(avg_log2FC) & !is.nan(avg_log2FC) & !is.na(p_val),
      .SD[order(-avg_log2FC)][seq_len(min(.N, n_top_markers))],
      by = cluster
    ]
  }

  # ── Save ──────────────────────────────────────────────────────────────────

  saveRDS(markers_all, file.path(dir_markers, "markers_all.rds"))
  data.table::fwrite(markers_top, file.path(dir_markers, "markers_top.csv"))
  saveRDS(
    list(
      annotation_col = ann_col,
      cluster_col    = cluster_col,
      active_res     = active_resolution,
      n_groups       = n_groups,
      n_top_markers  = n_top_markers,
      marker_ranking = marker_ranking,
      genes_tested   = length(hvg),
      ncores         = ncores_markers
    ),
    file.path(dir_markers, "markers_params.rds")
  )
  message("  Saved: ", ann_col, "/markers/  (",
          nrow(markers_all), " rows, ", n_groups, " groups)")
}

message("\n═══════════════════════════════════════════════════════════════")
message("Step 05 complete")
message("═══════════════════════════════════════════════════════════════")
