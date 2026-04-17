# =============================================================================
# Step 04 — Clustering and UMAP
# =============================================================================
#
# Purpose:
#   Build a kNN similarity graph from PCA embeddings using uwot, compute
#   a UMAP layout for visualization, and run unsupervised clustering via
#   Seurat::FindClusters at one or more resolutions.
#
#   The similarity graph is computed with uwot::similarity_graph() (geometric
#   kNN, cosine distance by default) and converted to a Seurat Graph object
#   for clustering. This avoids manual igraph conversion and allows Seurat to
#   handle multiple resolutions in a single call.
#
# Inputs:
#   pca/embeddings.rds              — cells × npcs PCA embeddings
#   stats/preprocessing_log.csv     — slide_id per cell (for cell metadata)
#
# Outputs:
#   clustering/clustering.rds       — data.table: cell_id, slide_id,
#                                     umap_1, umap_2, cluster_res{r} ...
#   clustering/umap.rds             — matrix cells × 2 (same as umap_1/2)
#   clustering/clustering_params.rds — parameters used
#
# Usage:
#   Rscript 04_cluster.R       (standalone)
#   source("04_cluster.R")     (from run_pipeline.R)
# =============================================================================

suppressPackageStartupMessages({
  library(uwot)
  library(Seurat)
  library(data.table)
})

if (!exists("out_dir")) source("config.R")

stopifnot(cluster_algorithm %in% c("louvain", "leiden"))

# Seurat algorithm codes: 1 = Louvain, 4 = Leiden
cluster_algorithm_id <- if (cluster_algorithm == "louvain") 1L else 4L

# ── Load inputs ───────────────────────────────────────────────────────────────

embeddings      <- readRDS(file.path(out_dir, "pca", "embeddings.rds"))
cell_slide_csv  <- file.path(out_dir, "stats", "cell_slide_lookup.csv")

if (file.exists(cell_slide_csv)) {
  slide_lookup <- data.table::fread(cell_slide_csv)
} else {
  # Fallback for outputs produced before cell_slide_lookup.csv was added.
  # Build the lookup from the count matrices and save it for future runs.
  message("cell_slide_lookup.csv not found — building from count matrices...")
  global_stats <- readRDS(file.path(out_dir, "stats", "global_stats.rds"))
  preproc_log  <- data.table::fread(
    file.path(out_dir, "stats", "preprocessing_log.csv")
  )
  slide_rows <- vector("list", nrow(preproc_log))
  for (i in seq_len(nrow(preproc_log))) {
    cells_i <- colnames(readRDS(global_stats$count_files[[i]]))
    slide_rows[[i]] <- data.table::data.table(
      cell_id  = cells_i,
      slide_id = preproc_log$slide_id[[i]]
    )
  }
  slide_lookup <- data.table::rbindlist(slide_rows)
  data.table::fwrite(slide_lookup, cell_slide_csv)
  message("Saved: stats/cell_slide_lookup.csv (", nrow(slide_lookup), " cells)")
}

n_cells <- nrow(embeddings)

message("═══════════════════════════════════════════════════════════════")
message("Step 04 — Clustering and UMAP")
message("═══════════════════════════════════════════════════════════════")
message("Cells       : ", n_cells)
message("n_neighbors : ", n_neighbors)
message("metric      : ", metric)
message("algorithm   : ", cluster_algorithm)
message("resolutions : ", paste(resolutions, collapse = ", "))
message("═══════════════════════════════════════════════════════════════")

# ── Output directory ──────────────────────────────────────────────────────────

dir_clust <- file.path(out_dir, "clustering")
dir.create(dir_clust, showWarnings = FALSE, recursive = TRUE)

# ── Build kNN similarity graph ────────────────────────────────────────────────

message("\nBuilding kNN similarity graph...")
simgrph <- uwot::similarity_graph(
  X           = embeddings,
  n_neighbors = n_neighbors,
  metric      = metric,
  nn_method   = nn_method
)
# uwot does not set dimnames; Seurat::as.Graph() requires rownames
rownames(simgrph) <- rownames(embeddings)
colnames(simgrph) <- rownames(embeddings)

# ── Compute UMAP ─────────────────────────────────────────────────────────────

message("Computing UMAP...")
set.seed(umap_seed)
umap_coords <- uwot::umap(
  X           = embeddings,
  n_neighbors = n_neighbors,
  metric      = metric,
  min_dist    = min_dist,
  nn_method   = nn_method,
  verbose     = FALSE
)
rownames(umap_coords) <- rownames(embeddings)
colnames(umap_coords) <- c("umap_1", "umap_2")

# ── Clustering ────────────────────────────────────────────────────────────────

message("Running ", cluster_algorithm, " clustering at ",
        length(resolutions), " resolution(s)...")

seurat_graph   <- Seurat::as.Graph(simgrph)
cluster_result <- Seurat::FindClusters(
  object    = seurat_graph,
  resolution = resolutions,
  algorithm  = cluster_algorithm_id,
  verbose    = FALSE
)
# cluster_result is a data.frame; rownames = cell IDs; columns = "res.X.X"

# ── Build output data.table ───────────────────────────────────────────────────

cell_ids <- rownames(embeddings)

# slide_lookup loaded above from stats/cell_slide_lookup.csv

# Cluster columns: rename from Seurat's "res.X.X" → "cluster_res{r}"
# Pattern substitution on whatever res.* columns Seurat produced — robust to
# numeric formatting differences (e.g. "res.1" vs "res.1.0").
clust_dt  <- data.table::as.data.table(
  cluster_result, keep.rownames = "cell_id"
)
res_found <- grep("^res\\.", colnames(clust_dt), value = TRUE)
data.table::setnames(
  clust_dt, res_found, sub("^res\\.", "cluster_res", res_found)
)

# Warn (not stop) if a resolution is missing — clustering was expensive
expected_cols <- paste0("cluster_res", resolutions)
missing_cols  <- setdiff(expected_cols, colnames(clust_dt))
if (length(missing_cols) > 0)
  warning(
    "FindClusters did not produce expected columns: ",
    paste(missing_cols, collapse = ", "),
    "\nFound: ", paste(colnames(clust_dt), collapse = ", ")
  )

umap_dt <- data.table::as.data.table(umap_coords, keep.rownames = "cell_id")

clustering_dt <- Reduce(
  function(a, b) merge(a, b, by = "cell_id", all = FALSE),
  list(
    data.table(cell_id = cell_ids),
    slide_lookup,
    umap_dt,
    clust_dt
  )
)

# Report cluster sizes per resolution
for (res in resolutions) {
  col <- paste0("cluster_res", res)
  if (col %in% colnames(clustering_dt)) {
    n_clust <- length(unique(clustering_dt[[col]]))
    message("  res=", res, " → ", n_clust, " clusters")
  }
}

# ── Save outputs ──────────────────────────────────────────────────────────────

saveRDS(clustering_dt,  file.path(dir_clust, "clustering.rds"))
message("\nSaved: clustering/clustering.rds (", nrow(clustering_dt), " cells)")

saveRDS(umap_coords, file.path(dir_clust, "umap.rds"))
message("Saved: clustering/umap.rds")

saveRDS(
  list(
    n_neighbors       = n_neighbors,
    metric            = metric,
    min_dist          = min_dist,
    nn_method         = nn_method,
    cluster_algorithm = cluster_algorithm,
    resolutions       = resolutions
  ),
  file.path(dir_clust, "clustering_params.rds")
)
message("Saved: clustering/clustering_params.rds")

# ── Initialise res{active_resolution}/cluster_labels.csv ─────────────────────
# Created once; never overwritten. Contains:
#   cluster_id  — integer cluster numbers at active_resolution
# Users add annotation columns to this file for iterative annotation.
# cluster_id itself is treated as the first annotation pass (identity mapping).

labels_file <- file.path(out_dir, "annotation",
                         paste0("res", active_resolution),
                         "cluster_labels.csv")

if (!file.exists(labels_file)) {
  dir.create(dirname(labels_file), showWarnings = FALSE, recursive = TRUE)
  cluster_col_active <- paste0("cluster_res", active_resolution)
  cluster_ids_vec    <- sort(unique(as.integer(
    as.character(clustering_dt[[cluster_col_active]])
  )))
  labels_template    <- data.table::data.table(
    cluster_id = cluster_ids_vec
  )
  data.table::fwrite(labels_template, labels_file)
  message("Saved: res", active_resolution, "/cluster_labels.csv  (",
          length(cluster_ids_vec), " clusters — add annotation columns to annotate)")
} else {
  message("Kept : res", active_resolution, "/cluster_labels.csv  (already exists)")
}

message("\n═══════════════════════════════════════════════════════════════")
message("Step 04 complete")
message("═══════════════════════════════════════════════════════════════")
