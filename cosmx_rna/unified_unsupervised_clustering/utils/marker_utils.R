# =============================================================================
# Marker Gene Utility Functions — Batched Welch t-test Markers
# =============================================================================
#
# Source: Approach1_scpearson/code/batched_marker_stats_sketch.R
#
# Compute cluster-vs-rest Welch t-test marker gene statistics for a large
# gene-expression dataset split across many .rds batch files, without loading
# all cells into memory at once.
#
# Algorithm:
#   Pass 1 (optional) — compute per-cell library sizes. Skipped when
#                        totalcounts is supplied.
#   Pass 2            — accumulate sum_expr, sum_sq, n_pos per gene per cluster.
#   Finalise          — derive mean / sd / pct, then Welch t-test (cluster vs. rest).
#
# Normalisation: y_{gc} = raw_{gc} * (meanLibsize / libsize_c)
#   Supply totalcounts (named per-cell library sizes from ALL genes) to avoid
#   an extra pass; otherwise estimated from colSums of each batch file (only
#   correct when batch files contain ALL genes).
#
# Output format matches Seurat::FindAllMarkers() with extensions:
#   p_val, p_val_adj, avg_log2FC, pct.1, pct.2, cluster_ncells,
#   cluster_expr, clusterprime_expr, cluster, feats
#
# Dependencies: Matrix, data.table, parallel
# =============================================================================

library(Matrix)
library(data.table)


# =============================================================================
#  Pass 1 — compute per-cell library sizes (skip when totalcounts supplied)
# =============================================================================
#
# Returns a named numeric vector of per-cell colSums across ALL genes in each
# batch file. Only valid when batch files contain the complete gene set.
#
batched_marker_libsizes <- function(batch_files) {
  tc_list <- vector("list", length(batch_files))
  for (i in seq_along(batch_files)) {
    message("Pass 1 – reading ", basename(batch_files[[i]]))
    x_batch      <- readRDS(batch_files[[i]])
    tc_list[[i]] <- Matrix::colSums(x_batch)
  }
  unlist(tc_list)
}


# =============================================================================
#  Pass 2 — accumulate per-cluster, per-gene summary statistics
# =============================================================================
#
# For each batch:
#   1. Normalise columns:  x_norm = x_batch %*% diag(meanLibsize / libsize)
#   2. For each cluster present in this batch, accumulate:
#        sum_expr[g, k]  +=  rowSums(x_norm for cells in k)
#        sum_sq[g, k]    +=  rowSums(x_norm^2 for cells in k)
#        n_pos[g, k]     +=  number of cells in k with x_norm[g] > 0
#
# Returns list(n_cells, sum_expr, sum_sq, n_pos)
#
batched_marker_accumulate <- function(batch_files, genes, metadata,
                                      cluster_column, cell_id_column,
                                      totalcounts, meanLibsize) {

  meta     <- data.table::as.data.table(metadata)
  meta     <- meta[!is.na(meta[[cluster_column]])]
  clusters <- sort(unique(as.character(meta[[cluster_column]])))
  n_genes  <- length(genes)
  n_clust  <- length(clusters)

  # n_cells per cluster from metadata (no I/O needed)
  n_cells <- setNames(
    as.integer(tabulate(match(as.character(meta[[cluster_column]]), clusters))),
    clusters
  )

  # Accumulators: genes x clusters (dense, usually small)
  sum_expr <- matrix(0,  nrow = n_genes, ncol = n_clust, dimnames = list(genes, clusters))
  sum_sq   <- matrix(0,  nrow = n_genes, ncol = n_clust, dimnames = list(genes, clusters))
  n_pos    <- matrix(0L, nrow = n_genes, ncol = n_clust, dimnames = list(genes, clusters))

  for (f in batch_files) {
    message("Pass 2 – reading ", basename(f))
    x_batch     <- readRDS(f)[genes, , drop = FALSE]

    # Only process cells that have cluster assignments
    batch_cells <- intersect(colnames(x_batch), meta[[cell_id_column]])
    if (length(batch_cells) == 0L) next

    x_batch     <- x_batch[, batch_cells, drop = FALSE]
    tc          <- totalcounts[batch_cells]
    tc[tc == 0] <- 1L                               # guard zero-libsize cells

    # Normalise: x_norm = X %*% diag(meanLibsize / libsize)
    scale_fac   <- meanLibsize / tc
    x_norm      <- x_batch %*% Matrix::Diagonal(x = scale_fac)
    dimnames(x_norm) <- dimnames(x_batch)

    # Squared values (only stored non-zeros change)
    x_sq        <- x_norm
    x_sq@x      <- x_sq@x^2

    # Cluster-wise accumulation
    batch_meta <- meta[meta[[cell_id_column]] %in% batch_cells]
    batch_meta[[cluster_column]] <- as.character(batch_meta[[cluster_column]])

    for (ct in unique(batch_meta[[cluster_column]])) {
      cells_ct <- batch_meta[[cell_id_column]][batch_meta[[cluster_column]] == ct]

      sum_expr[, ct] <- sum_expr[, ct] + Matrix::rowSums(x_norm[,   cells_ct, drop = FALSE])
      sum_sq[, ct]   <- sum_sq[, ct]   + Matrix::rowSums(x_sq[,     cells_ct, drop = FALSE])
      # pct: fraction > 0 is the same pre- and post-normalisation (scale_fac > 0)
      n_pos[, ct]    <- n_pos[, ct]    +
        as.integer(Matrix::rowSums(x_batch[, cells_ct, drop = FALSE] > 0))
    }
  }

  list(n_cells  = n_cells,
       sum_expr = sum_expr,
       sum_sq   = sum_sq,
       n_pos    = n_pos)
}

# =============================================================================
#  Finalise — derive mean / sd / pct, then Welch t-test
# =============================================================================
#
# The Welch step is vectorised across genes; only the outer loop is over
# clusters. Set ncores > 1 to parallelise over clusters via parallel::mclapply.
#
batched_marker_welch <- function(acc, ncores = 1) {
  n_cells  <- acc$n_cells
  sum_expr <- acc$sum_expr
  sum_sq   <- acc$sum_sq
  n_pos    <- acc$n_pos
  clusters <- names(n_cells)
  genes    <- rownames(sum_expr)

  # ---- Derive per-cluster statistics -----------------------------------------
  mean_expr <- sweep(sum_expr, 2, n_cells, "/")

  sd_expr <- matrix(0, nrow = nrow(mean_expr), ncol = ncol(mean_expr),
                    dimnames = dimnames(mean_expr))
  for (ct in clusters) {
    n <- n_cells[ct]
    if (n > 1L) {
      sd_expr[, ct] <- sqrt(
        pmax((sum_sq[, ct] - n * mean_expr[, ct]^2) / (n - 1L), 0)
      )
    }
  }

  pct_expr <- sweep(n_pos, 2, n_cells, "/")

  # ---- Welch t-test: cluster ct vs. all others -------------------------------
  outl <- parallel::mclapply(clusters, function(ct) {

    n_ii     <- n_cells[ct]
    others   <- clusters[clusters != ct]
    n_others <- n_cells[others]
    n_prime  <- sum(n_others)
    w        <- n_others / n_prime            # normalised weights (sum to 1)

    mean_ii  <- mean_expr[, ct]
    sd_ii    <- sd_expr[,  ct]
    pct_ii   <- pct_expr[, ct]

    M_others <- mean_expr[, others, drop = FALSE]   # genes x |others|
    S_others <- sd_expr[,  others, drop = FALSE]
    P_others <- pct_expr[, others, drop = FALSE]

    mean_prime <- as.numeric(M_others %*% w)
    pct_prime  <- as.numeric(P_others %*% w)

    # grand.sd vectorised across genes:
    # grand.sd = sqrt( weighted.mean(S^2 + M^2, N) - weighted.mean(M, N)^2 )
    sd_prime_sq <- as.numeric((S_others^2 + M_others^2) %*% w) - mean_prime^2
    sd_prime    <- sqrt(pmax(sd_prime_sq, 0))

    # Welch t-test
    welch_var  <- sd_ii^2    / n_ii    + sd_prime^2 / n_prime
    welch_df   <- welch_var^2 / (
      (sd_ii^2    / n_ii)^2    / (n_ii    - 1L) +
      (sd_prime^2 / n_prime)^2 / (n_prime - 1L)
    )
    welch_t <- (mean_ii - mean_prime) / sqrt(welch_var)
    welch_p <- 2 * pt(-abs(welch_t), df = welch_df)

    data.table::data.table(
      p_val             = welch_p,
      p_val_adj         = p.adjust(welch_p, method = "BH"),
      avg_log2FC        = log2(mean_ii / mean_prime),
      pct.1             = pct_ii,
      pct.2             = pct_prime,
      cluster_ncells    = n_ii,
      cluster_expr      = mean_ii,
      clusterprime_expr = mean_prime,
      cluster           = ct,
      feats             = genes
    )
  }, mc.cores = ncores)

  data.table::rbindlist(outl)
}


# =============================================================================
#  Top-level wrapper
# =============================================================================
#
# batch_files    : character vector of .rds paths, each a sparse genes x cells
#                  raw counts matrix.
# genes          : genes to test (must be present in all batch files).
# metadata       : data.frame / data.table with cell_id_column and
#                  cluster_column. Cells with NA cluster are excluded.
# cluster_column : column in metadata giving cluster assignment.
# cell_id_column : column in metadata giving cell ID (default "cell_ID").
# totalcounts    : optional named numeric vector of per-cell library sizes
#                  (colSums over ALL genes). If NULL, estimated via colSums()
#                  of each batch file — only correct when files contain all genes.
# meanLibsize    : optional scalar normalisation target. If NULL, set to
#                  mean(totalcounts). Does not affect p-values or log2FC
#                  (it cancels in the t-statistic).
# ncores         : parallel threads for the Welch t-test step.
#
# Returns a data.table with columns:
#   p_val, avg_log2FC, pct.1, pct.2, cluster_ncells, cluster_expr,
#   cluster, feats
# (matches Seurat::FindAllMarkers() output format)
#
batched_find_markers <- function(batch_files,
                                 genes,
                                 metadata,
                                 cluster_column,
                                 cell_id_column = "cell_id",
                                 totalcounts    = NULL,
                                 meanLibsize    = NULL,
                                 ncores         = 1) {

  # ---- Pass 1: library sizes (skip if totalcounts supplied) ----
  if (is.null(totalcounts)) {
    message("=== Pass 1: computing per-cell library sizes ===")
    totalcounts <- batched_marker_libsizes(batch_files)
  }
  if (is.null(meanLibsize)) {
    meanLibsize <- mean(totalcounts)
  }

  # ---- Pass 2: per-cluster, per-gene accumulators ----
  message("=== Pass 2: accumulating per-cluster expression stats ===")
  acc <- batched_marker_accumulate(
    batch_files    = batch_files,
    genes          = genes,
    metadata       = metadata,
    cluster_column = cluster_column,
    cell_id_column = cell_id_column,
    totalcounts    = totalcounts,
    meanLibsize    = meanLibsize
  )

  # ---- Welch t-test ----
  message("=== Computing Welch t-test marker stats ===")
  batched_marker_welch(acc, ncores = ncores)
}
