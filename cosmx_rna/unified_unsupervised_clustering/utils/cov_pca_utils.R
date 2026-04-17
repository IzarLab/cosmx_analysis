# =============================================================================
# Covariance PCA Utility Functions
# =============================================================================
#
# Source: Approach1_scpearson/code/batched_cov_pca_sketch.R
#
# Compute PCA on a large gene-expression dataset split across many batch files
# using the covariance approach, without ever loading all cells into memory.
#
# Mathematical background:
#   For scaled/centered data X_sc = (X_obs - mu) / sigma, the unnormalised
#   covariance matrix is:
#
#     X_sc X_sc^T  =  D^{-1} * [X_obs X_obs^T
#                                - n * mu_obs * mu^T
#                                - n * mu * mu_obs^T
#                                + n * mu * mu^T] * D^{-1}
#
#   All quantities on the RHS are sums across batches, enabling streaming PCA.
#
# Input format: each batch is a pre-normalised sparse Matrix (genes x cells).
# For PCA on raw counts, use cov_pca_from_raw_batches() which applies
# library-size normalisation (and optionally log1p) on-the-fly.
#
# Dependencies: Matrix, RSpectra
# =============================================================================

library(Matrix)
library(RSpectra)


# =============================================================================
#  Helper: row-wise standard deviations for a sparse matrix
# =============================================================================
row_sds_sparse <- function(x) {
  n       <- ncol(x)
  means   <- Matrix::rowMeans(x)
  x2      <- x; x2@x <- x2@x^2
  mean_sq <- Matrix::rowSums(x2) / n
  sqrt(pmax(mean_sq - means^2, 0))
}


# =============================================================================
#  Helper: clip expression to scale_max SDs above the (pre-clip) gene mean
# =============================================================================
threshold_by_gene <- function(x, max_val_thresh) {
  vals       <- x@x
  gene_names <- rownames(x)[(x@i + 1L)]
  cap        <- max_val_thresh[gene_names]
  vals[vals > cap] <- cap[vals > cap]
  x@x        <- vals
  x
}


# =============================================================================
#  PASS 1 — accumulate per-gene sums and sums-of-squares to get mu and sigma
# =============================================================================
#
# Returns: list(n, mu, sigma)
#   n     : total cell count
#   mu    : global per-gene means  (length = #genes)
#   sigma : global per-gene SDs    (length = #genes)
#
batched_summary_stats <- function(batch_files, genes) {
  n_total <- 0L
  sum_x   <- setNames(numeric(length(genes)), genes)
  sum_x2  <- setNames(numeric(length(genes)), genes)

  for (f in batch_files) {
    message("Pass 1 – reading ", basename(f))
    x_batch  <- readRDS(f)[genes, , drop = FALSE]

    n_total  <- n_total + ncol(x_batch)
    sum_x    <- sum_x  + Matrix::rowSums(x_batch)

    x2_batch <- x_batch; x2_batch@x <- x2_batch@x^2
    sum_x2   <- sum_x2 + Matrix::rowSums(x2_batch)
  }

  mu    <- sum_x  / n_total
  sigma <- sqrt(pmax(sum_x2 / n_total - mu^2, 0))

  list(n = n_total, mu = mu, sigma = sigma)
}


# =============================================================================
#  PASS 2 — accumulate X_obs X_obs^T and mu_obs using clipped expression
# =============================================================================
batched_crossprod <- function(batch_files, genes, mu, sigma, scale_max = 10) {
  max_val_thresh <- mu + scale_max * sigma

  n_total <- 0L
  sum_obs <- setNames(numeric(length(genes)), genes)
  XtX     <- Matrix::Matrix(0, nrow = length(genes), ncol = length(genes),
                             dimnames = list(genes, genes))

  for (f in batch_files) {
    message("Pass 2 – reading ", basename(f))
    x_batch  <- readRDS(f)[genes, , drop = FALSE]
    x_obs    <- threshold_by_gene(x_batch, max_val_thresh)

    n_total  <- n_total + ncol(x_obs)
    sum_obs  <- sum_obs + Matrix::rowSums(x_obs)
    XtX      <- XtX    + Matrix::tcrossprod(x_obs)
  }

  mu_obs <- sum_obs / n_total
  list(n = n_total, mu_obs = mu_obs, XtX = XtX, max_val_thresh = max_val_thresh)
}


# =============================================================================
#  Build covariance matrix and run truncated SVD
# =============================================================================
compute_pca_loadings <- function(stats1, stats2, npcs = 50) {
  n      <- stats1$n
  mu     <- stats1$mu
  sigma  <- stats1$sigma
  mu_obs <- stats2$mu_obs
  XtX    <- stats2$XtX

  inner_prod <- n * (mu_obs %*% t(mu))
  cov_raw    <- XtX - inner_prod - t(inner_prod) + n * tcrossprod(mu)

  D_inv  <- Matrix::Diagonal(x = 1 / sigma)
  cov_sc <- D_inv %*% cov_raw %*% D_inv

  message("Computing truncated SVD (", npcs, " components)…")
  svdd <- RSpectra::svds(cov_sc, k = npcs)

  feature_loadings <- svdd$u
  rownames(feature_loadings) <- names(mu)
  colnames(feature_loadings) <- paste0("PC_", seq_len(npcs))
  sdev <- sqrt(svdd$d / (n - 1))

  list(loadings = feature_loadings, sdev = sdev,
       mu = mu, sigma = sigma, max_val_thresh = stats2$max_val_thresh)
}


# =============================================================================
#  PASS 3 — project each batch into PCA space
# =============================================================================
#
# chunk_size : max cells to transpose/scale at once in Pass 3.
#              Peak memory: chunk_size * n_genes * 8 bytes.
#              Default Inf = entire batch at once.
#
batched_embeddings <- function(batch_files, genes, pca_result,
                               chunk_size = Inf) {
  loadings       <- pca_result$loadings
  mu             <- pca_result$mu
  sigma          <- pca_result$sigma
  max_val_thresh <- pca_result$max_val_thresh

  embed_list <- vector("list", length(batch_files))

  for (i in seq_along(batch_files)) {
    message("Pass 3 – reading ", basename(batch_files[[i]]))
    x_batch  <- readRDS(batch_files[[i]])[genes, , drop = FALSE]
    x_obs    <- threshold_by_gene(x_batch, max_val_thresh)

    n_cells      <- ncol(x_obs)
    chunk_size_i <- min(n_cells, chunk_size)
    chunk_starts <- seq(1, n_cells, by = chunk_size_i)
    chunk_embeds <- vector("list", length(chunk_starts))

    for (j in seq_along(chunk_starts)) {
      idx            <- chunk_starts[j]:min(chunk_starts[j] + chunk_size_i - 1, n_cells)
      x_sc           <- scale(Matrix::t(x_obs[, idx, drop = FALSE]),
                               center = mu, scale = sigma)
      chunk_embeds[[j]] <- x_sc %*% loadings
    }

    embed_list[[i]] <- do.call(rbind, chunk_embeds)
  }

  do.call(rbind, embed_list)
}


# =============================================================================
#  Top-level wrapper — pre-normalised data
# =============================================================================
#
# batch_files : .rds paths, each a (sparse) genes x cells NORMALISED matrix.
# genes       : character vector of genes to use.
# npcs        : number of principal components.
# scale_max   : clipping threshold in SDs.
# chunk_size  : max cells per chunk in Pass 3.
#
# Returns list(embeddings, loadings, sdev, mu, sigma)
#
pca_from_batches <- function(batch_files,
                             genes,
                             npcs       = 50,
                             scale_max  = 10,
                             chunk_size = Inf) {

  message("=== Pass 1: global mean and SD ===")
  stats1 <- batched_summary_stats(batch_files, genes)

  # Drop genes with zero variance (cannot be scaled)
  keep         <- stats1$sigma > 0
  if (!all(keep)) {
    message("  Dropping ", sum(!keep), " zero-variance genes")
    genes        <- genes[keep]
    stats1$mu    <- stats1$mu[keep]
    stats1$sigma <- stats1$sigma[keep]
  }

  message("=== Pass 2: accumulate cross-products ===")
  stats2 <- batched_crossprod(batch_files, genes,
                              mu = stats1$mu, sigma = stats1$sigma,
                              scale_max = scale_max)

  message("=== Computing PCA loadings ===")
  pca_result <- compute_pca_loadings(stats1, stats2, npcs = npcs)

  message("=== Pass 3: computing cell embeddings ===")
  embeddings <- batched_embeddings(batch_files, genes, pca_result,
                                   chunk_size = chunk_size)

  list(
    embeddings = embeddings,
    loadings   = pca_result$loadings,
    sdev       = pca_result$sdev,
    mu         = pca_result$mu,
    sigma      = pca_result$sigma
  )
}


# =============================================================================
#  Wrapper for raw count data: normalise on-the-fly, then run covariance PCA
# =============================================================================
#
# Use this when batch files contain raw counts (as produced by 01_extract_counts.R).
# Normalisation is applied in-memory for each batch during each pass — no
# additional shard files are written to disk.
#
# batch_files      : .rds paths, each a sparse genes x cells RAW COUNTS matrix.
# genes            : genes to use for PCA.
# mean_totalcounts : global scaling target (floor of mean per-cell total counts).
# apply_log1p      : apply log1p after library-size scaling (default TRUE).
# npcs, scale_max, chunk_size : as in pca_from_batches().
#
cov_pca_from_raw_batches <- function(batch_files,
                                     genes,
                                     mean_totalcounts,
                                     apply_log1p = TRUE,
                                     npcs        = 50,
                                     scale_max   = 10,
                                     chunk_size  = Inf) {

  # Inline normalisation function applied inside each pass
  normalise_batch <- function(x_raw) {
    tc     <- Matrix::colSums(x_raw)
    tc[tc == 0] <- 1L                       # guard zero-libsize cells
    x_norm <- x_raw %*% Matrix::Diagonal(x = mean_totalcounts / tc)
    if (apply_log1p) x_norm@x <- log1p(x_norm@x)
    x_norm
  }

  # Wrap the three internal functions to apply normalisation before accumulation

  # ---- Pass 1: global mean and SD (on normalised data) ----
  message("=== Pass 1: global mean and SD (normalising on-the-fly) ===")
  n_total <- 0L
  sum_x   <- setNames(numeric(length(genes)), genes)
  sum_x2  <- setNames(numeric(length(genes)), genes)

  for (f in batch_files) {
    message("Pass 1 – reading ", basename(f))
    x_batch  <- normalise_batch(readRDS(f))[genes, , drop = FALSE]
    n_total  <- n_total + ncol(x_batch)
    sum_x    <- sum_x  + Matrix::rowSums(x_batch)
    x2       <- x_batch; x2@x <- x2@x^2
    sum_x2   <- sum_x2 + Matrix::rowSums(x2)
  }
  mu    <- sum_x  / n_total
  sigma <- sqrt(pmax(sum_x2 / n_total - mu^2, 0))

  # Drop zero-variance genes
  keep <- sigma > 0
  if (!all(keep)) {
    message("  Dropping ", sum(!keep), " zero-variance genes")
    genes <- genes[keep]
    mu    <- mu[keep]
    sigma <- sigma[keep]
  }

  stats1 <- list(n = n_total, mu = mu, sigma = sigma)

  # ---- Pass 2: cross-products (on normalised data) ----
  message("=== Pass 2: accumulate cross-products (normalising on-the-fly) ===")
  max_val_thresh <- mu + scale_max * sigma
  sum_obs <- setNames(numeric(length(genes)), genes)
  XtX     <- Matrix::Matrix(0, nrow = length(genes), ncol = length(genes),
                             dimnames = list(genes, genes))

  for (f in batch_files) {
    message("Pass 2 – reading ", basename(f))
    x_batch  <- normalise_batch(readRDS(f))[genes, , drop = FALSE]
    x_obs    <- threshold_by_gene(x_batch, max_val_thresh)
    sum_obs  <- sum_obs + Matrix::rowSums(x_obs)
    XtX      <- XtX    + Matrix::tcrossprod(x_obs)
  }
  mu_obs <- sum_obs / n_total
  stats2 <- list(n = n_total, mu_obs = mu_obs, XtX = XtX,
                 max_val_thresh = max_val_thresh)

  # ---- Loadings ----
  message("=== Computing PCA loadings ===")
  pca_result <- compute_pca_loadings(stats1, stats2, npcs = npcs)

  # ---- Pass 3: embeddings (on normalised data) ----
  message("=== Pass 3: computing cell embeddings (normalising on-the-fly) ===")
  embed_list <- vector("list", length(batch_files))

  for (i in seq_along(batch_files)) {
    message("Pass 3 – reading ", basename(batch_files[[i]]))
    x_batch  <- normalise_batch(readRDS(batch_files[[i]]))[genes, , drop = FALSE]
    x_obs    <- threshold_by_gene(x_batch, pca_result$max_val_thresh)

    n_cells      <- ncol(x_obs)
    chunk_size_i <- min(n_cells, chunk_size)
    chunk_starts <- seq(1, n_cells, by = chunk_size_i)
    chunk_embeds <- vector("list", length(chunk_starts))

    for (j in seq_along(chunk_starts)) {
      idx            <- chunk_starts[j]:min(chunk_starts[j] + chunk_size_i - 1, n_cells)
      x_sc           <- scale(Matrix::t(x_obs[, idx, drop = FALSE]),
                               center = pca_result$mu, scale = pca_result$sigma)
      chunk_embeds[[j]] <- x_sc %*% pca_result$loadings
    }

    embed_list[[i]] <- do.call(rbind, chunk_embeds)
  }

  embeddings <- do.call(rbind, embed_list)

  list(
    embeddings = embeddings,
    loadings   = pca_result$loadings,
    sdev       = pca_result$sdev,
    mu         = pca_result$mu,
    sigma      = pca_result$sigma
  )
}
