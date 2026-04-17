# =============================================================================
# Pearson Residual PCA Utility Functions
# =============================================================================
#
# Source: Approach1_scpearson/code/batched_pearson_pca_sketch.R
#
# Compute PCA via quasi-Poisson Pearson residuals (Lause 2021) on large
# gene-expression data split across many batch files, without ever loading
# all cells into memory at once.
#
# Mathematical background:
#   For gene g, cell c with raw count y_{gc}, total counts n_c, global gene
#   frequency p_g (per batch), and quasi-Poisson overdispersion phi:
#
#     ytilde_{gc}  =  y_{gc} / sqrt(p_g * n_c * phi)          [variance-stabilised]
#     r_{gc}       =  ytilde_{gc} - sqrt(p_g/phi) * sqrt(n_c) [Pearson residual]
#
#   The genes x genes cross-product of Pearson residuals decomposes additively
#   across batches, enabling streaming PCA without a full cells x genes matrix.
#   The genes x genes cross-product of Pearson residuals is:
#
#     Q  =  ytilde ytilde^T
#            - ytilde_muhat
#            - ytilde_muhat^T
#            + sum_c(n_c) * sqrt(p/phi) sqrt(p/phi)^T
#
#   where  ytilde_muhat = [sum_c ytilde_c * sqrt(n_c)] * sqrt(p/phi)^T
#                                  (genes x 1)              (1 x genes)
#
# Passes:
#   Pass 1  — per-file gene frequencies (grate, genes x n_files) and totalcounts.
#             Each batch file is treated as one biological batch for the purpose
#             of gene-frequency estimation, matching the behaviour of
#             scPearsonPCA::gene_frequency(batch_variable=...).
#   Pass 2a — accumulate ytilde cross-products WITHOUT clipping → derive mean and SD of pearson residuals.
#   Pass 2b — (only when scale_max < Inf) re-accumulate WITH clipping.
#   SVD     — decompose centered/scaled Q to get loadings.
#   Pass 3  — project each batch into PCA space in user-specified cell chunks.
#
# Assumed input format:
#   Each batch is a sparse Matrix (genes x cells) stored as an .rds file.
#   All batches share the same row set (genes).
#   totalcounts should be supplied when batch files contain only a subset of
#   genes; if NULL, colSums() of the batch file is used (correct only when
#   ALL genes are present in every file).
#   grate, if supplied, must be either:
#     - a named numeric vector (genes) for a single global gene frequency, OR
#     - a genes x n_files matrix where column b gives gene frequencies for
#       batch file b (matching scPearsonPCA's per-batch grate behaviour).
#
# Dependencies: Matrix, RSpectra
# =============================================================================

library(Matrix)
library(RSpectra)


# =============================================================================
#  PASS 1 — compute per-file gene frequencies (grate) and per-cell total counts
# =============================================================================
#
# grate is a genes x n_files matrix. Column b contains gene frequencies for
# the cells in batch file b:
#   grate[g, b] = (sum of gene g counts in file b) /
#                 (total counts across all genes in file b)
# Each column sums to 1.
#
# Returns: list(n, grate, totalcounts)
#
pearson_batched_summary_stats <- function(batch_files,
                                          totalcounts = NULL) {
  n_total    <- 0L
  tc_list    <- if (is.null(totalcounts)) list() else NULL
  grate_cols <- vector("list", length(batch_files))

  for (i in seq_along(batch_files)) {
    message("Pass 1 – reading ", basename(batch_files[[i]]))
    x_batch <- readRDS(batch_files[[i]])

    if (is.null(totalcounts)) {
      tc_batch <- Matrix::colSums(x_batch)
      tc_list  <- c(tc_list, list(tc_batch))
    }

    n_total <- n_total + ncol(x_batch)

    gene_sums_b  <- Matrix::rowSums(x_batch)
    total_cts_b  <- sum(Matrix::colSums(x_batch))
    grate_b      <- gene_sums_b / total_cts_b
    grate_cols[[i]] <- grate_b / sum(grate_b)       # normalise to sum 1
  }

  if (is.null(totalcounts)) totalcounts <- unlist(tc_list)

  grate <- do.call(cbind, grate_cols)
  rownames(grate) <- rownames(x_batch)
  colnames(grate) <- basename(batch_files)

  list(n = n_total, grate = grate, totalcounts = totalcounts)
}


# =============================================================================
#  PASS 2 — accumulate ytilde cross-products and summary statistics
# =============================================================================
#
# Called twice when scale_max < Inf:
#   Pass 2a: clip_vals = NULL  → used to derive mean_pr and sd_pr
#   Pass 2b: clip_vals = named vector per gene → re-accumulate after clipping
#
# Returns list with accumulators for building Q and computing mean Pearson residuals.
#
pearson_batched_crossprod <- function(batch_files, genes, grate, totalcounts,
                                      phi, clip_vals = NULL) {
  n_files         <- length(batch_files)
  n_total         <- 0L
  sum_tc          <- numeric(n_files)
  sum_sqrt_tc     <- numeric(n_files)
  sum_ytilde      <- setNames(numeric(length(genes)), genes)
  ytilde_root_tc  <- matrix(0, nrow = length(genes), ncol = n_files,
                             dimnames = list(genes, basename(batch_files)))
  ytilde_ytilde_t <- Matrix::Matrix(0,
                                    nrow = length(genes), ncol = length(genes),
                                    dimnames = list(genes, genes))

  for (i in seq_along(batch_files)) {
    message(if (is.null(clip_vals)) "Pass 2a" else "Pass 2b",
            " – reading ", basename(batch_files[[i]]))
    x_batch  <- readRDS(batch_files[[i]])[genes, , drop = FALSE]
    tc       <- totalcounts[colnames(x_batch)]
    sqrt_tc  <- sqrt(tc)
    grate_i  <- grate[genes, i] # gene frequencies for this file
    
    # ytilde = D_g^{-1} %*% X %*% D_c^{-1}
    #   D_g = diag(sqrt(grate_i * phi)),  D_c = diag(sqrt(tc))
    D_g_inv  <- Matrix::Diagonal(x = sqrt(1 / (grate_i * phi)))
    D_c_inv  <- Matrix::Diagonal(x = 1 / sqrt_tc)
    ytilde   <- D_g_inv %*% x_batch %*% D_c_inv
    dimnames(ytilde) <- dimnames(x_batch)

    if (!is.null(clip_vals)) {
      ytilde <- clip_ytilde(ytilde, clip_vals, grate_i, phi, tc)
    }

    n_total             <- n_total + ncol(x_batch)
    sum_tc[i]           <- sum(tc)
    sum_sqrt_tc[i]      <- sum(sqrt_tc)
    sum_ytilde          <- sum_ytilde + Matrix::rowSums(ytilde)
    ytilde_root_tc[, i] <- as.numeric(ytilde %*% sqrt_tc)
    ytilde_ytilde_t     <- ytilde_ytilde_t + Matrix::tcrossprod(ytilde)
  }

  list(ytilde_ytilde_t = ytilde_ytilde_t,
       ytilde_root_tc  = ytilde_root_tc,
       sum_ytilde      = sum_ytilde,
       sum_tc          = sum_tc,
       sum_sqrt_tc     = sum_sqrt_tc,
       n_total         = n_total)
}


# =============================================================================
#  Helper: clip ytilde values so that Pearson residuals do not exceed clip_vals
# =============================================================================
# ...do not exceed: clip_vals[gene] above the mean
#
# The Pearson residual is:   r_{gc} = ytilde_{gc} - sqrt(p_g/phi) * sqrt(n_c)   
# Clipping r_{gc} <= clip_vals[g] is equivalent to:
#   ytilde_{gc} <= clip_vals[g] + sqrt(p_g/phi) * sqrt(n_c)
#
clip_ytilde <- function(ytilde, clip_vals, grate, phi, tc) {
  gene_names   <- rownames(ytilde)[(ytilde@i + 1L)]
  cell_idx     <- rep.int(seq_len(ncol(ytilde)), diff(ytilde@p))
  cell_names   <- colnames(ytilde)[cell_idx]
  ceiling_vals <- clip_vals[gene_names] +
                  sqrt(grate[gene_names] / phi) * sqrt(tc[cell_names])
  ytilde@x     <- pmin(ytilde@x, ceiling_vals)
  ytilde
}


# =============================================================================
#  Helper: build Q from accumulated cross-product statistics
# =============================================================================
#
# grate is a genes x n_files matrix; cp$ytilde_root_tc is also genes x n_files.
#
# Q  =  ytilde ytilde^T
#        - ytilde_root_tc %*% t(sqrt_gp_mat)
#        - sqrt_gp_mat %*% t(ytilde_root_tc)
#        + sqrt_gp_mat %*% diag(sum_tc) %*% t(sqrt_gp_mat)
#
# This is the per-batch generalisation of the single-batch formula and matches
# scPearsonPCA::sparse_quasipoisson_pca_seurat_batch() exactly.
#
build_qp <- function(cp, grate, phi) {
  sqrt_gp_mat  <- sqrt(grate / phi)                       # genes x n_files
  ytilde_muhat <- Matrix::Matrix(cp$ytilde_root_tc) %*%
                  Matrix::t(Matrix::Matrix(sqrt_gp_mat))  # genes x genes
  correction   <- Matrix::Matrix(sqrt_gp_mat) %*%
                  Matrix::Diagonal(x = cp$sum_tc) %*%
                  Matrix::t(Matrix::Matrix(sqrt_gp_mat))  # genes x genes
  cp$ytilde_ytilde_t -
    ytilde_muhat - Matrix::t(ytilde_muhat) +
    correction
}


# =============================================================================
#  Helper: mean Pearson residual from accumulated stats
# =============================================================================
#  mean(r_g) = mean_c(ytilde_{gc})
#              - sum_b( sqrt(grate[g,b] / phi) * sum_sqrt_tc_b ) / N          #
#                                                                              #
# With per-batch grate, the expected mean of ytilde_{gc} for a cell in batch b
# is sqrt(grate[g,b]/phi) * mean_sqrt_tc_b, so the global mean is the weighted
# sum across batches (weighted by n_b/N).  In accumulated form:
#   expected = (1/N) * sum_b( sqrt(grate[:,b]/phi) * sum_sqrt_tc_b )
#            = (1/N) * sqrt_gp_mat %*% sum_sqrt_tc_b
# --------------------------------------------------------------------------- #
compute_mean_pr <- function(cp, grate, phi) {
  N           <- cp$n_total
  sqrt_gp_mat <- sqrt(grate / phi)
  expected    <- as.numeric(Matrix::Matrix(sqrt_gp_mat) %*% cp$sum_sqrt_tc) / N
  cp$sum_ytilde / N - expected
}


# =============================================================================
#  Decompose Q: center, scale, SVD
# =============================================================================
compute_pearson_pca_loadings <- function(qp, N, mean_pr_obs, mean_pr_center,
                                         sd_pr, do.center, do.scale, npcs) {
  # mean_pr_obs    : mean of the (possibly clipped) Pearson residuals
  #                  used as the "observed" side of the centering correction
  # mean_pr_center : mean of the UNCLIPPED residuals
  #                  used as the centering target (matches the projection step)
  # These are the same when scale_max = Inf.

  if (do.center) {
    message("Centering Pearson residuals")
    inner_prod <- N * Matrix::Matrix(mean_pr_obs,    ncol = 1) %*%
                      Matrix::Matrix(mean_pr_center, nrow = 1)
    qp <- qp - inner_prod - Matrix::t(inner_prod) +
          N  * Matrix::tcrossprod(Matrix::Matrix(mean_pr_center, ncol = 1))
  }

  if (do.scale) {
    message("Scaling Pearson residuals")
    D  <- Matrix::Diagonal(x = 1 / sd_pr)
    qp <- D %*% qp %*% D
  }

  message("Computing truncated SVD (", npcs, " components)…")
  svdd <- RSpectra::svds(qp, k = npcs)

  feature_loadings <- svdd$u
  rownames(feature_loadings) <- names(sd_pr)
  colnames(feature_loadings) <- paste0("PC_", seq_len(npcs))

  list(loadings = feature_loadings,
       sdev     = sqrt(svdd$d / (N - 1)))
}


# =============================================================================
#  PASS 3 — project each batch into PCA space in manageable cell chunks
# =============================================================================
#
# chunk_size : max cells to densify at once. Peak memory per chunk:
#              chunk_size * n_genes * 8 bytes. Default Inf = whole batch.
#
pearson_batched_embeddings <- function(batch_files, genes, totalcounts, grate,
                                       phi, clip_vals,
                                       feature_loadings, mean_pr, sd_pr,
                                       do.center, do.scale,
                                       chunk_size = Inf) {
  embed_list <- vector("list", length(batch_files))

  for (i in seq_along(batch_files)) {
    message("Pass 3 – reading ", basename(batch_files[[i]]))
    x_batch  <- readRDS(batch_files[[i]])[genes, , drop = FALSE]
    tc       <- totalcounts[colnames(x_batch)]
    sqrt_tc  <- sqrt(tc)
    grate_i  <- grate[, i]                    # per-file gene frequencies
    sqrt_gp_i <- sqrt(grate_i / phi)

    D_g_inv  <- Matrix::Diagonal(x = sqrt(1 / (grate_i * phi)))
    D_c_inv  <- Matrix::Diagonal(x = 1 / sqrt_tc)
    ytilde   <- D_g_inv %*% x_batch %*% D_c_inv
    dimnames(ytilde) <- dimnames(x_batch)

    if (!is.null(clip_vals)) {
      ytilde <- clip_ytilde(ytilde, clip_vals, grate_i, phi, tc)
    }

    n_cells      <- ncol(ytilde)
    chunk_size_i <- min(n_cells, chunk_size)
    chunk_starts <- seq(1, n_cells, by = chunk_size_i)
    chunk_embeds <- vector("list", length(chunk_starts))

    for (j in seq_along(chunk_starts)) {
      idx <- chunk_starts[j]:min(chunk_starts[j] + chunk_size_i - 1L, n_cells)

      # Pearson residual for chunk (cells x genes, dense)
      scx <- as.matrix(Matrix::t(ytilde[, idx, drop = FALSE])) -
             matrix(sqrt_tc[idx], ncol = 1) %*%
             matrix(sqrt_gp_i,    nrow = 1)

      if (do.center && do.scale) {
        scx <- scale(scx, center = mean_pr, scale = sd_pr)
      } else if (do.center) {
        scx <- scale(scx, center = mean_pr, scale = FALSE)
      } else if (do.scale) {
        scx <- scale(scx, center = FALSE,   scale = sd_pr)
      }

      chunk_embeds[[j]] <- scx %*% feature_loadings
    }

    embed_list[[i]] <- do.call(rbind, chunk_embeds)
  }

  do.call(rbind, embed_list)
}


# =============================================================================
#  Top-level wrapper
# =============================================================================
#
# batch_files  : character vector of .rds file paths, each a sparse genes x
#                cells RAW COUNTS matrix.
# genes        : genes to use for PCA (must be present in all batch files).
# totalcounts  : optional named numeric vector of per-cell total UMI counts
#                computed from ALL genes. If NULL, estimated via colSums()
#                during Pass 1 — only correct when files contain all genes.
# grate        : optional gene frequency matrix. If NULL, computed in Pass 1.
# phi          : quasi-Poisson variance inflation factor (default 1.01).
# npcs         : number of principal components.
# scale_max    : clip Pearson residuals at this many SDs above mean (Inf = none).
# do.center    : subtract mean Pearson residual before PCA.
# do.scale     : divide by SD of Pearson residual before PCA.
# chunk_size   : max cells per chunk in Pass 3 embedding projection.
#
# Returns list(embeddings, loadings, sdev, mean_pr, sd_pr, grate, phi)
#
pearson_pca_from_batches <- function(batch_files,
                                     genes,
                                     totalcounts = NULL,
                                     grate       = NULL,
                                     phi         = 1.01,
                                     npcs        = 50,
                                     scale_max   = 10,
                                     do.center   = TRUE,
                                     do.scale    = TRUE,
                                     chunk_size  = Inf) {

  # ---- Pass 1: grate and totalcounts (skip if both supplied) ----
  if (is.null(grate) || is.null(totalcounts)) {
    message("=== Pass 1: gene frequencies and total counts ===")
    stats1 <- pearson_batched_summary_stats(batch_files,
                                            totalcounts = totalcounts)
    if (is.null(grate))       grate       <- stats1$grate
    if (is.null(totalcounts)) totalcounts <- stats1$totalcounts
  }
  # Accept either a global vector or a per-file matrix.
  # Coerce vector → matrix (same frequencies replicated for each file).
  # Coerce global vector → per-file matrix if needed
  if (is.null(dim(grate))) {
    grate <- matrix(grate[genes], nrow = length(genes), ncol = length(batch_files),
                    dimnames = list(genes, basename(batch_files)))
  } else {
    grate <- grate[genes, , drop = FALSE]
  }
  grate[grate == 0] <- .Machine$double.eps    # guard against p_g = 0

  # ---- Pass 2a: initial cross-products (no clipping) ----
  message("=== Pass 2a: accumulating cross-products (unclipped) ===")
  cp_init  <- pearson_batched_crossprod(batch_files, genes, grate,
                                        totalcounts, phi, clip_vals = NULL)
  N        <- cp_init$n_total
  qp_init  <- build_qp(cp_init, grate, phi)
  mean_pr  <- compute_mean_pr(cp_init, grate, phi)

  # SD of Pearson residuals: Var(r_g) = (diag(Q) - N * mean^2) / (N-1)
  sd_pr    <- sqrt(pmax((Matrix::diag(qp_init) - N * mean_pr^2) / (N - 1), 0))

  # ---- Pass 2b: re-accumulate with clipping (only if scale_max < Inf) ----
  clip_vals   <- NULL
  qp          <- qp_init
  mean_pr_obs <- mean_pr # "observed" (clipped) mean, used for centering

  if (scale_max < Inf) {
    message("=== Pass 2b: re-accumulating with clipping (scale_max=",
            scale_max, ") ===")
    clip_vals   <- scale_max * sd_pr + mean_pr          # per-gene clip threshold
    cp_clip     <- pearson_batched_crossprod(batch_files, genes, grate,
                                             totalcounts, phi,
                                             clip_vals = clip_vals)
    qp          <- build_qp(cp_clip, grate, phi)
    mean_pr_obs <- compute_mean_pr(cp_clip, grate, phi)       # clipped mean
  }

  # ---- Compute loadings via SVD ----
  message("=== Computing PCA loadings ===")
  svd_result <- compute_pearson_pca_loadings(
    qp             = qp,
    N              = N,
    mean_pr_obs    = mean_pr_obs,         # clipped mean (= unclipped if no clipping)
    mean_pr_center = mean_pr,             # unclipped mean (centering target)
    sd_pr          = sd_pr,
    do.center      = do.center,
    do.scale       = do.scale,
    npcs           = npcs
  )

  # ---- Pass 3: cell embeddings ----
  # Center embeddings by the PRE-CLIP mean, matching scPearsonPCA convention:
  # the clipped mean (mean_pr_obs) is used only for the asymmetric QP centering
  # correction; the unclipped mean (mean_pr) is the centering target for
  # projecting cells, both here and when projecting new data later.
  message("=== Pass 3: computing cell embeddings ===")
  embeddings <- pearson_batched_embeddings(
    batch_files      = batch_files,
    genes            = genes,
    totalcounts      = totalcounts,
    grate            = grate,
    phi              = phi,
    clip_vals        = clip_vals,
    feature_loadings = svd_result$loadings,
    mean_pr          = mean_pr,         # pre-clip mean — centering target
    sd_pr            = sd_pr,
    do.center        = do.center,
    do.scale         = do.scale,
    chunk_size       = chunk_size
  )

  list(
    embeddings = embeddings,            # all cells x npcs  
    loadings   = svd_result$loadings,   # genes x npcs
    sdev       = svd_result$sdev,       # length npcs
    mean_pr    = mean_pr,        # per-gene mean Pearson residual, pre-clip
    sd_pr      = sd_pr,                 # per-gene SD Pearson residual
    grate      = grate,
    phi        = phi
  )
}


# =============================================================================
#  Example usage (not run)                                                      #
# =============================================================================
if (FALSE) {

  batch_files <- list.files("path/to/batches", pattern = "\\.rds$",
                             full.names = TRUE)

  # Highly variable genes to use for PCA
  hvg <- readLines("path/to/hvg.txt")

  # Total UMI per cell from ALL genes (recommended to supply pre-computed)
  totalcounts <- readRDS("path/to/totalcounts.rds")   # named numeric vector

  result <- pearson_pca_from_batches(
    batch_files = batch_files,
    genes       = hvg,
    totalcounts = totalcounts,   # omit to estimate from batch files
    phi         = 1.01,
    npcs        = 50,
    scale_max   = 10,
    do.center   = TRUE,
    do.scale    = TRUE,
    chunk_size  = 5000           # cells per chunk in Pass 3
  )

  # result$embeddings : all cells x 50 PCs
  # result$loadings   : genes x 50 PCs
  # result$sdev       : per-PC standard deviations
  # result$mean_pr    : per-gene mean Pearson residual (for projecting new data)
  # result$sd_pr      : per-gene SD Pearson residual
}
