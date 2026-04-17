# =============================================================================
# HVG Utility Functions — Batched VST Highly Variable Gene Selection
# =============================================================================
#
# Source: Approach1_scpearson/code/batched_hvg_sketch.R
#
# Identify highly variable genes (HVGs) from raw count data split across many
# batch files using the Seurat VST approach: fit a loess curve of
# log10(variance) ~ log10(mean) on raw counts, then rank genes by the ratio
# of observed to loess-predicted variance.
#
# Algorithm (one pass):
#   Pass 1  — accumulate per-gene sum_x and sum_x2 across all batch files.
#   Finalise — derive global mean and sample variance per gene.
#              Fit loess: log10(var) ~ log10(mean) for genes with var > 0.
#              variance.standardized = observed_var / 10^(loess fitted value).
#              Return top n_hvg genes ranked by variance.standardized.
#
# All batch files must contain raw (un-normalised) counts.
# All batch files must share the same gene set (rows), in the same order.
#
# Dependencies: Matrix
# =============================================================================

library(Matrix)


# =============================================================================
#  Pass 1 — accumulate per-gene sums and sums-of-squares across batches
# =============================================================================
#
# Returns list(n, sum_x, sum_x2)
#   n      : total cell count (scalar integer)
#   sum_x  : per-gene sum of raw expression   (length = #genes)
#   sum_x2 : per-gene sum of squared raw expression (length = #genes)
#
# genes : optional character vector of genes to accumulate. Default NULL
#         uses all rows present in the first batch file. Pass a subset to
#         restrict HVG selection (e.g. exclude NegProbe / FalseCode targets).
#
batched_hvg_summary_stats <- function(batch_files, genes = NULL) {
  n_total <- 0L
  sum_x   <- NULL
  sum_x2  <- NULL

  for (f in batch_files) {
    message("Pass 1 – reading ", basename(f))
    x_batch <- readRDS(f)

    if (!is.null(genes)) {
      x_batch <- x_batch[genes, , drop = FALSE]
    }

    if (is.null(sum_x)) {
      g      <- rownames(x_batch)
      sum_x  <- setNames(numeric(length(g)), g)
      sum_x2 <- setNames(numeric(length(g)), g)
    }

    n_total <- n_total + ncol(x_batch)
    sum_x   <- sum_x  + Matrix::rowSums(x_batch)

    x2_batch   <- x_batch
    x2_batch@x <- x2_batch@x^2
    sum_x2  <- sum_x2 + Matrix::rowSums(x2_batch)
  }

  list(n = n_total, sum_x = sum_x, sum_x2 = sum_x2)
}


# =============================================================================
#  Derive mean and sample variance from accumulated statistics
# =============================================================================
derive_mean_variance <- function(stats) {
  n      <- stats$n
  mean_x <- stats$sum_x / n
  # Sample variance: (sum_x2 - n * mean^2) / (n - 1)
  var_x  <- pmax((stats$sum_x2 - n * mean_x^2) / (n - 1L), 0)
  list(mean = mean_x, variance = var_x)
}


# =============================================================================
#  VST: fit loess, compute standardised variance, rank genes
# =============================================================================
#
# gene_stats : list(mean, variance) as returned by derive_mean_variance()
# n_hvg      : number of highly variable genes to return
# loess_span : smoothing span for the loess fit (Seurat default 0.3)
#
# Returns a data.frame with one row per gene, columns:
#   gene                  : gene name
#   mean                  : global per-gene mean raw expression
#   variance              : global per-gene sample variance
#   variance.expected     : loess-predicted variance from the mean-variance fit
#   variance.standardized : observed / expected variance (ranking criterion)
#   highly.variable       : logical, TRUE for the top n_hvg genes
#
vst_select_hvg <- function(gene_stats, n_hvg = 2000, loess_span = 0.3) {
  mean_x <- gene_stats$mean
  var_x  <- gene_stats$variance
  genes  <- names(mean_x)

  # Only fit loess on genes with positive variance and positive mean
  fit_idx <- which(var_x > 0 & mean_x > 0)
  if (length(fit_idx) < 2L)
    stop("Fewer than 2 genes with positive variance and positive mean; ",
         "cannot fit loess.")

  loess_fit <- loess(
    log10(var_x[fit_idx]) ~ log10(mean_x[fit_idx]),
    span = loess_span
  )

  var_expected <- rep(NA_real_, length(genes))
  names(var_expected) <- genes
  var_expected[fit_idx] <- 10^loess_fit$fitted

  # Standardised variance: observed / expected
  # Genes excluded from the fit get NA — never selected as HVGs.
  var_std <- var_x / var_expected

  n_eligible <- sum(!is.na(var_std) & var_std > 0)
  n_keep     <- min(n_hvg, n_eligible)
  hvg_rank   <- rank(-var_std, na.last = "keep", ties.method = "first")
  is_hvg     <- !is.na(hvg_rank) & hvg_rank <= n_keep

  data.frame(
    gene                  = genes,
    mean                  = mean_x,
    variance              = var_x,
    variance.expected     = var_expected,
    variance.standardized = var_std,
    highly.variable       = is_hvg,
    row.names             = genes,
    stringsAsFactors      = FALSE
  )
}


# =============================================================================
#  Top-level wrapper
# =============================================================================
#
# batch_files : character vector of .rds paths, each a sparse genes x cells
#               RAW COUNTS matrix.
# genes       : optional character vector of genes to consider for HVG
#               selection. Default NULL uses all rows in the batch files.
#               Pass a pre-filtered list to exclude negative probes, etc.
# n_hvg       : number of highly variable genes to return (default 2000).
# loess_span  : loess smoothing parameter (default 0.3, matching Seurat).
#
# Returns a data.frame (one row per gene) with columns:
#   gene, mean, variance, variance.expected, variance.standardized,
#   highly.variable
#
# The character vector of selected HVG names is:
#   result$gene[result$highly.variable]
#
batched_find_hvg <- function(batch_files,
                             genes      = NULL,
                             n_hvg      = 2000,
                             loess_span = 0.3) {

  message("=== Pass 1: accumulating per-gene mean and variance ===")
  raw_stats  <- batched_hvg_summary_stats(batch_files, genes)

  message("=== Deriving mean and variance ===")
  gene_stats <- derive_mean_variance(raw_stats)

  message("=== Fitting loess and selecting top ", n_hvg, " HVGs ===")
  vst_select_hvg(gene_stats, n_hvg = n_hvg, loess_span = loess_span)
}
