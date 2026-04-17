# =============================================================================
# Pipeline Configuration
# =============================================================================
#
# All user-facing parameters live here. Every step script sources this file
# at startup. Run an individual step with:
#   Rscript 01_extract_counts.R
# or run the full pipeline with:
#   Rscript run_pipeline.R
# =============================================================================

# ── Pipeline control ──────────────────────────────────────────────────────────

# Annotation column names to force-rebuild (ignores folder-sentinel caching).
# Phases 1–3 (stats, hvg, pca, clustering, context plots) are never affected by
# force_rerun; delete the relevant folder manually to rebuild those.
#
# Examples:
#   force_rerun <- c("annotations_initial")
#   force_rerun <- c("annotations_initial", "annotations_merged")
force_rerun  <- c()

# ── Study ─────────────────────────────────────────────────────────────────────

study_name   <- "TMAs_16_17_18_v2"

# Paths (relative to pipeline/ directory)
base_dir     <- ".."                                           # project root
seurat_dir   <- file.path(base_dir, "rawdata")
out_dir      <- file.path(base_dir, "outputs",  study_name)

# Regex to match Seurat .rds files in seurat_dir (case-insensitive)
file_pattern <- "\\.RDS$"

# ── Step 01: Count extraction ─────────────────────────────────────────────────
assay_name   <- "RNA"
# (no extra parameters — file discovery uses seurat_dir and file_pattern above)

# ── Step 02: HVG selection ────────────────────────────────────────────────────

hvg_n          <- 2000     # number of highly variable genes to select
hvg_loess_span <- 0.3      # LOESS span for mean-variance fit (Seurat default)

# ── Step 03: PCA ──────────────────────────────────────────────────────────────

pca_method   <- "pearson"  # "pearson" (default) or "covariance"
npcs         <- 50         # number of principal components
scale_max    <- 10         # clip residuals / expression at this many SDs above mean
chunk_size   <- 5000       # cells per chunk in Pass 3 (tune to available RAM)

# Pearson PCA only:
phi          <- 1.01       # quasi-Poisson overdispersion factor

# Covariance PCA only (when pca_method = "covariance"):
apply_log1p  <- TRUE       # log1p after library-size normalisation

# ── Step 04: Clustering ───────────────────────────────────────────────────────

n_neighbors        <- 30
metric             <- "cosine"
min_dist           <- 0.01
nn_method          <- "annoy"
umap_seed          <- 42               # random seed for reproducible UMAP
cluster_algorithm  <- "louvain"        # "louvain" (default) or "leiden"
resolutions        <- c(1.2, 0.8, 1.0)

# Which resolution drives all downstream annotation steps.
# Must be a member of `resolutions` above.
# Switching this value (within pre-computed resolutions) requires no
# recomputation — all heavy outputs are shared at the top level.
active_resolution  <- resolutions[1]

# ── Step 06: LLM prompt / annotation ─────────────────────────────────────────

tissue_type         <- "FFPE human Breast Cancer (6k plex cosmx data)"   # description used in the LLM prompt

# ── Step 05: Marker genes ─────────────────────────────────────────────────────

n_top_markers      <- 10    # top N marker genes per cluster in CSV summary
ncores_markers     <- 1     # parallel cores for Welch t-test (mclapply)
marker_ranking     <- "priority"  # "priority" (default) or "logfc"
# priority: expression ratio vs other clusters; rewards specific + detected
# logfc:    avg_log2FC only; may surface lowly-detected genes

# ── Step 07: Annotate Seurat objects ─────────────────────────────────────────

pca_slot   <- "study_pca"    # DimReduc slot name written into Seurat objects
umap_slot  <- "study_umap"

# ── Step 08: Visualization ────────────────────────────────────────────────────

# Set TRUE to regenerate context plots (plots/) even if the directory exists.
# Annotation plots are controlled per-column by force_rerun.
force_overview_plots <- FALSE

scatter_sample_n <- 50000    # cells to sample for scatter / UMAP plots

# Spatial coordinate columns in Seurat meta.data (spatial section)
# "auto" = try known CosMX column patterns; or set explicitly, e.g. "x_slide_mm"
spatial_coord_x  <- "auto"
spatial_coord_y  <- "auto"

# Top marker genes per cluster to display as UMAP feature plots
feature_plot_n   <- 2

# Seurat meta.data column for patient/study grouping in composition plots
# Set to NULL to skip patient composition plots entirely.
study_id_col     <- "study_id"

# Seurat meta.data columns to extract from rawdata objects (QC violin plots,
# per-cluster means in step 09). Add annotation column names here to include
# them in spatial/QC plots.
meta_vis_vars <- c("Mean.PanCK", "Mean.CD45", "Mean.DAPI",
                   "nCount_RNA", "nFeature_RNA", "Area.um2")

# ── Step 10: Cell manifest export ─────────────────────────────────────────────

# Single flat file combining all per-cell non-gene-expression info.
# Useful for interactive visualization (Shiny, etc.).
# Output: {out_dir}/cell_manifest_res{X.X}.rds  (+ optional .csv)
force_manifest         <- TRUE   # TRUE to regenerate even if file already exists
export_manifest_csv    <- TRUE   # TRUE to also write .csv (large at scale)
export_manifest_seurat <- TRUE   # TRUE to also write a Seurat RDS (metadata only,
                                  #   no assay counts; UMAP injected as DimReduc)
