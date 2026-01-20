# These are helper functions for CosMX data analysis
# import them to your analysis using source("helper_functions.R")
# Author: Frederik Stihler

## IO UTILS

# Function: Load all Seurat objects from a folder into a named list
# Args:
#   folder: directory containing .RDS Seurat objects
# Returns:
#   Named list of Seurat objects
load_seurat_objects <- function(folder, assays_to_remove = c("QC_Normalization.RNA.1_1")) {
  files <- list.files(folder, pattern = "^seuratObject_.*\\.RDS$", full.names = TRUE)
  print(files)
  # Extract slide names from filenames
  slide_names <- sub("^seuratObject_(.*)\\.RDS$", "\\1", basename(files))
  print(slide_names)
  
  # Load into list
  seurat_list <- setNames(lapply(files, readRDS), slide_names)
  
  # Remove specified assays
  seurat_list <- lapply(seurat_list, function(obj) {
    for (assay_name in assays_to_remove) {
      print(paste("Removing assay:", assay_name))
      if (assay_name %in% names(obj@assays)) {
        obj[[assay_name]] <- NULL
      }
    }
    return(obj)
  })

  return(seurat_list)
}

# Function: Merge all Seurat objects with slide-specific cell ID prefixes
# Args:
#   seurat_list: named list of Seurat objects
#   project_name: name for merged object
# Returns:
#   Merged Seurat object
merge_seurat_objects <- function(seurat_list, project_name = "CosMx") {
  merged <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = project_name
  )
  return(merged)
}

## ANNOTATIONS

make_casewhen <- function(borders) {
  borders <- sort(borders)
  starts <- c(1, borders[-length(borders)] + 1)
  ends <- borders
  regions <- seq_along(borders)
  
  lines <- paste0(
    "fov %in% ", starts, ":", ends, " ~ \"", regions, "\""
  )
  
  case_block <- paste(lines, collapse = ",\n    ")
  full_code <- paste0(
    "case_when(\n    ",
    case_block,
    ",\n    TRUE ~ \"Unknown\"\n  )"
  )
  
  cat(full_code)
  invisible(full_code)
}

process_tma_csv <- function(file_path, sep = ",", skip_rows = 1) {
  # Load CSV
  metadata_csv <- read.csv(file_path, stringsAsFactors = FALSE, sep = sep, skip = skip_rows)
  
  # Rename columns
  colnames(metadata_csv) <- c(
    "study_id",
    "accession_id",
    "accession_year",
    "part_id",
    "block_id",
    "tma_number",
    "tma_map_tumor",
    "tma_map_normal",
    "batch"
  )
  
  # Keep only relevant columns
  metadata_csv <- metadata_csv[, c("study_id", "tma_number", "tma_map_tumor", "tma_map_normal", "batch")]

  # Remove leading and trailing whitespace
  metadata_csv <- metadata_csv %>%
    mutate(across(everything(), ~trimws(.)))
  
  # Convert to long format
  metadata_long <- metadata_csv %>%
    mutate(
      tma_map_tumor = ifelse(tma_map_tumor == "No Tumor", NA, tma_map_tumor),
      tma_map_normal = ifelse(tma_map_normal == "No Normal", NA, tma_map_normal)
    ) %>%
    pivot_longer(
      cols = c(tma_map_tumor, tma_map_normal),
      names_to = "type",
      values_to = "location"
    ) %>%
    filter(!is.na(location)) %>%
    mutate(type = ifelse(type == "tma_map_tumor", "T", "N")) %>%
    separate_rows(location, sep = ",\\s*") %>%
    select(study_id, tma_number, location, batch, type)
  
  return(metadata_long)
}


## VISUALIZATION

xyplot <- function(cluster_column,
                   x_column = "x_slide_mm",
                   y_column = "y_slide_mm",
                   cls = NULL,
                   clusters = NULL,
                   metadata,
                   ptsize = 0.25,
                   plotfirst = NULL,
                   plotfirst_on_top = FALSE,
                   alphasize = 1,
                   show_legend = TRUE,
                   coord_equal = TRUE,
                   continuous_palette = function(n) viridis::viridis(n, option = "plasma"),
                   aes_mappings = list(size = NULL, shape = NULL, alpha = NULL),
                   order = NULL,
                   na_color = "black",
                   theme = ggplot2::theme_bw(),
                   show_labels = FALSE,         
                   label_size = 4,               
                   label_color = "black",        
                   label_fontface = "bold",     
                   label_method = "median",      
                   label_repel = TRUE,           
                   label_max_overlaps = 50       
) {

  pd <- data.table::copy(data.table::data.table(metadata))

  if (!is.null(order)) {
    order <- as.character(order)
  } else if (is.factor(pd[[cluster_column]])) {
    order <- levels(pd[[cluster_column]])
  }

  if (is.null(clusters)) clusters <- unique(pd[[cluster_column]])
  clusters <- unique(clusters)

  mask <- pd[[cluster_column]] %in% clusters
  mask[is.na(mask)] <- FALSE
  if (any(is.na(pd[[cluster_column]]))) mask <- mask | is.na(pd[[cluster_column]])
  pd_plot <- pd[mask, , drop = FALSE]

  if (!is.null(order)) {
    present_levels <- order
    pd_plot[[cluster_column]] <- factor(as.character(pd_plot[[cluster_column]]), levels = present_levels)
    clusters_plot <- present_levels
  } else {
    if (all(!is.na(suppressWarnings(as.numeric(as.character(clusters)))))) {
      clusters_plot <- sort(as.numeric(as.character(clusters)))
      clusters_plot <- as.character(clusters_plot)
    } else {
      clusters_plot <- sort(as.character(clusters))
    }
    if (any(is.na(pd_plot[[cluster_column]]))) clusters_plot <- unique(c(clusters_plot, NA))
  }

  mapping_args <- list(x = x_column, y = y_column, color = cluster_column)
  if (!is.null(aes_mappings$size)) mapping_args$size <- aes_mappings$size
  if (!is.null(aes_mappings$shape)) mapping_args$shape <- aes_mappings$shape
  if (!is.null(aes_mappings$alpha)) mapping_args$alpha <- aes_mappings$alpha
  mapping <- do.call(ggplot2::aes_string, mapping_args)

  size_is_mapped <- !is.null(aes_mappings$size)

  if (!is.null(plotfirst)) {
    plotfirst <- intersect(clusters_plot, plotfirst)
    notplotfirst <- setdiff(clusters_plot, plotfirst)

    df_first <- pd_plot[as.character(pd_plot[[cluster_column]]) %in% as.character(plotfirst), , drop = FALSE]
    df_rest  <- pd_plot[as.character(pd_plot[[cluster_column]]) %in% as.character(notplotfirst), , drop = FALSE]

    layer_first <- ggplot2::geom_point(
      data = df_first, mapping = mapping,
      size = if (!size_is_mapped) ptsize else NULL,
      alpha = if (is.null(aes_mappings$alpha)) 1 else NULL
    )
    layer_rest <- ggplot2::geom_point(
      data = df_rest, mapping = mapping,
      size = if (!size_is_mapped) ptsize else NULL,
      alpha = if (is.null(aes_mappings$alpha)) alphasize else NULL
    )

    p <- ggplot2::ggplot() + theme
    if (plotfirst_on_top) {
      p <- p + layer_rest + layer_first
    } else {
      p <- p + layer_first + layer_rest
    }
  } else {
    p <- ggplot2::ggplot(pd_plot, mapping = mapping) +
      ggplot2::geom_point(
        size = if (!size_is_mapped) ptsize else NULL,
        alpha = if (is.null(aes_mappings$alpha)) alphasize else NULL
      ) +
      theme
  }

  if (is.numeric(pd[[cluster_column]])) {
    cols <- try(continuous_palette(256), silent = TRUE)
    if (inherits(cols, "try-error") || !is.character(cols)) {
      p <- p + ggplot2::scale_color_viridis_c(
        guide = if (show_legend) ggplot2::guide_colorbar() else "none",
        na.value = na_color
      )
    } else {
      p <- p + ggplot2::scale_color_gradientn(
        colors = cols,
        guide = if (show_legend) ggplot2::guide_colorbar() else "none",
        na.value = na_color
      )
    }
  } else {
    n_clusters <- length(clusters_plot)
    if (is.null(cls)) {
      base_pal <- unname(pals::alphabet())
      if (n_clusters <= length(base_pal)) {
        palette_used <- base_pal[seq_len(n_clusters)]
      } else {
        palette_used <- grDevices::colorRampPalette(base_pal)(n_clusters)
      }
      names(palette_used) <- as.character(clusters_plot)
      cls_final <- palette_used
    } else {
      if (is.character(cls) && length(cls) < n_clusters) {
        cls_final <- grDevices::colorRampPalette(cls)(n_clusters)
        names(cls_final) <- as.character(clusters_plot)
      } else {
        if (is.null(names(cls))) {
          cls_final <- cls[seq_len(min(length(cls), n_clusters))]
          names(cls_final) <- as.character(clusters_plot)[seq_len(length(cls_final))]
          if (length(cls) < n_clusters) {
            extra <- grDevices::colorRampPalette(cls)(n_clusters - length(cls))
            cls_final <- c(
              cls_final,
              stats::setNames(extra, as.character(clusters_plot)[(length(cls) + 1):n_clusters])
            )
          }
        } else {
          cls_final <- cls[as.character(clusters_plot)]
          missing_idx <- which(is.na(cls_final))
          if (length(missing_idx) > 0) {
            fill_colors <- grDevices::colorRampPalette(unname(cls))(length(missing_idx))
            cls_final[missing_idx] <- fill_colors
          }
        }
      }
    }

    p <- p + ggplot2::scale_color_manual(
      values = cls_final,
      guide = if (show_legend)
        ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1))
      else "none",
      na.value = na_color
    )
  }

  if (isTRUE(coord_equal)) p <- p + ggplot2::coord_fixed()

  # New label layer
  if (isTRUE(show_labels) && !is.numeric(pd_plot[[cluster_column]])) {
    summary_fun <- switch(label_method, "mean" = mean, "median" = median, median)
    centroids <- pd_plot[, .(
      x = summary_fun(get(x_column), na.rm = TRUE),
      y = summary_fun(get(y_column), na.rm = TRUE)
    ), by = cluster_column]

    if (requireNamespace("ggrepel", quietly = TRUE) && isTRUE(label_repel)) {
      p <- p + ggrepel::geom_label_repel(
        data = centroids,
        aes(x = x, y = y, label = .data[[cluster_column]]),
        inherit.aes = FALSE,
        color = label_color,
        size = label_size,
        fontface = label_fontface,
        show.legend = FALSE,
        max.overlaps = label_max_overlaps
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = centroids,
        aes(x = x, y = y, label = .data[[cluster_column]]),
        inherit.aes = FALSE,
        color = label_color,
        size = label_size,
        fontface = label_fontface,
        show.legend = FALSE
      )
    }
  }

  p
}

plot_embedding <- function(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label_size = 3.5,
  alpha = 0.6,
  point_size = 0.3,
  raster = TRUE,
  shuffle = TRUE,
  seed = 123,
  palette = NULL,
  na_color = "grey80",
  palette_continuous = function(n) viridis::viridis(n, option = "plasma"),
  legend = TRUE
) {
  # Validate inputs
  if (!reduction %in% names(seu@reductions))
    stop("Reduction not found: ", reduction)
  if (!group.by %in% colnames(seu@meta.data))
    stop("Column not found in meta.data: ", group.by)

  # Extract embeddings and metadata
  emb <- as.data.frame(Embeddings(seu, reduction))
  colnames(emb)[1:2] <- c("Dim1", "Dim2")
  emb$group <- seu@meta.data[[group.by]]

  # Optional shuffle for even plotting
  if (shuffle) {
    set.seed(seed)
    emb <- emb[sample(nrow(emb)), , drop = FALSE]
  }

  # Define plotting geometry
  geom_fun <- if (raster) ggrastr::geom_point_rast else geom_point

  # Handle numeric vs categorical
  if (is.numeric(emb$group)) {
    # Continuous variable (e.g., expression, score)
    p <- ggplot(emb, aes(Dim1, Dim2, color = group)) +
      geom_fun(size = point_size, alpha = alpha) +
      scale_color_gradientn(
        colors = palette_continuous(256),
        na.value = na_color
      ) +
      labs(color = group.by)
  } else {
    emb$group <- as.factor(emb$group)
    levels_vec <- levels(emb$group)
    # palette logic
    if (is.null(palette)) {
      pal_vec <- unname(pals::alphabet())
      pal_use <- setNames(pal_vec[seq_along(levels_vec)], levels_vec)
    } else {
      # Check if the user provided a named vector
      if (!is.null(names(palette))) {
        # Filter and reorder the palette to match levels exactly
        # This ensures cluster "a" always gets the color assigned to "a"
        pal_use <- palette[levels_vec]
        
        # Handle cases where some levels might be missing from the palette
        pal_use[is.na(pal_use)] <- na_color
      } else {
        # Fallback for unnamed palettes: match by position
        pal_vec <- palette
        if (length(pal_vec) < length(levels_vec))
          pal_vec <- rep(pal_vec, length.out = length(levels_vec))
        pal_use <- setNames(pal_vec[seq_along(levels_vec)], levels_vec)
      }
    }

    p <- ggplot(emb, aes(Dim1, Dim2, color = group)) +
      geom_fun(size = point_size, alpha = alpha) +
      scale_color_manual(values = pal_use, na.value = na_color) +
      labs(color = group.by)

    # optional labeling
    if (label) {
      centers <- aggregate(cbind(Dim1, Dim2) ~ group, data = emb, FUN = median)
      p <- p +
        ggrepel::geom_text_repel(
          data = centers,
          aes(x = Dim1, y = Dim2, label = group),
          size = label_size,
          color = "black",
          fontface = "bold",
          max.overlaps = Inf,
          show.legend = FALSE
        )
    }
  }

  # styling
  p <- p +
    theme_bw() +
    coord_fixed() +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = if (legend) "right" else "none"
    )

  return(p)
}

umap_plots <- function( seu,
                        slide_name,
                        outdir = ".",
                        file_prefix = "Results_Report_",
                        reduction = "umapharmony",
                        cluster_col = "clusters_unsup_harmony",
                        slide_id_col = "slide_id",
                        region_col = "region",
                        condition_col = "condition",
                        slide_cls = brewer.pal(8, "Set2"),
                        region_cls = NULL,
                        condition_cls = c("T" = "red2", "N" = "thistle"),
                        shuffle = TRUE) {

    pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")),
        width = 12, height = 10)
    
    md_cols <- colnames(seu@meta.data)

    # Normalize cluster_col to character vector
    if (is.list(cluster_col)) {
        cluster_cols <- unlist(cluster_col)
    } else {
        cluster_cols <- as.character(cluster_col)
    }

    # DimPlots cluster columns
    for (cluster in cluster_cols) {
        if (cluster %in% md_cols) {
            p <- plot_embedding(seu, reduction = reduction, group.by = cluster, shuffle = shuffle) + 
                ggtitle(paste0("UMAP: ", reduction, ": ", cluster))
            print(p)
        }
        else {
            message("Skipping DimPlot for missing metadata column: ", cluster)
        }
    }

    # DimPlots standard columns
    if (slide_id_col %in% md_cols) {
        p <- plot_embedding(seu, reduction = reduction, group.by = slide_id_col, shuffle = shuffle, label = FALSE, palette = slide_cls) + 
            ggtitle(paste0("UMAP: ", reduction, ": ", slide_id_col))
        print(p)
    }

    if (region_col %in% md_cols) {
        p <- plot_embedding(seu, reduction = reduction, group.by = region_col, shuffle = shuffle, label = FALSE, legend = FALSE, palette = region_cls) + 
            ggtitle(paste0("UMAP: ", reduction, ": ", region_col))
        print(p)
    }
    if (condition_col %in% md_cols) {
        seu@meta.data[[condition_col]] <- factor(seu@meta.data[[condition_col]],
                                         levels = names(condition_cls))
        p <- plot_embedding(seu, reduction = reduction, group.by = condition_col, shuffle = shuffle, label = FALSE, palette = condition_cls) + 
            ggtitle(paste0("UMAP: ", reduction, ": ", condition_col))
        print(p)
    }

    dev.off()
    message("Summary plots saved for slide: ", slide_name)
}

  # Function to get cluster colors
  get_cluster_colors <- function(clusters, palette_spec = NULL) {
    n_clusters <- length(clusters)
    base_pal <- unname(pals::alphabet())

    if (is.null(palette_spec)) {
      # default alphabet with extension if needed
      if (n_clusters <= length(base_pal)) {
        cols <- base_pal[seq_len(n_clusters)]
      } else {
        cols <- grDevices::colorRampPalette(base_pal)(n_clusters)
      }
    } else if (is.character(palette_spec)) {
      # user-supplied palette vector
      if (length(palette_spec) < n_clusters) {
        cols <- c(
          palette_spec,
          grDevices::colorRampPalette(palette_spec)(n_clusters - length(palette_spec))
        )
      } else {
        cols <- palette_spec[seq_len(n_clusters)]
      }
    } else if (is.function(palette_spec)) {
      cols <- palette_spec(n_clusters)
    } else {
      cols <- grDevices::rainbow(n_clusters)
    }

    stats::setNames(cols, as.character(clusters))
  }

xy_plots <- function(seu,
                          slide_name,
                          outdir = ".",
                          file_prefix = "Spatial_Plots_",
                          cluster_col = "clusters_unsup_harmony",
                          shuffle = TRUE,
                          palettes = NULL,
                          show_legend = TRUE) {

  # Ensure cluster_col is a character vector
  cluster_cols <- as.character(unlist(cluster_col))

  # ensure palettes align with cluster_cols
  if (!is.null(palettes)) {
    if (length(palettes) != length(cluster_cols)) {
      stop("Length of 'palettes' does not match number of cluster columns. Stopping...")
    }
  }

  pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")),
      width = 12, height = 10)

  md_cols <- colnames(seu@meta.data)

  # ------------------------
  # Full tissue plots (one page per clustering)
  # ------------------------
    for (i in seq_along(cluster_cols)) {
    grp <- cluster_cols[i]
    palette_spec <- if (!is.null(palettes)) palettes[[i]] else NULL

    if (all(c("x_slide_mm", "y_slide_mm") %in% md_cols) && grp %in% md_cols) {
      all_clusters <- sort(unique(seu@meta.data[[grp]]))

      # Use global get_cluster_colors function
      cols <- get_cluster_colors(all_clusters, palette_spec)

      print(
        xyplot(
          grp,
          metadata = seu@meta.data,
          ptsize = 0.01,
          cls = cols,
          show_legend = show_legend
        ) +
          coord_fixed() +
          labs(title = paste("SpatialPlot:", grp))
      )
    } else {
      message("Skipping FeatureScatter for ", grp,
              ": required columns missing (x_slide_mm, y_slide_mm, or cluster_col)")
    }
  }

  dev.off()
  message("Summary plots saved for slide: ", slide_name)
}

xy_plots_by_region <- function(
  seu,
  slide_name,
  outdir = ".",
  file_prefix = "XY_Zoom_",
  cluster_col = "clusters_unsup_harmony",
  region_col = "region",
  regions_per_page = 2,
  shuffle = TRUE,
  palettes = NULL,
  show_legend = TRUE
) {
  # Ensure cluster_col is a character vector
  cluster_cols <- as.character(unlist(cluster_col))

  # ensure palettes align with cluster_cols
  if (!is.null(palettes) && length(palettes) != length(cluster_cols)) {
    stop("Length of 'palettes' must match number of cluster columns")
  }
  
  pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")), width = 12, height = 10)

  md_cols <- colnames(seu@meta.data)
  if (!region_col %in% md_cols)
    stop("Specified region_col '", region_col, "' not found in metadata.")
  
  region_vals <- seu@meta.data[[region_col]]
  if (is.factor(region_vals)) {
    regions <- levels(region_vals)
  } else {
    regions <- mixedsort(unique(region_vals))
  }

  # helper: plot a single region for a given clustering column
  plot_region_cluster <- function(obj_region, clust, palette_spec) {
    region_md_cols <- colnames(obj_region@meta.data)
    if (!all(c("x_slide_mm", "y_slide_mm") %in% region_md_cols) || !(clust %in% region_md_cols)) {
      return(ggplot() + theme_void() + labs(title = paste(clust, "missing for region")))
    }

    clusters_plot <- sort(unique(obj_region@meta.data[[clust]]))
    cols <- get_cluster_colors(clusters_plot, palette_spec)

    xyplot(
      clust,
      metadata = obj_region@meta.data,
      ptsize = 0.01,
      cls = cols,
      show_legend = show_legend
    ) +
      coord_fixed() +
      labs(title = paste(clust, "\nRegion", unique(obj_region@meta.data[[region_col]])))
  }

  if (length(cluster_cols) == 1) {
    # single cluster column: possibly multiple regions per page
    clust <- cluster_cols[1]
    palette_spec <- if (!is.null(palettes)) palettes[[1]] else NULL

    region_plots <- lapply(regions, function(r) {
      obj_region <- subset(seu, cells = rownames(seu@meta.data[seu@meta.data[[region_col]] == r, ]))
      plot_region_cluster(obj_region, clust, palette_spec)
    })

    # Remove NULLs
    region_plots <- Filter(Negate(is.null), region_plots)

    # arrange multiple regions per page
    for (i in seq(1, length(region_plots), by = regions_per_page)) {
      idx <- i:min(i + regions_per_page - 1, length(region_plots))
      print(wrap_plots(region_plots[idx], ncol = 2))
    }

  } else {
    # multiple cluster columns: one region per page, all clusters on same page
    for (r in regions) {
      cells <- rownames(seu@meta.data[seu@meta.data[[region_col]] == r, ])
      if (length(cells) == 0) next
      obj_region <- subset(seu, cells = cells)

      # if palettes is NULL, create a list of NULLs matching cluster columns
      region_palettes <- palettes
      if (is.null(region_palettes)) {
        region_palettes <- vector("list", length(cluster_cols))
      }

      region_plots <- mapply(function(clust, pal) {
        plot_region_cluster(obj_region, clust, pal)
      }, cluster_cols, region_palettes, SIMPLIFY = FALSE)

      # Remove NULLs if any
      region_plots <- Filter(Negate(is.null), region_plots)
      if (length(region_plots) == 0) next

      print(wrap_plots(region_plots, ncol = 2))
    }
  }

  dev.off()
  message("Summary plots saved for slide: ", slide_name)
}

spatial_plots <- function(seu,
                          slide_name,
                          outdir = ".",
                          file_prefix = "XY_Zoom_",
                          cluster_col = "clusters_unsup_harmony",
                          shuffle = TRUE,
                          show_legend = TRUE) {

  # Normalize cluster_col to character vector
  if (is.list(cluster_col)) {
    cluster_cols <- unlist(cluster_col)
  } else {
    cluster_cols <- as.character(cluster_col)
  }

  pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")),
      width = 12, height = 10)

  md_cols <- colnames(seu@meta.data)
  regions <- mixedsort(unique(seu$region))

  # ------------------------
  # 2) Per-region plots (one page per region, all clusterings side by side)
  # ------------------------
  for (r in regions) {
    obj_region <- subset(seu, subset = region == r)
    region_md_cols <- colnames(obj_region@meta.data)

    region_plots <- lapply(cluster_cols, function(clust) {

      if (!all(c("x_slide_mm", "y_slide_mm") %in% region_md_cols) || !(clust %in% region_md_cols)) {
        return(ggplot() + theme_void() +
                 labs(title = paste(clust, "missing for region", r)))
      }

      all_clusters <- sort(unique(obj_region@meta.data[[clust]]))
      cols <- setNames(alphabet()[1:length(all_clusters)], all_clusters)

      xyplot(clust,
             metadata = obj_region@meta.data,
             ptsize = 0.01,
             cls = cols,
             show_legend = show_legend) +
        coord_fixed() +
        labs(title = paste(clust, "\n- Region", r))
    })

    # Combine all clusterings side by side for this region
    combined_region <- wrap_plots(region_plots, ncol = 2)
    print(combined_region)
  }

  dev.off()
  message("Summary plots saved for slide: ", slide_name)
}


## QC

summarize_qc_by <- function(metadata, group_col,
                            qc_vars = c("nCount_RNA", "nFeature_RNA", "Area.um2", "nCell", "nCount_negprobes", "nCount_falsecode"),
                                        digits = 0) {
  stopifnot(group_col %in% colnames(metadata))
  
  metadata %>%
    group_by(.data[[group_col]]) %>%
    summarize(across(all_of(qc_vars),
                     list(
                       mean = ~mean(.x, na.rm = TRUE),
                       median = ~median(.x, na.rm = TRUE),
                       min = ~min(.x, na.rm = TRUE),
                       max = ~max(.x, na.rm = TRUE),
                       range = ~diff(range(.x, na.rm = TRUE)),
                       sd = ~sd(.x, na.rm = TRUE),
                       variance = ~var(.x, na.rm = TRUE),
                       IQR = ~IQR(.x, na.rm = TRUE),
                       count = ~sum(!is.na(.x)),
                       total = ~sum(.x, na.rm = TRUE)
                     ),
                     .names = "{.col}_{.fn}"
    ),
    .groups = "drop") %>%
    pivot_longer(
      cols = -all_of(group_col),
      names_to = c("Variable", "Statistic"),
      names_sep = "_(?=[^_]+$)",  # split only at the LAST underscore
      values_to = "Value"
    ) %>%
    dplyr::mutate(Value = round(Value, digits))
}


calculate_qc_summaries <- function(metadata,
                                   slide_col = "slide_name",
                                   region_col = "region",
                                   fov_col = "fov",
                                   qc_vars = c("nCount_RNA", "nFeature_RNA", "Area.um2", "nCell", "nCount_negprobes", "nCount_falsecode"),
                                   digits = 2) {
  stats_list <- list(
    slide_level  = summarize_qc_by(metadata, slide_col, qc_vars, digits),
    region_level = summarize_qc_by(metadata, region_col, qc_vars, digits),
    fov_level    = summarize_qc_by(metadata, fov_col, qc_vars, digits)
  )

  # # compute SNR per group
  # if ("nCount_negprobes" %in% qc_vars && "nCount_RNA" %in% qc_vars) {
    
  #   stats_list <- lapply(stats_list, function(df) {
  #     df %>%
  #       group_by(.data[[names(df)[1]]]) %>%
  #       mutate(SNR = Value[Variable == "nCount_RNA" & Statistic == "mean"] /
  #                   Value[Variable == "nCount_negprobes" & Statistic == "mean"])
  #   })  
  # }
  return(stats_list)
}


create_qc_plots <- function(seu, batch_col = "slide_id", fov_col = "fov", 
                             count_col = "nCount_RNA", feature_col = "nFeature_RNA",
                             neg_col = "nCount_negprobes", fc_col = "nCount_falsecode",
                             pal = NULL) {

    md <- seu@meta.data
  
    # Determine number of negative and system control probes
    nNeg <- length(unique(seu[["negprobes"]]@counts@Dimnames[[1]]))  # total negative probes
    nFC  <- length(unique(seu[["falsecode"]]@counts@Dimnames[[1]]))  # total system control probes
    
    # Compute FOV-level QC summary
    qc_stats <- md %>%
        group_by(.data[[batch_col]], .data[[fov_col]]) %>%
        summarise(
        Number_Cells_Per_FOV = mean(nCell),
        Mean_Transcripts_Per_Cell_Per_FOV = mean(.data[[count_col]]),
        Mean_Unique_Transcripts_Per_Cell_Per_FOV = mean(.data[[feature_col]]),
        Total_Transcripts_Per_FOV = sum(.data[[count_col]]),
        Mean_Negative_Probe_Per_Plex_Per_Cell_Per_FOV = sum(.data[[neg_col]]) / (unique(nCell) * nNeg) ,
        Mean_SystemControl_Per_Plex_Per_Cell_Per_FOV = sum(.data[[fc_col]]) / (unique(nCell) * nFC) ,
        .groups = "drop"
        )
    
    # Default palette if not provided
    if (is.null(pal)) {
        n_batches <- length(unique(md[[batch_col]]))
        pal <- scales::hue_pal()(n_batches)
    }
    
    # Metrics to plot
    plot_metrics <- c(
        "Mean_Transcripts_Per_Cell_Per_FOV",
        "Mean_Unique_Transcripts_Per_Cell_Per_FOV",
        "Total_Transcripts_Per_FOV",
        "Number_Cells_Per_FOV",
        "Mean_Negative_Probe_Per_Plex_Per_Cell_Per_FOV",
        "Mean_SystemControl_Per_Plex_Per_Cell_Per_FOV"
    )
    
    # Loop over metrics and plot
    for (metric in plot_metrics) {
        # skip if metric not present
        if (!metric %in% colnames(qc_stats)) next
        
        p <- ggplot(qc_stats, aes_string(x = batch_col, y = metric, fill = batch_col)) +
        geom_violin(alpha = 0.3) +
        geom_jitter(width = 0.25, height = 0, size = 0.5) +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_fill_manual(values = pal) +
        ylab(metric) +
        xlab("Batch / Slide") +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5)
        )
        
        print(p)
    }
    }

# Function: QC analysis and PDF report per Seurat object
# Args:
#   seu: Seurat object
#   slide_name: name of the slide for labeling outputs
#   outdir: directory to save PDF reports
# Returns:
#   Nothing, saves PDF report
qc_report <- function(seu, slide_name, outdir = ".", file_prefix = "QC_Report_", filter_log = NULL, count_column = "nCount_RNA", region_col = "region", flag_col = "qcCellsFlagged", boxplots = TRUE, fov_plots = TRUE) {
  pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")), width = 14, height = 10)

  # Tissue plot from metadata
  md <- seu@meta.data

  # If filter_log is provided, print the row for this slide
  if (!is.null(filter_log) && slide_name %in% filter_log$Slide) {
    slide_log <- filter_log[filter_log$Slide == slide_name, ]
    log_text <- paste(
      paste(names(slide_log), slide_log, sep = ": "), collapse = "\n"
    )
    grid::grid.newpage()
    grid::grid.text(
      log_text,
      x = 0.5, y = 0.5,
      gp = grid::gpar(fontsize = 12), just = "center"
    )
  }

  if (all(c("x_slide_mm", "y_slide_mm") %in% colnames(md))) {
    tissue_plot <- xyplot(region_col, x_column = "x_slide_mm", y_column = "y_slide_mm", 
                      metadata = md, 
                      alphasize = 1, show_legend = FALSE, show_labels = TRUE,
                      label_color = "black",
                      label_repel = TRUE) +
                      ggtitle(paste("Tissue Layout -", slide_name))
    print(tissue_plot)
    region_colors <- ggplot2::ggplot_build(tissue_plot)$plot$scales$get_scales("colour")$palette.cache
  }

  if ("nCount_RNA" %in% colnames(md) & "nFeature_RNA" %in% colnames(md)) {
    md$log10_nCount_RNA <- log10(md$nCount_RNA + 1)
    md$log10_nFeature_RNA <- log10(md$nFeature_RNA + 1)
    
    tissue_plot_nCount <- xyplot("log10_nCount_RNA", x_column = "x_slide_mm", y_column = "y_slide_mm", 
                      metadata = md, 
                      alphasize = 1, show_legend = TRUE, show_labels = FALSE,
                      continuous_palette = function(n) viridis::viridis(n, option = "inferno")) +
                      ggtitle(paste("log10(nCountRNA + 1) -", slide_name))
    print(tissue_plot_nCount)

    # tissue_plot_nFeature <- xyplot("log10_nFeature_RNA", x_column = "x_slide_mm", y_column = "y_slide_mm", 
    #                   metadata = md, 
    #                   alphasize = 1, show_legend = TRUE, show_labels = FALSE,
    #                   continuous_palette = function(n) viridis::viridis(n, option = "inferno")) +
    #                   ggtitle(paste("log10(FeatureRNA + 1) -", slide_name))
    # print(tissue_plot_nFeature)
  }
  
  if (flag_col %in% colnames(md)) {
    flag_cls <- c("TRUE" = "red4", "FALSE" = "gray")

    tissue_plot_qc <- xyplot(flag_col, x_column = "x_slide_mm", y_column = "y_slide_mm", 
                   cls = flag_cls, metadata = md, ptsize = 0.25,
                   alphasize = 1, show_legend = TRUE, show_labels = FALSE) +
                   ggtitle(paste(flag_cls,": AtoMX QC Flags"))
    print(tissue_plot_qc)
  }

  # Tissue plot with condition
  condition_col <- "condition"

  if (condition_col %in% colnames(md)) {
    condition_cls <- c("T" = "red2", "N" = "thistle")

    tissue_plot_cond <- xyplot(condition_col, x_column = "x_slide_mm", y_column = "y_slide_mm", 
                   cls = condition_cls,metadata = md, ptsize = 0.25,
                   alphasize = 1, show_legend = TRUE, show_labels = FALSE) +
                   ggtitle(paste("Tissue plot ", condition_col))
    print(tissue_plot_cond)
  }
  
  # Violin plots with median bars
  vp1 <- VlnPlot(seu, features = count_column, pt.size = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
    ggtitle(paste(count_column," per Cell"))

  count_column_2 <- "nFeature_RNA"
  vp2 <- VlnPlot(seu, features = count_column_2, pt.size = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
    ggtitle(paste(count_column_2," per Cell"))
  
  print(vp1 + vp2 + plot_layout(ncol = 2) + plot_annotation(title = paste("QC Metrics Violin Plots -", slide_name)))

  # Summary stats 
  qc_stats <- calculate_qc_summaries(
      metadata = md,
      slide_col = "slide_id",
      region_col = region_col,
      fov_col = "fov",
      digits = 2
    )

  ## Slide Level
  slide_table <- qc_stats$slide_level %>%
    pivot_wider(names_from = Statistic, values_from = Value)

  grid::grid.newpage()
  grid::grid.text("Slide-Level QC Summary", y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))

  qc_table_grob <- gridExtra::tableGrob(
      slide_table,
      rows = NULL,
      theme = gridExtra::ttheme_default(
        core = list(fg_params = list(cex = 0.7)),
        colhead = list(fg_params = list(cex = 0.8, fontface = "bold"))
      )
    )

  grid::grid.draw(qc_table_grob)

  # Signal-to-noise ratio
  if ("negprobes" %in% names(seu)) {
    snr <- mean(Matrix::colMeans(seu[["RNA"]]@counts)) / mean(Matrix::colMeans(seu[["negprobes"]]@counts))
  } else {
    snr <- NA
  }

  # Flagged cells at the bottom
  if (flag_col %in% colnames(md)) {
    flagged_count <- sum(md[[flag_col]], na.rm = TRUE)
    total_cells <- ncol(seu)
    flagged_text <- paste0("SNR: ", round(snr, 2),
                          "\nTotal cells: ", total_cells,
                          "\nFlagged cells: ", flagged_count,
                          " (", round(100 * flagged_count / total_cells, 2), "%)")
    grid::grid.text(
      flagged_text,
      x = 0.5, y = 0.3,  # lower vertical position
      gp = grid::gpar(fontsize = 14)
    )
  }

  # Violin plots split by condition
  if (condition_col %in% colnames(seu@meta.data)) {
    p_condition <- VlnPlot(seu, features = count_column, pt.size = 0, group.by = condition_col, cols = condition_cls) +
      ggtitle(paste(count_column," per Cell by Condition")) +
      stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none")
    print(p_condition)
  }

  # Violin plots split by region
  if (region_col %in% colnames(seu@meta.data)) {
    p_region <- VlnPlot(seu, features = count_column, pt.size = 0, group.by = region_col, cols = region_colors) +
      ggtitle(paste(count_column," per Cell by Region")) +
      stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none")
    print(p_region)
  }

  if (boxplots) {
    if (region_col %in% colnames(seu@meta.data)) {
      p_region <- ggplot(seu@meta.data, aes_string(x = region_col, y = count_column, fill = region_col)) +
        geom_boxplot(outlier.size = 0.5) +
        scale_fill_manual(values = region_colors) +
        labs(title = paste(count_column," per Region"), x = "Region", y = count_column) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none"
        )
      print(p_region)
    }
  }

  if (region_col %in% colnames(seu@meta.data)) {
    p_region <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0, group.by = region_col, cols = region_colors) +
      ggtitle("nFeature_RNA per Cell by Region") +
      stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none")
    print(p_region)
  }
  
  if (boxplots) {
    if (region_col %in% colnames(seu@meta.data)) {
      p_region <- ggplot(seu@meta.data, aes_string(x = region_col, y = "nFeature_RNA", fill = region_col)) +
        geom_boxplot(outlier.size = 0.5) +
        scale_fill_manual(values = region_colors) +
        labs(title = "nFeature_RNA per Region", x = "Region", y = "nFeature_RNA") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none"
        )
      print(p_region)
    }
  }
    
    # Bar plot of cell counts per region
    cell_counts_df <- as.data.frame(table(seu@meta.data[[region_col]]))
    colnames(cell_counts_df) <- c("Region", "Cell_Count")
    
    p_counts <- ggplot(cell_counts_df, aes(x = Region, y = Cell_Count, fill = Region)) +
      geom_bar(stat = "identity", color = "black") +
      labs(title = "Number of Cells per Region", x = "Region", y = "Cell Count") +
      scale_fill_manual(values = region_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none")
    print(p_counts)
  
  # Cell size distribution
  if ("Area.um2" %in% colnames(md)) {
    area_plot <- ggplot(md, aes(x = Area.um2)) +
      geom_histogram(bins = 100, fill = "grey70", color = "black") +
      geom_vline(xintercept = c(40, 600), color = "red", linetype = "dashed", linewidth = 0.5) +
      theme_minimal() +
      ggtitle(paste("Distribution of Cell Sizes (µm²) -", slide_name)) +
      xlab("Area (µm²)") + ylab("Cell count")
    print(area_plot)
  }

if (boxplots) {
  if (region_col %in% colnames(seu@meta.data)) {
  p_area <- ggplot(seu@meta.data, aes_string(x = region_col, y = "Area.um2", fill = region_col)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = region_colors) +
    labs(title = "Cell Area per Region", x = "Region", y = "Area (µm²)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "none"
    )
  print(p_area)
}
}

## FOV Level
if (fov_plots) {
create_qc_plots(seu, batch_col = "slide_id", fov_col = "fov", count_col = count_column, feature_col = "nFeature_RNA", 
                  neg_col = "nCount_negprobes", fc_col = "nCount_falsecode")
}

  
  dev.off()
}

# Function: Run QC for list of Seurat objects
qc_all_slides <- function(seurat_list, outdir = ".", file_prefix = "QC_Report_", filter_log = NULL, count_column = "nCount_RNA", region_col = "region", flag_col = "qcCellsFlagged", boxplots = TRUE, fov_plots = TRUE) {
  for (nm in names(seurat_list)) {
    message("Running QC for: ", nm)
    qc_report(seurat_list[[nm]], nm, outdir, file_prefix = file_prefix, filter_log = filter_log, count_column = count_column, region_col = region_col, flag_col = flag_col, boxplots = boxplots, fov_plots = fov_plots)
  }
}


filter_seurat_objects <- function(seurat_list, thresholds_list = NULL, low_thres = 0.025, high_thres = 0.975, area_upper = 600, area_lower = 40, flag_col = "remove_flagged_cells", region_col = "region", min_cells_per_region = 500) {
  filtered_list <- list()
  log_list <- list()
  
  for (slide_name in names(seurat_list)) {
    message("Filtering: ", slide_name)
    seu <- seurat_list[[slide_name]]
    # seu@images <- list()
    # seu <- suppressWarnings(suppressMessages(UpdateSeuratObject(seu)))
    # seu <- UpdateSeuratObject(seu)
    total_cells <- ncol(seu)

    # Count flagged cells before any filtering
    message("Counting flagged cells before filtering ", slide_name)
    flagged_before_total <- if (flag_col %in% colnames(seu@meta.data)) {
      sum(seu@meta.data[[flag_col]], na.rm = TRUE)
    } else 0
    
    # Determine thresholds
    message("Determining thresholds ", slide_name)
    counts_vector <- seu$nCount_RNA
    if (!is.null(thresholds_list) && slide_name %in% names(thresholds_list)) {
      lower <- thresholds_list[[slide_name]]$lower
      upper <- thresholds_list[[slide_name]]$upper
    } else {
      lower <- quantile(counts_vector, low_thres)
      upper <- quantile(counts_vector, high_thres)
    }
    
    # Quantile / threshold filtering
    bottom_removed <- sum(counts_vector <= lower)
    top_removed <- sum(counts_vector >= upper)
    seu_filtered <- subset(seu, subset = nCount_RNA > lower & nCount_RNA < upper)
    
    remaining_after_counts <- ncol(seu_filtered)
    
    bottom_pct <- bottom_removed / total_cells * 100
    top_pct <- top_removed / total_cells * 100

    # Count flagged cells after threshold filtering
    message("Filtering flagged cells for ", slide_name)
    flagged_after_threshold <- if (flag_col %in% colnames(seu_filtered@meta.data)) {
      sum(seu_filtered@meta.data[[flag_col]], na.rm = TRUE)
    } else 0
    
    # Filter flagged cells
    # seu_filtered <- subset(seu_filtered, subset = remove_flagged_cells == FALSE)
    keep_cells <- rownames(seu_filtered@meta.data)[!seu_filtered@meta.data[[flag_col]]]
    seu_filtered <- subset(seu_filtered, cells = keep_cells)
    flagged_removed <- flagged_after_threshold
    flagged_pct <- flagged_removed / total_cells * 100

    # Filter by Area.um2
    area_removed <- 0
    if ("Area.um2" %in% colnames(seu_filtered@meta.data)) {
      area_removed <- sum(seu_filtered$Area.um2 > area_upper | seu_filtered$Area.um2 < area_lower, na.rm = TRUE)
      if (area_removed > 0) {
        seu_filtered <- subset(seu_filtered, subset = Area.um2 <= area_upper & Area.um2 >= area_lower)
      }
    }

    # Region-based Filtering
    regions_removed_count <- 0
    if (region_col %in% colnames(seu_filtered@meta.data)) {
      # Count cells currently in each region
      region_counts <- table(seu_filtered@meta.data[[region_col]])
      # Identify which regions meet the minimum threshold
      keep_regions <- names(region_counts[region_counts >= min_cells_per_region])
      # Calculate how many cells are in the 'dropped' regions for the log
      regions_removed_count <- sum(region_counts[!(names(region_counts) %in% keep_regions)])
      # Perform the subsetting
      seu_filtered <- subset(seu_filtered, 
                             cells = rownames(seu_filtered@meta.data)[seu_filtered@meta.data[[region_col]] %in% keep_regions])
      seu_filtered@meta.data[[region_col]] <- droplevels(as.factor(seu_filtered@meta.data[[region_col]]))

      message(paste0("  Regions dropped: ", length(region_counts) - length(keep_regions), 
                     " (Total cells removed from these regions: ", regions_removed_count, ")"))
    }

    remaining_cells <- ncol(seu_filtered)
    
    # Store filtered object
    filtered_list[[slide_name]] <- seu_filtered
    
    #LOG
    log_list[[slide_name]] <- data.frame(
      Slide = slide_name,
      Total_Cells = total_cells,
      Bottom_Threshold = lower,
      Bottom_Removed = bottom_removed,
      Bottom_Percent = round(bottom_pct, 2),
      Top_Threshold = upper,
      Top_Removed = top_removed,
      Top_Percent = round(top_pct, 2),
      Flagged_Before_Threshold = flagged_before_total,
      Flagged_After_Threshold = flagged_after_threshold,
      Flagged_Removed = flagged_removed,
      Flagged_Percent = round(flagged_pct, 2),
      Area_Upper_Threshold = area_upper,
      Area_Lower_Threshold = area_lower,
      Area_Removed = area_removed,
      Remaining_Cells = remaining_cells,
      Remaining_Percent = round(remaining_cells / total_cells * 100, 2),
      Cells_Removed = bottom_removed + top_removed + flagged_removed + area_removed,
      Cells_Removed_Percent = round((bottom_removed + top_removed + flagged_removed + area_removed) / total_cells * 100, 2),
      stringsAsFactors = FALSE
    )
    
    message(paste0("Filtered slide: ", slide_name,
                   "  | Total cells: ", total_cells,
                   "  | Cells removed: ", bottom_removed + top_removed + flagged_removed + area_removed,
                   "  | Cells removed (calc): ", total_cells - remaining_cells,
                   " | Remaining cells: ", remaining_cells))
  }
  
  filter_log <- do.call(rbind, log_list)
  rownames(filter_log) <- NULL
  
  return(list(filtered_objects = filtered_list, filter_log = filter_log))
}


## BATCH CORRECTION
apply_scPearsonPCA <- function(seu, 
                              nfeatures = 3000,
                              slot_names = list(
                                                pca = "pearsonpca",
                                                umap = "pearsonumap",
                                                graph = "pearsongraph",
                                                clusters = "pearson_clusters"
                                              )
  ) {
    message("Calculating total counts and gene frequencies...")
    tc <- Matrix::colSums(seu[["RNA"]]@counts) ## total counts per cell (across all genes)
    genefreq <- scPearsonPCA::gene_frequency(seu[["RNA"]]@counts) ## gene frequency (across all cells)

    if (!sum(genefreq)==1) {
      message("Gene frequencies do not sum to 1.")
      return(NULL)
        } else {
        message("Gene frequencies sum to 1. Check passed.")
        }

    message("Finding highly variable genes...")    
    seu <- Seurat::FindVariableFeatures(seu, nfeatures = nfeatures)
    hvgs <- seu@assays$RNA@var.features 

    message("Performing scPearsonPCA...")   
    pcaobj <- sparse_quasipoisson_pca_seurat(seu[["RNA"]]@counts[hvgs,]
                               ,totalcounts = tc
                               ,grate = genefreq[hvgs]
                               ,scale.max = 10 ## PCs reflect clipping pearson residuals > 10 SDs above the mean pearson residual
                               ,do.scale = TRUE ## PCs reflect as if pearson residuals for each gene were scaled to have standard deviation=1
                               ,do.center = TRUE ## PCs reflect as if pearson residuals for each gene were centered to have mean=0
                               )

    message("Making umap...")
    umapobj <- scPearsonPCA::make_umap(pcaobj)
    seu[[slot_names$pca]] <- pcaobj$reduction.data
    seu[[slot_names$umap]] <- umapobj$ump  ## umap
    seu[[slot_names$graph]] <- Seurat::as.Graph(umapobj$grph) ## nearest neighbors / adjacency matrix used for unsupervised clustering
    
    message("Finding clusters...")
    seu <- Seurat::FindClusters(seu, graph = slot_names$graph)
    seu@meta.data[[slot_names$clusters]] <- seu@meta.data$seurat_clusters

    list(
    tc = tc,
    hvgs = hvgs,
    seu = seu
  )
}

## (minor side-note), this helper function `make_umap` calls the same package 
## and function as Seurat::RunUMAP, but also returns the 
## 'nearest-neighbor graph' used in the UMAP algorithm.
## Using the same nearest neighbor graph for clustering and UMAP, 
## can give better concordance between UMAP and unsupervised clusters.

make_umap <- function(pcaobj, min_dist=0.01, n_neighbors=30, metric="cosine",key ="UMAP_" ){
  ump <- 
    uwot::umap(pcaobj$reduction.data@cell.embeddings
               ,n_neighbors = n_neighbors
               ,nn_method = "annoy"
               ,metric = metric
               ,min_dist = min_dist
               ,ret_extra = c("fgraph","nn")
               ,verbose = TRUE)
  
  umpgraph <- ump$fgraph
  dimnames(umpgraph) <- list(rownames(ump$nn[[1]]$idx), rownames(ump$nn[[1]]$idx))
  colnames(ump$embedding) <- paste0(key, c(1,2)) 
  ump <- Seurat::CreateDimReducObject(embeddings = ump$embedding, key = key)
  return(list(grph = umpgraph
              ,ump = ump))
}


# From: https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/blob/Main/_code/HieraType/R/utils.R
#' Normalize cell expression in a raw counts matrix by their totalcounts
#' param counts_matrix a cells x genes matrix of raw counts to be normalized
#' param tc optional vector of totalcounts.  Useful if providing a counts_matrix based on a subset of genes.  
#' If not provided, totalcounts per cells is taken to be the rowSums of the counts matrix.
totalcount_norm <- function(counts_matrix, tc = NULL){
  if(is.null(tc)) tc <- Matrix::rowSums(counts_matrix)
  scale.factor <- mean(tc)
  tc[tc==0] <- 1
  return(Matrix::Diagonal(x = scale.factor/tc, names = TRUE) %*% counts_matrix)
}

## Functions relevant for Cell Typing with InSituType

plot_composition <- function(seu, cluster_col, split_by, type = "relative", palette = NULL) {
  
  # 1. Extract and aggregate data
  df <- seu@meta.data %>%
    dplyr::group_by(!!rlang::sym(cluster_col), !!rlang::sym(split_by)) %>%
    dplyr::tally() %>%
    dplyr::ungroup()
  
  colnames(df) <- c("Cluster", "SplitVar", "Count")
  
  # 2. Determine scaling and labels
  # 'fill' scales bars to 1 (100%), 'stack' keeps raw counts
  pos <- if (type == "relative") "fill" else "stack"
  label_y <- if (type == "relative") "Proportion of Cells" else "Number of Cells"

  # 3. Setup Palette (Default to alphabet if NULL)
  if (is.null(palette)) {
    unique_items <- as.character(unique(df$SplitVar))
    n_items <- length(unique_items)
    
    if (n_items <= 26) {
      # Use alphabet directly if small enough
      pal_vec <- as.character(pals::alphabet(n = n_items))
    } else {
      # Interpolate alphabet colors to handle 26+ items
      pal_vec <- colorRampPalette(as.character(pals::alphabet()))(n_items)
    }
    palette <- setNames(pal_vec, unique_items)
  }

  p <- ggplot(df, aes(x = Cluster, y = Count, fill = SplitVar)) +
    geom_bar(stat = "identity", position = pos) +
    scale_fill_manual(values = palette) + # Use the alphabet palette
    theme_bw() +
    labs(
      title = paste("Composition of", cluster_col, "by", split_by),
      subtitle = paste("Mode:", type),
      x = "Cluster",
      y = label_y,
      fill = split_by
    ) +
    theme(
      # Updated to 90 degrees and adjusted justification to align with ticks
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank() # Clean up the X-axis grid for bar charts
    )
  
  return(p)
}


# CT_QC_plot:
# This function generates a comprehensive PDF report of cluster-level quality control (QC) plots for cell typing validation.
#
# The PDF includes:
#   - Barplot of cell counts per cluster.
#   - Violin plots for QC features across clusters (e.g., marker expression, gene/cell counts).
#   - Optional heatmap and flightpath/trajectory plots if InSituType results provided.
#   - Embedding plots of clusters, annotations, and study IDs for visual inspection.

CT_QC_plot <- function(seu, cluster_col, cluster_pal=NULL, annotation_col = NULL, annotation_pal = NULL, IST_obj = NULL, out_dir = ".", reduction = NULL, split_by_col = "study_id") {
    
    pdf(file.path(out_dir, paste0("Cluster_QC_", cluster_col, ".pdf")), width = 10, height = 10)

    print(count_per_ct(seu, cluster_col=cluster_col, cluster_pal = cluster_pal))
    print(feature_vlns(seu, cluster_col=cluster_col, features = c("Mean.PanCK", "Mean.CD45", "nCount_RNA", "nFeature_RNA"), cluster_pal = cluster_pal))


    if (!is.null(IST_obj)) {
        heatmap(sweep(IST_obj$profiles, 1, pmax(apply(IST_obj$profiles, 1, max), .2), "/"), scale = "none",
            main = "Cluster mean expression profiles")
        fp_layout(IST_obj, cluster_pal = cluster_pal)
        print(flightpath_plot(flightpath_result = NULL, insitutype_result = IST_obj, col = cluster_pal[IST_obj$clust]))
    }

    if (reduction %in% names(seu@reductions)) {
        p1 <- plot_embedding(
            seu,
            reduction = reduction,
            group.by = cluster_col,
            label = TRUE,
            palette = cluster_pal,
            legend = TRUE
            ) 
        print(p1)

        if (!is.null(annotation_col)) {
            p2 <- plot_embedding(
                seu,
                reduction = reduction,
                group.by = annotation_col,
                label = TRUE,
                palette = annotation_pal,
                legend = TRUE
                ) 
            print(p2)
        }

        p3 <- plot_embedding(
            seu,
            reduction = reduction,
            group.by = split_by_col,
            label = TRUE,
            palette = NULL,
            legend = TRUE
            ) 
        print(p3)
    }

    print(plot_composition(seu, cluster_col = cluster_col, split_by = split_by_col, type = "relative"))
    print(plot_composition(seu, cluster_col = cluster_col, split_by = split_by_col, type = "absolute"))
    print(plot_composition(seu, cluster_col = split_by_col, split_by = cluster_col, type = "relative", palette = cluster_pal))
    print(plot_composition(seu, cluster_col = split_by_col, split_by = cluster_col, type = "absolute", palette = cluster_pal))

    dev.off()
    }

count_per_ct <- function(seu, cluster_col, xlab = "Cell Type",
                         ylab = "Cell Count", title = "Cell Count per Cell Type",
                         cluster_pal = NULL) {

    df <- seu@meta.data %>%
        dplyr::count(celltype = .data[[cluster_col]]) %>%
        dplyr::arrange(desc(n))

    p <- ggplot(df, aes(x = reorder(celltype, -n), y = n, fill = celltype)) +
        geom_bar(stat = "identity") +
        labs(x = xlab, y = ylab, title = title) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")

    if (!is.null(cluster_pal)) {
        p <- p + scale_fill_manual(values = cluster_pal)}
    p
}

feature_vlns <- function(seu, cluster_col, features = c("Mean.PanCK", "Mean.CD45"), cluster_pal = NULL) {
    levels_ct <- levels(seu[[cluster_col]][,1])
    cols_vec <- cluster_pal[levels_ct]

    plots <- lapply(features, function(feature) {
        VlnPlot(
            seu,
            features = feature,
            group.by = cluster_col,
            pt.size = 0,
            cols = cols_vec
        ) +
            ggtitle(paste0(feature, " by ", cluster_col)) +
            theme_bw() +
            theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5)
            ) +
            coord_cartesian(clip = "off")
        })

    plots
}

fp_layout <- function(IST_obj, cluster_pal = NULL) {
    if (is.null(cluster_pal)) {
        cols <- InSituType::colorCellTypes(freqs = table(IST_obj$clust), palette = "brewers")
    } else {
        cols <- cluster_pal
    }

    flightpath <- InSituType::flightpath_layout(logliks = IST_obj$logliks, profiles = IST_obj$profiles)
    par(mar = c(0,0,0,0))
    plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[IST_obj$clust])
    text(flightpath$clustpos[, 1], flightpath$clustpos[, 2], rownames(flightpath$clustpos), cex = 0.7)
    par(mar = c(5, 4, 4, 2) + 0.1)
}


plot_all_cts <- function(seu, cluster_col, cluster_pal = NULL, out_dir = ".", save_pngs = FALSE, folder_name = "CT_PLOTS") {
    # Extract the coordinates from metadata
    xy <- as.matrix(seu@meta.data[, c("x_slide_mm", "y_slide_mm")])

    if (save_pngs) {
        # Create the output directory if it doesn't exist
        out_dir <- file.path(out_dir, folder_name)
        if (!dir.exists(out_dir)) {
            dir.create(out_dir)
        }
    }

    for (ct in unique(seu@meta.data[[cluster_col]])) {
        if (save_pngs) {
            png(file.path(out_dir, paste0("celltype_", ct, "_spread.png")), 
            width = diff(range(xy[,1]))*.7, height = diff(range(xy[,2]))*.7, units = "in", 
            res = 400)  # res of 400 is pretty good; 600 is publication-quality
            par(mar = c(0,0,0,0))
        }
        plot(xy, pch = 16, col = scales::alpha(cluster_pal[seu@meta.data[[cluster_col]]], 0.3), cex = 0.1,
            xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        points(xy[semisup$clust == ct, ], pch = 16, cex = 0.1, col = "black")
        legend("top", legend = ct)
        if (save_pngs) {
            dev.off()
        }
    }

    if (save_pngs) {
        par(mar = c(5, 4, 4, 2) + 0.1)
    }
}

plot_hm_pdf <- function(IST_obj, out_dir = ".", file_name = "Cluster_mean_expression_profiles.pdf") {
    pdf(file.path(out_dir, paste0(file_name)), width = 6, height = 20)
    
    heatmap(sweep(IST_obj$profiles, 1, pmax(apply(IST_obj$profiles, 1, max), .2), "/"), scale = "none",
        main = "Cluster mean expression profiles")

    dev.off()
}



## DEG specific functions
plot_volcano <- function(df, 
                         x_col, 
                         y_col, 
                         x_trans = NULL, 
                         y_trans = "-log10", 
                         x_cut = 0.5, 
                         p_cut = 0.05, 
                         genes_of_interest = NULL,
                         sig_color = "red", 
                         insig_color = "grey50", 
                         poi_color = "blue") {
  
  # 1. Prepare data and apply transformations
  plt_df <- df
  plt_df$X <- plt_df[[x_col]]
  plt_df$Y <- plt_df[[y_col]]
  
  if (!is.null(x_trans) && x_trans == "log2") plt_df$X <- log2(plt_df$X)
  if (!is.null(y_trans) && y_trans == "-log10") plt_df$Y_plot <- -log10(plt_df$Y) else plt_df$Y_plot <- plt_df$Y
  
  # 2. Define Significance categories
  # Note: logic uses the transformed X (Effect) but the RAW p-value threshold (p_cut)
  plt_df$status <- "Not Significant"
  plt_df$status[plt_df$Y < p_cut & abs(plt_df$X) > x_cut] <- "Significant"
  
  # 3. Handle specific genes of interest (POI)
  plt_df$is_poi <- FALSE
  if (!is.null(genes_of_interest)) {
    plt_df$is_poi <- rownames(plt_df) %in% genes_of_interest
  }
  
  # 4. Determine plotting order (so blue points are on top)
  plt_df <- plt_df[order(plt_df$status, plt_df$is_poi), ]

  # 5. Build Plot
  g <- ggplot(plt_df, aes(x = X, y = Y_plot)) +
    # Background points
    geom_point(aes(color = status), alpha = 0.6, size = 1.5) +
    # Highlight Genes of Interest
    geom_point(data = subset(plt_df, is_poi), color = poi_color, size = 2.5) +
    # Vertical lines for effect size
    geom_vline(xintercept = c(-x_cut, x_cut), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
    # Horizontal line for p-value
    geom_hline(yintercept = ifelse(y_trans == "-log10", -log10(p_cut), p_cut), 
               linetype = "dashed", color = "darkgrey") +
    # Labels for significant OR points of interest
    geom_text_repel(
      data = subset(plt_df, status == "Significant" | is_poi),
      aes(label = rownames(subset(plt_df, status == "Significant" | is_poi))),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5,
      # Color labels blue if they are POIs, otherwise black
      color = ifelse(subset(plt_df, status == "Significant" | is_poi)$is_poi, poi_color, "black")
    ) +
    scale_color_manual(values = c("Significant" = sig_color, "Not Significant" = insig_color)) +
    labs(
      x = paste("Effect Size (", x_col, ")"),
      y = paste(y_trans, "(", y_col, ")"),
      color = "Status"
    ) +
    theme_classic() +
    theme(legend.position = "top")

  return(g)
}

get_gene_bin_plot_list <- function(seu_obj, 
                                   genes, 
                                   bin_col, 
                                   bins_to_include = NULL, 
                                   bin_order = NULL,
                                   bin_colors = NULL) {
  
  # 1. Verify genes
  genes <- intersect(genes, rownames(seu_obj))
  if(length(genes) == 0) stop("None of the provided genes found.")
  
  # 2. Extract and Prepare Data
  # We use data.table for efficiency
  meta <- as.data.table(seu_obj@meta.data)
  expr_data <- Seurat::GetAssayData(seu_obj, slot = "data")[genes, , drop = FALSE]
  
  plot_list <- list()
  
  # 3. Loop through genes and create individual plots
  for (goi in genes) {
    
    # Construct DT for this specific gene
    dt <- data.table(
      expr = as.numeric(expr_data[goi, ]),
      bin = meta[[bin_col]]
    )
    
    # Filter and order bins
    if (!is.null(bins_to_include)) dt <- dt[bin %in% bins_to_include]
    if (!is.null(bin_order)) dt[, bin := factor(bin, levels = bin_order)]
    
    # Calculate Mean and Standard Error
    summary_dt <- dt[, .(
      mean_expr = mean(expr, na.rm = TRUE),
      se = sd(expr, na.rm = TRUE) / sqrt(.N)
    ), by = bin]
    
    # 4. Create the Plot
    p <- ggplot(summary_dt, aes(x = bin, y = mean_expr, fill = bin)) +
      geom_col(color = "black", width = 0.7) +
      geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se), 
                    width = 0.2, color = "black") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Adds 10% space at top
      theme_classic() +
      labs(
        title = paste("Gene:", goi),
        subtitle = paste("Cell type:", celltype_oi),
        x = "Spatial Zone",
        y = "Mean Normalized Expression (+/- SE)",
        fill = "Zone"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none" # Hide legend to maximize space on the PDF page
      )
    
    if (!is.null(bin_colors)) p <- p + scale_fill_manual(values = bin_colors)
    
    plot_list[[goi]] <- p
  }
  
  return(plot_list)
}

generate_spatial_dotplot <- function(seu_obj, 
                                     genes, 
                                     bin_col,
                                     scale = TRUE,
                                     bins_to_include = NULL, 
                                     bin_order = NULL) {
  
  # 1. Filter genes to ensure they exist
  genes <- intersect(genes, rownames(seu_obj))
  
  # 2. Handle subsetting and ordering
  plot_obj <- seu_obj
  
  # Filter bins if requested
  if (!is.null(bins_to_include)) {
    plot_obj <- subset(plot_obj, cells = colnames(plot_obj)[plot_obj[[bin_col, drop=TRUE]] %in% bins_to_include])
  }
  
  # Enforce factor levels for the bin column to control plot order
  if (!is.null(bin_order)) {
    plot_obj@meta.data[[bin_col]] <- factor(plot_obj@meta.data[[bin_col]], levels = bin_order)
  }

  # 3. Generate DotPlot using group.by
  p <- DotPlot(
    object = plot_obj, 
    features = genes, 
    group.by = bin_col,      # Use group.by instead of changing Idents
    cols = c("lightgrey", "firebrick"), 
    scale = scale, 
    col.min = -2.5, 
    col.max = 2.5
  ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Spatial Gene Expression",
      x = "Genes",
      y = "Spatial Zone",
      size = "Percent Expressed",
      color = "Average Expression (Z-score)"
    )
  
  return(p)
}