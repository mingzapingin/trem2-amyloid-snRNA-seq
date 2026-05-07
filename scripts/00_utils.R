# ============================================================
# 00_utils.R
# Shared helpers sourced by downstream pipeline scripts.
# Source this at the top of any script that needs these functions.
# ============================================================

library(ggplot2)

# ------------------------------------------------------------------
# Marker-list helpers
# ------------------------------------------------------------------

#' Keep only genes present in a Seurat object or character vector
filter_marker_list <- function(marker_list, ref) {
  genes_present <- if (is.character(ref)) ref else rownames(ref)
  out <- lapply(marker_list, function(x) x[x %in% genes_present])
  out[lengths(out) > 0]
}

#' Table showing which markers are / are not found in the object
marker_presence_table <- function(marker_list, ref) {
  genes_present <- if (is.character(ref)) ref else rownames(ref)
  dplyr::bind_rows(lapply(names(marker_list), function(nm) {
    tibble::tibble(
      marker_group       = nm,
      gene               = marker_list[[nm]],
      present_in_object  = marker_list[[nm]] %in% genes_present
    )
  }))
}

# ------------------------------------------------------------------
# Safe wrappers
# ------------------------------------------------------------------

safe_join_layers <- function(obj, assay = "RNA") {
  tryCatch(JoinLayers(obj, assay = assay), error = function(e) obj)
}

safe_save_plot <- function(plot_obj, filename,
                           width = 12, height = 8, dpi = 300) {
  tryCatch(
    ggsave(filename, plot = plot_obj,
           width = width, height = height, dpi = dpi),
    error = function(e) {
      png(filename, width = width, height = height, units = "in", res = dpi)
      print(plot_obj)
      dev.off()
    }
  )
}

# ------------------------------------------------------------------
# matrixplot()
# Heatmap of mean expression per group, inspired by sc.pl.matrixplot.
# ------------------------------------------------------------------

matrixplot <- function(
    seurat_obj,
    features,
    group.by,
    assay            = NULL,
    layer            = "data",
    standard_scale   = "var",
    cluster_rows     = FALSE,
    cluster_cols     = FALSE,
    categories_order = NULL,
    swap_axes        = FALSE,
    col_palette      = c("white", "#08519c"),
    vmin             = 0,
    vmax             = 1,
    vcenter          = NULL,
    border_color     = "grey80",
    border_width     = 0.1,
    colorbar_title   = "Mean expression\nin group",
    title            = NULL,
    x_label          = NULL,
    y_label          = NULL,
    x_text_angle     = 90,
    base_size        = 11,
    return_data      = FALSE
) {
  # --- Handle named-list features (gene grouping) ---
  feature_groups <- NULL
  if (is.list(features) && !is.null(names(features))) {
    feature_groups <- stack(features)
    colnames(feature_groups) <- c("feature", "feature_group")
    feature_groups$feature       <- as.character(feature_groups$feature)
    feature_groups$feature_group <- factor(feature_groups$feature_group,
                                           levels = names(features))
    # Deduplicate: keep the first group assignment for each gene
    feature_groups <- feature_groups[!duplicated(feature_groups$feature), ]
    features <- feature_groups$feature
  }
  
  # --- Extract expression matrix ---
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
  
  expr_mat <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer),
    error = function(e)
      Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
  )
  
  missing <- setdiff(features, rownames(expr_mat))
  if (length(missing) > 0) {
    warning("Features not found and will be dropped: ",
            paste(missing, collapse = ", "))
    features <- intersect(features, rownames(expr_mat))
  }
  # Ensure no duplicates (intersect already deduplicates, but be safe)
  features <- unique(features)
  if (length(features) == 0) stop("No valid features remaining.")
  
  # Keep feature_groups in sync with surviving features
  if (!is.null(feature_groups)) {
    feature_groups <- feature_groups[feature_groups$feature %in% features, ]
  }
  
  expr_mat <- expr_mat[features, , drop = FALSE]
  
  # --- Grouping vector ---
  groups <- seurat_obj@meta.data[[group.by]]
  
  # --- Compute group means ---
  if (inherits(expr_mat, "dgCMatrix") || inherits(expr_mat, "sparseMatrix")) {
    expr_mat <- as.matrix(expr_mat)
  }
  
  group_levels <- if (!is.null(categories_order)) {
    categories_order
  } else if (is.factor(groups)) {
    levels(groups)
  } else {
    sort(unique(groups))
  }
  
  mean_df <- do.call(rbind, lapply(group_levels, function(g) {
    cells <- which(groups == g)
    if (length(cells) == 0) return(rep(NA_real_, length(features)))
    if (length(cells) == 1) return(expr_mat[features, cells])
    rowMeans(expr_mat[features, cells, drop = FALSE])
  }))
  rownames(mean_df) <- group_levels
  colnames(mean_df) <- features
  mean_df <- as.data.frame(mean_df)
  
  # --- Standard scaling ---
  if (!is.null(standard_scale)) {
    if (standard_scale == "group") {
      row_min <- apply(mean_df, 1, min, na.rm = TRUE)
      mean_df <- mean_df - row_min
      row_max <- apply(mean_df, 1, max, na.rm = TRUE)
      row_max[row_max == 0] <- 1
      mean_df <- mean_df / row_max
    } else if (standard_scale == "var") {
      col_min <- apply(mean_df, 2, min, na.rm = TRUE)
      mean_df <- sweep(mean_df, 2, col_min)
      col_max <- apply(mean_df, 2, max, na.rm = TRUE)
      col_max[col_max == 0] <- 1
      mean_df <- sweep(mean_df, 2, col_max, FUN = "/")
    }
  }
  
  if (return_data) return(mean_df)
  
  # --- Hierarchical clustering ---
  row_order <- group_levels
  if (cluster_rows && nrow(mean_df) > 2) {
    hc_row   <- hclust(dist(as.matrix(mean_df)))
    row_order <- rownames(mean_df)[hc_row$order]
  }
  col_order <- features
  if (cluster_cols && ncol(mean_df) > 2) {
    hc_col   <- hclust(dist(t(as.matrix(mean_df))))
    col_order <- colnames(mean_df)[hc_col$order]
  }
  
  # --- Long format ---
  mean_df$group <- rownames(mean_df)
  plot_df <- .mp_reshape_long(mean_df, id_col = "group")
  plot_df$group   <- factor(plot_df$group,   levels = rev(row_order))
  plot_df$feature <- factor(plot_df$feature, levels = col_order)
  
  # --- ggplot ---
  if (swap_axes) {
    p <- ggplot(plot_df, aes(x = group, y = feature, fill = value))
  } else {
    p <- ggplot(plot_df, aes(x = feature, y = group, fill = value))
  }
  
  p <- p + geom_tile(color = border_color, linewidth = border_width)
  p <- p + .mp_fill_scale(col_palette, vmin, vmax, vcenter, colorbar_title)
  
  p <- p +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x      = element_text(angle = x_text_angle,
                                      hjust = if (x_text_angle == 90) 1 else 0.5,
                                      vjust = 0.5),
      panel.grid        = element_blank(),
      panel.border      = element_rect(fill = NA, color = "black",
                                       linewidth = 0.8),
      axis.ticks        = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      legend.position   = "right"
    ) +
    labs(title = title, x = x_label, y = y_label) +
    coord_cartesian(expand = FALSE)
  
  if (!is.null(feature_groups) && !swap_axes) {
    p <- p + .mp_group_strip(feature_groups, col_order, base_size,
                             n_groups = length(row_order))
  }
  
  if (cluster_rows && nrow(mean_df) > 2) {
    p <- p + labs(caption = paste0("Rows clustered (", hc_row$method, " linkage)"))
  }
  if (cluster_cols && ncol(mean_df) > 2) {
    existing <- p$labels$caption
    note <- paste0("Columns clustered (", hc_col$method, " linkage)")
    p <- p + labs(caption = if (is.null(existing)) note
                  else paste(existing, note, sep = " | "))
  }
  
  p
}

# --- internal matrixplot helpers (prefixed with .mp_) ---

.mp_reshape_long <- function(df, id_col = "group") {
  features <- setdiff(colnames(df), id_col)
  data.frame(
    group   = rep(df[[id_col]], each = length(features)),
    feature = rep(features, times = nrow(df)),
    value   = as.vector(t(df[, features, drop = FALSE])),
    stringsAsFactors = FALSE
  )
}

.mp_fill_scale <- function(col_palette, vmin, vmax, vcenter, legend_title) {
  limits <- if (!is.null(vmin) || !is.null(vmax)) c(vmin, vmax) else NULL
  
  if (!is.null(vcenter)) {
    return(scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = vcenter, limits = limits,
      name = legend_title, na.value = "grey90", oob = scales::squish
    ))
  }
  if (length(col_palette) == 1 && !grepl("^#", col_palette) &&
      !col_palette %in% grDevices::colours()) {
    if (col_palette == "viridis") {
      return(scale_fill_viridis_c(
        limits = limits, name = legend_title,
        na.value = "grey90", oob = scales::squish
      ))
    }
    return(scale_fill_distiller(
      palette = col_palette, limits = limits, name = legend_title,
      na.value = "grey90", oob = scales::squish, direction = 1
    ))
  }
  scale_fill_gradientn(
    colours = col_palette, limits = limits, name = legend_title,
    na.value = "grey90", oob = scales::squish
  )
}

.mp_group_strip <- function(feature_groups, col_order, base_size, n_groups) {
  fg <- feature_groups
  fg$feature <- factor(fg$feature, levels = col_order)
  fg <- fg[order(fg$feature), ]
  
  rle_g  <- rle(as.character(fg$feature_group))
  ends   <- cumsum(rle_g$lengths) + 0.5
  starts <- c(0.5, ends[-length(ends)])
  mids   <- (starts + ends) / 2
  labels <- rle_g$values
  spans  <- ends - starts
  
  grp_colors <- scales::hue_pal()(length(unique(labels)))
  color_map  <- setNames(grp_colors, unique(labels))
  
  # --- Decide whether labels need rotation ---
  # Heuristic: ~3.5 characters fit per column unit at typical plot widths
  chars_per_col <- 3.5
  use_rotate <- any(nchar(labels) > spans * chars_per_col)
  
  # Scale factor: more y-groups → each data-unit is physically smaller,
  # so the strip needs more data-units for the same visual depth
  yscale <- n_groups / 10
  
  if (use_rotate) {
    max_nchar <- max(nchar(labels))
    # Strip depth grows with the longest label
    strip_depth      <- yscale * (0.5 + max_nchar * 0.06)
    bottom_margin_pt <- 15 + max_nchar * 2.5
  } else {
    strip_depth      <- yscale * 0.7
    bottom_margin_pt <- 25
  }
  
  gap        <- yscale * 0.1
  strip_ymax <- -gap
  strip_ymin <- -(strip_depth + gap)
  text_y     <- (strip_ymin + strip_ymax) / 2
  text_size  <- (base_size - 3) / ggplot2::.pt
  
  annotations <- list()
  for (i in seq_along(labels)) {
    annotations <- c(annotations, list(
      annotate("rect",
               xmin = starts[i], xmax = ends[i],
               ymin = strip_ymin, ymax = strip_ymax,
               fill = color_map[labels[i]], color = NA),
      annotate("text",
               x = mids[i], y = text_y,
               label = labels[i],
               size  = text_size,
               angle = if (use_rotate) 90 else 0,
               hjust = 0.5, vjust = 0.5,
               fontface = "bold")
    ))
  }
  
  c(annotations, list(
    coord_cartesian(clip = "off", expand = FALSE),
    theme(plot.margin = margin(5.5, 5.5, bottom_margin_pt, 5.5, unit = "pt"))
  ))
}