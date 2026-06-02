#!/usr/bin/env Rscript

# Generate mouse-level figures directly from an analysis-ready save_list object.
# This script does not require raw KPTracer data.

suppressPackageStartupMessages({
  library(ape)
  library(dendextend)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_file) == 0L) getwd() else dirname(normalizePath(script_file))

save_list_file <- Sys.getenv("SAVE_LIST_FILE", unset = file.path(script_dir, "outputs", "save_list_all_mice_linemap.rds"))
mouse_id <- Sys.getenv("MOUSE_ID", unset = "3522")
base_out_dir <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = file.path(script_dir, "outputs", "figures"))

save_list <- readRDS(save_list_file)

plot_one_mouse <- function(mouse_id, x, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cell2T <- setNames(x$table$tumor, x$table$cellID)
  tumor_color_map <- tapply(
    x$colors$tip_col[x$table$cellID],
    x$table$tumor,
    function(z) names(sort(table(z), decreasing = TRUE))[1]
  )
  tumor_color_map <- tumor_color_map[order(names(tumor_color_map))]

  add_tumor_legend <- function() {
    legend(
      "topright",
      legend = names(tumor_color_map),
      fill = unname(tumor_color_map),
      border = NA,
      title = "Tumor",
      bty = "n",
      cex = 0.85
    )
  }

  plot_dend_with_bars <- function(dend, main, file) {
    pdf(file, width = 12, height = 8)
    plot(dend, main = main, ylab = "height")
    add_tumor_legend()
    colored_bars(
      cbind(
        tumor = x$colors$tip_col,
        fitness = x$colors$fitness,
        plas = x$colors$plasticity,
        celltype = x$colors$celltype,
        lentiBC = x$colors$lentiBC
      ),
      dend,
      rowLabels = c("tumor", "fitnesses", "plas", "cell type", "lentiBC")
    )
    dev.off()
  }

  make_display_dend <- function(hclust_obj) {
    dend <- as.dendrogram(hclust_obj)
    old_labels <- labels(dend)
    clustercut <- unname(x$colors$tip_col[x$cellBC])
    labels(dend) <- unname(cell2T[old_labels])
    labels_colors(dend) <- clustercut[order.dendrogram(dend)]
    hang.dendrogram(dend, hang = 0.01)
  }

  dend1 <- make_display_dend(x$hclust1)
  plot_dend_with_bars(
    dend1,
    paste(mouse_id, "dist(distance matrix): clone seperation"),
    file.path(out_dir, paste0(mouse_id, "_dist_distance_matrix_clone_seperation.pdf"))
  )

  dend2 <- make_display_dend(x$hclust2)
  plot_dend_with_bars(
    dend2,
    paste(mouse_id, "reconstructed tree from distance matrix by UPGMA"),
    file.path(out_dir, paste0(mouse_id, "_reconstructed_tree_UPGMA.pdf"))
  )

  tree_nj <- x$tree_nj
  original_tip_labels <- tree_nj$tip.label
  tip_color <- x$colors$tip_col[original_tip_labels]
  tree_nj$tip.label <- cell2T[original_tip_labels]

  pdf(file.path(out_dir, paste0(mouse_id, "_reconstructed_tree_NJ.pdf")), width = 12, height = 9)
  par(mar = c(7, 4, 4, 2) + 0.1, xpd = NA)
  plot(
    tree_nj,
    direction = "downwards",
    show.tip.label = TRUE,
    tip.color = tip_color,
    cex = 0.6,
    no.margin = FALSE,
    main = paste(mouse_id, "reconstructed tree from distance matrix by NJ")
  )
  add_tumor_legend()
  add.scale.bar(length = 0.1, lwd = 2, cex = 0.9)
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_x <- pp$xx[seq_along(original_tip_labels)]
  tip_y <- pp$yy[seq_along(original_tip_labels)]
  annotation_cols <- rbind(
    tumor = x$colors$tip_col[original_tip_labels],
    fitness = x$colors$fitness[match(original_tip_labels, x$cellBC)],
    plas = x$colors$plasticity[match(original_tip_labels, x$cellBC)],
    celltype = x$colors$celltype[match(original_tip_labels, x$cellBC)],
    lentiBC = x$colors$lentiBC[match(original_tip_labels, x$cellBC)]
  )
  y_step <- diff(range(pp$yy, na.rm = TRUE)) * 0.035
  y_start <- min(tip_y, na.rm = TRUE) - y_step * 1.5
  x_width <- median(diff(sort(unique(tip_x)))) * 0.85
  if (!is.finite(x_width) || x_width <= 0) x_width <- 0.35
  for (i in seq_len(nrow(annotation_cols))) {
    y <- y_start - (i - 1L) * y_step
    rect(
      tip_x - x_width / 2,
      y - y_step * 0.32,
      tip_x + x_width / 2,
      y + y_step * 0.32,
      col = annotation_cols[i, ],
      border = NA
    )
    text(min(tip_x, na.rm = TRUE) - x_width * 4, y, rownames(annotation_cols)[i], adj = 1, cex = 0.7)
  }
  dev.off()

  invisible(TRUE)
}

if (toupper(mouse_id) == "ALL") {
  failed <- character()
  for (id in names(save_list)) {
    x <- save_list[[id]]
    out_dir <- file.path(base_out_dir, id)
    ok <- tryCatch(
      {
        plot_one_mouse(id, x, out_dir)
        message("Wrote figures to: ", out_dir)
        TRUE
      },
      error = function(e) {
        message("Skipping mouse ", id, " because plotting failed: ", conditionMessage(e))
        FALSE
      }
    )
    if (!ok) failed <- c(failed, id)
  }
  if (length(failed) > 0L) {
    message("Plotting failed for: ", paste(failed, collapse = ", "))
  }
  quit(save = "no", status = 0)
}

x <- save_list[[mouse_id]]
if (is.null(x)) stop("Mouse `", mouse_id, "` was not found in: ", save_list_file)
out_dir <- if (basename(base_out_dir) == mouse_id) base_out_dir else file.path(base_out_dir, mouse_id)
plot_one_mouse(mouse_id, x, out_dir)
message("Wrote figures to: ", out_dir)
