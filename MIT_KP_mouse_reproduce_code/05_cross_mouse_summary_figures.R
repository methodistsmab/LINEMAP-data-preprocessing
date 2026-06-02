#!/usr/bin/env Rscript

# Cross-mouse summary figures generated directly from the all-mouse save_list.
# This script does not require raw KPTracer allele tables.
#
# Figures generated here:
# 1. celltype_height_boxplot.pdf: UPGMA height by cell state.
# 2. UMAP_height2_all_cells.pdf: UMAP colored by UPGMA height.
# 3. Mesenchymal_zoom_UMAP.pdf: Mesenchymal cell-state UMAP zoom.
# 4. Early_EMT1_zoom_UMAP.pdf: Early EMT-1 UMAP zoom.
# 5. EMT_distribution_bymonth.pdf: Early EMT-1 height density by age.
# 6. Mesenchymal_distribution_bymonth.pdf: Early EMT-1 height histogram by age.
# 7. gene_height_correlation_examples.pdf: genes correlated with UPGMA height.

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_file) == 0L) getwd() else dirname(normalizePath(script_file))

save_list_file <- Sys.getenv("SAVE_LIST_FILE", unset = file.path(script_dir, "outputs", "save_list_all_mice_linemap.rds"))
out_dir <- Sys.getenv("SUMMARY_FIGURE_OUTPUT_DIR", unset = file.path(script_dir, "outputs", "cross_mouse_figures"))
kptracer_dir <- Sys.getenv("KPTRACER_DATA_DIR", unset = "")
umap_file <- Sys.getenv("UMAP_COORD_FILE", unset = "")
mouse_metadata_file <- Sys.getenv(
  "MOUSE_METADATA_FILE",
  unset = if (nzchar(kptracer_dir)) file.path(kptracer_dir, "trcr_master.txt") else ""
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_list <- readRDS(save_list_file)

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    as.data.frame(readRDS(path))
  } else if (ext == "csv") {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported table format: ", path)
  }
}

build_all_cells <- function(save_list) {
  all_cells <- do.call(rbind, lapply(names(save_list), function(mouse_id) {
    x <- save_list[[mouse_id]]
    tab <- x$table
    if (is.null(tab) || nrow(tab) == 0L) return(NULL)

    required <- c("cellID", "distancetoroot", "celltype", "height1", "height2", "height.hj")
    missing_cols <- setdiff(required, colnames(tab))
    if (length(missing_cols) > 0L) {
      stop("Mouse ", mouse_id, " table is missing columns: ", paste(missing_cols, collapse = ", "))
    }

    n_intBC <- x$meta$n_intBC
    data.frame(
      sample_key = mouse_id,
      mouse_id = mouse_id,
      MouseID = mouse_id,
      cellID = tab$cellID,
      tumor = if ("tumor" %in% colnames(tab)) as.character(tab$tumor) else NA_character_,
      clonal = if ("clonal" %in% colnames(tab)) tab$clonal else NA,
      celltype = as.character(tab$celltype),
      n_intBC = n_intBC,
      distance_to_root_raw = as.numeric(tab$distancetoroot),
      height1_raw = as.numeric(tab$height1),
      height2_raw = as.numeric(tab$height2),
      height_hj_raw = as.numeric(tab$height.hj),
      distance_to_root = as.numeric(tab$distancetoroot) / n_intBC,
      height1 = as.numeric(tab$height1) / n_intBC,
      height2 = as.numeric(tab$height2) / n_intBC,
      height.hj = as.numeric(tab$height.hj) / n_intBC,
      stringsAsFactors = FALSE
    )
  }))

  all_cells[!is.na(all_cells$celltype) & !is.na(all_cells$distance_to_root), , drop = FALSE]
}

merge_mouse_metadata <- function(all_cells, metadata_file) {
  if (!nzchar(metadata_file)) return(all_cells)
  if (!file.exists(metadata_file)) {
    message("Skipping age/genotype metadata merge because file was not found: ", metadata_file)
    return(all_cells)
  }
  metadata <- read_table_auto(metadata_file)

  mouse_col <- intersect(c("mouse_id", "MouseID", "mouseID", "Mouse", "mouse"), colnames(metadata))
  if (length(mouse_col) == 0L) {
    stop("Mouse metadata file must contain a mouse ID column, for example MouseID or mouse_id.")
  }
  colnames(metadata)[match(mouse_col[1], colnames(metadata))] <- "mouse_id"
  metadata$mouse_id <- as.character(metadata$mouse_id)
  metadata <- metadata[!duplicated(metadata$mouse_id), , drop = FALSE]

  keep <- unique(c("mouse_id", intersect(c("Aging_Month", "Aging_day", "Genotype", "ES_clone"), colnames(metadata))))
  merge(all_cells, metadata[, keep, drop = FALSE], by = "mouse_id", all.x = TRUE, sort = FALSE)
}

load_umap_coords <- function(umap_file) {
  if (!nzchar(umap_file)) return(NULL)
  umap <- read_table_auto(umap_file)
  id_col <- intersect(c("cellID", "Cells", "cell_id", "id"), colnames(umap))
  x_col <- intersect(c("X", "UMAP_1", "UMAP1", "x"), colnames(umap))
  y_col <- intersect(c("Y", "UMAP_2", "UMAP2", "y"), colnames(umap))

  if (length(id_col) == 0L || length(x_col) == 0L || length(y_col) == 0L) {
    stop("UMAP coordinate file must contain cell ID, X, and Y columns.")
  }

  data.frame(
    cellID = as.character(umap[[id_col[1]]]),
    X = as.numeric(umap[[x_col[1]]]),
    Y = as.numeric(umap[[y_col[1]]]),
    stringsAsFactors = FALSE
  )
}

load_umap_from_main_inputs <- function(kptracer_dir) {
  expression_index <- file.path(kptracer_dir, "Expression_data_index.xlsx")
  if (!file.exists(expression_index)) {
    message("Skipping UMAP figures because Expression_data_index.xlsx was not found in: ", kptracer_dir)
    return(NULL)
  }

  kp_gene_ids <- as.data.frame(read_excel(expression_index, sheet = "KP_gene_IDs"))
  umap <- as.data.frame(read_excel(expression_index, sheet = "UMAP"))
  colnames(umap) <- c("id", "X", "Y")

  cell_col <- if ("Cells" %in% colnames(kp_gene_ids)) "Cells" else colnames(kp_gene_ids)[8]
  if (nrow(kp_gene_ids) != nrow(umap)) {
    stop("KP_gene_IDs and UMAP sheets have different row counts in: ", expression_index)
  }
  data.frame(
    cellID = as.character(kp_gene_ids[[cell_col]]),
    X = as.numeric(umap$X),
    Y = as.numeric(umap$Y),
    stringsAsFactors = FALSE
  )
}

plot_celltype_height_boxplot <- function(all_cells, outfile) {
  colors <- c(
    "#F8766D", "#E9842C", "#D69100", "#BC9D00", "#9CA700",
    "#6FB000", "#00B813", "#00BD61", "#00C08E", "#00C0B4",
    "#00BDD4", "#00B5EE", "#00A7FF", "#7F96FF", "#BC81FF",
    "#E26EF7", "#F763DF", "#FF62BF", "#FF6A9A"
  )
  cellstate <- c(
    "Apc Early", "Apc Mesenchymal-1 (Met)", "Apc Mesenchymal-2",
    "AT1-like", "AT2-like", "Early EMT-1", "Early EMT-2",
    "Early gastric", "Endoderm-like", "Gastric-like", "High plasticity",
    "Late Gastric", "Lkb1 subtype", "Lung progenitor-like",
    "Mesenchymal-1", "Mesenchymal-1 (Met)", "Mesenchymal-2",
    "Mesenchymal-2 (Met)", "Pre-EMT"
  )

  all_cells$celltype <- as.factor(all_cells$celltype)
  fill_values <- colors[match(levels(all_cells$celltype), cellstate)]
  names(fill_values) <- levels(all_cells$celltype)

  p <- ggplot(all_cells, aes(x = reorder(celltype, -height2, FUN = mean), y = height2, fill = celltype)) +
    geom_boxplot() +
    scale_fill_manual(values = fill_values, na.value = "gray70") +
    geom_jitter(color = "black", size = 0.01, alpha = 0.3) +
    scale_y_continuous(limits = c(0, 20)) +
    labs(x = "cell stage (ordered by mean)", y = "Height") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank()
    )

  ggsave(outfile, p, width = 8, height = 6)
  p
}

height_color_creat <- function(height) {
  breaks <- quantile(height, seq(0, 1, 0.001), na.rm = TRUE)
  ii <- cut(height, breaks = unique(breaks), include.lowest = TRUE)
  ii <- droplevels(ii)
  colorRampPalette(c("red", "yellow", "green"))(length(levels(ii)))[ii]
}

plot_umap_highlight <- function(all_cells, umap, highlight.cell, heights,
                                outfile, x_range = NULL, y_range = NULL,
                                cex = 2, title = "Heights Plot") {
  xlim_val <- if (!is.null(x_range)) x_range else range(umap$X, na.rm = TRUE)
  ylim_val <- if (!is.null(y_range)) y_range else range(umap$Y, na.rm = TRUE)

  plot_core <- function() {
    col_height <- rep(NA_character_, nrow(umap))
    col_height1 <- height_color_creat(as.numeric(heights))
    idx <- match(highlight.cell, umap$cellID)
    keep <- !is.na(idx)
    col_height[idx[keep]] <- col_height1[keep]
    col_height[is.na(col_height)] <- "gray90"

    plot(
      umap$X, umap$Y,
      pch = ".",
      col = col_height,
      cex = c(cex, 0.1)[as.factor(col_height == "gray90")],
      xlim = xlim_val,
      ylim = ylim_val,
      main = title
    )
  }

  pdf(outfile, width = 7, height = 7)
  plot_core()
  dev.off()

  invisible(TRUE)
}

plot_selected_umap_figures <- function(all_cells, umap, out_dir) {
  if (is.null(umap)) {
    message("Skipping UMAP figures because UMAP coordinates are not available.")
    return(invisible(NULL))
  }

  plot_umap_highlight(
    all_cells = all_cells,
    umap = umap,
    highlight.cell = all_cells$cellID,
    heights = all_cells$height2,
    outfile = file.path(out_dir, "UMAP_height2_all_cells.pdf")
  )

  tmp <- all_cells[all_cells$celltype %in% c("Mesenchymal-1", "Mesenchymal-1 (Met)"), , drop = FALSE]
  plot_umap_highlight(
    all_cells = all_cells,
    umap = umap,
    highlight.cell = tmp$cellID,
    heights = tmp$height2,
    outfile = file.path(out_dir, "Mesenchymal_zoom_UMAP.pdf"),
    x_range = c(0, 3.7),
    y_range = c(9.7, 14),
    cex = 3
  )

  tmp <- all_cells[all_cells$celltype %in% "Early EMT-1", , drop = FALSE]
  plot_umap_highlight(
    all_cells = all_cells,
    umap = umap,
    highlight.cell = tmp$cellID,
    heights = tmp$height2,
    outfile = file.path(out_dir, "Early_EMT1_zoom_UMAP.pdf"),
    x_range = c(-2, 2.5),
    y_range = c(4.5, 7.5),
    cex = 3
  )
}

plot_age_genotype_figures <- function(all_cells, out_dir) {
  required <- c("Genotype", "Aging_Month")
  if (!all(required %in% colnames(all_cells))) {
    message("Skipping age/genotype figures because Genotype and Aging_Month are not available.")
    return(invisible(NULL))
  }

  tmp <- all_cells[all_cells$celltype %in% "Early EMT-1", , drop = FALSE]
  tmp1 <- tmp[tmp$Genotype == "sgNT", , drop = FALSE]
  tmp1 <- tmp1[!is.na(tmp1$height2) & !is.na(tmp1$Aging_Month), , drop = FALSE]
  if (nrow(tmp1) == 0L) {
    message("Skipping age/genotype figures because no Early EMT-1 sgNT cells were found.")
    return(invisible(NULL))
  }

  tmp1$Aging_Month <- as.factor(tmp1$Aging_Month)
  p_emt <- ggplot(tmp1, aes(x = height2, color = Aging_Month, fill = Aging_Month)) +
    geom_density(alpha = 0.2) +
    coord_cartesian(xlim = c(0, 15)) +
    theme_minimal() +
    labs(
      title = "Density of height by Age(months)",
      x = "height",
      y = "Density",
      color = "Age (Months)",
      fill = "Age (Months)"
    )
  ggsave(file.path(out_dir, "EMT_distribution_bymonth.pdf"), p_emt, width = 8, height = 8)

  p_hist <- ggplot(tmp1, aes(x = height2, fill = Aging_Month)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.4, position = "identity") +
    theme_minimal() +
    labs(
      title = "Histogram Density of height by Age(months)",
      x = "height",
      y = "Density",
      fill = "Age (Months)"
    )
  ggsave(file.path(out_dir, "Mesenchymal_distribution_bymonth.pdf"), p_hist, width = 8, height = 8)
}

plot_gene_height_correlation_examples <- function(all_cells, kptracer_dir, out_dir) {
  if (!nzchar(kptracer_dir)) {
    message("Skipping gene-height correlation figure because KPTRACER_DATA_DIR was not set.")
    return(invisible(NULL))
  }

  expression_file <- file.path(kptracer_dir, "expression", "adata_processed.nt.h5ad")
  if (!file.exists(expression_file)) {
    message("Skipping gene-height correlation figure because file was not found: ", expression_file)
    return(invisible(NULL))
  }
  if (!requireNamespace("zellkonverter", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE) ||
      !requireNamespace("mgcv", quietly = TRUE)) {
    message("Skipping gene-height correlation figure because zellkonverter/SummarizedExperiment/Matrix/mgcv is not available.")
    return(invisible(NULL))
  }

  meta <- all_cells[all_cells$celltype %in% c("Mesenchymal-1", "Mesenchymal-1 (Met)"), , drop = FALSE]
  if ("Genotype" %in% colnames(meta)) {
    meta <- meta[meta$Genotype == "sgNT", , drop = FALSE]
  }
  meta <- meta[!is.na(meta$height2) & meta$height2 < 10, , drop = FALSE]
  if (nrow(meta) < 10L) {
    message("Skipping gene-height correlation figure because too few cells passed filtering.")
    return(invisible(NULL))
  }

  adata <- zellkonverter::readH5AD(expression_file)
  assay_name <- if ("X" %in% SummarizedExperiment::assayNames(adata)) "X" else SummarizedExperiment::assayNames(adata)[1]
  expr <- SummarizedExperiment::assay(adata, assay_name)

  common_cells <- intersect(meta$cellID, colnames(expr))
  if (length(common_cells) < 10L) {
    message("Skipping gene-height correlation figure because too few cells matched the expression matrix.")
    return(invisible(NULL))
  }
  meta <- meta[match(common_cells, meta$cellID), , drop = FALSE]
  expr <- expr[, common_cells, drop = FALSE]

  expressed_counts <- Matrix::rowSums(expr > 0)
  candidate_genes <- names(sort(expressed_counts, decreasing = TRUE))[seq_len(min(2000L, length(expressed_counts)))]
  candidate_genes <- candidate_genes[!is.na(candidate_genes)]
  if (length(candidate_genes) == 0L) {
    message("Skipping gene-height correlation figure because no candidate genes were found.")
    return(invisible(NULL))
  }

  cor_res <- apply(as.matrix(expr[candidate_genes, , drop = FALSE]), 1, function(x) {
    test <- suppressWarnings(cor.test(as.numeric(x), meta$height2, method = "spearman", exact = FALSE))
    c(cor = unname(test$estimate), pval = test$p.value)
  })
  cor_df <- as.data.frame(t(cor_res))
  cor_df$gene <- rownames(cor_df)
  cor_df$padj <- p.adjust(cor_df$pval, method = "BH")

  sig <- cor_df[cor_df$padj < 0.05 & abs(cor_df$cor) > 0.2, , drop = FALSE]
  if (nrow(sig) > 0L) {
    top_pos <- head(sig[order(sig$cor, decreasing = TRUE), "gene"], 3L)
    top_neg <- head(sig[order(sig$cor), "gene"], 3L)
    selected_genes <- unique(c(top_pos, top_neg))
  } else {
    selected_genes <- head(cor_df[order(abs(cor_df$cor), decreasing = TRUE), "gene"], 6L)
  }
  selected_genes <- selected_genes[selected_genes %in% rownames(expr)]
  if (length(selected_genes) == 0L) {
    message("Skipping gene-height correlation figure because no genes were selected.")
    return(invisible(NULL))
  }

  plot_df <- do.call(rbind, lapply(selected_genes, function(gene) {
    data.frame(
      gene = gene,
      height2 = meta$height2,
      gene_exp = as.numeric(expr[gene, ]),
      stringsAsFactors = FALSE
    )
  }))

  label_df <- cor_df[match(selected_genes, cor_df$gene), , drop = FALSE]
  plot_df$gene_label <- factor(
    plot_df$gene,
    levels = selected_genes,
    labels = sprintf("%s (rho=%.2f)", label_df$gene, label_df$cor)
  )

  p <- ggplot(plot_df, aes(x = height2, y = gene_exp)) +
    geom_point(alpha = 0.2, size = 0.3, color = "black") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red") +
    scale_x_reverse() +
    facet_wrap(~gene_label, scales = "free_y", ncol = 3) +
    theme_minimal() +
    labs(
      title = "Correlation between cell height and gene expression",
      x = "Cell Height (reverse)",
      y = "Gene Expression"
    )

  ggsave(file.path(out_dir, "gene_height_correlation_examples.pdf"), p, width = 9, height = 6)
  invisible(TRUE)
}

all_cells <- build_all_cells(save_list)
all_cells <- merge_mouse_metadata(all_cells, mouse_metadata_file)

plot_celltype_height_boxplot(
  all_cells,
  file.path(out_dir, "celltype_height_boxplot.pdf")
)

umap <- if (nzchar(umap_file)) load_umap_coords(umap_file) else load_umap_from_main_inputs(kptracer_dir)
plot_selected_umap_figures(all_cells, umap, out_dir)
plot_age_genotype_figures(all_cells, out_dir)
plot_gene_height_correlation_examples(all_cells, kptracer_dir, out_dir)

message("Wrote cross-mouse summary figures to: ", out_dir)
