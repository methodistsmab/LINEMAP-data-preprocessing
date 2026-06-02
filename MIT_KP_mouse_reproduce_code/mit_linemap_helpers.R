suppressPackageStartupMessages({
  library(ape)
  library(dendextend)
  library(readxl)
  library(LINEMAP)
})

default_kptracer_dir <- Sys.getenv("KPTRACER_DATA_DIR", unset = "")
helper_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
helper_dir <- if (!is.na(helper_file)) dirname(helper_file) else getwd()
mouse_config_static_candidates <- c(
  file.path(helper_dir, "mouse_configs_static.R"),
  file.path(getwd(), "updated_linemap_pipeline", "mouse_configs_static.R"),
  file.path(getwd(), "mouse_configs_static.R")
)
mouse_config_static_file <- mouse_config_static_candidates[file.exists(mouse_config_static_candidates)][1]
if (!is.na(mouse_config_static_file)) source(mouse_config_static_file)

load_kptracer_tables <- function(kptracer_dir = default_kptracer_dir) {
  if (!nzchar(kptracer_dir)) {
    stop("Please set KPTRACER_DATA_DIR to the extracted KPTracer-Data directory.")
  }
  if (!dir.exists(kptracer_dir)) {
    stop("KPTRACER_DATA_DIR does not exist: ", kptracer_dir)
  }
  data <- read.csv(file.path(kptracer_dir, "KPTracer_meta.csv"))
  data.mutation <- read.delim(file.path(kptracer_dir, "KPTracer.alleleTable.unfiltered.txt"))
  data.plasticity <- read.delim(file.path(kptracer_dir, "plasticity_scores.tsv"))
  data.celltype <- as.data.frame(read_excel(file.path(kptracer_dir, "KP_mice_metadata.xlsx"), sheet = "KP_mice_metadata"))
  colnames(data.celltype)[1] <- "X"

  fitness_dir <- file.path(kptracer_dir, "fitnesses")
  fitness_files <- list.files(fitness_dir, full.names = TRUE)
  data.fitness <- do.call(rbind, lapply(fitness_files, read.delim))

  list(
    data = data,
    data.mutation = data.mutation,
    data.plasticity = data.plasticity,
    data.celltype = data.celltype,
    data.fitness = data.fitness
  )
}

mouse_config <- function(mouse_id, data.celltype = NULL, n_intBC = NULL) {
  if (exists("mouse_configs_static", inherits = TRUE) && mouse_id %in% names(mouse_configs_static)) {
    return(mouse_configs_static[[mouse_id]])
  }

  stop(
    "No mouse-specific config found for mouse `", mouse_id, "`. ",
    "Add this mouse to `updated_linemap_pipeline/mouse_configs_static.R` ",
    "with explicit sub_tumors, colors, n_intBC, clonal_threshold, ",
    "clonal_height_threshold, and clonal_plot_index before rebuilding."
  )
}

create_sequence_matrix <- function(sample.cellBC, data.mutation.filter, intBC) {
  sequence <- matrix(NA_character_, nrow = length(sample.cellBC), ncol = 3L * length(intBC))
  rownames(sequence) <- sample.cellBC
  colnames(sequence) <- as.vector(rbind(
    paste0(intBC, "_r1"),
    paste0(intBC, "_r2"),
    paste0(intBC, "_r3")
  ))

  for (i in seq_along(sample.cellBC)) {
    cell_mut <- data.mutation.filter[data.mutation.filter$cellBC == sample.cellBC[i], ]
    if (nrow(cell_mut) == 0L) next
    for (j in seq_len(nrow(cell_mut))) {
      int_idx <- match(cell_mut$intBC[j], intBC)
      if (is.na(int_idx)) next
      sequence[i, 3L * int_idx - 2L] <- as.character(cell_mut$r1[j])
      sequence[i, 3L * int_idx - 1L] <- as.character(cell_mut$r2[j])
      sequence[i, 3L * int_idx] <- as.character(cell_mut$r3[j])
    }
  }
  sequence
}

build_tabs_for_intBC <- function(transition_inputs, n_intBC) {
  base_tabs <- transition_inputs$dictionary
  if (is.null(base_tabs)) {
    base_tabs <- transition_inputs$tabs
  }
  unlist(replicate(n_intBC, list(base_tabs$r1, base_tabs$r2, base_tabs$r3), simplify = FALSE), recursive = FALSE)
}

lentiBC_color_create <- function(sample.cellBC, lentiBC.list, d1) {
  lentiBC <- rep(NA_character_, length(sample.cellBC))
  for (i in seq_along(lentiBC.list)) {
    if (!is.na(lentiBC.list[i]) && nchar(lentiBC.list[i]) == 16L) {
      is_x <- which(sample.cellBC %in% d1$X[d1$lentiBC == lentiBC.list[i]])
      lentiBC[is_x] <- lentiBC.list[i]
    }
  }
  lentiBC <- factor(lentiBC)
  cols_l <- c("purple", "pink", "cyan", "yellow", "brown", "gray", "black", "magenta", "navy", "gold", "maroon")
  cols_l[seq_len(length(levels(lentiBC)))][lentiBC]
}

plasticity_color_create <- function(sample.cellBC, data.plasticity) {
  tmp <- data.plasticity[data.plasticity$X %in% sample.cellBC, ]
  if (nrow(tmp) == 0L) return(rep(NA_character_, length(sample.cellBC)))
  values <- tmp$scPlasticity[match(sample.cellBC, tmp$X)]
  ii <- cut(values, breaks = seq(min(tmp$scPlasticity), max(tmp$scPlasticity), length.out = 100L), include.lowest = TRUE)
  colorRampPalette(c("yellow", "red"))(length(levels(ii)))[ii]
}

celltype_color_create <- function(sample.cellBC, data.celltype) {
  tmp <- data.celltype[data.celltype$X %in% sample.cellBC, ]
  values <- factor(tmp$Cluster.Name[match(sample.cellBC, tmp$X)])
  colorRampPalette(c("blue", "red"))(length(levels(values)))[values]
}

fitness_color_create <- function(sample.cellBC, data.fitness) {
  values <- data.fitness$mean_fitness[match(sample.cellBC, data.fitness$X)]
  ii <- cut(values, breaks = seq(min(data.fitness$mean_fitness), max(data.fitness$mean_fitness), length.out = 100L), include.lowest = TRUE)
  colorRampPalette(c("red", "yellow"))(length(levels(ii)))[ii]
}

branch_node_find <- function(clusters1) {
  node.list <- list()
  merge <- clusters1$merge
  for (i in seq_len(nrow(merge))) {
    tmp <- c()
    stack <- c(merge[i, 1], merge[i, 2])
    while (length(stack) > 0L) {
      if (stack[1] < 0L) {
        tmp <- c(tmp, stack[1])
        stack <- stack[-1]
      } else {
        stack <- c(stack, merge[stack[1], 1], merge[stack[1], 2])
        stack <- stack[-1]
      }
    }
    node.list[[length(node.list) + 1L]] <- -tmp
  }
  node.list
}

leaf_height_get <- function(clusters1) {
  leaf.height <- rep(NA_real_, length(clusters1$order))
  for (i in seq_along(clusters1$height)) {
    left_samples <- clusters1$merge[i, 1]
    right_samples <- clusters1$merge[i, 2]
    if (left_samples < 0L) leaf.height[-left_samples] <- clusters1$height[i]
    if (right_samples < 0L) leaf.height[-right_samples] <- clusters1$height[i]
  }
  leaf.height
}

find_two_branch <- function(height, branch.node, l) {
  daughter.height <- c()
  daughter.node <- list()
  index <- c()
  if (l > 2L) {
    for (ll in (l - 1L):1L) {
      if (all(branch.node[[ll]] %in% branch.node[[l]])) {
        daughter.height[1] <- height[ll]
        daughter.node[[1]] <- branch.node[[ll]]
        index <- c(index, ll)
        break
      }
    }
    if (length(setdiff(branch.node[[l]], branch.node[[ll]])) == 1L) {
      daughter.height[2] <- height[l]
      daughter.node[[2]] <- setdiff(branch.node[[l]], branch.node[[ll]])
      index <- c(index, ll)
    } else {
      for (lll in (l - 1L):1L) {
        if (setequal(branch.node[[lll]], setdiff(branch.node[[l]], branch.node[[ll]]))) {
          daughter.height[2] <- height[lll]
          daughter.node[[2]] <- branch.node[[lll]]
          index <- c(index, lll)
          break
        }
      }
    }
  }
  list(daughter.height = daughter.height, daughter.node = daughter.node, branch.index = index)
}

get_clonal_branch <- function(height, branch.node, thres, thres.height = 5) {
  l <- length(height)
  branch.index <- c()
  for (i in which(height > thres.height)) {
    res <- find_two_branch(height, branch.node, i)
    if (length(intersect(unlist(branch.node[branch.index]), unlist(branch.node[[i]]))) == 0L) {
      branch.index <- c(branch.index, res$branch.index[which((height[i] - res$daughter.height) > thres)])
    }
  }
  remain <- setdiff(branch.node[[l]], unlist(branch.node[branch.index]))
  for (i in l:max(1L, l - 100L)) {
    if (all(branch.node[[i]] %in% remain)) {
      branch.index <- c(branch.index, i)
      remain <- setdiff(remain, branch.node[[i]])
      if (length(remain) == 0L) break
    }
  }
  branch.index
}

build_mouse_save_entry <- function(mouse_id, tables, transition_inputs, cfg = NULL) {
  if (is.null(cfg)) cfg <- mouse_config(mouse_id, tables$data.celltype)
  data.mutation <- tables$data.mutation
  data.celltype <- tables$data.celltype

  intBC <- names(sort(table(data.mutation$intBC[grepl(mouse_id, data.mutation$Tumor)]), decreasing = TRUE))[seq_len(cfg$n_intBC)]
  data.mutation.filter <- data.mutation[data.mutation$intBC %in% intBC, ]
  d1 <- data.celltype[grepl(mouse_id, data.celltype$Tumor), ]
  cellBC <- unique(data.mutation.filter$cellBC[grepl(mouse_id, data.mutation.filter$Tumor)])

  subgroup_cells <- lapply(cfg$sub_tumors, function(sub_tumor) {
    intersect(cellBC, d1[d1$SubTumor == sub_tumor, "X"])
  })
  if (!is.null(cfg$max_cells_per_subgroup)) {
    subgroup_cells <- lapply(subgroup_cells, function(cells) {
      if (length(cells) <= cfg$max_cells_per_subgroup) {
        cells
      } else {
        cells[sample.int(length(cells), cfg$max_cells_per_subgroup)]
      }
    })
  }
  sample.cellBC <- unlist(subgroup_cells, use.names = FALSE)
  sequence <- create_sequence_matrix(sample.cellBC, data.mutation.filter, intBC)
  tabs_list <- build_tabs_for_intBC(transition_inputs, length(intBC))

  D <- compute_prob_distance(
    sequence = sequence,
    tabs_list = tabs_list,
    gRNA.num = ncol(sequence),
    barcode.label = "NC",
    add_C0 = TRUE,
    missing_label = "MISSING",
    fast = TRUE,
    show_progress = FALSE
  )

  distance_probability_matrix <- D[-1, -1, drop = FALSE]
  distance_to_root <- D[1, -1]
  colnames(distance_probability_matrix) <- rownames(distance_probability_matrix) <- sample.cellBC
  names(distance_to_root) <- sample.cellBC

  subgroup_cell_ids <- unlist(subgroup_cells, use.names = FALSE)
  if (!identical(sample.cellBC, subgroup_cell_ids)) {
    stop("Internal cell order mismatch while building tumor labels for mouse `", mouse_id, "`.")
  }
  cell2T <- setNames(rep(names(subgroup_cells), lengths(subgroup_cells)), subgroup_cell_ids)
  tumor_color <- setNames(unname(rep(cfg$colors, lengths(subgroup_cells))), subgroup_cell_ids)
  clustercut <- unname(tumor_color)

  dist_object <- dist(distance_probability_matrix)
  hclust1 <- hclust(dist_object, method = "average")
  dend1 <- as.dendrogram(hclust1)
  old_labels <- labels(dend1)
  labels(dend1) <- unname(cell2T[old_labels])
  labels_colors(dend1) <- clustercut[order.dendrogram(dend1)]
  dend1 <- hang.dendrogram(dend1, hang = 0.01)

  distance_for_asdist <- distance_probability_matrix
  diag(distance_for_asdist) <- 0
  asdist_matrix <- as.dist(distance_for_asdist)
  hclust2 <- hclust(asdist_matrix, method = "average")
  dend2 <- as.dendrogram(hclust2)
  old_labels <- labels(dend2)
  labels(dend2) <- unname(cell2T[old_labels])
  labels_colors(dend2) <- clustercut[order.dendrogram(dend2)]
  dend2 <- hang.dendrogram(dend2, hang = 0.01)

  diag(D) <- 0
  colnames(D) <- rownames(D) <- c("C0", sample.cellBC)
  tree_nj <- build_tree_from_distance(D, method = "NJ")
  tip_height <- node.depth.edgelength(tree_nj)[seq_len(Ntip(tree_nj))]
  names(tip_height) <- tree_nj$tip.label

  branch.node <- branch_node_find(hclust1)
  clonal.index <- get_clonal_branch(hclust1$height, branch.node, cfg$clonal_threshold, cfg$clonal_height_threshold)
  clonal <- rep(NA_integer_, length(sample.cellBC))
  for (i in seq_along(clonal.index)) {
    clonal[branch.node[clonal.index][[i]]] <- i
  }

  table <- data.frame(
    cellID = sample.cellBC,
    height1 = leaf_height_get(hclust1),
    height2 = leaf_height_get(hclust2),
    height.hj = tip_height,
    distancetoroot = distance_to_root,
    tumor = unname(cell2T[sample.cellBC]),
    celltype = data.celltype$Cluster.Name[match(sample.cellBC, data.celltype$X)],
    clonal = clonal
  )

  lentiBC.list <- unique(d1$lentiBC)

  list(
    meta = list(
      intBC = intBC,
      n_cells = length(sample.cellBC),
      n_intBC = length(intBC),
      subgroups = subgroup_cells
    ),
    sample = sample.cellBC,
    sequence = sequence,
    cellBC = sample.cellBC,
    distance_probability_matrix = distance_probability_matrix,
    distance_to_root = distance_to_root,
    dist_object = dist_object,
    hclust1 = hclust1,
    dend1 = dend1,
    asdist_matrix = asdist_matrix,
    hclust2 = hclust2,
    dend2 = dend2,
    D = D,
    tree_nj = tree_nj,
    table = table,
    colors = list(
      tip_col = tumor_color,
      fitness = fitness_color_create(sample.cellBC, tables$data.fitness),
      plasticity = plasticity_color_create(sample.cellBC, tables$data.plasticity),
      celltype = celltype_color_create(sample.cellBC, tables$data.celltype),
      lentiBC = lentiBC_color_create(sample.cellBC, lentiBC.list, d1)
    )
  )
}
