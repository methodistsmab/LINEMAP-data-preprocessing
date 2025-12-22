#########################################
## ------------------------------------------------------------
## Purpose:
## For each read (row), compare its annotation results against
## two candidate hgRNA targets (A21 and B25) and determine the
## most likely target of origin.
##
## Multiple quality and sequence features are integrated to
## assign each read to A21, B25, or mark it as unassigned when
## the evidence is insufficient or conflicting.
##
## Read-level assignments are then consolidated at the cell level
## to derive a unique barcode (mutation allele) per cell and per
## target.
##
## Outputs are saved as cell-level barcode tables under:
##   result/cell_barcode/<target>/
## These tables serve as the direct input for downstream
## mutation evolution network construction.
## ------------------------------------------------------------


# library(dplyr)

# res = run_step4_cell_level_barcodes(len_gap = 4)



#############################################
## ------------------------------------------------------------
##   We construct a structured mutation description matrix by 
##   collapsing cell-level alignment results into unique mutation 
##   alleles, each defined by its read sequence and mutation 
##   report (encoding mutation type, position, and affected bases).
##
##   Based on this structured allele representation and 
##   guide-specific mutation constraints, we then infer directed 
##   relationships between mutation alleles using 
##   `make.network.direction`, yielding a directed edge list that 
##   defines the mutation evolution network.
##
##   The resulting node table (alleles with support counts) and 
##   directed edges jointly constitute the mutation evolution 
##   network used in downstream analyses.
## ------------------------------------------------------------

library('igraph')
library(dplyr)

mutation_hgRNA_identify <- function(targets, GuideRNAs, enable_plot = TRUE) {
  # Your function implementation here

  # ---- 0) Read Step 4 result from disk ----
  base_dir <- file.path("result", "cell_barcode")
  # targets  <- c("A21", "B25")

  ## hgRNA spacer sequences (named by target ID)
  # GuideRNAs <- c(
  #   A21 = "GTTCCCGTCCAGTAATCGTG",
  #   B25 = "GTCGTTGTAGCAACCTATCGGGTG"
  # )

  mutation.network.plot <- enable_plot

  out_dir <- file.path("result", "mutation_network")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


  message("===== Step 5: Mutation evolution network =====")
  message("Constructing the mutation evolution network by collapsing cell-level mutations into alleles and inferring directed ancestral relationships.")


  for ( target in targets){
    
    message("\n===== Step5: ", target, " =====")
    

    cell_level_barcode <- readRDS(file.path(base_dir, target, "cell_level_barcode.rds"))
    stopifnot(is.data.frame(cell_level_barcode))
    
    data <- cell_level_barcode
    
    # ---- (1) filtering / QC ----
    data <- subset(data, number.of.notinread == 0)
    data <- subset(data, !is.na(align.short))
    
    message("[", target, "] Unique cells: ", nrow(data))
    
    # ---- (2) collapse to allele-level ----
    allele_tbl <-
      data %>%
      group_by(read.short, align.short, mutation.report) %>%
      summarise(n = dplyr::n(), .groups = "drop") %>%
      as.data.frame()
    
    allele_tbl <- allele_tbl[order(allele_tbl$n, decreasing = TRUE), ]
    message("[", target, "] Collapsed alleles: ", nrow(allele_tbl))
    
    # ---- (3) build node table ----
    # stable unique id
    allele_tbl$allele_id <- allele_tbl$read.short

    # Your current preference: node = read.short; node1 = mutation.report
    node <- allele_tbl[, c("allele_id", "read.short", "align.short", "mutation.report", "n")]
    colnames(node) = c("node", "read.short", "align.short", "mutation.report", "n")

    # ---- (4) choose top nodes covering 90% of total reads ----
    cum_prop <- cumsum(node$cell.number) / sum(node$cell.number)
    idx <- which(cum_prop >= 0.9)
    top_n <- if (length(idx) == 0) nrow(node) else min(idx)
    
    node_top <- node[seq_len(top_n), , drop = FALSE]
    
    # message("[", target, "] Top nodes covering >=90% counts: ", nrow(node_top))
    
    
    # ---- (5) infer directed edges among alleles (mutation.report-level) ----
    # make.network.direction is assumed to take a node table containing mutation.report in column "node1"
    # and return edges with from/to referring to mutation.report strings.
    network.direction <- make.network.direction(node_top, GuideRNAs[[target]])
    edges <- as.data.frame(network.direction)
    
    # sanity: require endpoints
    stopifnot(all(c("from", "to") %in% names(edges)))
    
    # Preserve original endpoints (mutation.report)
    edges$mutation.report.from <- edges$from
    edges$mutation.report.to   <- edges$to
    
    # Map mutation.report -> node (read.short) for visualization
    map_rep_to_node <- setNames(node_top$node, node_top$mutation.report)
    edges$from <- unname(map_rep_to_node[edges$mutation.report.from])
    edges$to   <- unname(map_rep_to_node[edges$mutation.report.to])
    
    # Drop edges that couldn't be mapped (shouldn't happen if consistent, but safe)
    edges <- edges[!is.na(edges$from) & !is.na(edges$to), , drop = FALSE]
    
    # Keep a clean edge table for downstream usage
    if (!("score" %in% names(edges))) edges$score <- NA_real_
    edges <- edges[, c("from", "mutation.report.from", "to", "mutation.report.to", "score"), drop = FALSE]


    # ---- (6) save outputs ----
    # ensure per-target output directory exists
    target_dir <- file.path(out_dir, target)
    dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)
    
    saveRDS(node, file.path(target_dir, paste0("node", target, ".rds")))
    saveRDS(edges, file.path(target_dir, paste0("edges", target, ".rds")))
    
    # optional: also save as TSV for inspection
    write.table(node,
                file = file.path(target_dir, paste0("node", target, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(edges,
                file = file.path(target_dir,  paste0("edges", target, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # ---- (7) optional plot ----
    if (isTRUE(mutation.network.plot)) {
      
      plot_file_pdf  <- file.path(target_dir, paste0("mutation_evolution_network_", target, ".pdf"))
      plot_file_tiff <- file.path(target_dir, paste0("mutation_evolution_network_", target, ".tiff"))
      
      pdf(plot_file_pdf, width = 8, height = 8)
      plot_mutation_network(
        node = node_top,
        edges = edges,
        node_id_col = "node",
        node_size_col = "n"
      )
      dev.off()
      
      ## ---- TIFF (raster, journal-ready) ----
      tiff(plot_file_tiff,
          width = 8,
          height = 8,
          units = "in",
          res = 600)
      plot_mutation_network(
        node = node_top,
        edges = edges,
        node_id_col = "node",
        node_size_col = "n"
      )
      dev.off()
      
    }
    
    
  }

}
  

