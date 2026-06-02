
library("stringr")
library('RecordLinkage')
source("utility_functions.R")

#########################################
## ------------------------------------------------------------
## Step: Per-read mutation annotation and feature extraction
##
## In this step, we convert raw amplicon-level alignments into
## read-level mutation representations. Specifically, for each
## aligned read we:
##   (i)  parse the spacer and scaffold regions,
##   (ii) correct alignment artifacts (e.g., leading-base shifts),
##   (iii) annotate all mutation events on the spacer sequence, and
##   (iv) compute quantitative mutation features (mismatch, indel
##        counts, lengths, and complexity scores).
##   (v) we construct a structured mutation description
##       matrix that explicitly encodes each mutation allele (including
##       its type, position, and affected bases).
##
## This matrix provides a unified and machine-readable representation of mutation states
## and serves as the core input for downstream mutation evolution
## network construction.
##
##
## Outputs are saved separately for each hgRNA target (A21, B25) under:
##   result/mutation_annotation/<target>/
## ------------------------------------------------------------

# targets <- c("A21", "B25")

# ## hgRNA spacer sequences (named by target ID)
# GuideRNAs <- c(
#   A21 = "GTTCCCGTCCAGTAATCGTG",
#   B25 = "GTCGTTGTAGCAACCTATCGGGTG"
# )

## hgRNA spacer sequences (named by target ID)
mutation_annotation <- function(targets, GuideRNAs) {
  message("===== Step 3: Mutation annotation =====")
  message("Annotating mutation types, positions, and affected bases from aligned reads.")

  ## Create a new folder under result/ to store Step 3 outputs
  step3_dir <- file.path("result", "mutation_annotation")
  dir.create(step3_dir, showWarnings = FALSE, recursive = TRUE)
  message("Step 3 output folder: ", step3_dir)

  for (target in targets) {
    
    message("Processing target: ", target)
    
    ## Input: alignment summary file produced in Step 2
    aln_file <- file.path("result", "alignment", target, "alignment_all.tsv")
    if (!file.exists(aln_file)) stop("Missing alignment file: ", aln_file)
    
    ## Output folder for this target
    out_dir <- file.path(step3_dir, target)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    message("[", target, "] Output folder: ", out_dir)
    
    ## Load alignment table (no header by design)
    cell.read.align <- read.delim(aln_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    ## Downstream processing
    cell.read.matrix <- make.align.matrix(cell.read.align)
    
    ## Annotate (split + optional correction + mutation.report)
    cell.read.matrix <- annotate_mutations(
      cell.read.matrix = cell.read.matrix,
      target_id = target,
      GuideRNA = GuideRNAs[target],
      verbose_every = 1000
    )
    
    ## Compute additional per-read metrics
    cell.read.matrix <- calculate_mutation_metrics(cell.read.matrix)
    
    cell.read.matrix <- add_spacer_similarity_features(cell.read.matrix, spacer_seq = GuideRNAs[target])
    
    # ## Build structured mutation representation (list of matrices)
    # mutation.report.matrix <- make.mutation.report.matrix(
    #   mutation.report = cell.read.matrix$mutation.report,
    #   GuideRNA = GuideRNAs[target]
    # )
    
    ## -------------------- Save outputs --------------------
    
    ## 1) Save cell.read.matrix (TSV + RDS)
    write.table(
      cell.read.matrix,
      file = file.path(out_dir, "cell_read_matrix.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    saveRDS(cell.read.matrix, file = file.path(out_dir, "cell_read_matrix.rds"))
    
    
    ## Optional: save minimal metadata
    meta <- list(
      target = target,
      guideRNA = GuideRNAs[target],
      input_alignment_file = aln_file
    )
    saveRDS(meta, file = file.path(out_dir, "meta.rds"))
  }
}


############################################
