
#########################################
## ------------------------------------------------------------
## Amplicon configuration for amplican
##
## This run uses the example FASTQ data subsampled from the D5 time point.
## We perform the amplican alignment/indel calling separately for the two
## hgRNA targets (A21 and B25). The corresponding guide spacers and amplicon
## reference sequences are specified below.
##
## NOTE:
## - To analyze additional time points (e.g., D14), update `Group` and use
##   the matching FASTQ inputs.
## ------------------------------------------------------------

amplican_alignment <- function(targets, GuideRNAs, Forward_Primer, Amplicon, Amplicon_random){

Group = 'd5'   # time point label for this run (can also be set to 'd14')

# targets <- c("A21", "B25")
# ######################xiaohui
# GuideRNAs <- c(
#   A21 = "GTTCCCGTCCAGTAATCGTG",            # hgRNA spacer: A21
#   B25 = "GTCGTTGTAGCAACCTATCGGGTG"         # hgRNA spacer: B25
# )

# Forward_Primer = 'AAGCAGTGGTATCAACGCAGAGTACATGGG'  # forward primer sequence

# Amplicon = c('AAGCAGTGGTATCAACGCAGAGTACATGGGGTTCCCGTCCAGTAATCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTC', # A21
#              'AAGCAGTGGTATCAACGCAGAGTACATGGGGTCGTTGTAGCAACCTATCGGGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCT') # B25



#########################################
## ------------------------------------------------------------
## Load filtered FASTQ files and extract per-read information
##
## R1 reads are used to extract cell barcodes (first 16 bp),
## which define cell identities and enable per-cell grouping.
## R2 reads contain the hgRNA sequences and are used for
## downstream alignment and mutation (indel) analysis.
## ------------------------------------------------------------

message("===== Step 2: Amplican processing =====")
message("Running amplican to align reads to the reference amplicon and detect indels and mismatches.")


step1_dir <- file.path( "result", "read_filtered_data")

fqr1 <- readFastq(file.path(step1_dir, "R1.filtered.fastq.gz"))
fqr2 <- readFastq(file.path(step1_dir, "R2.filtered.fastq.gz"))


barcodelist = substring(as.character(sread(fqr1)),1,16)

umilist     <- substring(as.character(sread(fqr1)),17,28)

tb = table(barcodelist)

names.cell.barcode <- names(tb)                  
names.cell.barcode <- names.cell.barcode[names.cell.barcode != ""]   # 防御性处理
length(names.cell.barcode)

indexlist = match(barcodelist,names.cell.barcode)

read.string = as.character(sread(fqr2))

#########################################
## ------------------------------------------------------------
## Main alignment loop
##
## The pipeline runs amplican separately for each hgRNA target
## (e.g., A21 and B25). For each target, all outputs are written
## to a dedicated subdirectory under result/alignment/<target>/.
##
## Within each target, alignment and indel calling are performed
## on a per-cell basis to preserve single-cell resolution.
## ------------------------------------------------------------

## Unified output root directory
output_dir <- "result"

## Alignment output root
alignment_root <- file.path(output_dir, "alignment")

print (alignment_root)

dir.create(alignment_root, showWarnings = FALSE, recursive = TRUE)



for (t in seq(length(targets))){
  
  message("Processing target: ", targets[t])
  
  target_id <- targets[t]
  
  print (paste0("Processing target: ", target_id))

  ## Current target label (e.g., "A21" or "B25")
  ## (Assumes you set `target_id` inside the loop)
  alignment_dir <- file.path(alignment_root, target_id)
  dir.create(alignment_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## Combined alignment events TSV (recommended for downstream analysis)
  out_alignment <- file.path(alignment_dir, "alignment_all.tsv")
  if (file.exists(out_alignment)) file.remove(out_alignment)  # start fresh for this target
  
  ## Event-level table (one row per indel/mutation event)
  out_events <- file.path(alignment_dir, "events_all.tsv")
  if (file.exists(out_events)) file.remove(out_events)  # start fresh for this target
  
  ## (Optional) Save per-cell AlignmentsExperimentSet objects for traceability
  save_rds_per_cell <- TRUE
  rds_dir <- file.path(alignment_dir, "aes_by_cell")
  if (save_rds_per_cell) dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  ## Temporary files/folders used during alignment for this target
  tmp_fastq <- file.path(alignment_dir, "tmp.read.fastq.gz")
  tmp_cfg   <- file.path(alignment_dir, "tmp.config.csv")
  tmp_res   <- file.path(alignment_dir, "results_folder_tmp")
  
  
  ## ----------------------------------------------------------
  ## Per-cell alignment loop
  ##
  ## Each cell is processed independently:
  ## 1) extract hgRNA reads belonging to the cell
  ## 2) run amplican on the cell-specific FASTQ
  ## 3) extract and append alignment / indel events
  ## ----------------------------------------------------------
  
  cell.number <- length(names.cell.barcode)
  print ("Total cell number:")
  print (cell.number)
  for (i in seq(cell.number)) {
    
    index.tmp <- which(indexlist == i)
    all.sequence <- read.string[index.tmp]

    if (length(all.sequence) == 0) next
   
    
    ## Require the presence of the forward primer sequence
    if (!any(grepl(Forward_Primer, all.sequence, fixed = TRUE))) next
    
    ## Guard against cells with only a single unique read
    ## (required by amplican to perform alignment)
    unique.all.sequence <- unique(
      all.sequence[grepl(Forward_Primer, all.sequence, fixed = TRUE)]
    )
    
    if (length(unique.all.sequence) == 1) {
      if (unique.all.sequence == Amplicon[t]) {
        all.sequence <- c(all.sequence, Amplicon_random[t])
      }
      if (unique.all.sequence == Amplicon[t]) {
        all.sequence <- c(all.sequence, Amplicon_random[t])
      }
    }
    
    
    ## Construct temporary FASTQ records for the current cell
    d <- length(all.sequence)
    id <- paste0("@", names.cell.barcode[i], "_", seq_len(d))
    
    tmp.read = read.construct(id,all.sequence)
    
    
    if (file.exists(tmp_fastq)) file.remove(tmp_fastq)
    write.table(tmp.read, file = tmp_fastq, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    
    ## Cell identifier
    ID <- names.cell.barcode[i]
    Forward_Reads <- basename(tmp_fastq)
    
    ## Generate amplican configuration for this cell and target
    tmp.config <- make.config(ID, Forward_Reads, Group, GuideRNAs[t],
                              Forward_Primer, Amplicon[t])
    
    if (file.exists(tmp_cfg)) file.remove(tmp_cfg)
    write.csv(tmp.config, file = tmp_cfg, quote = FALSE, row.names = FALSE)
    
    if (dir.exists(tmp_res)) unlink(tmp_res, recursive = TRUE)
    dir.create(tmp_res, showWarnings = FALSE)
    
    ## --------------------------------------------------------
    ## Run amplican (single cell, single target)
    ## Console output is suppressed to keep logs clean.
    ## --------------------------------------------------------
    zz <- file(tempfile(), open = "wt")
    sink(zz)
    sink(zz, type = "message")
    
    tryCatch({
      suppressWarnings(
        suppressMessages(
          
          amplicanPipeline(
            config = tmp_cfg,
            fastq_folder = alignment_dir,
            results_folder = tmp_res,
            fastqfiles = 1,
            knit_reports = FALSE
          )
          
        )
      )
    }, error = function(e) {
      message("amplicaN error after alignments, ignored")
      print (e)
    })
    
    sink(type = "message")
    sink()
    close(zz)
    
    ## --------------------------------------------------------
    ## Collect alignment results and export events
    ## --------------------------------------------------------
    
    aes_path <- file.path(tmp_res, "alignments", "AlignmentsExperimentSet.rds")
    if (!file.exists(aes_path)) {
      message("Missing AlignmentsExperimentSet.rds for cell ", ID, " group ", Group)
      next
    }
    
    aln <- readRDS(aes_path)
    
    append_alignment_blocks(
      aes = aln,
      cell_id = ID,
      out_file = out_alignment,
      group_id = Group,
      max_show = Inf,          
      include_score = TRUE    
    )
    
    ## (Optional) Save per-cell RDS for traceability
    if (save_rds_per_cell) {
      saveRDS(aln, file = file.path(rds_dir, paste0(ID, "_", Group, ".rds")))
    }
    
    ## Extract structured indel events and append to global table
    ev <- extractEvents(aln)
    ev$CellID <- ID
    ev$Group  <- Group
    
    append_tsv(ev, out_events)
    
    if (i %% 100 == 0) message("Processed cells: ", i, " / ", cell.number)
  }
}
}


