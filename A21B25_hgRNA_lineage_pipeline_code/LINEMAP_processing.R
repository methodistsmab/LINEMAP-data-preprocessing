## ------------------------------------------------------------
## Set working directory
##
## Users should set the working directory to the root folder
## of the project before running the pipeline.
## All file paths in the script are defined relative to this
## project root.
##
## NOTE:
## Do NOT hard-code absolute paths. Replace the path below
## with the local project directory on your system.
## ------------------------------------------------------------
# setwd("/home/yxh/data/LINEMAP")

## ------------------------------------------------------------
## Load required libraries before printing pipeline messages
##
suppressWarnings(suppressPackageStartupMessages({
  library(ggplot2)
  library(ShortRead)
  library("dplyr")
  library(Biostrings)
  library(amplican)
  library("stringr")
  library('RecordLinkage')
}))

suppressWarnings(suppressPackageStartupMessages({
  source("utility_functions.R")
  source("amplican_alignment.R")
  source("mutation_annotation.R")
  source("mutation_hgRNA_identify.R")
}))

message("Current working directory:")
message(getwd())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide exactly two arguments: <R1_FASTQ> <R2_FASTQ>")
}
fq_r1 <- args[1]
fq_r2 <- args[2]

if (file.exists(fq_r1)){
  message("R1 FASTQ file found: ", fq_r1)
} else {
  stop (paste("R1 FASTQ file not found:", fq_r1))
}

if (file.exists(fq_r2)){
  message("R2 FASTQ file found: ", fq_r2)
} else {
  stop (paste("R2 FASTQ file not found:", fq_r2))
}

message("Reading input FASTQ files...")


## ------------------------------------------------------------
## Read paired-end FASTQ files for CRISPR barcode sequencing
##
## R1: cell barcode (CB) + UMI
## R2: hgRNA spacer + PAM + scaffold
##
## For reproducibility and quick verification:
## - The full FASTQ files (Day 5 / Day 14) are used for the
##   actual analyses reported in the manuscript.
## - A small example dataset (sample_R1/sample_R2) is provided
##   to allow users to quickly verify that the pipeline runs
##   end-to-end on their system (a "dry run"/"sanity check").
##   This example dataset is NOT intended for biological inference.
## ------------------------------------------------------------

## =========================
## Full dataset (Day 5)
## =========================
# Paired-end sequencing reads collected at Day 5
# These files contain tens of millions of reads and are
# computationally expensive to process.


# <FASTQ_DIR>//D5_CRISPR_S5_L003_R1_001.fastq.gz
# <FASTQ_DIR>//D5_CRISPR_S5_L003_R2_001.fastq.gz


## =========================
## Full dataset (Day 14)
## =========================
# Paired-end sequencing reads collected at Day 14
# Used for mutation frequency estimation at the later time point.

# <FASTQ_DIR>//D14_CRISPR_S6_L003_R1_001.fastq.gz
# <FASTQ_DIR>//D14_CRISPR_S6_L003_R2_001.fastq.gz

## =========================
## Subsampled dataset (for testing)
## =========================
# A small subset of reads extracted from the original FASTQ files.
# This sample preserves the exact read structure (R1/R2 format),
# but dramatically reduces file size and runtime.


# <FASTQ_DIR>//sample_R1.fastq.gz
# <FASTQ_DIR>//sample_R2.fastq.gz



## ------------------------------------------------------------
## FASTQ filtering and preprocessing
##
## IMPORTANT:
## - fq1_in corresponds to sequencing Read 2 (R2),
##   which contains the hgRNA spacer, PAM, and scaffold.
## - fq2_in corresponds to sequencing Read 1 (R1),
##   which contains the cell barcode (CB) and UMI.
##
## The naming (fq1 / fq2) reflects the internal design of
## filter_fastq_stream(), not the Illumina R1/R2 convention.
## ------------------------------------------------------------

## Input FASTQ files
# fq_r2 <- file.path("samples", "sample_R2.fastq.gz")  # sequencing R2
# fq_r1 <- file.path("samples", "sample_R1.fastq.gz")  # sequencing R1
message("Input R2 FASTQ: ", fq_r2)
message("Input R1 FASTQ: ", fq_r1)

## Output directory (relative to project root)
output_dir <- file.path(getwd(), "result")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 1 output directory
filtered_dir <- file.path( output_dir,"read_filtered_data")
dir.create(filtered_dir, showWarnings = FALSE, recursive = TRUE)

## Output FASTQ files (filtered)
fq_r2_filt <- file.path(filtered_dir, "R2.filtered.fastq.gz")
fq_r1_filt <- file.path(filtered_dir, "R1.filtered.fastq.gz")


res <- filter_fastq_stream(
  fq1_in  = fq_r2,        # input FASTQ: sequencing R2 (hgRNA spacer + PAM + scaffold)
  fq1_out = fq_r2_filt,   # output FASTQ: filtered R2 reads
  fq2_in  = fq_r1,        # input FASTQ: sequencing R1 (cell barcode + UMI)
  fq2_out = fq_r1_filt,   # output FASTQ: filtered R1 reads (paired with R2)
  scaffold_seq = "TTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTC", # hgRNA scaffold reference sequence
  primer_seq   = "AAGCAGTGGTATCAACGCAGAGTACATGGG",        # sequencing primer reference sequence
  k = 10,                 # minimum scaffold match length (bp)
  meanQ_min = 20,         # minimum mean Phred quality score per read
  chunk = 5e5,            # number of reads processed per streaming chunk
  verbose_every = 1       # progress reporting frequency (in chunks)
)

## targets
targets <- c("A21", "B25")

## hgRNA spacer sequences (named by target ID)
GuideRNAs <- c(
  A21 = "GTTCCCGTCCAGTAATCGTG",
  B25 = "GTCGTTGTAGCAACCTATCGGGTG"
)

Forward_Primer = 'AAGCAGTGGTATCAACGCAGAGTACATGGG'  # forward primer sequence

Amplicon = c('AAGCAGTGGTATCAACGCAGAGTACATGGGGTTCCCGTCCAGTAATCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTC', # A21
             'AAGCAGTGGTATCAACGCAGAGTACATGGGGTCGTTGTAGCAACCTATCGGGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCT') # B25
Amplicon_random = c('AAGCAGTGGTATCAACGCAGAGTACATGGGGTTCCCGTCCAGTAATCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTT', # A21
             'AAGCAGTGGTATCAACGCAGAGTACATGGGGTCGTTGTAGCAACCTATCGGGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCC') # B25

amplican_alignment(targets, GuideRNAs, Forward_Primer, Amplicon, Amplicon_random)
mutation_annotation(targets, GuideRNAs)

run_step4_cell_level_barcodes(len_gap = 4)

mutation_hgRNA_identify(targets, GuideRNAs, enable_plot =TRUE)









