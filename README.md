# LINEMAP data preprocessing

### A Temporal Dynamics Framework for Single-Cell Lineage Reconstruction: Integrating Molecular Barcoding and Advanced Machine Learning

## Introduction

This pipeline provides an end-to-end workflow for processing CRISPR-based lineage tracing sequencing data generated from the A21B25 PyMT-Cas9 tracing cell line. The goal is to identify cell-level hgRNA barcode states and reconstruct mutation evolution networks for two homing guide RNA (hgRNA) barcodes, A21 and B25.

Starting from raw paired-end FASTQ files, the pipeline performs read filtering and quality control, followed by amplicon-aware alignment and mutation calling using `amplican`. Mutations are annotated at the barcode level and aggregated into cell-level barcode identities. These cell-resolved mutation profiles are then used to construct directed mutation evolution networks that model irreversible and accumulative hgRNA mutation dynamics.

The preprocessing output can also be used for downstream LINEMAP modeling, including two-time-point estimation of the hgRNA transition matrix `P`.

## 1. Installation

Install the required R packages:

```r
install.packages("devtools")
devtools::install_github("hrbrmstr/waffle")

install.packages("ggplot2")
install.packages("dplyr")
install.packages("RecordLinkage")
install.packages("stringr")
install.packages("igraph")
install.packages("pracma")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("amplican")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
```

## 2. Download The Code

Download or clone this repository, then enter the source-code folder:

```bash
git clone https://github.com/methodistsmab/LINEMAP-data-preprocessing.git

```


