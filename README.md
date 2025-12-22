# LINEMAP: 
### A Temporal Dynamics Framework for Single-Cell Lineage Reconstruction: Integrating Molecular Barcoding and Advanced Machine Learning

### Introduction:

This pipeline provides an end-to-end workflow for processing CRISPR-based lineage tracing sequencing data generated from the A21B25 PyMT-Cas9 tracing cell line, with the goal of reconstructing mutation evolution networks for two homing guide RNA (hgRNA) barcodes.
Starting from raw paired-end FASTQ files, the pipeline performs sequential read filtering and quality control, followed by amplicon-aware alignment and mutation calling using amplican. Identified mutations are then annotated at the barcode level and aggregated to infer cell-level barcode identities, enabling robust tracking of hgRNA mutation states across individual cells. Based on these cell-resolved mutation profiles, the pipeline constructs directed mutation evolution networks that model irreversible and accumulative mutation dynamics of the two homing guide barcodes.
This workflow is specifically designed to support downstream probabilistic modeling of mutation transitions and time-scaled lineage reconstruction, serving as the data preprocessing backbone for the LINEMAP framework.

### 1. Required R packages:

```
install.packages("devtools")
devtools::install_github("hrbrmstr/waffle")

install.packages("ggplot2")
install.packages("dplyr")
install.packages("RecordLinkage")
install.packages("stringr")
install.packages('igraph')

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("amplican")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
```


### Pipeline command:
```
Rscript LINEMAP_processing.R <R1_FASTQ> <R2_FASTQ>
```

### For sample cases:
```
Rscript LINEMAP_processing.R samples/sample_R1.fastq.gz samples/sample_R2.fastq.gz
```
#### Sample data download links:
```
https://ibrisk.houstonmethodist.org/LINEMAP/samples/sample_R1.fastq.gz
https://ibrisk.houstonmethodist.org/LINEMAP/samples/sample_R2.fastq.gz
```
