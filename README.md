# LINEMAP: 
### A Temporal Dynamics Framework for Single-Cell Lineage Reconstruction: Integrating Molecular Barcoding and Advanced Machine Learning

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
#### sample data download links:
```
https://ibrisk.houstonmethodist.org/LINEMAP/samples/sample_R1.fastq.gz
https://ibrisk.houstonmethodist.org/LINEMAP/samples/sample_R2.fastq.gz
```
