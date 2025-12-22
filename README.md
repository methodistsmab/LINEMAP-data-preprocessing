# LINEMAP: 
### A Temporal Dynamics Framework for Single-Cell Lineage Reconstruction: Integrating Molecular Barcoding and Advanced Machine Learning

### Instruction:

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


### The command:
Rscript data_preprocessing.R <R1_FASTQ> <R2_FASTQ>

### For sample cases:
Rscript data_preprocessing.R samples/sample_R1.fastq.gz samples/sample_R2.fastq.gz
