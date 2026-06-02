
# A Temporal Dynamics Framework for Single-Cell Lineage Reconstruction: Integrating Molecular Barcoding and Advanced Machine Learning

## Introduction

This pipeline provides an end-to-end workflow for processing CRISPR-based lineage tracing sequencing data generated from the A21B25 PyMT-Cas9 tracing cell line. The goal is to identify cell-level hgRNA barcode states and reconstruct mutation evolution networks for two homing guide RNA (hgRNA) barcodes, A21 and B25.

Starting from raw paired-end FASTQ files, the pipeline performs read filtering and quality control, followed by amplicon-aware alignment and mutation calling using `amplican`. Mutations are annotated at the barcode level and aggregated into cell-level barcode identities. These cell-resolved mutation profiles are then used to construct directed mutation evolution networks that model irreversible and accumulative hgRNA mutation dynamics.

The preprocessing output can also be used for downstream LINEMAP modeling, including two-time-point estimation of the hgRNA transition matrix `P`.

## 0. Use pipepline code

Download or clone this repository, then enter the source-code folder:

```bash
git clone https://github.com/methodistsmab/LINEMAP-data-preprocessing.git
cd LINEMAP-data-preprocessing/A21B25_hgRNA_lineage_pipeline_code
```

All commands below should be run from `A21B25_hgRNA_lineage_pipeline_code`.


## 1. R package installation

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

## 2. Run The FASTQ Preprocessing Pipeline

Use paired-end FASTQ files as input:

```bash
Rscript LINEMAP_processing.R <R1_FASTQ> <R2_FASTQ>
```

Example:

```bash
Rscript LINEMAP_processing.R samples/sample_R1.fastq.gz samples/sample_R2.fastq.gz
```

Sample FASTQ files:

```text
https://ibrisk.houstonmethodist.org/LINEMAP/samples/sample_R1.fastq.gz
https://ibrisk.houstonmethodist.org/LINEMAP/samples/sample_R2.fastq.gz
```

## 3. Main Outputs

After running `LINEMAP_processing.R`, results are saved in the `result/` folder.

Important output folders:

```text
result/cell_barcode/
result/mutation_network/
result/mutation_annotation/
result/alignment/
```

The main final outputs are:

```text
result/cell_barcode/A21/cell_level_barcode.rds
result/cell_barcode/B25/cell_level_barcode.rds
result/mutation_network/A21/nodeA21.tsv
result/mutation_network/A21/edgesA21.tsv
result/mutation_network/B25/nodeB25.tsv
result/mutation_network/B25/edgesB25.tsv
```

The `node` files contain mutation states, and the `edges` files contain directed mutation evolution relationships.

## 4. Estimate Transition Matrix P From Two Time Points

This step is used when two time-point LINEMAP preprocessing results are available, for example day 5 and day 14.

To generate two-time-point inputs from FASTQ files, run the preprocessing
pipeline separately for each time point. `LINEMAP_processing.R` writes results
to a folder named `result/`, so rename the result folder after each run:

```bash
Rscript LINEMAP_processing.R d5_R1.fastq.gz d5_R2.fastq.gz
mv result result_d5

Rscript LINEMAP_processing.R d14_R1.fastq.gz d14_R2.fastq.gz
mv result result_d14
```

Then estimate the transition matrix for one target using the generated result
folders:

```bash
Rscript estimate_transition_P.R \
  --early_result result_d5 \
  --late_result result_d14 \
  --target A21 \
  --selection min_cells \
  --min_cells 10 \
  --output result_transition_P/A21
```

Preprocessed day 5 and day 14 result folders may also be provided with this
repository. In that case, users do not need to rerun the FASTQ pipeline and can
directly run:

```bash
Rscript estimate_transition_P.R \
  --early_file data/cell.read.d5.RData \
  --early_object cell.read.d5A21 \
  --late_file data/cell.read.d14.RData \
  --late_object cell.read.d14A21 \
  --target A21 \
  --selection min_cells \
  --min_cells 10 \
  --output result_transition_P/A21
```

Main outputs:

```text
result_transition_P/A21/transition_node_A21.tsv
result_transition_P/A21/transition_edges_A21.tsv
result_transition_P/A21/transition_network_pie_A21.pdf
result_transition_P/A21/P_transition_fsolve_A21.tsv
result_transition_P/A21/P_transition_fsolve_A21.rds
result_transition_P/A21/P_fsolve_fit_check_A21.tsv
```

`P_transition_fsolve_A21.tsv` is the final estimated transition matrix.

## 5. Notes

- This pipeline is designed for A21B25 PyMT-Cas9 hgRNA lineage tracing data.
- The main preprocessing workflow starts from paired-end FASTQ files and produces cell-level barcodes and mutation evolution networks.
- The transition matrix step is an advanced downstream step for two-time-point data.
- By default, mutation states are retained if they contain at least `10` cells in at least one time point.
