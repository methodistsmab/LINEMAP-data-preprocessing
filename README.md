# LINEMAP Data Preprocessing

This repository contains data preprocessing and reproducibility code for
lineage tracing analyses used with LINEMAP. The code is organized into two
independent modules, each with its own README and workflow instructions.

## Code Modules

### `A21B25_hgRNA_lineage_pipeline`

This module processes CRISPR-based lineage tracing sequencing data generated
from the A21B25 PyMT-Cas9 tracing cell line. Starting from paired-end FASTQ
files, it performs read filtering, amplicon-aware mutation calling, barcode
annotation, and cell-level hgRNA barcode assignment for the A21 and B25 hgRNA
barcodes. It also constructs mutation evolution networks and estimates the
transition matrix `P` for the A21 hgRNA barcode.

### `MIT_KP_mouse_reprodce_code`

This module reproduces the MIT/KPTracer mouse lineage tracing analysis using
public KPTracer data. It estimates sgRNA transition matrices, builds distance
dictionaries, generates one-mouse analysis objects, reconstructs lineage trees,
and produces downstream mouse-level and cross-mouse figures from either a
user-generated combined all-mouse `save_list` or the provided
`LINEMAP_MIT_processed_analysis_objects.RDS` object downloaded from Zenodo.

## Data Availability

The A21B25 module uses paired-end FASTQ data from the A21B25 PyMT-Cas9 tracing
cell line collected at two time points, d5 and d14. The full FASTQ data are
available from GEO:

```text
GSE335072
```

The Zenodo record below provides example FASTQ files for testing the A21B25
pipeline and the MIT/KPTracer processed all-mouse analysis object,
`LINEMAP_MIT_processed_analysis_objects.RDS`:

```text
https://doi.org/10.5281/zenodo.20576933
```

The MIT/KPTracer raw data are not redistributed in this repository. Users should
download the public processed KPTracer data from Zenodo:

```text
https://doi.org/10.5281/zenodo.5847461
```

## Usage

Please see the README file inside each module for detailed installation,
required inputs, command-line examples, and expected outputs.
