# iSCORE-PD: scRNA-seq of Irradiated WIBR3 hESCs

Analysis code for single-cell RNA sequencing of irradiated wild-type WIBR3 human embryonic stem cells (hESCs), part of the [iSCORE-PD](https://doi.org/10.1101/2024.02.12.579917) isogenic stem cell collection for Parkinson's disease research.

## Overview

WIBR3 hESCs were exposed to five radiation doses (0, 0.5, 2, 5, and 10 Gy) and multiplexed using MULTI-seq barcoding for pooled single-cell RNA sequencing on 10x Chromium 3' v3.1 (NovaSeq 6000). This repository contains downstream analysis scripts for Seurat preprocessing, differential expression, and visualization of p53 DNA damage response gene activation across radiation doses.

## Upstream Processing (Not Included)

MULTI-seq demultiplexing and sample assignment were performed prior to the analysis in this repository. Raw FASTQ files for both the gene expression (GEX) and MULTI-seq barcode libraries, along with the processed Seurat object and barcode classification table, are available on the **ASAP CRN Cloud** under dataset `team_rio_hesc_sc_rnaseq_irradiated`.

MULTI-seq barcode-to-dose assignments:
| Barcode | Dose |
|---------|------|
| GGAGAAGA | 0 Gy (CTRL) |
| CCACAATG | 0.5 Gy |
| TGAGACCT | 2 Gy |
| GCACACGC | 5 Gy |
| AGAGAGAG | 10 Gy |

## Repository Structure

```
scripts/
  R/
    1-create_seurat.R           # QC filtering, normalization, PCA, clustering, UMAP
    2-DEAnalysis.Rmd            # Differential expression: CTRL vs each dose (Wilcoxon/Bonferroni)
    3-BoxplotGenerator.Rmd      # Boxplots of 9 p53 target genes across radiation doses
```

## Analysis Pipeline

1. **Preprocessing** (`1-create_seurat.R`) — Seurat v4 object creation from 10x GEX count matrix, QC filtering (percent.mt < 15%), LogNormalize (scale factor 10,000), 2,000 variable features (VST), PCA, clustering (10 dims, resolution 0.5), UMAP
2. **Differential expression** (`2-DEAnalysis.Rmd`) — `FindMarkers()` (Wilcoxon rank sum test, Bonferroni correction) comparing CTRL to each radiation dose (0.5, 2, 5, 10 Gy), outputting per-comparison CSV tables
3. **Visualization** (`3-BoxplotGenerator.Rmd`) — Boxplots of 9 p53 DNA damage response target genes (GDF15, FDXR, DDB2, RPS27L, PHPT1, BAX, CDKN1A, PHLDA3, GADD45A) showing dose-dependent expression, saved as PNG and EPS

## Dependencies

- R 4.x
- Seurat v4
- ggplot2, cowplot, dplyr

## Data Availability

Raw sequencing data (GEX and MULTI-seq FASTQs), processed Seurat objects, count matrices, differential expression statistics, and visualization outputs are deposited in the ASAP CRN Cloud platform under dataset `team_rio_hesc_sc_rnaseq_irradiated`.

## Citation

Busquets O, Li H, Syed KM, et al. "iSCORE-PD: an isogenic stem cell collection to research Parkinson Disease." *bioRxiv* (2024). DOI: [10.1101/2024.02.12.579917](https://doi.org/10.1101/2024.02.12.579917)

## Funding

ASAP-024409, ASAP-000486 — Aligning Science Across Parkinson's (ASAP) Collaborative Research Network

