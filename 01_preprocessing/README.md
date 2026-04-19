# Pre-processing of 10X Single-Cell RNA Datasets

[![Tool](https://img.shields.io/badge/Tool-Scanpy%201.9+-blue?style=flat-square)](https://scanpy.readthedocs.io/)
[![QC](https://img.shields.io/badge/QC-Scrublet%20%7C%20Mitochondrial-green?style=flat-square)](https://github.com/swolock/scrublet)
[![Data](https://img.shields.io/badge/Input-10X%20Cell%20Ranger%20MTX-orange?style=flat-square)](https://www.10xgenomics.com/)
[![Output](https://img.shields.io/badge/Output-.h5ad%20AnnData-purple?style=flat-square)](https://anndata.readthedocs.io/)
[![License](https://img.shields.io/badge/License-CC--BY--4.0-lightgrey?style=flat-square)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Faiqa Zarar
**Institution:** NUST — National University of Sciences and Technology
**Course:** Bioinformatics
**Date:** April 2026

---

## Abstract

This section documents the complete pre-processing pipeline applied to raw 10X Genomics Chromium single-cell RNA sequencing data. Starting from Cell Ranger output (sparse matrix format: `.mtx`, `barcodes.tsv`, `features.tsv`), the pipeline performs six sequential steps: data loading, quality control metric calculation, low-quality cell filtering, doublet detection and removal using Scrublet, library-size normalization with log1p transformation, and highly variable gene selection. The output is a clean, normalized AnnData `.h5ad` object ready for downstream dimensionality reduction and clustering analysis. All steps are implemented in Python using Scanpy and follow current best practices as outlined by Luecken & Theis (2019).

**Keywords:** 10X Genomics, Cell Ranger, quality control, normalization, doublet detection, Scrublet, highly variable genes, AnnData, Scanpy

---

## Table of Contents

1. [Overview](#1-overview)
2. [Folder Structure](#2-folder-structure)
3. [Data Acquisition](#3-data-acquisition)
4. [Pipeline Steps](#4-pipeline-steps)
   - 4.1 [Load 10X Data](#41-load-10x-data)
   - 4.2 [Quality Control](#42-quality-control)
   - 4.3 [Cell Filtering](#43-cell-filtering)
   - 4.4 [Doublet Detection](#44-doublet-detection)
   - 4.5 [Normalization](#45-normalization)
   - 4.6 [Feature Selection](#46-feature-selection)
   - 4.7 [Save Output](#47-save-output)
5. [Results](#5-results)
6. [References](#6-references)

---

## 1. Overview

Raw 10X Genomics scRNA-seq data contains thousands of droplet barcodes, many of which do not correspond to real cells (empty droplets, doublets, or low-quality cells). Pre-processing removes these artifacts and normalizes expression values so that cells are comparable to each other. This pipeline transforms raw Cell Ranger output into a clean, analysis-ready AnnData object.

```
Cell Ranger Output (MTX)
        │
        ▼
┌───────────────────┐
│  Load Data        │  sc.read_10x_mtx()
└────────┬──────────┘
         │
         ▼
┌───────────────────┐
│  Quality Control  │  MT%, n_genes, total_counts
└────────┬──────────┘
         │
         ▼
┌───────────────────┐
│  Cell Filtering   │  Remove low-quality cells
└────────┬──────────┘
         │
         ▼
┌───────────────────┐
│  Doublet Removal  │  Scrublet
└────────┬──────────┘
         │
         ▼
┌───────────────────┐
│  Normalization    │  CPM + log1p
└────────┬──────────┘
         │
         ▼
┌───────────────────┐
│  Feature Select   │  Highly Variable Genes
└────────┬──────────┘
         │
         ▼
  preprocessed_adata.h5ad
```

---

## 2. Folder Structure

```
01_preprocessing/
├── README.md                    ← This document
├── scripts/
│   ├── quality_control.py       # QC metric calculation and violin plots
│   ├── normalization.py         # Library-size normalization + log1p
│   ├── feature_selection.py     # Highly variable gene identification
│   └── doublet_detection.py     # Scrublet-based doublet removal
└── outputs/
    ├── qc_metrics.csv           # Per-cell QC statistics table
    ├── filtered_matrix/         # Post-QC sparse count matrix
    └── preprocessed_adata.h5ad  # Final analysis-ready AnnData object
```

---

## 3. Data Acquisition

For this tutorial, we use the publicly available PBMC 3K dataset from 10X Genomics:

```
Dataset : PBMC 3K (Peripheral Blood Mononuclear Cells)
Source  : 10X Genomics
URL     : https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
Cells   : ~2,700 single cells
Reads   : ~70,000 reads per cell
Genome  : GRCh38 (human)
```

```bash
# Download and extract
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

**Input files (Cell Ranger output):**

| File | Description |
|---|---|
| `matrix.mtx` | Sparse count matrix (genes × cells) |
| `barcodes.tsv` | Cell barcode sequences |
| `features.tsv` | Gene IDs and gene symbols |

---

## 4. Pipeline Steps

### 4.1 Load 10X Data

```python
import scanpy as sc
import numpy as np
import pandas as pd

sc.settings.verbosity = 3
sc.settings.figdir = 'outputs/'

# Load from Cell Ranger output directory
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',
    cache=True
)

# Make gene names unique (some may be duplicated)
adata.var_names_make_unique()

print(adata)
# AnnData object with n_obs × n_vars = 2700 × 32738
```

---

### 4.2 Quality Control

Three standard QC metrics are computed per cell:

| Metric | Description | Biological Meaning |
|---|---|---|
| `n_genes_by_counts` | Number of genes detected per cell | Low = empty droplet or dead cell |
| `total_counts` | Total UMI counts per cell | Very high = potential doublet |
| `pct_counts_mt` | % reads mapping to mitochondrial genes | High = dying/damaged cell |

```python
# Flag mitochondrial genes (human: prefix "MT-")
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt'],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Visualize distributions
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save='_qc_before_filtering.png'
)

# Scatter plots to identify outliers
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_mt_scatter.png')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_genes_scatter.png')
```

---

### 4.3 Cell Filtering

```python
# Filter low-quality cells
sc.pp.filter_cells(adata, min_genes=200)    # Remove empty droplets
sc.pp.filter_genes(adata, min_cells=3)      # Remove unexpressed genes

# Remove cells with too many genes (likely doublets)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]

# Remove cells with high mitochondrial content (dead/dying cells)
adata = adata[adata.obs.pct_counts_mt < 5, :]

print(f"Cells remaining after filtering: {adata.n_obs}")
print(f"Genes remaining after filtering: {adata.n_vars}")
```

---

### 4.4 Doublet Detection

```python
import scrublet as scr
import scipy.sparse as sp

# Convert to dense for Scrublet if needed
counts_matrix = adata.X if not sp.issparse(adata.X) else adata.X.toarray()

# Run Scrublet
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2,
    min_cells=3,
    n_prin_comps=30
)

# Add results to AnnData
adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets

# Remove doublets
n_before = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']]
print(f"Doublets removed: {n_before - adata.n_obs}")
print(f"Cells remaining: {adata.n_obs}")
```

---

### 4.5 Normalization

```python
# Normalize total counts per cell to 10,000 (library-size normalization)
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform (stabilizes variance)
sc.pp.log1p(adata)

# Store raw (normalized) counts before further processing
adata.raw = adata

print("Normalization complete.")
```

---

### 4.6 Feature Selection

```python
# Identify highly variable genes (most informative for downstream analysis)
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

# Visualize HVG selection
sc.pl.highly_variable_genes(adata, save='_hvg.png')

print(f"Highly variable genes selected: {adata.var.highly_variable.sum()}")

# Filter to HVGs only
adata = adata[:, adata.var.highly_variable]

# Scale data (zero mean, unit variance) — clip at max_value=10
sc.pp.scale(adata, max_value=10)
```

---

### 4.7 Save Output

```python
# Export QC metrics
adata.obs.to_csv('outputs/qc_metrics.csv')

# Save final preprocessed AnnData
adata.write('outputs/preprocessed_adata.h5ad')

print("Pre-processing complete!")
print(adata)
```

---

## 5. Results

### 5.1 QC Summary

| Metric | Before Filtering | After Filtering |
|---|---|---|
| Total cells | 2,700 | ~2,638 |
| Total genes | 32,738 | ~13,714 |
| Median genes/cell | ~1,100 | ~1,200 |
| Cells removed (MT%) | — | ~30 |
| Doublets removed | — | ~32 |
| Highly variable genes | — | ~1,838 |

### 5.2 QC Plots

QC violin plots and scatter plots are saved to `outputs/`:

| Figure | Description |
|---|---|
| `qc_before_filtering.png` | Violin plots of QC metrics before filtering |
| `mt_scatter.png` | MT% vs total counts scatter |
| `genes_scatter.png` | n_genes vs total counts scatter |
| `hvg.png` | Highly variable gene dispersion plot |

---

## 6. References

1. Luecken, M. D., & Theis, F. J. (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*, 15(6), e8746. https://doi.org/10.15252/msb.20188746
2. Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15. https://doi.org/10.1186/s13059-017-1382-0
3. Wolock, S. L., Lopez, R., & Klein, A. M. (2019). Scrublet: computational identification of cell doublets. *Cell Systems*, 8(4), 281–291. https://doi.org/10.1016/j.cels.2018.11.005
4. 10X Genomics PBMC 3K Dataset. https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-v-1-2-1-2-standard-1-2-0

---

⬅️ [Back to Main README](../README.md) &nbsp;|&nbsp; ➡️ [Next: Basic Tutorial](../02_basic_scrna_tutorial/README.md)
