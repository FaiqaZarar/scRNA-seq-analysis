# AnnData Workflows

[![Tool](https://img.shields.io/badge/Tool-AnnData%200.10+-blue?style=flat-square)](https://anndata.readthedocs.io/)
[![Source](https://img.shields.io/badge/Source-anndata.readthedocs.io-green?style=flat-square)](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
[![Source](https://img.shields.io/badge/Source-scverse--tutorials-orange?style=flat-square)](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)
[![Format](https://img.shields.io/badge/Format-.h5ad%20HDF5-purple?style=flat-square)](https://anndata.readthedocs.io/)
[![License](https://img.shields.io/badge/License-CC--BY--4.0-lightgrey?style=flat-square)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Faiqa Zarar
**Institution:** NUST — National University of Sciences and Technology
**Course:** Bioinformatics
**Date:** April 2026

---

## Abstract

This section provides hands-on tutorials for working with the AnnData (Annotated Data) object — the standard data structure for single-cell genomics in Python. Two tutorials are presented, based on the official documentation from anndata.readthedocs.io and scverse-tutorials.readthedocs.io. Tutorial 1 introduces the AnnData object from scratch: constructing an object from a simulated count matrix, adding cell-level and gene-level metadata via `.obs` and `.var`, storing multidimensional embeddings in `.obsm` and `.varm`, adding unstructured metadata via `.uns`, working with data layers in `.layers`, and reading and writing to disk in HDF5 format (`.h5ad`). Tutorial 2 covers more advanced operations including shared indexing via `obs_names` and `var_names`, pairwise cell matrices in `.obsp`, memory-efficient views versus in-memory copies, and multi-sample concatenation using `ad.concat()`. These tutorials establish foundational competency with AnnData required for all scRNA-seq analyses in the scverse ecosystem.

**Keywords:** AnnData, annotated data, obs, var, obsm, obsp, uns, layers, h5ad, HDF5, views, concatenation, scverse, single-cell

---

## Table of Contents

1. [What is AnnData?](#1-what-is-anndata)
2. [Folder Structure](#2-folder-structure)
3. [Tutorial 1 — AnnData Getting Started](#3-tutorial-1--anndata-getting-started)
   - 3.1 [Installation](#31-installation)
   - 3.2 [Creating an AnnData Object](#32-creating-an-anndata-object)
   - 3.3 [Adding Metadata — obs and var](#33-adding-metadata--obs-and-var)
   - 3.4 [Multidimensional Metadata — obsm and varm](#34-multidimensional-metadata--obsm-and-varm)
   - 3.5 [Unstructured Metadata — uns](#35-unstructured-metadata--uns)
   - 3.6 [Layers](#36-layers)
   - 3.7 [Subsetting AnnData](#37-subsetting-anndata)
   - 3.8 [Reading and Writing h5ad](#38-reading-and-writing-h5ad)
4. [Tutorial 2 — scverse AnnData Tutorial](#4-tutorial-2--scverse-anndata-tutorial)
   - 4.1 [Shared Indexing](#41-shared-indexing)
   - 4.2 [Pairwise Matrices — obsp](#42-pairwise-matrices--obsp)
   - 4.3 [Views vs Copies](#43-views-vs-copies)
   - 4.4 [Concatenating AnnData Objects](#44-concatenating-anndata-objects)
5. [AnnData Quick Reference](#5-anndata-quick-reference)
6. [Results](#6-results)
7. [References](#7-references)

---

## 1. What is AnnData?

AnnData is specifically designed for matrix-like data — n observations, each of which can be represented as a d-dimensional vector, where each dimension corresponds to a variable or feature. Both the rows and columns of this matrix are indexed, and additional metadata can be stored at both the observation and variable levels.

```
AnnData object with n_obs × n_vars

┌──────────────────────────────────────────────────────────┐
│                  .X  (n_obs × n_vars)                    │
│              Main count matrix (cells × genes)           │
├───────────────────────┬──────────────────────────────────┤
│  .obs  (n_obs rows)   │  Cell-level metadata             │
│                       │  e.g. cell_type, batch, n_genes  │
├───────────────────────┼──────────────────────────────────┤
│  .var  (n_vars rows)  │  Gene-level metadata             │
│                       │  e.g. gene_name, highly_variable │
├───────────────────────┴──────────────────────────────────┤
│  .obsm   → Multi-dim cell data (PCA, UMAP embeddings)    │
│  .varm   → Multi-dim gene data (PCA loadings)            │
│  .obsp   → Pairwise cell matrices (KNN graph, distances) │
│  .uns    → Unstructured metadata (dict: colors, params)  │
│  .layers → Additional matrices (raw, spliced, unspliced) │
│  .raw    → Stored raw count backup                       │
└──────────────────────────────────────────────────────────┘
```

---

## 2. Folder Structure

```
03_anndata_workflows/
├── README.md                            ← This document
├── notebooks/
│   ├── anndata_getting_started.ipynb    # Tutorial 1 (anndata.readthedocs.io)
│   └── anndata_scverse_tutorial.ipynb   # Tutorial 2 (scverse-tutorials)
└── outputs/
    ├── sample_anndata.h5ad              # AnnData created in Tutorial 1
    └── tutorial_outputs/
        ├── subset_adata.h5ad            # Subsetted AnnData
        └── concatenated_adata.h5ad      # Multi-batch concatenated AnnData
```

---

## 3. Tutorial 1 — AnnData Getting Started

> **Source:** https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html

### 3.1 Installation

```bash
pip install anndata numpy pandas
```

```python
import anndata as ad
import numpy as np
import pandas as pd

print(ad.__version__)
```

---

### 3.2 Creating an AnnData Object

```python
# Simulate a count matrix: 100 cells × 2000 genes
np.random.seed(42)
n_obs, n_vars = 100, 2000

# Poisson-distributed counts (realistic for scRNA-seq)
X = np.random.poisson(lam=5, size=(n_obs, n_vars)).astype(float)

# Create the AnnData object
adata = ad.AnnData(X=X)

# Assign cell names (obs_names) and gene names (var_names)
adata.obs_names = [f'Cell_{i}' for i in range(n_obs)]
adata.var_names = [f'Gene_{i}' for i in range(n_vars)]

print(adata)
# AnnData object with n_obs × n_vars = 100 × 2000
```

---

### 3.3 Adding Metadata — obs and var

`adata.obs` stores per-cell metadata (a Pandas DataFrame with one row per cell).
`adata.var` stores per-gene metadata (a Pandas DataFrame with one row per gene).

```python
# Add cell-level metadata to .obs
adata.obs['cell_type'] = pd.Categorical(
    np.random.choice(['B', 'T', 'Monocyte'], size=n_obs)
)
adata.obs['batch'] = pd.Categorical(
    np.random.choice(['batch1', 'batch2'], size=n_obs)
)
adata.obs['n_counts'] = X.sum(axis=1)

# Add gene-level metadata to .var
adata.var['gene_name'] = [f'Gene_{i}' for i in range(n_vars)]
adata.var['highly_variable'] = np.random.choice([True, False], size=n_vars)

print(adata.obs.head())
print(adata.var.head())
print(adata)
# AnnData object with n_obs × n_vars = 100 × 2000
# obs: 'cell_type', 'batch', 'n_counts'
# var: 'gene_name', 'highly_variable'
```

---

### 3.4 Multidimensional Metadata — obsm and varm

For embeddings (PCA, UMAP) and gene loadings that have more than one dimension per cell/gene, AnnData uses `.obsm` and `.varm`.

```python
# Add a 2D UMAP embedding (randomly simulated here)
adata.obsm['X_umap'] = np.random.normal(0, 1, size=(n_obs, 2))

# Add 50D PCA embedding
adata.obsm['X_pca'] = np.random.normal(0, 1, size=(n_obs, 50))

# Add gene-level metadata: PCA loadings (50 components × n_vars)
adata.varm['gene_stuff'] = np.random.normal(0, 1, size=(n_vars, 5))

print(adata.obsm)
# OverloadedDict with keys: ['X_umap', 'X_pca']

print(adata)
# obs: 'cell_type', 'batch', 'n_counts'
# var: 'gene_name', 'highly_variable'
# obsm: 'X_umap', 'X_pca'
# varm: 'gene_stuff'
```

---

### 3.5 Unstructured Metadata — uns

`.uns` stores any additional unstructured metadata as a Python dictionary — useful for storing analysis parameters, color palettes, or experimental information.

```python
# Add unstructured metadata
adata.uns['experiment'] = {
    'organism': 'Homo sapiens',
    'tissue': 'PBMC',
    'date': '2026-04-01',
    'instrument': '10X Chromium v3'
}

adata.uns['random'] = [1, 2, 3]

print(adata.uns['experiment'])
```

---

### 3.6 Layers

Layers store alternative versions of the count matrix (e.g., raw counts, log-normalized counts, spliced/unspliced reads). All layers share the same dimensions as `.X`.

```python
# Store raw counts before normalization
adata.layers['raw_counts'] = adata.X.copy()

# Store log-transformed counts
adata.layers['log_transformed'] = np.log1p(adata.X)

print(adata.layers['log_transformed'])

# AnnData now shows:
# layers: 'raw_counts', 'log_transformed'
print(adata)
```

---

### 3.7 Subsetting AnnData

Subsetting AnnData works like NumPy/pandas — rows are cells, columns are genes.

```python
# Subset by cell name
subset_cells = adata[['Cell_1', 'Cell_10', 'Cell_50']]
print(subset_cells)   # 3 × 2000

# Subset by boolean mask (only B cells)
b_cells = adata[adata.obs['cell_type'] == 'B']
print(b_cells)

# Subset by highly variable genes
hvg = adata[:, adata.var['highly_variable']]
print(hvg)

# Combined: B cells + HVGs only
subset = adata[adata.obs['cell_type'] == 'B', adata.var['highly_variable']]
print(subset)

# Save subset
subset.write('outputs/tutorial_outputs/subset_adata.h5ad')
```

---

### 3.8 Reading and Writing h5ad

AnnData uses the HDF5-based `.h5ad` format for efficient on-disk storage of all components.

```python
# Save to disk
adata.write('outputs/sample_anndata.h5ad')
print("Saved to outputs/sample_anndata.h5ad")

# Load from disk
adata_loaded = ad.read_h5ad('outputs/sample_anndata.h5ad')
print(adata_loaded)

# Verify all components preserved
print(adata_loaded.obs.columns.tolist())
print(list(adata_loaded.obsm.keys()))
print(list(adata_loaded.layers.keys()))
```

---

## 4. Tutorial 2 — scverse AnnData Tutorial

> **Source:** https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html

### 4.1 Shared Indexing

One important aspect of AnnData objects is the shared indexing. No matter if we look at a cell's raw counts, its cluster assignments, or its k-nearest neighbors in PCA space, we always want to be able to refer to the same cell with the same name. That is why all cell annotations share the same index: the obs_names. Similarly, var_names indexes genes and their annotations.

```python
# obs_names and var_names are the shared index
print(adata.obs_names[:5].tolist())
# ['Cell_0', 'Cell_1', 'Cell_2', 'Cell_3', 'Cell_4']

print(adata.var_names[:5].tolist())
# ['Gene_0', 'Gene_1', 'Gene_2', 'Gene_3', 'Gene_4']

# Index must be unique — check:
assert adata.obs_names.is_unique, "Cell names must be unique!"
assert adata.var_names.is_unique, "Gene names must be unique!"
```

---

### 4.2 Pairwise Matrices — obsp

`.obsp` stores pairwise relationships between cells, such as distance matrices or KNN connectivity graphs. These are always square matrices of shape `(n_obs, n_obs)`.

```python
import scipy.sparse as sp

# Simulate a KNN connectivity matrix (sparse)
n = adata.n_obs
knn_matrix = sp.random(n, n, density=0.05, format='csr')
knn_matrix = (knn_matrix + knn_matrix.T) / 2   # Make symmetric

adata.obsp['connectivities'] = knn_matrix
adata.obsp['distances'] = knn_matrix.copy()

print(adata.obsp['connectivities'].shape)   # (100, 100)
print(adata)
# obsp: 'connectivities', 'distances'
```

---

### 4.3 Views vs Copies

Technically, subsetting AnnData objects returns a view of the AnnData you are subsetting from, instead of making a copy. That avoids having the same data in memory twice.

```python
# Subsetting creates a VIEW (memory-efficient, read-only)
view = adata[:50]
print(view.is_view)   # True
print(view.n_obs)     # 50

# Views are read-only — modifying them converts to a copy automatically
# To safely modify, explicitly copy first:
copy = adata[:50].copy()
print(copy.is_view)           # False
copy.obs['new_col'] = 'value' # Safe to modify
```

---

### 4.4 Concatenating AnnData Objects

```python
# Simulate two batches
adata_b1 = ad.read_h5ad('outputs/sample_anndata.h5ad')
adata_b2 = ad.read_h5ad('outputs/sample_anndata.h5ad')

# Rename obs_names to avoid conflicts
adata_b2.obs_names = [f'Cell_b2_{i}' for i in range(adata_b2.n_obs)]

# Concatenate (outer join preserves all genes from both datasets)
adata_combined = ad.concat(
    [adata_b1, adata_b2],
    keys=['batch1', 'batch2'],
    label='batch',
    join='outer',
    fill_value=0
)

print(adata_combined)
# AnnData object with n_obs × n_vars = 200 × 2000
# obs: 'cell_type', 'batch_label', 'n_counts', 'batch'

adata_combined.write('outputs/tutorial_outputs/concatenated_adata.h5ad')
print("Concatenated AnnData saved.")
```

---

## 5. AnnData Quick Reference

| Attribute | Type | Shape | Description |
|---|---|---|---|
| `.X` | ndarray / sparse | (n_obs, n_vars) | Main count matrix |
| `.obs` | DataFrame | (n_obs, n_attrs) | Cell-level metadata |
| `.var` | DataFrame | (n_vars, n_attrs) | Gene-level metadata |
| `.obs_names` | Index | (n_obs,) | Unique cell identifiers |
| `.var_names` | Index | (n_vars,) | Unique gene identifiers |
| `.obsm` | dict of arrays | (n_obs, k) | Cell embeddings (PCA, UMAP) |
| `.varm` | dict of arrays | (n_vars, k) | Gene loadings |
| `.obsp` | dict of sparse | (n_obs, n_obs) | Pairwise cell matrices |
| `.uns` | dict | any | Unstructured metadata |
| `.layers` | dict of matrices | (n_obs, n_vars) | Alternative count matrices |
| `.raw` | Raw object | — | Backup of raw counts |
| `.n_obs` | int | — | Number of cells |
| `.n_vars` | int | — | Number of genes |

---

## 6. Results

### Output Files

| File | Format | Description |
|---|---|---|
| `sample_anndata.h5ad` | H5AD | Full AnnData from Tutorial 1 (100 cells × 2000 genes) |
| `subset_adata.h5ad` | H5AD | B cell subset with highly variable genes |
| `concatenated_adata.h5ad` | H5AD | Two-batch concatenated AnnData (200 cells × 2000 genes) |

### AnnData Components Summary

| Component | Content Added in Tutorial |
|---|---|
| `.X` | Simulated Poisson count matrix (100 × 2000) |
| `.obs` | `cell_type`, `batch`, `n_counts` |
| `.var` | `gene_name`, `highly_variable` |
| `.obsm` | `X_umap` (2D), `X_pca` (50D) |
| `.varm` | `gene_stuff` (5D) |
| `.obsp` | `connectivities`, `distances` |
| `.uns` | `experiment` dict, `random` list |
| `.layers` | `raw_counts`, `log_transformed` |

---

## 7. References

1. Virshup, I., et al. (2023). The scverse project provides a computational ecosystem for single-cell omics data analysis. *Nature Biotechnology*, 41, 604–606. https://doi.org/10.1038/s41587-023-01733-8
2. Virshup, I., et al. (2024). anndata: Annotated data. *Journal of Open Source Software*, 9(101), 4371. https://doi.org/10.21105/joss.04371
3. AnnData Documentation — Getting Started. https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html
4. scverse Tutorials — Getting Started with AnnData. https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html
5. Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15. https://doi.org/10.1186/s13059-017-1382-0

---

⬅️ [Back: Basic Tutorial](../02_basic_scrna_tutorial/README.md) &nbsp;|&nbsp; [Back to Main README](../README.md)
