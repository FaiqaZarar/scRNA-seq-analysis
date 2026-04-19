# Single-Cell RNA Sequencing Analysis

[![Platform](https://img.shields.io/badge/Platform-Python%203.10-blue?style=flat-square)](https://www.python.org/)
[![Framework](https://img.shields.io/badge/Framework-Scanpy%20%7C%20AnnData-4CAF50?style=flat-square)](https://scanpy.readthedocs.io/)
[![Data](https://img.shields.io/badge/Data-10X%20Genomics-orange?style=flat-square)](https://www.10xgenomics.com/)
[![Tutorials](https://img.shields.io/badge/Tutorials-3%20Sections-purple?style=flat-square)](https://github.com/)
[![License](https://img.shields.io/badge/License-CC--BY--4.0-lightgrey?style=flat-square)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Faiqa Zarar
**Institution:** NUST — National University of Sciences and Technology
**Course:** Bioinformatics
**Date:** April 2026

---

## Abstract

This repository documents a complete single-cell RNA sequencing (scRNA-seq) analysis workflow structured across three independent sections. Section 1 covers end-to-end pre-processing of 10X Genomics Chromium scRNA-seq data, including quality control, doublet detection, normalization, and highly variable gene selection using Scanpy. Section 2 presents an updated scRNA-seq analysis tutorial encompassing dimensionality reduction via PCA, neighborhood graph construction, UMAP visualization, and Leiden-based cell clustering with marker gene identification using the PBMC 3K dataset. Section 3 provides hands-on tutorials for the AnnData data format following the official anndata.readthedocs.io and scverse-tutorials documentation, covering object construction, metadata handling, subsetting, views, layers, and disk-based I/O. Together, these sections represent a reproducible, well-documented reference for scRNA-seq data analysis using the Python scverse ecosystem.

**Keywords:** single-cell RNA sequencing, scRNA-seq, 10X Genomics, Scanpy, AnnData, UMAP, Leiden clustering, quality control, bioinformatics, scverse

---

## Table of Contents

1. [Repository Structure](#1-repository-structure)
2. [Environment Setup](#2-environment-setup)
3. [Section 1 — Pre-processing of 10X scRNA Datasets](#3-section-1--pre-processing-of-10x-scrna-datasets)
4. [Section 2 — Basic scRNA-seq Tutorial (Updated)](#4-section-2--basic-scrna-seq-tutorial-updated)
5. [Section 3 — AnnData Workflows](#5-section-3--anndata-workflows)
6. [References](#6-references)

---

## 1. Repository Structure

```
scRNA-seq-analysis/
│
├── README.md                              ← This document
│
├── 01_preprocessing/
│   ├── README.md                          ← Pre-processing pipeline documentation
│   ├── scripts/
│   │   ├── quality_control.py             # QC metrics and filtering
│   │   ├── normalization.py               # Library-size normalization
│   │   ├── feature_selection.py           # Highly variable gene selection
│   │   └── doublet_detection.py           # Scrublet doublet removal
│   └── outputs/
│       ├── qc_metrics.csv                 # Per-cell QC statistics
│       ├── filtered_matrix/               # Post-QC sparse matrix
│       └── preprocessed_adata.h5ad        # Final preprocessed AnnData
│
├── 02_basic_scrna_tutorial/
│   ├── README.md                          ← Tutorial documentation
│   ├── notebooks/
│   │   ├── 01_load_and_explore.ipynb      # Data loading and exploration
│   │   ├── 02_clustering.ipynb            # PCA, UMAP, Leiden clustering
│   │   └── 03_marker_genes.ipynb          # Marker gene identification
│   └── outputs/
│       ├── clustered_adata.h5ad           # AnnData with cluster labels
│       ├── umap_clusters.png              # UMAP plot colored by cluster
│       └── marker_genes.csv              # Top marker genes per cluster
│
└── 03_anndata_workflows/
    ├── README.md                          ← AnnData tutorial documentation
    ├── notebooks/
    │   ├── anndata_getting_started.ipynb  # Based on anndata.readthedocs.io
    │   └── anndata_scverse_tutorial.ipynb # Based on scverse-tutorials
    └── outputs/
        ├── sample_anndata.h5ad            # Example AnnData object
        └── tutorial_outputs/
            ├── subset_adata.h5ad
            └── concatenated_adata.h5ad
```

---

## 2. Environment Setup

```bash
# Step 1: Create and activate conda environment
conda create -n scrna python=3.10
conda activate scrna

# Step 2: Install all required packages
pip install scanpy anndata pandas numpy matplotlib seaborn
pip install scrublet leidenalg jupyterlab

# Step 3: Verify installation
python -c "import scanpy; print('Scanpy:', scanpy.__version__)"
python -c "import anndata; print('AnnData:', anndata.__version__)"

# Step 4: Launch Jupyter
jupyter lab
```

---

## 3. Section 1 — Pre-processing of 10X scRNA Datasets

Complete pre-processing pipeline from raw Cell Ranger output to analysis-ready AnnData objects. Covers quality control filtering (mitochondrial %, n_genes, total_counts), Scrublet-based doublet detection, library-size normalization, log1p transformation, and highly variable gene selection.

➡️ [View Full Section 1 README](./01_preprocessing/README.md)

---

## 4. Section 2 — Basic scRNA-seq Tutorial (Updated)

Updated end-to-end scRNA-seq analysis tutorial using the PBMC 3K dataset. Covers PCA, neighborhood graph construction, UMAP embedding, Leiden clustering (replacing deprecated Louvain), Wilcoxon-based marker gene ranking, and visualization via dot plots and violin plots.

➡️ [View Full Section 2 README](./02_basic_scrna_tutorial/README.md)

---

## 5. Section 3 — AnnData Workflows

Hands-on tutorials covering the AnnData data format based on official documentation:
- https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html
- https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html

Covers AnnData object construction, `.obs`/`.var`/`.obsm`/`.uns`/`.layers` attributes, subsetting, views vs copies, sparse matrix handling, concatenation, and HDF5 (`.h5ad`) I/O.

➡️ [View Full Section 3 README](./03_anndata_workflows/README.md)

---

## 6. References

1. Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15. https://doi.org/10.1186/s13059-017-1382-0
2. Virshup, I., et al. (2023). The scverse project provides a computational ecosystem for single-cell omics data analysis. *Nature Biotechnology*, 41, 604–606. https://doi.org/10.1038/s41587-023-01733-8
3. Luecken, M. D., & Theis, F. J. (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*, 15(6), e8746. https://doi.org/10.15252/msb.20188746
4. Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233. https://doi.org/10.1038/s41598-019-41695-z
5. Wolock, S. L., Lopez, R., & Klein, A. M. (2019). Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Systems*, 8(4), 281–291. https://doi.org/10.1016/j.cels.2018.11.005

---

*This repository was prepared as part of a Bioinformatics course assignment at NUST.*
*All analyses were performed using publicly available datasets and open-source tools.*
