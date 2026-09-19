# Connectoids snRNA-seq analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22847191.svg)](https://doi.org/10.5281/zenodo.22847191)

R scripts for the single-nucleus RNA-seq analyses in:

> Corsi, Sozzi, et al. "Modeling long-range human nigrostriatal connectivity *in vitro* reveals striatal target-dependent maturation of dopamine neurons." 2025.

The study compares two human stem cell-derived models of the nigrostriatal
pathway: fused **assembloids** (ventral midbrain + striatal organoids fused
together) and **connectoids** (the two regions kept spatially separated and
connected via a microfluidic device, preserving the physical separation
between dopaminergic cell bodies and their striatal targets seen *in vivo*).

## Scripts

- **`Assembloid_Connectoid.R`** — Preprocessing, QC, and Harmony integration
  of vMB-STR assembloid vs. connectoid samples (d60/d90). Clusters and
  annotates cell types, then compares dopaminergic (DA) neurons between the
  two models (differential expression, GO enrichment, gene set scoring,
  marker dot plots).
- **`Connectoids.R`** — Same pipeline applied within connectoids only,
  comparing the vMB side of vMB-vMB vs. vMB-STR connectoids to assess how the
  striatal target influences DA neuron maturation.

Both scripts take Cell Ranger output (`filtered_feature_bc_matrix`) plus
`demuxlet` genetic demultiplexing results as input, and produce the UMAP,
dot plot, volcano plot, and heatmap figures used in the manuscript.

## Requirements

R (≥4.6) with: `tidyverse`, `Seurat`, `harmony`, `writexl`, `pheatmap`,
`EnhancedVolcano`, `org.Hs.eg.db`, `clusterProfiler`, `escape`, `msigdbr`,
`msigdb`.

## Data availability

Raw sequencing data are not included in this repository. File paths in the
scripts point to the original analysis environment and should be updated to
local data locations before rerunning.
