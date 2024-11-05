# KdeggeR

## Brief Description

The KdeggeR package is designed to perform peptide and protein turnover rates estimation from dynamic SILAC labeling (pulse-SILAC, or pSILAC) proteomic experiments analyzed using multiplex DIA-MS. The package is optimized to handle DIA-MS data processed using commong raw MS data processing tools such as Spectronaut, DIA-NN, or Fragpipe, but can handle data measured using DDA-MS and analyzed with software tools like MaxQuant. The package offers optimized input data filtering, several kloss estimation methods for precursor/peptide level data, and functions for protein kloss and kdeg aggregation. The curve fits can be also inspected and exported using provided functions.

See the preprint for a more detailed description of the packaga functions: doi: https://doi.org/10.1101/2024.10.28.620709
A detailed tutorial coming soon. :) 

## How to Install 

To install all dependencies. 

```{r}
# Required packages
if(!require(pacman)) install.packages("pacman")
pacman::p_load(dplyr, outliers, parallel, purrr, stats, stringr)

# Optional R package for robust linear model fitting
install.packages("MASS")
```

To install the package. 

```{r}
library("devtools")
install_github("yslproteomics/KdeggeR", build_vignettes = TRUE)
```

## How to Cite

If you use the KdeggeR package in your work, please cite the following preprint. 

A Comprehensive and Robust Multiplex-DIA Workflow Profiles Protein Turnover Regulations Associated with Cisplatin Resistance
Barbora Salovska, Wenxue Li, Oliver M. Bernhardt, Pierre-Luc Germain, Tejas Gandhi, Lukas Reiter, Yansheng Liu
bioRxiv 2024.10.28.620709; doi: https://doi.org/10.1101/2024.10.28.620709

