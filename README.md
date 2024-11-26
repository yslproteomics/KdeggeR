# KdeggeR

## Overview

The KdeggeR package is designed to perform peptide and protein turnover rates estimation from dynamic SILAC labeling (pulse-SILAC, or pSILAC) proteomic experiments analyzed using multiplex DIA-MS. The package is optimized to handle DIA-MS data processed using commong raw MS data processing tools such as Spectronaut, DIA-NN, or Fragpipe, but can handle data measured using DDA-MS and analyzed with software tools like MaxQuant. The package offers optimized input data filtering, several kloss estimation methods for precursor/peptide level data, and functions for protein kloss and kdeg aggregation. The curve fits can be also inspected and exported using provided functions.

## Documentation

See the vignette for detailed instructions and example code.

The .pdf version of the vignette is provided [here](https://github.com/yslproteomics/KdeggeR/blob/main/inst/manual/KdeggerUserManual.pdf). 

See the preprint for a more detailed description of the package's functions and applications: doi: https://doi.org/10.1101/2024.10.28.620709

## Hardware Requirements

This package requires only a standard computer with enough RAM to support the operations defined by a user.

The package was developed and tested in a computer with the following specs: RAM 32 GB, CPU: 6 cores, 3.2 GHz/core.

## Software requirements

This package was tested on macOS and Windows operating systems. The development version of the package has been tested on the following systems:

- macOS Sonoma, version 14.3.1
- Microsoft Windows 11, version 10.0.22631

Before setting up this package, users should have R version 4.3.0 or higher, and install the dependencies as specified below. 


## Installation Guide

To install all dependencies: 

```{r}
# Required packages
if(!require(pacman)) install.packages("pacman")
pacman::p_load(dplyr, purrr, stringr, tibble, outliers)

# Optional R package for robust linear model fitting
install.packages("MASS")
```

To install the package: 

```{r}
library("devtools")
install_github("yslproteomics/KdeggeR", build_vignettes = TRUE)
```

To open the vignette with detailed instructions and example code: 

```{r}
vignette("KdeggerUserManual", package = "KdeggeR")
```

The installation takes about 30 seconds in a computer with the following specs: RAM 32 GB, CPU: 6 cores, 3.2 GHz/core. 

## How to Cite

If you use the KdeggeR package in your work, please cite the following preprint. 

A Comprehensive and Robust Multiplex-DIA Workflow Profiles Protein Turnover Regulations Associated with Cisplatin Resistance
Barbora Salovska, Wenxue Li, Oliver M. Bernhardt, Pierre-Luc Germain, Tejas Gandhi, Lukas Reiter, Yansheng Liu
bioRxiv 2024.10.28.620709; doi: https://doi.org/10.1101/2024.10.28.620709


## Quick demo analysis guide

This guide provides a general workflow how to run the demo analysis using provided datasets. For more details and options please see the vignette and documentation. 

The example dataset contains a data.frame containing the first 20,000 unique precursors from the A2780Cis and A2780 parental cell line dataset analyzed using the labeled workflow in Spectronaut v19 using the Group Q-value filtering. The full dataset will be available through the PRIDE repository with the dataset identifier PXD057632. 

### Example data

To list the example data with description, run the following code.  

```{r}
data(package = "KdeggeR")
```

See the example input from the analysis in Spectronaut v19.

```{r}
KdeggeR::example_spectronaut %>%
  dplyr::glimpse()
```

See the example design tables.

```{r}
# Design table without replicate design
KdeggeR::example_spectronaut_design %>%
  dplyr::glimpse()

# Design table with replicate design
KdeggeR::example_spectronaut_design_replicates %>%
  dplyr::glimpse()
```

See the example output pSILAC object. 

```{r}
KdeggeR::example_spectronaut_pSILAC_object %>% 
  View()
```

### 1. Generate pSILAC class object 

```{r}
# Analysis without replicate design
input_data <- KdeggeR::example_spectronaut
input_design <- KdeggeR::example_spectronaut_design

pSILAC_object <- KdeggeR::generatepSILACObject(dataset = input_data, 
                                   design = input_design, 
                                   inputDataType = "spectronaut",
                                   aggregate.replicates = NA, # if NA, replicates are not aggregated
                                   filterPeptides = T, 
                                   ncores = NULL,
                                   noiseCutoff = 8)
                                   

# Analysis with replicate design
input_data <- KdeggeR::example_spectronaut
input_design <- KdeggeR::example_spectronaut_design_replicates

pSILAC_object_replicates <- KdeggeR::generatepSILACObject(dataset = input_data, 
                                   design = input_design, 
                                   inputDataType = "spectronaut",
                                   aggregate.replicates = "mean", # can be "mean" or "median"
                                   filterPeptides = T, 
                                   ncores = NULL,
                                   noiseCutoff = 8)
```

This step takes about 6 seconds and 13 seconds (with replicate aggregation) in a computer with the following specs: RAM 32 GB, CPU: 6 cores, 3.2 GHz/core. 

### 2. Data quality filtering

```{r}
# Filter based on valid values
pSILAC_object <- KdeggeR::filterValidValues(pSILAC_object, values_cutoff = 2, skip_time_point = 1)

# Filter based on monotone trend
pSILAC_object <- KdeggeR::filterMonotone(pSILAC_object, skip_time_point = 1)
pSILAC_object <- KdeggeR::filterMonotoneTimePoint1(pSILAC_object)

# Filter based on linear regression
pSILAC_object <- KdeggeR::filterLinearRegression(pSILAC_object, skip_time_point = 1, R2_cutoff = 0.9, p_cutoff = 0.05)
```

This step takes about 1.5 min in a computer with the following specs: RAM 32 GB, CPU: 6 cores, 3.2 GHz/core. 

### 3. Calculate k_loss

This is a convenient wrapper function to perform precursor-level k_loss estimation using all three models and protein-level k_loss aggregation using weighted average of precursor-level k_loss values using default parameters. 

```{r}
pSILAC_object <- KdeggeR::calcAllRates(pSILAC_object, method = "RIA", 
                  ag.metric = "mean", 
                  ag.weights = "both")
```

This step takes about 1.6 min in a computer with the following specs: RAM 32 GB, CPU: 6 cores, 3.2 GHz/core. 

### 4. Calculate protein degradation rates (k_deg) and halflives

See the example k_cd table (experimentally-determined cell division rates). 

```{r}
KdeggeR::example_kcd %>%
  dplyr::glimpse()
```

Calculate k_deg: 

```{r}
# Using experimentally-derived k_cd
input_kcd <- KdeggeR::example_kcd
pSILAC_object <- KdeggeR::calcKdeg(pSILAC_object, rate_df = input_kcd, type = "kcd")

# Using estimated k_cd (k_perc)
pSILAC_object <- KdeggeR::calcKdeg(pSILAC_object, rate_df = NULL, type = "kperc", perc_neg = 0.01)
```

Calculate t(1/2): 

```{r}
pSILAC_object <- KdeggeR::calcHalflife(pSILAC_object)
```

These steps take less than 1 s in a computer with the following specs: RAM 32 GB, CPU: 6 cores, 3.2 GHz/core. 

### 5. Export protein degradation rates

```{r}
# Protein degradation rates are stored in the protein.kdeg dataframe
protein_table <- pSILAC_object$protein.kdeg %>%
  tibble::rownames_to_column("Protein_ID") %>%
  dplyr::glimpse()
```

### 6. Visualize results

Plot precursor RIA model: 

```{r}
KdeggeR::plotPeptideRIA(pSILAC_object, peptide = "_MLIPYIEHWPR_.3")
```

Plot precursor ln(H/L +1) model: 

```{r}
KdeggeR::plotPeptideHoL(pSILAC_object, peptide = "_MLIPYIEHWPR_.3")
```

Plot protein RIA:

```{r}
KdeggeR::plotProteinRIA(pSILAC_object, protein = "O75208")
```

Plot protein HoL:

```{r}
KdeggeR::plotProteinHol(pSILAC_object, protein = "O75208")
```

Plot protein summary:

```{r}
KdeggeR::plotProtein(pSILAC_object, protein = "O75208")
```
