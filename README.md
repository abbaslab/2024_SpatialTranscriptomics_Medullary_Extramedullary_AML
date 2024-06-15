# 2024_SpatialTranscriptomics_Medullary_Extramedullary_AML
 
# Spatial Transcriptomic Profiling of Medullary and Extramedullary Acute Myeloid Leukemia (AML)

This repository contains the code and analysis for the paper titled "Spatial Transcriptomic Profiling Reveals Inflammation and Trans-differentiation States of Acute Myeloid Leukemia in Extramedullary and Medullary Tissues."

## Overview

This study provides a comprehensive spatial transcriptomic analysis of both medullary and extramedullary AML tissues. Using Visium and GeoMx spatial transcriptomics technologies, we explored the spatial distribution of cell types, identified key inflammatory pathways, and mapped the microenvironment of AML. The repository includes the R scripts used for data processing, visualization, and analysis as described in the manuscript.

## Installation

To run the analysis, you will need to have R and the following packages installed:

- Seurat
- ggplot2
- dplyr
- SCP
- SpaCET
- ggthemes
- RColorBrewer
- harmony
- ggpubr
- AUCell
- CellChat
- 

You can install these packages using the following commands:

```R
install.packages(c("ggplot2", "dplyr", "ggthemes", "RColorBrewer", "ggpubr", "harmony", ))
install.packages("Seurat")
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("AUCell")


# Install GitHub packages
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("jinworks/CellChat")
devtools::install_github("zhanghao-njmu/SCP")
devtools::install_github("data2intelligence/SpaCET")
