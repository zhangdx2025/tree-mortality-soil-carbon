# Code for tree mortality and soil carbon study

This code repository is intended solely to reproduce the analyses and figures reported in the accompanying manuscript. 
Because of repository size limits, only a subset of the input data is stored here; 
The full metagenomic gene count matrices and functional annotation tables are archived in Zenodo (https://doi.org/10.5281/zenodo.17606854). 
After downloading those files, place them in the Data/microbiome_profiles directory before running the scripts.

## Contents

- `Code/` – R scripts implementing all data processing, statistical analyses and figure generation described in the manuscript.  
- `Data/` – input data tables required by the R scripts (soil properties and microbial summary data).  
  The full metagenomic gene count matrices and gene annotation catalogues are archived in Zenodo (https://doi.org/10.5281/zenodo.17606854) and are not stored in this repository; after download, they should be placed in `Data/microbiome_profiles`.  
- `LICENSE` – license terms governing reuse of the code.  
- `.gitignore` – patterns specifying files and folders that are not tracked by Git.  
- `README.md` – project overview and instructions for reproducing the results.


## Software requirements

- R (version ≥ 4.5.0)
- R packages: lme4, vegan, smatr, edgeR, xgboost, data.table, ggplot2 and others listed at the top of each script.

## How to use

1. Download or clone this repository.  
2. Open R (or RStudio), set the working directory to the repository root.  
3. Run the scripts in the order indicated in the comments to reproduce the main and supplementary results.
