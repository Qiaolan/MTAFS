# MTAFS: An Adaptive and Robust Method for Multi-trait Analysis of Genome-wide Association Studies Using Summary Statistics

MTAFS is an efficient and robust method for multi-trait analysis of GWAS. It uses GWAS summary statistics and construct its test statistic adaptively. Compared with many existing methods, MTAFS has three compelling features making it a useful method in practice:

- MTAFS can control type 1 error well even when there are hundreds of traits and significance levels are small.
- MTAFS has robust power performance given various settings.
- MTAFS is an efficient method by avoiding permutations.


## Installation

`devtools::install_github("Qiaolan/MTAFS")` or `remotes::install_github("Qiaolan/MTAFS")`

`library(MTAFS)`

## Vignette

The vignette of MTAFS is [here](http://htmlpreview.github.io/?https://github.com/Qiaolan/MTAFS/blob/main/index.html)
