# MTAFS: An Adaptive and Robust Method for Multi-trait Analysis of Genome-wide Association Studies Using Summary Statistics

MTAFS is an efficient and robust method for multi-trait analysis of GWAS. It uses GWAS summary statistics and construct its test statistic adaptively. Compared with many existing methods, MTAFS has three compelling features making it a useful method in practice:

- MTAFS can control type 1 error well even when there are hundreds of traits and significance levels are small, while some exisiting methods are inflated at small significance levels.
- MTAFS has robust power performance given various settings. In constrast, many methods are not robust.
- MTAFS is an efficient method by avoiding permutations. In the real data application, it can analyze 593,416 SNPs and 212 traits in 10 minutes by using 60 cores of 4GB memory.


## Installation

`devtools::install_github("Qiaolan/MTAFS")` or `remotes::install_github("Qiaolan/MTAFS")`
`library(MTAFS)`

