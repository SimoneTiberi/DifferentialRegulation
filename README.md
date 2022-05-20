# DifferentialRegulation: a method for differential analyses via hierarchical permutation tests

`DifferentialRegulation` is a method for detecting differentially regulated genes between two groups of samples (e.g., healthy vs. disease, or treated vs. untreated samples), by targeting differences in the balance of spliced and unspliced mRNA abundances, obtained from single-cell RNA-sequencing (scRNA-seq) data.

Check the vignettes for a description of the main conceptual and mathematical aspects, as well as usage guidelines.

A pre-print will follow shortly (~summer 2022).

## Bioconductor installation 
`DifferentialRegulation` is available on [Bioconductor](https://bioconductor.org/packages/DifferentialRegulation) and can be installed with the command:
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("DifferentialRegulation")
```

## Vignette
The vignette illustrating how to use the package can be accessed on 
[Bioconductor](https://bioconductor.org/packages/DifferentialRegulation)
or from R via:
``` r
vignette("DifferentialRegulation")
```
or
``` r
browseVignettes("DifferentialRegulation")
```
