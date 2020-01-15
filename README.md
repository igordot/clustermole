# clustermole: blindly digging for cell types in scRNA-seq clusters

[![Travis Build Status](https://travis-ci.org/igordot/clustermole.svg?branch=master)](https://travis-ci.org/igordot/clustermole)
[![codecov](https://codecov.io/gh/igordot/clustermole/branch/master/graph/badge.svg)](https://codecov.io/gh/igordot/clustermole)

> See, children, the misguided Mole.  
> He lives down in a deep, dark hole;  
> Sweetness, and light, and good fresh air  
> Are things for which he does not care.  
> He has not even that makeshift  
> Of feeble minds - the social gift.  
> But say not that he has no soul,  
> Lest haply we misjudge the Mole;  
> Nay, if we measure him by men,  
> No doubt he sits in his dark den  
> Instructing others blind as he  
> Exactly how the world should be.  
>
> -- Oliver Herford

A typical computational pipeline to process single-cell RNA sequencing (scRNA-seq) data  involves clustering of cells. Assignment of cell type labels to those clusters is often a time-consuming process that involves manual inspection of the cluster marker genes complemented with a detailed literature search. This is especially challenging if you are not familiar with all the captured subpopulations or have unexpected contaminants. `clustermole` is an R package that provides a comprehensive meta collection of cell identity markers for thousands of human and mouse cell types sourced from a variety of databases as well as methods to query them.

Install clustermole (development version):

```r
BiocManager::install("igordot/clustermole", update = FALSE)
```

Load clustermole:

```r
library(clustermole)
```

Retrieve a table of all cell type markers:

```r
clustermole_markers(genes, species)
```

Perform cell type overrepresentation analysis for a given set of genes:

```r
clustermole_overlaps(expr_mat, species)
```

Perform cell type enrichment for a given full gene expression matrix:

```r
clustermole_enrichment(expr_mat, species)
```
