# clustermole: blindly digging for cell types in scRNA-seq clusters

[![CRAN](https://www.r-pkg.org/badges/version/clustermole)](https://cran.r-project.org/package=clustermole)
[![Travis Build Status](https://travis-ci.org/igordot/clustermole.svg?branch=master)](https://travis-ci.org/igordot/clustermole)
[![codecov](https://codecov.io/gh/igordot/clustermole/branch/master/graph/badge.svg)](https://codecov.io/gh/igordot/clustermole)

![clustermole-book](https://user-images.githubusercontent.com/6363505/72761156-12414280-3ba9-11ea-87de-57ff6cd690bb.png)

A typical computational pipeline to process single-cell RNA sequencing (scRNA-seq) data includes clustering of cells as one of the steps. Assignment of cell type labels to those clusters is often a time-consuming process that involves manual inspection of the cluster marker genes complemented with a detailed literature search. This is especially challenging for those who are not familiar with all the captured subpopulations or have unexpected contaminants. The clustermole R package provides a comprehensive meta collection of cell identity markers for thousands of human and mouse cell types sourced from a variety of databases as well as methods to query them.

The clustermole package provides three primary features:

* cell type prediction based on marker genes
* cell type prediction based on a full expression matrix
* a database of cell type markers

---

Install clustermole from CRAN:

```r
install.packages("clustermole")
```

Alternatively, you can install the development version from GitHub (not recommended):

```r
BiocManager::install("igordot/clustermole", update = FALSE)
```

Load clustermole:

```r
library(clustermole)
```

Perform cell type overrepresentation analysis for a given set of genes:

```r
clustermole_overlaps(genes, species = "hs")
```

Perform cell type enrichment for a given full gene expression matrix:

```r
clustermole_enrichment(expr_mat, species = "hs")
```

Retrieve a table of all cell type markers:

```r
clustermole_markers(species = "hs")
```

---

*Image credit: "A Child's Primer Of Natural History" by Oliver Herford*
