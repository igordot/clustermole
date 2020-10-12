# clustermole: blindly digging for cell types in scRNA-seq clusters

[![CRAN](https://www.r-pkg.org/badges/version/clustermole)](https://cran.r-project.org/package=clustermole)
[![R build status](https://github.com/igordot/clustermole/workflows/R-CMD-check/badge.svg)](https://github.com/igordot/clustermole/actions)
[![Travis Build Status](https://travis-ci.com/igordot/clustermole.svg?branch=master)](https://travis-ci.com/igordot/clustermole)
[![codecov](https://codecov.io/gh/igordot/clustermole/branch/master/graph/badge.svg)](https://codecov.io/gh/igordot/clustermole)

![clustermole-book](https://user-images.githubusercontent.com/6363505/72761156-12414280-3ba9-11ea-87de-57ff6cd690bb.png)

## About

Assignment of cell type labels to single-cell RNA sequencing (scRNA-seq) clusters is often a time-consuming process that involves manual inspection of the cluster marker genes complemented with a detailed literature search. This is especially challenging when unexpected or poorly described populations are present. The clustermole R package provides methods to query thousands of human and mouse cell identity markers sourced from a variety of databases.

The clustermole package provides three primary features:

* cell type prediction based on marker genes
* cell type prediction based on a full expression matrix
* a database of cell type markers

## Usage

A [vignette](https://CRAN.R-project.org/package=clustermole/vignettes/clustermole-intro.html) is available with usage examples.

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
