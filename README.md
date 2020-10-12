# clustermole: exploratory scRNA-seq cell type analysis

[![CRAN](https://www.r-pkg.org/badges/version/clustermole)](https://cran.r-project.org/package=clustermole)
[![R build status](https://github.com/igordot/clustermole/workflows/R-CMD-check/badge.svg)](https://github.com/igordot/clustermole/actions)
[![Travis Build Status](https://travis-ci.com/igordot/clustermole.svg?branch=master)](https://travis-ci.com/igordot/clustermole)
[![codecov](https://codecov.io/gh/igordot/clustermole/branch/master/graph/badge.svg)](https://codecov.io/gh/igordot/clustermole)

![clustermole-book](https://user-images.githubusercontent.com/6363505/72761156-12414280-3ba9-11ea-87de-57ff6cd690bb.png)

## Overview

Assignment of cell type labels to single-cell RNA sequencing (scRNA-seq) clusters is often a time-consuming process that involves manual inspection of the cluster marker genes complemented with a detailed literature search.
This can be especially challenging when unexpected or poorly described populations are present.
The clustermole R package provides methods to query thousands of human and mouse cell identity markers sourced from a variety of databases.

The clustermole package provides three primary features:

* a database of markers for hundreds of cell types
* cell type prediction based on marker genes
* cell type prediction based on the full expression matrix

Check the [documentation website](https://igordot.github.io/clustermole) for more information.

---

*Image credit: "A Child's Primer Of Natural History" by Oliver Herford*
