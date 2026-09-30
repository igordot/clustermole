# Introduction to clustermole

## Overview

The clustermole R package is designed to simplify the assignment of cell
type labels to unknown cell populations, such as scRNA-seq clusters. It
provides methods to query cell identity markers sourced from a variety
of databases. The package includes three primary features:

- a meta-database of human and mouse markers for thousands of cell types
  ([`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md))
- cell type prediction based on a set of marker genes
  ([`clustermole_overlaps()`](https://igordot.github.io/clustermole/reference/clustermole_overlaps.md))
- cell type prediction based on a table of expression values
  ([`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md))

## Setup

You can install clustermole from
[CRAN](https://cran.r-project.org/package=clustermole).

``` r

install.packages("clustermole")
```

Load clustermole.

``` r

library(clustermole)
```

## Cell type markers

You can use clustermole as a simple database and get a data frame of all
cell type markers.

``` r

markers <- clustermole_markers(species = "hs")
markers
#> # A tibble: 521,262 × 8
#>    celltype_full         db    species organ celltype gene_origi…¹ gene  n_genes
#>    <chr>                 <chr> <chr>   <chr> <chr>    <chr>        <chr>   <int>
#>  1 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … ADAMTS2      ADAM…      10
#>  2 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … FN1          FN1        10
#>  3 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … KLHL1        KLHL1      10
#>  4 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … LIPM         LIPM       10
#>  5 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … NPSR1        NPSR1      10
#>  6 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … NTS          NTS        10
#>  7 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … RAB38        RAB38      10
#>  8 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … RXFP1        RXFP1      10
#>  9 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … STAC         STAC       10
#> 10 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … TLE4         TLE4       10
#> # ℹ 521,252 more rows
#> # ℹ abbreviated name: ¹​gene_original
```

Each row contains a gene and a cell type associated with it. The `gene`
column is the gene symbol (human or mouse), the `gene_original` column
is the gene symbol from the source database, and the `celltype_full`
column contains the full cell type string, including the species and the
original database.

Many tools that work with gene sets require input as a list. To convert
the markers from a data frame to a list, you can use `gene` as the
values and `celltype_full` as the grouping variable.

``` r

markers_list <- split(x = markers$gene, f = markers$celltype_full)
```

## Cell types based on marker genes

If you have a character vector of genes, such as cluster markers, you
can compare them to known cell type markers to see if they overlap any
of the known cell type markers (overrepresentation analysis).

``` r

my_overlaps <- clustermole_overlaps(genes = my_genes_vec, species = "hs")
```

## Cell types based on an expression matrix

If you have expression values, such as average expression for each
cluster, you can perform cell type enrichment based on the full gene
expression matrix (log-transformed CPM/TPM/FPKM values). The matrix
should have genes as rows and clusters/samples as columns. The
underlying enrichment method can be changed using the `method`
parameter.

``` r

my_enrichment <- clustermole_enrichment(expr_mat = my_expr_mat, species = "hs")
```
