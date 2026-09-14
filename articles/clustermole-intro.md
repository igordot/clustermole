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
#> # A tibble: 422,292 × 8
#>    celltype_full         db    species organ celltype n_genes gene_origi…¹ gene 
#>    <chr>                 <chr> <chr>   <chr> <chr>      <int> <chr>        <chr>
#>  1 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 ACCSL        ACCSL
#>  2 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 ACVR1B       ACVR…
#>  3 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 ASF1B        ASF1B
#>  4 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 BCL2L10      BCL2…
#>  5 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 BLCAP        BLCAP
#>  6 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 CASC3        CASC3
#>  7 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 CLEC10A      CLEC…
#>  8 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 CNOT11       CNOT…
#>  9 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 DCLK2        DCLK2
#> 10 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 DHCR7        DHCR7
#> # ℹ 422,282 more rows
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
