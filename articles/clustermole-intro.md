# Introduction to clustermole

## Overview

The clustermole R package was developed to assist with identifying
candidate cell type labels for single-cell RNA-seq (scRNA-seq) clusters.
The package includes three primary features:

- Retrieve human and mouse cell type signatures with
  [`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md).
- Compare a set of marker genes with the signatures using
  [`clustermole_overlaps()`](https://igordot.github.io/clustermole/reference/clustermole_overlaps.md).
- Score signature enrichment across an expression matrix without
  selecting marker genes using
  [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md).

All three functions return tables that can be easily filtered,
summarized, and used in other workflows.

The reference meta-database combines and cleans cell type markers from
multiple sources. See [database
details](https://igordot.github.io/clustermole/articles/db.md) for an
overview of the contents and the curation process. Sources differ in
their studies, tissues, curation methods, and selection criteria. A
marker may identify a broad cell class, a subtype, or a cell state.
Comparing results across sources can reveal both recurring candidates
and alternative labels worth investigating that a single source may
miss.

For examples of applications across species, tissues, and technologies,
see [publications using
clustermole](https://igordot.github.io/clustermole/articles/pubs.md).

## Setup

clustermole is available from
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
cell type markers. The `species` argument selects human or mouse gene
symbols. It does not restrict signatures to the selected species.

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

The output data frame is in an R-friendly tidy/long format with one
gene-to-signature mapping per row. Each row contains a gene and a cell
type associated with it. The `gene` column is the gene symbol and the
`celltype_full` column contains the full cell type string (cell type,
organ/tissue, species, and source). The `species` column contains the
source species which can differ from the selected gene symbol species.

Many tools that work with gene sets require input as a list. You can use
[`split()`](https://rdrr.io/r/base/split.html) to convert the data frame
to a list using `gene` as the value and `celltype_full` as the grouping
variable.

``` r

markers <- dplyr::distinct(markers, celltype_full, gene)
markers_list <- split(x = markers$gene, f = markers$celltype_full)
```

## Cell types based on marker genes

If you have a vector of genes, such as cluster markers, you can compare
them to known cell type markers to see if they overlap any of the known
cell type markers (over-representation analysis).

``` r

my_overlaps <- clustermole_overlaps(genes = my_genes, species = "hs")
```

## Cell types based on an expression matrix

If you have expression values, such as average expression for each
cluster, you can perform cell type enrichment based on the full gene
expression matrix (log-transformed CPM/TPM/FPKM values). The underlying
enrichment method can be changed using the `method` parameter.

``` r

my_enrichment <- clustermole_enrichment(expr_mat = my_expr_mat, species = "hs")
```

Bioconductor packages GSVA, GSEABase, and singscore need to be installed
to use
[`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md).

``` r

BiocManager::install(c("GSVA", "GSEABase", "singscore"))
```

See the [bone marrow
tutorial](https://igordot.github.io/clustermole/articles/example-bm-seurat.md)
for examples based on a scRNA-seq dataset.
