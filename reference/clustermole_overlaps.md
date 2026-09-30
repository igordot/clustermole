# Cell types based on overlap of marker genes

Perform overrepresentation analysis for a set of genes compared to all
cell type signatures.

## Usage

``` r
clustermole_overlaps(genes, species)
```

## Arguments

- genes:

  A character vector of gene symbols.

- species:

  Gene symbol species: `hs` for human or `mm` for mouse.

## Value

A data frame with one row per returned signature:

- `overlap`: Unique gene count shared by the input and signature.

- `p_value`: Hypergeometric test p-value.

- `fdr`: Benjamini-Hochberg adjusted p-value across all tested
  signatures.

- `n_genes`: Unique gene count in the signature for the requested
  species.

- Signature metadata (see
  [`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md)).

## Examples

``` r
my_genes <- c("CD2", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "LTB")
my_overlaps <- clustermole_overlaps(genes = my_genes, species = "hs")
head(my_overlaps)
#> # A tibble: 6 × 9
#>   celltype_full   db    species organ celltype n_genes overlap  p_value      fdr
#>   <chr>           <chr> <chr>   <chr> <chr>      <int>   <dbl>    <dbl>    <dbl>
#> 1 CD4+ T cell (G… Cell… "HS"    Stom… CD4+ T …      23       7 4.13e-22 4.93e-18
#> 2 Effector CD4+ … ScTy… ""      Immu… Effecto…      27       7 1.50e-21 4.93e-18
#> 3 Naive CD4+ T c… ScTy… ""      Immu… Naive C…      27       7 1.50e-21 4.93e-18
#> 4 Memory CD4+ T … ScTy… ""      Immu… Memory …      28       7 1.99e-21 4.93e-18
#> 5 Effector CD8+ … ScTy… ""      Immu… Effecto…      29       7 2.63e-21 4.93e-18
#> 6 Naive CD8+ T c… ScTy… ""      Immu… Naive C…      29       7 2.63e-21 4.93e-18
```
