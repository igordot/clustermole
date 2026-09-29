# Cell types based on overlap of marker genes

Perform overrepresentation analysis for a set of genes compared to all
cell type signatures.

## Usage

``` r
clustermole_overlaps(genes, species)
```

## Arguments

- genes:

  A vector of genes.

- species:

  Species: `hs` for human or `mm` for mouse.

## Value

A data frame of enrichment results with hypergeometric test p-values.

## Examples

``` r
my_genes <- c("CD2", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "LTB")
my_overlaps <- clustermole_overlaps(genes = my_genes, species = "hs")
head(my_overlaps)
#> # A tibble: 6 × 10
#>   celltype_full    db    species_original species organ celltype n_genes overlap
#>   <chr>            <chr> <chr>            <chr>   <chr> <chr>      <int>   <dbl>
#> 1 CD4+ T cell (Ga… Cell… "Human"          "HS"    Stom… CD4+ T …      23       7
#> 2 Effector CD4+ T… ScTy… ""               ""      Immu… Effecto…      28       7
#> 3 Naive CD4+ T ce… ScTy… ""               ""      Immu… Naive C…      28       7
#> 4 Memory CD4+ T c… ScTy… ""               ""      Immu… Memory …      29       7
#> 5 Effector CD8+ T… ScTy… ""               ""      Immu… Effecto…      31       7
#> 6 Naive CD8+ T ce… ScTy… ""               ""      Immu… Naive C…      31       7
#> # ℹ 2 more variables: p_value <dbl>, fdr <dbl>
```
