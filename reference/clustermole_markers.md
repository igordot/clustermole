# Available cell type markers

Retrieve the full list of cell type markers in the `clustermole`
database.

## Usage

``` r
clustermole_markers(species = c("hs", "mm"))
```

## Arguments

- species:

  Species: `hs` for human or `mm` for mouse.

## Value

A data frame of cell type markers (one gene per row).

## Examples

``` r
markers <- clustermole_markers()
head(markers)
#> # A tibble: 6 × 9
#>   celltype_full            db    species_original species organ celltype n_genes
#>   <chr>                    <chr> <chr>            <chr>   <chr> <chr>      <int>
#> 1 (Pro-) Subiculum | Hipp… ScTy… ""               ""      Hipp… (Pro-) …      10
#> 2 (Pro-) Subiculum | Hipp… ScTy… ""               ""      Hipp… (Pro-) …      10
#> 3 (Pro-) Subiculum | Hipp… ScTy… ""               ""      Hipp… (Pro-) …      10
#> 4 (Pro-) Subiculum | Hipp… ScTy… ""               ""      Hipp… (Pro-) …      10
#> 5 (Pro-) Subiculum | Hipp… ScTy… ""               ""      Hipp… (Pro-) …      10
#> 6 (Pro-) Subiculum | Hipp… ScTy… ""               ""      Hipp… (Pro-) …      10
#> # ℹ 2 more variables: gene_original <chr>, gene <chr>
```
