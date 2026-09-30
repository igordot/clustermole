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
#> # A tibble: 6 × 8
#>   celltype_full         db    species organ celltype gene_original gene  n_genes
#>   <chr>                 <chr> <chr>   <chr> <chr>    <chr>         <chr>   <int>
#> 1 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … ADAMTS2       ADAM…      10
#> 2 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … FN1           FN1        10
#> 3 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … KLHL1         KLHL1      10
#> 4 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … LIPM          LIPM       10
#> 5 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … NPSR1         NPSR1      10
#> 6 (Pro-) Subiculum | H… ScTy… ""      Hipp… (Pro-) … NTS           NTS        10
```
