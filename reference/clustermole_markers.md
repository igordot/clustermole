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
#>   celltype_full         db    species organ celltype n_genes gene_original gene 
#>   <chr>                 <chr> <chr>   <chr> <chr>      <int> <chr>         <chr>
#> 1 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 ACCSL         ACCSL
#> 2 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 ACVR1B        ACVR…
#> 3 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 ASF1B         ASF1B
#> 4 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 BCL2L10       BCL2…
#> 5 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 BLCAP         BLCAP
#> 6 1-cell stage cell (B… Cell… Human   Embr… 1-cell …      32 CASC3         CASC3
```
