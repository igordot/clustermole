# Database details

We will load clustermole along with dplyr to help with summarizing the
data.

``` r

library(clustermole)
library(dplyr)
```

You can use clustermole as a simple database and get a table of all cell
type markers.

``` r

markers <- clustermole_markers(species = "hs")
markers
#> # A tibble: 522,188 × 9
#>    celltype_full           db    species_original species organ celltype n_genes
#>    <chr>                   <chr> <chr>            <chr>   <chr> <chr>      <int>
#>  1 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#>  2 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#>  3 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#>  4 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#>  5 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#>  6 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#>  7 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#>  8 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#>  9 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#> 10 (Pro-) Subiculum | Hip… ScTy… ""               ""      Hipp… (Pro-) …      10
#> # ℹ 522,178 more rows
#> # ℹ 2 more variables: gene_original <chr>, gene <chr>
```

Each row contains a gene and a cell type associated with it. The `gene`
column is the gene symbol (human or mouse), the `gene_original` column
is the gene symbol from the source database, and the `celltype_full`
column contains the detailed cell type string including the species and
the original database.

## Cell types

Check the total number of available cell types.

``` r

length(unique(markers$celltype_full))
#> [1] 13350
```

## Cell types by source database

Check the source databases and the number of cell types from each.

``` r

distinct(markers, celltype_full, db) |> count(db)
#> # A tibble: 11 × 2
#>    db               n
#>    <chr>        <int>
#>  1 CellMarker    8097
#>  2 CellMatch      836
#>  3 CellTaxonomy   667
#>  4 DISCO          388
#>  5 HPA             70
#>  6 MSigDB        1084
#>  7 PanglaoDB      336
#>  8 SaVanT         619
#>  9 ScType         247
#> 10 TISSUES        517
#> 11 xCell          489
```

## Cell types by species

Check the number of cell types per species (not available for all cell
types).

``` r

distinct(markers, celltype_full, species) |> count(species)
#> # A tibble: 3 × 2
#>   species     n
#>   <chr>   <int>
#> 1 ""       1262
#> 2 "HS"     8951
#> 3 "MM"     3137
```

## Cell types by organ

Check the number of available cell types per organ (not available for
all cell types).

``` r

distinct(markers, celltype_full, organ) |> count(organ, sort = TRUE)
#> # A tibble: 386 × 2
#>    organ                  n
#>    <chr>              <int>
#>  1 ""                  3698
#>  2 "Brain"              809
#>  3 "Lung"               681
#>  4 "Liver"              546
#>  5 "Peripheral blood"   524
#>  6 "Skin"               510
#>  7 "Bone marrow"        440
#>  8 "Kidney"             434
#>  9 "Pancreas"           318
#> 10 "Breast"             295
#> # ℹ 376 more rows
```

## Package version

Check the package version since the database contents may change.

``` r

packageVersion("clustermole")
#> [1] '1.1.1.9000'
```
