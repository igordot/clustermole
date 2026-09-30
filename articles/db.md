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
#> # A tibble: 521,262 × 8
#>    celltype_full        db    species organ celltype gene_original gene  n_genes
#>    <chr>                <chr> <chr>   <chr> <chr>    <chr>         <chr>   <int>
#>  1 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … ADAMTS2       ADAM…      10
#>  2 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … FN1           FN1        10
#>  3 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … KLHL1         KLHL1      10
#>  4 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … LIPM          LIPM       10
#>  5 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … NPSR1         NPSR1      10
#>  6 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … NTS           NTS        10
#>  7 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … RAB38         RAB38      10
#>  8 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … RXFP1         RXFP1      10
#>  9 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … STAC          STAC       10
#> 10 (Pro-) Subiculum | … ScTy… ""      Hipp… (Pro-) … TLE4          TLE4       10
#> # ℹ 521,252 more rows
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
#> [1] 13130
```

## Cell types by source database

Check the source databases and the number of cell types from each.

``` r

distinct(markers, celltype_full, db) |> count(db)
#> # A tibble: 11 × 2
#>    db               n
#>    <chr>        <int>
#>  1 CellMarker    7922
#>  2 CellMatch      818
#>  3 CellTaxonomy   647
#>  4 DISCO          388
#>  5 HPA             70
#>  6 MSigDB        1082
#>  7 PanglaoDB      336
#>  8 SaVanT         619
#>  9 ScType         247
#> 10 TISSUES        512
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
#> 1 ""       1216
#> 2 "HS"     8946
#> 3 "MM"     2968
```

## Cell types by organ

Check the number of available cell types per organ (not available for
all cell types).

``` r

distinct(markers, celltype_full, organ) |> count(organ, sort = TRUE)
#> # A tibble: 383 × 2
#>    organ                  n
#>    <chr>              <int>
#>  1 ""                  3686
#>  2 "Brain"              798
#>  3 "Lung"               659
#>  4 "Liver"              534
#>  5 "Peripheral blood"   521
#>  6 "Skin"               497
#>  7 "Kidney"             428
#>  8 "Bone marrow"        427
#>  9 "Pancreas"           315
#> 10 "Breast"             292
#> # ℹ 373 more rows
```

## Package version

Check the package version since the database contents may change.

``` r

packageVersion("clustermole")
#> [1] '1.1.1.9000'
```
