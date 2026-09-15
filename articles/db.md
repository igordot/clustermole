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
#> # A tibble: 417,251 × 8
#>    celltype_full        db    species organ celltype n_genes gene_original gene 
#>    <chr>                <chr> <chr>   <chr> <chr>      <int> <chr>         <chr>
#>  1 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 ACCSL         ACCSL
#>  2 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 ACVR1B        ACVR…
#>  3 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 ASF1B         ASF1B
#>  4 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 BCL2L10       BCL2…
#>  5 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 BLCAP         BLCAP
#>  6 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 CASC3         CASC3
#>  7 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 CLEC10A       CLEC…
#>  8 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 CNOT11        CNOT…
#>  9 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 DCLK2         DCLK2
#> 10 1-cell stage cell (… Cell… Human   Embr… 1-cell …      32 DHCR7         DHCR7
#> # ℹ 417,241 more rows
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
#> [1] 3039
```

## Cell types by source database

Check the source databases and the number of cell types from each.

``` r

distinct(markers, celltype_full, db) |> count(db)
#> # A tibble: 7 × 2
#>   db             n
#>   <chr>      <int>
#> 1 ARCHS4       108
#> 2 CellMarker   692
#> 3 MSigDB       295
#> 4 PanglaoDB    322
#> 5 SaVanT       619
#> 6 TISSUES      537
#> 7 xCell        466
```

## Cell types by species

Check the number of cell types per species (not available for all cell
types).

``` r

distinct(markers, celltype_full, species) |> count(species)
#> # A tibble: 3 × 2
#>   species     n
#>   <chr>   <int>
#> 1 ""        323
#> 2 "Human"  1866
#> 3 "Mouse"   850
```

## Cell types by organ

Check the number of available cell types per organ (not available for
all cell types).

``` r

distinct(markers, celltype_full, organ) |> count(organ, sort = TRUE)
#> # A tibble: 93 × 2
#>    organ                  n
#>    <chr>              <int>
#>  1 ""                  2160
#>  2 "Brain"              122
#>  3 "Immune system"       50
#>  4 "Lung"                47
#>  5 "Kidney"              43
#>  6 "Bone marrow"         42
#>  7 "Liver"               38
#>  8 "Blood"               33
#>  9 "Embryo"              30
#> 10 "Peripheral blood"    29
#> # ℹ 83 more rows
```

## Package version

Check the package version since the database contents may change.

``` r

packageVersion("clustermole")
#> [1] '1.1.1.9000'
```
