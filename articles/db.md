# Database details

The clustermole meta-database collects gene signatures from public
marker databases and cell atlases. The same cell type can have several
signatures that can differ by organ, species, and study.

The database includes data from the following sources:

- CellMarker - Hu et al. *Nucleic Acids Research* (2023)
  [10.1093/nar/gkac947](https://doi.org/10.1093/nar/gkac947)
- CellMatch - Shao et al. *iScience* (2020)
  [10.1016/j.isci.2020.100882](https://doi.org/10.1016/j.isci.2020.100882)
- Cell Taxonomy - Jiang et al. *Nucleic Acids Research* (2023)
  [10.1093/nar/gkac816](https://doi.org/10.1093/nar/gkac816)
- DISCO - Li et al. *Nucleic Acids Research* (2022)
  [10.1093/nar/gkab1020](https://doi.org/10.1093/nar/gkab1020)
- HPA - Karlsson et al. *Science Advances* (2021)
  [10.1126/sciadv.abh2169](https://doi.org/10.1126/sciadv.abh2169)
- MSigDB - Liberzon et al. *Bioinformatics* (2011)
  [10.1093/bioinformatics/btr260](https://doi.org/10.1093/bioinformatics/btr260)
- PanglaoDB - Franzen et al. *Database* (2019)
  [10.1093/database/baz046](https://doi.org/10.1093/database/baz046)
- SaVanT - Lopez et al. *BMC Genomics* (2017)
  [10.1186/s12864-017-4167-7](https://doi.org/10.1186/s12864-017-4167-7)
- ScType - Ianevski et al. *Nature Communications* (2022)
  [10.1038/s41467-022-28803-w](https://doi.org/10.1038/s41467-022-28803-w)
- TISSUES - Palasca et al. *Database* (2018)
  [10.1093/database/bay003](https://doi.org/10.1093/database/bay003)
- xCell - Aran et al. *Genome Biology* (2017)
  [10.1186/s13059-017-1349-1](https://doi.org/10.1186/s13059-017-1349-1)

Many of these are static, but some are updated periodically.

The raw source data is processed to make the structure and content
consistent across databases. Database preparation steps include:

- Converting recognized human and mouse gene aliases to canonical
  symbols.
- Checking that the species labels are supported by the provided genes.
- Removing signatures with fewer than five or more than 1,000 genes.
- Mapping genes to orthologs in the other species.

## Database data frame

The summaries below describe the database based on the actual package
contents.

``` r

library(clustermole)
library(dplyr)
```

Retrieve a data frame of all cell type markers in the database. The
`species` argument selects the appropriate gene symbols and does not
restrict signatures to the selected species.

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
gene-to-signature mapping per row.

## Total cell types

Check the total number of available cell type signatures. Each
`celltype_full` value uniquely identifies a signature (cell type,
organ/tissue, species, and source).

``` r

length(unique(markers$celltype_full))
#> [1] 13130
```

## Cell types per source

Check the number of cell type signatures from each source database.

``` r

distinct(markers, db, celltype_full) |> count(db)
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

## Cell types per species

Check the number of cell type signatures per original source species.
Not all sources provide information about the species, so it can be
unknown or blank.

``` r

distinct(markers, species, celltype_full) |> count(species)
#> # A tibble: 3 × 2
#>   species     n
#>   <chr>   <int>
#> 1 ""       1216
#> 2 "HS"     8946
#> 3 "MM"     2968
```

## Cell types per organ

Check the number of available cell type signatures per organ or tissue.
This label is not standardized and is not always available.

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

Check the package version since the database contents can change.

``` r

packageVersion("clustermole")
#> [1] '1.1.1.9000'
```
