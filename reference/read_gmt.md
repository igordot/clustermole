# Read a GMT file into a data frame

Read a GMT file into a data frame

## Usage

``` r
read_gmt(file, geneset_label = "celltype", gene_label = "gene")
```

## Arguments

- file:

  A connection object or a character string (can be a URL).

- geneset_label:

  Column name for gene sets (first column of the GMT file) in the output
  data frame.

- gene_label:

  Column name for genes (variable columns of the GMT file) in the output
  data frame.

## Value

A data frame with gene sets as the first column and genes as the second
column (one gene per row).

## Examples

``` r
gmt <- "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"
gmt_tbl <- read_gmt(gmt)
head(gmt_tbl)
#> # A tibble: 6 × 2
#>   celltype                                              gene    
#>   <chr>                                                 <chr>   
#> 1 Zheng_Cord_Blood_C1_Putative_Megakaryocyte_Progenitor ABCC3   
#> 2 Zheng_Cord_Blood_C1_Putative_Megakaryocyte_Progenitor ABCC4   
#> 3 Zheng_Cord_Blood_C1_Putative_Megakaryocyte_Progenitor ACTN1   
#> 4 Zheng_Cord_Blood_C1_Putative_Megakaryocyte_Progenitor ARHGAP18
#> 5 Zheng_Cord_Blood_C1_Putative_Megakaryocyte_Progenitor ARHGAP6 
#> 6 Zheng_Cord_Blood_C1_Putative_Megakaryocyte_Progenitor BANK1   
```
