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
if (FALSE) { # \dontrun{
gmt <- "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"
gmt_tbl <- read_gmt(gmt)
head(gmt_tbl)
} # }
```
