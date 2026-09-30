# Read a GMT file into a data frame

Read a GMT file into a data frame

## Usage

``` r
read_gmt(file, geneset_label = "celltype", gene_label = "gene")
```

## Arguments

- file:

  A file path, URL, or connection.

- geneset_label:

  Output column name for gene sets (GMT column 1).

- gene_label:

  Output column name for genes (GMT columns 3 onward).

## Value

A data frame with gene sets and genes, one gene per row.

## Examples

``` r
if (FALSE) { # \dontrun{
gmt <- "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"
gmt_tbl <- read_gmt(gmt)
head(gmt_tbl)
} # }
```
