# Available cell type markers

Retrieve cell type markers from the clustermole database.

## Usage

``` r
clustermole_markers(species = c("hs", "mm"))
```

## Arguments

- species:

  Gene symbol species: `hs` for human or `mm` for mouse.

## Value

A data frame of cell type markers with these columns:

- `gene`: Canonical gene symbol for the requested species.

- `gene_original`: Original source gene symbol.

- `celltype_full`: Full cell type signature identifier.

- `db`: Source database.

- `celltype`: Cell type label.

- `organ`: Organ label.

- `species`: Source signature species, if known.

- `n_genes`: Gene count per signature.

## Details

The `gene` column uses the official NCBI symbols for the requested
species. The package maps the aliases in each source database to these
symbols. The `gene_original` column preserves the symbol from the
source. The source databases include both human and mouse cell type
markers. The package maps them across species with ortholog data from
the Alliance of Genome Resources and keeps only reciprocal best matches.
The output does not include genes that have no ortholog in the requested
species and does not include signatures with fewer than five genes.

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
