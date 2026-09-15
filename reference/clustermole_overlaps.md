# Cell types based on overlap of marker genes

Perform overrepresentation analysis for a set of genes compared to all
cell type signatures.

## Usage

``` r
clustermole_overlaps(genes, species)
```

## Arguments

- genes:

  A vector of genes.

- species:

  Species: `hs` for human or `mm` for mouse.

## Value

A data frame of enrichment results with hypergeometric test p-values.

## Examples

``` r
my_genes <- c("CD2", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "LTB")
my_overlaps <- clustermole_overlaps(genes = my_genes, species = "hs")
head(my_overlaps)
#> # A tibble: 6 × 9
#>   celltype_full   db    species organ celltype n_genes overlap  p_value      fdr
#>   <chr>           <chr> <chr>   <chr> <chr>      <int>   <dbl>    <dbl>    <dbl>
#> 1 DURANTE_ADULT_… MSig… Human   ""    DURANTE…      22       6 4.36e-18 1.33e-14
#> 2 T memory cells… Pang… Human   "Imm… T memor…      38       6 1.61e-16 2.45e-13
#> 3 T memory cells… Pang… Mouse   "Imm… T memor…      40       6 3.57e-16 3.61e-13
#> 4 T cells | Immu… Pang… Human   "Imm… T cells       70       6 7.67e-15 5.55e-12
#> 5 T cells | Immu… Pang… Mouse   "Imm… T cells       69       6 9.14e-15 5.55e-12
#> 6 AIZARANI_LIVER… MSig… Human   ""    AIZARAN…     120       6 2.14e-13 1.08e-10
```
