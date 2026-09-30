# Cell types based on the expression of all genes

Score cell type signatures using the full gene expression matrix.

## Usage

``` r
clustermole_enrichment(expr_mat, species, method = "gsva")
```

## Arguments

- expr_mat:

  Expression matrix (logCPMs, logFPKMs, or logTPMs) with genes as rows
  and clusters/populations/samples as columns.

- species:

  Gene symbol species: `hs` for human or `mm` for mouse.

- method:

  Enrichment method: `gsva` (default), `ssgsea`, `singscore`, or `all`
  to combine ranks from all three methods. See references below.

## Value

A data frame with one row per returned signature and input column:

- `cluster`: Input column name.

- `score`: Enrichment score (higher means greater enrichment).

- `score_rank`: Signature rank (lower means greater enrichment).

- Signature metadata (see
  [`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md)).

With `method = "all"`, these columns replace `score`:

- `score_rank_{method}`: The ranks from each method.

- `score_ranks_{stat}`: Minimum, mean, and median ranks across methods.

## References

Barbie, D., Tamayo, P., Boehm, J. et al. Systematic RNA interference
reveals that oncogenic KRAS-driven cancers require TBK1. *Nature* 462,
108–112 (2009).
[doi:10.1038/nature08460](https://doi.org/10.1038/nature08460)

Hänzelmann, S., Castelo, R. & Guinney, J. GSVA: Gene set variation
analysis for microarray and RNA-Seq data. *BMC Bioinformatics* 14, 7
(2013).
[doi:10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7)

Foroutan, M., Bhuva, D.D., Lyu, R. et al. Single sample scoring of
molecular phenotypes. *BMC Bioinformatics* 19, 404 (2018).
[doi:10.1186/s12859-018-2435-4](https://doi.org/10.1186/s12859-018-2435-4)

## Examples

``` r
# my_enrichment <- clustermole_enrichment(
#   expr_mat = my_expr_mat, species = "hs"
# )
```
