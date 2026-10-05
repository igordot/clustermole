# Changelog

## clustermole (development version)

- Known gene aliases in the source databases are converted to canonical
  gene symbols.
- Adds HPA, ScType, CellTaxonomy, CellMatch, and DISCO databases.
- Updates CellMarker to version 3.
- Removes the ARCHS4 marker database because most of its signatures
  contain too many genes.
- Refreshes the markers from all source databases.
- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  gains `max_rank` to control the rank cutoff.
- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  accepts data frames and automatically converts them to matrices.
- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  removes duplicate marker genes within signatures before analysis.
- [`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md)
  computes species-specific signature gene counts.
- [`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md)
  excludes genes without a symbol for the requested species.
- [`clustermole_overlaps()`](https://igordot.github.io/clustermole/reference/clustermole_overlaps.md)
  gains `max_p` and `max_fdr` cutoffs.
- [`clustermole_overlaps()`](https://igordot.github.io/clustermole/reference/clustermole_overlaps.md)
  uses the proportion of mismatched genes for its species-mismatch
  check.
- [`clustermole_overlaps()`](https://igordot.github.io/clustermole/reference/clustermole_overlaps.md)
  removes duplicate marker genes within signatures before analysis.

## clustermole 1.1.1

CRAN release: 2024-01-08

- Fixes compatibility issues with GSVA and tidyselect.

## clustermole 1.1.0

CRAN release: 2021-01-26

- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  gains a `singscore` method.
- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  gains a combined enrichment method.
- Refreshes the markers from all source databases.

## clustermole 1.0.1

CRAN release: 2020-01-27

- Documentation is improved.

## clustermole 1.0.0

CRAN release: 2020-01-20

- Initial CRAN submission.
