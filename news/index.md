# Changelog

## clustermole (development version)

- Adds HPA, ScType, CellTaxonomy, CellMatch, and DISCO marker databases
  (11 total).
- Updates CellMarker to version 3.0.
- Converts known gene aliases to their canonical symbol in the source
  databases.
- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  accepts data frames and converts them to matrices automatically.
- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  removes duplicate markers before analysis.
- [`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md)
  computes species-specific signature gene counts.
- [`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md)
  excludes genes without a symbol for the requested species.
- [`clustermole_markers()`](https://igordot.github.io/clustermole/reference/clustermole_markers.md)
  removes duplicate rows caused by tied ortholog mappings.
- [`clustermole_overlaps()`](https://igordot.github.io/clustermole/reference/clustermole_overlaps.md)
  now uses a proportional threshold for its species-mismatch check.
- [`clustermole_overlaps()`](https://igordot.github.io/clustermole/reference/clustermole_overlaps.md)
  removes duplicate markers before analysis.

## clustermole 1.1.1

CRAN release: 2024-01-08

- Fixes internal GSVA and tidyselect function calls.

## clustermole 1.1.0

CRAN release: 2021-01-26

- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  gains a `singscore` method.
- [`clustermole_enrichment()`](https://igordot.github.io/clustermole/reference/clustermole_enrichment.md)
  gains a combined enrichment method.
- Updates cell type markers.

## clustermole 1.0.1

CRAN release: 2020-01-27

- Updates the documentation.

## clustermole 1.0.0

CRAN release: 2020-01-20

- Initial CRAN submission.
