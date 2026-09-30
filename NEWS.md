# clustermole (development version)

* Adds HPA, ScType, CellTaxonomy, CellMatch, and DISCO marker databases (11 total).
* Updates CellMarker to version 3.0.
* Converts known gene aliases to their canonical symbol in the source databases.
* `clustermole_enrichment()` accepts data frames and converts them to matrices automatically.
* `clustermole_enrichment()` removes duplicate markers before analysis.
* `clustermole_markers()` computes species-specific signature gene counts.
* `clustermole_markers()` excludes genes without a symbol for the requested species.
* `clustermole_markers()` removes duplicate rows caused by tied ortholog mappings.
* `clustermole_overlaps()` now uses a proportional threshold for its species-mismatch check.
* `clustermole_overlaps()` removes duplicate markers before analysis.

# clustermole 1.1.1

* Fixes internal GSVA and tidyselect function calls.

# clustermole 1.1.0

* `clustermole_enrichment()` gains a `singscore` method.
* `clustermole_enrichment()` gains a combined enrichment method.
* Updates cell type markers.

# clustermole 1.0.1

* Updates the documentation.

# clustermole 1.0.0

* Initial CRAN submission.
