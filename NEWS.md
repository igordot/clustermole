# clustermole (development version)

* Known gene aliases in the source databases are converted to canonical gene symbols.
* Adds HPA, ScType, CellTaxonomy, CellMatch, and DISCO databases.
* Updates CellMarker to version 3.
* Removes the ARCHS4 marker database because most of its signatures contain too many genes.
* Refreshes the markers from all source databases.
* `clustermole_enrichment()` gains `max_rank` to control the rank cutoff.
* `clustermole_enrichment()` accepts data frames and automatically converts them to matrices.
* `clustermole_enrichment()` removes duplicate marker genes within signatures before analysis.
* `clustermole_markers()` computes species-specific signature gene counts.
* `clustermole_markers()` excludes genes without a symbol for the requested species.
* `clustermole_overlaps()` gains `max_p` and `max_fdr` cutoffs.
* `clustermole_overlaps()` uses the proportion of mismatched genes for its species-mismatch check.
* `clustermole_overlaps()` removes duplicate marker genes within signatures before analysis.

# clustermole 1.1.1

* Fixes compatibility issues with GSVA and tidyselect.

# clustermole 1.1.0

* `clustermole_enrichment()` gains a `singscore` method.
* `clustermole_enrichment()` gains a combined enrichment method.
* Refreshes the markers from all source databases.

# clustermole 1.0.1

* Documentation is improved.

# clustermole 1.0.0

* Initial CRAN submission.
