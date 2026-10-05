#' Available cell type markers
#'
#' Retrieve cell type markers from the clustermole database.
#'
#' @details
#' The `gene` column uses the official NCBI symbols for the requested species.
#' The package maps the aliases in each source database to these symbols. The
#' `gene_original` column preserves the symbol from the source. The source
#' databases include both human and mouse cell type markers. The package maps
#' them across species with ortholog data from the Alliance of Genome Resources
#' and keeps only reciprocal best matches. The output does not include genes
#' that have no ortholog in the requested species and does not include
#' signatures with fewer than five genes.
#'
#' @param species Gene symbol species: `hs` for human or `mm` for mouse.
#'
#' @return A data frame of cell type markers with these columns:
#'
#'   - `gene`: Canonical gene symbol for the requested species.
#'   - `gene_original`: Original source gene symbol.
#'   - `celltype_full`: Full cell type signature identifier.
#'   - `db`: Source database.
#'   - `celltype`: Cell type label.
#'   - `organ`: Organ label.
#'   - `species`: Source signature species, if known.
#'   - `n_genes`: Gene count per signature.
#'
#' @importFrom dplyr distinct filter mutate n_distinct rename select
#' @importFrom tidyr drop_na
#' @export
#'
#' @examples
#' markers <- clustermole_markers()
#' head(markers)
clustermole_markers <- function(species = c("hs", "mm")) {
  species <- match.arg(species)
  tbl <- clustermole_markers_tbl
  if (species == "hs") {
    tbl <- rename(tbl, gene = "gene_hs")
    tbl <- select(tbl, !"gene_mm")
  } else if (species == "mm") {
    tbl <- rename(tbl, gene = "gene_mm")
    tbl <- select(tbl, !"gene_hs")
  }
  tbl <- drop_na(tbl, "gene")

  # a tied ortholog can repeat a row
  tbl <- distinct(tbl)

  # n_distinct() avoids overcounting genes repeated by different aliases
  tbl <- mutate(tbl, n_genes = n_distinct(.data$gene), .by = "celltype_full")
  filter(tbl, .data$n_genes >= 5)
}
