#' Available cell type markers
#'
#' Retrieve the full list of cell type markers in the `clustermole` database.
#'
#' @param species Species: `hs` for human or `mm` for mouse.
#'
#' @return A data frame of cell type markers (one gene per row).
#'
#' @export
#'
#' @examples
#' markers <- clustermole_markers()
#' head(markers)
clustermole_markers <- function(species = c("hs", "mm")) {
  species <- match.arg(species)
  m_tbl <- clustermole_markers_tbl
  if (species == "hs") {
    m_tbl <- dplyr::rename(m_tbl, gene = "gene_hs")
    m_tbl <- dplyr::select(m_tbl, !"gene_mm")
  } else if (species == "mm") {
    m_tbl <- dplyr::rename(m_tbl, gene = "gene_mm")
    m_tbl <- dplyr::select(m_tbl, !"gene_hs")
  }
  # a tied ortholog can repeat a row
  dplyr::distinct(m_tbl)
}
