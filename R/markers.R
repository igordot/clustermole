#' Available cell type markers
#'
#' Retrieve the full list of cell type markers in the \code{clustermole} database.
#'
#' @param species Species: \code{hs} for human or \code{mm} for mouse.
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
    names(m_tbl)[names(m_tbl) == "gene_hs"] <- "gene"
    m_tbl[, !(names(m_tbl) %in% "gene_mm")]
  } else if (species == "mm") {
    names(m_tbl)[names(m_tbl) == "gene_mm"] <- "gene"
    m_tbl[, !(names(m_tbl) %in% "gene_hs")]
  }
}
