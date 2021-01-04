#' Available cell type markers
#'
#' Retrieve the full list of cell type markers in the \code{clustermole} database.
#'
#' @param species Species: \code{hs} for human or \code{mm} for mouse.
#'
#' @return A data frame of cell type markers (one gene per row).
#'
#' @import dplyr
#' @export
#'
#' @examples
#' markers <- clustermole_markers()
#' head(markers)
clustermole_markers <- function(species = c("hs", "mm")) {
  species <- match.arg(species)
  m_tbl <- clustermole_markers_tbl
  if (species == "hs") {
    m_tbl %>%
      dplyr::select(-.data$gene_mm) %>%
      dplyr::rename(gene = .data$gene_hs)
  } else if (species == "mm") {
    m_tbl %>%
      dplyr::select(-.data$gene_hs) %>%
      dplyr::rename(gene = .data$gene_mm)
  }
}
