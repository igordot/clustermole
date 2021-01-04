#' Cell types based on overlap of marker genes
#'
#' Perform overrepresentation analysis for a set of genes compared to all cell type signatures.
#'
#' @param genes A vector of genes.
#' @param species Species: \code{hs} for human or \code{mm} for mouse.
#'
#' @return A data frame of enrichment results with hypergeometric test p-values.
#'
#' @import methods
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom stats phyper p.adjust
#' @export
#'
#' @examples
#' my_genes <- c("CD2", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "LTB")
#' my_overlaps <- clustermole_overlaps(genes = my_genes, species = "hs")
#' head(my_overlaps)
clustermole_overlaps <- function(genes, species) {

  # check that the genes vector seems reasonable
  if (!is(genes, "character")) {
    stop("genes is not a character vector")
  }
  genes <- unique(genes)
  if (length(genes) < 5) {
    stop("too few genes")
  }
  if (length(genes) > 5000) {
    stop("too many genes")
  }

  # retrieve markers
  markers_tbl <- clustermole_markers(species = species)
  markers_list <- split(x = markers_tbl$gene, f = markers_tbl$celltype_full)
  celltypes_tbl <-
    markers_tbl %>%
    dplyr::select(-dplyr::starts_with("gene")) %>%
    dplyr::distinct()

  # check that input genes overlap species genes
  all_genes <- unique(markers_tbl$gene)
  genes <- intersect(genes, all_genes)
  if (length(genes) < 5) {
    stop("the genes do not appear to correspond to the given species")
  }

  # run the enrichment analysis
  overlaps_mat <-
    sapply(markers_list, function(celltype_genes) {
      n_overlap <- length(intersect(genes, celltype_genes))
      n_query <- length(genes)
      n_celltype <- length(celltype_genes)
      n_all <- length(all_genes)
      # phyper(success-in-sample, success-in-bg, fail-in-bg, sample-size)
      p_val <- phyper(n_overlap - 1, n_celltype, n_all - n_celltype, n_query, lower.tail = FALSE)
      c("overlap" = n_overlap, "p_value" = p_val, "fdr" = 1)
    })
  overlaps_mat <- t(overlaps_mat)
  overlaps_mat[, "fdr"] <- p.adjust(overlaps_mat[, "p_value"], method = "fdr")

  # clean up the enrichment table
  overlaps_tbl <- tibble::as_tibble(overlaps_mat, rownames = "celltype_full")
  overlaps_tbl <- dplyr::filter(overlaps_tbl, .data$p_value < 0.05)
  overlaps_tbl <- dplyr::inner_join(celltypes_tbl, overlaps_tbl, by = "celltype_full")
  overlaps_tbl <- dplyr::arrange(overlaps_tbl, .data$fdr, .data$p_value, .data$celltype_full)
  overlaps_tbl
}
