#' Read a GMT file into a data frame
#'
#' @param file A connection object or a character string (can be a URL).
#' @param geneset_label Column name for gene sets (first column of the GMT file) in the output data frame.
#' @param gene_label Column name for genes (variable columns of the GMT file) in the output data frame.
#'
#' @return A data frame with gene sets as the first column and genes as the second column (one gene per row).
#'
#' @import utils
#' @importFrom tibble enframe
#' @importFrom tidyr unnest
#' @export
#'
#' @examples
#' gmt <- "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"
#' gmt_tbl <- read_gmt(gmt)
#' head(gmt_tbl)
read_gmt <- function(file, geneset_label = "celltype", gene_label = "gene") {
  gmt_split <- strsplit(readLines(file), "\t")
  gmt_list <- lapply(gmt_split, tail, -2)
  names(gmt_list) <- sapply(gmt_split, head, 1)
  gmt_df <- tibble::enframe(gmt_list, name = geneset_label, value = gene_label)
  gmt_df <- tidyr::unnest(gmt_df, gene_label)
  gmt_df
}
