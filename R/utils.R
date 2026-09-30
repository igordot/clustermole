#' Read a GMT file into a data frame
#'
#' @param file A file path, URL, or connection.
#' @param geneset_label Output column name for gene sets (GMT column 1).
#' @param gene_label Output column name for genes (GMT columns 3 onward).
#'
#' @return A data frame with gene sets and genes, one gene per row.
#'
#' @importFrom tibble enframe
#' @importFrom tidyr unnest
#' @export
#'
#' @examples
#' \dontrun{
#' gmt <- "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"
#' gmt_tbl <- read_gmt(gmt)
#' head(gmt_tbl)
#' }
read_gmt <- function(file, geneset_label = "celltype", gene_label = "gene") {
  gmt_split <- strsplit(readLines(file), "\t")
  gmt_list <- lapply(gmt_split, tail, -2)
  names(gmt_list) <- sapply(gmt_split, head, 1)
  gmt_df <- tibble::enframe(gmt_list, name = geneset_label, value = gene_label)
  gmt_df <- unnest(gmt_df, all_of(gene_label))
  gmt_df
}
