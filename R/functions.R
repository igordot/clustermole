
#' Retrieve the cell type markers data frame
#'
#' @return a data frame of markers with one gene per row
#' @export
clustermole_markers <- function() {
  markers
}

#' Retrieve the cell type markers data frame
#'
#' @param expr_mat normalized expression matrix (log-scale)
#' @param species species ("hs" or "mm")
#'
#' @return a data frame of enrichment results
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom GSVA gsva
#' @export
clustermole_enrichment <- function(expr_mat, species) {

  # check that the expression matrix seems reasonable
  if (class(expr_mat) != "matrix") {
    stop("expression matrix is not a matrix")
  }
  if (max(expr_mat) > 100) {
    stop("expression values do not appear to be log-scaled")
  }

  # remove genes without expression values
  expr_mat <- expr_mat[rowMeans(expr_mat) > 0, ]

  # retrieve markers
  markers_tbl <- clustermole_markers()
  celltypes_tbl <-
    markers_tbl %>%
    select(.data$db, .data$species, .data$organ, .data$celltype, .data$celltype_long, .data$n_genes) %>%
    distinct()
  markers_tbl <- clustermole_markers() %>% select(.data$celltype_long, .data$gene_h, .data$gene_m)
  if (species == "hs") {
    markers_tbl$gene <- markers_tbl$gene_h
  }
  if (species == "mm") {
    markers_tbl$gene <- markers_tbl$gene_m
  }

  # convert markers to list format
  markers_list <- split(x = markers_tbl$gene, f = markers_tbl$celltype_long)

  # run the actual enrichment analysis
  gsva_mat <- gsva(
    expr = expr_mat, gset.idx.list = markers_list,
    method = "gsva", kcdf = "Gaussian", verbose = FALSE
  )

  # clean up the enrichment table
  gsva_tbl <-
    gsva_mat %>%
    round(8) %>%
    as_tibble(rownames = "celltype_long") %>%
    gather(key = "cluster", value = "score", -.data$celltype_long) %>%
    select(.data$cluster, .data$celltype_long, .data$score) %>%
    group_by(.data$cluster) %>%
    top_n(n = 50) %>%
    ungroup() %>%
    inner_join(celltypes_tbl, by = "celltype_long") %>%
    arrange(.data$cluster, desc(.data$score))
  gsva_tbl
}
