
#' Retrieve a table of all cell type markers
#'
#' @return a data frame of markers with one gene per row
#' @export
#'
#' @examples
#' markers_tbl <- clustermole_markers()
clustermole_markers <- function() {
  markers
}

#' Perform cell type enrichment for a given gene expression matrix
#'
#' @param expr_mat expression matrix (logCPMs, logFPKMs, or logTPMs)
#' @param species species ("hs" for human or "mm" for mouse)
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
  if (nrow(expr_mat) < 5000) {
    stop("expression matrix does not appear to be complete (too few genes)")
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

  # create a markers list for gene set enrichment
  if (species == "hs") {
    markers_list <- split(x = markers_tbl$gene_h, f = markers_tbl$celltype_long)
  }
  if (species == "mm") {
    markers_list <- split(x = markers_tbl$gene_m, f = markers_tbl$celltype_long)
  }

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
