#' Cell types based on the expression of all genes
#'
#' Perform enrichment of cell type signatures based on the full gene expression matrix.
#'
#' @param expr_mat Expression matrix (logCPMs, logFPKMs, or logTPMs) with genes as rows and clusters/populations/samples as columns.
#' @param species Species: \code{hs} for human or \code{mm} for mouse.
#' @param method Enrichment method: \code{ssgsea}, \code{gsva}, \code{singscore}, or \code{all}. The method to use for the estimation of gene set enrichment scores. The options are ssGSEA (Barbie et al, 2009), GSVA (Hänzelmann et al, 2013), singscore (Foroutan et al, 2018), or a combination of all three methods.
#'
#' @return A data frame of enrichment results.
#'
#' @references
#' Barbie, D., Tamayo, P., Boehm, J. et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. \emph{Nature} 462, 108–112 (2009). \doi{10.1038/nature08460}
#'
#' Hänzelmann, S., Castelo, R. & Guinney, J. GSVA: Gene set variation analysis for microarray and RNA-Seq data. \emph{BMC Bioinformatics} 14, 7 (2013). \doi{10.1186/1471-2105-14-7}
#'
#' Foroutan, M., Bhuva, D.D., Lyu, R. et al. Single sample scoring of molecular phenotypes. \emph{BMC Bioinformatics} 19, 404 (2018). \doi{10.1186/s12859-018-2435-4}
#'
#' @import methods
#' @import dplyr
#' @export
#'
#' @examples
#' # my_enrichment <- clustermole_enrichment(expr_mat = my_expr_mat, species = "hs")
clustermole_enrichment <- function(expr_mat, species, method = "gsva") {

  # check that the expression matrix seems reasonable
  if (!is(expr_mat, "matrix")) {
    stop("expression matrix is not a matrix")
  }
  if (nrow(expr_mat) < 5000) {
    stop("expression matrix does not appear to be complete (too few rows)")
  }
  if (ncol(expr_mat) < 5) {
    stop("expression matrix does not appear to be complete (too few columns)")
  }
  if (max(expr_mat) > 100) {
    stop("expression values do not appear to be log-scaled")
  }

  # remove genes without low or not variable values
  expr_mat <- expr_mat[rowMeans(expr_mat) > min(expr_mat), ]

  # retrieve markers and filter for genes present in the expression table
  markers_tbl <- clustermole_markers(species = species)
  markers_tbl <- dplyr::filter(markers_tbl, .data$gene %in% rownames(expr_mat))
  markers_tbl <- dplyr::add_count(markers_tbl, .data$celltype_full, name = "n_genes")
  markers_tbl <- dplyr::filter(markers_tbl, .data$n_genes >= 5)

  # convert markers to a list
  markers_list <- dplyr::distinct(markers_tbl, .data$celltype_full, .data$gene)
  markers_list <- split(x = markers_list$gene, f = markers_list$celltype_full)

  # create a table of cell types (without genes)
  celltypes_tbl <-
    markers_tbl %>%
    dplyr::select(-dplyr::starts_with("gene")) %>%
    dplyr::distinct()

  # run the actual enrichment analysis
  scores_tbl <- get_scores(expr_mat = expr_mat, markers_list = markers_list, method = method)

  scores_tbl <- scores_tbl %>%
    dplyr::filter(.data$score_rank <= 100) %>%
    dplyr::inner_join(celltypes_tbl, by = "celltype_full") %>%
    dplyr::arrange(.data$cluster, .data$score_rank)
  scores_tbl
}

#' @import dplyr
#' @importFrom stats median
#' @importFrom GSVA gsva
#' @importFrom GSEABase GeneSet GeneSetCollection
#' @importFrom singscore rankGenes multiScore
get_scores <- function(expr_mat, markers_list, method = c("gsva", "ssgsea", "singscore", "all")) {
  method <- match.arg(method)

  if (method == "gsva" || method == "all") {
    scores_mat <- GSVA::gsva(
      expr = expr_mat, gset.idx.list = markers_list,
      method = "gsva", kcdf = "Gaussian", parallel.sz = 1, verbose = FALSE
    )
    scores_tbl <- lengthen_scores(scores_mat)
    scores_gsva_tbl <- dplyr::select(scores_tbl, .data$cluster, .data$celltype_full, score_rank_gsva = .data$score_rank)
  }

  if (method == "ssgsea" || method == "all") {
    scores_mat <- GSVA::gsva(
      expr = expr_mat, gset.idx.list = markers_list,
      method = "ssgsea", kcdf = "Gaussian", parallel.sz = 1, verbose = FALSE
    )
    scores_tbl <- lengthen_scores(scores_mat)
    scores_ssgsea_tbl <- dplyr::select(scores_tbl, .data$cluster, .data$celltype_full, score_rank_ssgsea = .data$score_rank)
  }

  if (method == "singscore" || method == "all") {
    markers_gsc <- Map(function(x, y) GSEABase::GeneSet(x, setName = y), markers_list, names(markers_list))
    markers_gsc <- GSEABase::GeneSetCollection(markers_gsc)
    scores_mat <- singscore::multiScore(rankData = rankGenes(expr_mat), upSetColc = markers_gsc)
    scores_mat <- scores_mat$Scores
    scores_tbl <- lengthen_scores(scores_mat)
    scores_singscore_tbl <- dplyr::select(scores_tbl, .data$cluster, .data$celltype_full, score_rank_singscore = .data$score_rank)
  }

  if (method == "all") {
    # combine all scores into a single table
    scores_tbl <- scores_gsva_tbl
    scores_tbl <- dplyr::full_join(scores_tbl, scores_ssgsea_tbl, by = c("cluster", "celltype_full"))
    scores_tbl <- dplyr::full_join(scores_tbl, scores_singscore_tbl, by = c("cluster", "celltype_full"))
    # create a score matrix for easier stats
    scores_mat <- dplyr::select(scores_tbl, dplyr::starts_with("score_rank_"))
    scores_mat <- as.matrix(scores_mat)
    # calculate stats
    scores_tbl$score_ranks_min <- apply(scores_mat, 1, min)
    scores_tbl$score_ranks_mean <- round(apply(scores_mat, 1, mean), 3)
    scores_tbl$score_ranks_median <- round(apply(scores_mat, 1, median), 3)
    # set the average rank as the default rank
    # not using the median as ssGSEA and singscore ranks tend to correlate well
    scores_tbl$score_rank <- scores_tbl$score_ranks_mean
  }

  scores_tbl
}

#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
lengthen_scores <- function(scores_mat) {
  scores_mat %>%
    round(10) %>%
    tibble::as_tibble(rownames = "celltype_full") %>%
    tidyr::gather(key = "cluster", value = "score", -.data$celltype_full) %>%
    dplyr::select(.data$cluster, .data$celltype_full, .data$score) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::mutate(score_rank = rank(desc(.data$score), ties.method = "first")) %>%
    dplyr::ungroup()
}
