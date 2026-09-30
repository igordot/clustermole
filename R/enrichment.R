#' Cell types based on the expression of all genes
#'
#' Score cell type signatures using the full gene expression matrix.
#'
#' @param expr_mat Expression matrix (logCPMs, logFPKMs, or logTPMs) with genes
#'   as rows and clusters/populations/samples as columns.
#' @inheritParams clustermole_markers species
#' @param method Enrichment method: `gsva` (default), `ssgsea`, `singscore`, or
#'   `all` to combine ranks from all three methods. See references below.
#'
#' @return A data frame with one row per returned signature and input column:
#'
#'   - `cluster`: Input column name.
#'   - `score`: Enrichment score (higher means greater enrichment).
#'   - `score_rank`: Signature rank (lower means greater enrichment).
#'   - Signature metadata (see [clustermole_markers()]).
#'
#'   With `method = "all"`, these columns replace `score`:
#'
#'   - `score_rank_{method}`: The ranks from each method.
#'   - `score_ranks_{stat}`: Minimum, mean, and median ranks across methods.
#'
#' @references
#' Barbie, D., Tamayo, P., Boehm, J. et al. Systematic RNA interference reveals
#' that oncogenic KRAS-driven cancers require TBK1. _Nature_ 462, 108–112
#' (2009). \doi{10.1038/nature08460}
#'
#' Hänzelmann, S., Castelo, R. & Guinney, J. GSVA: Gene set variation analysis
#' for microarray and RNA-Seq data. _BMC Bioinformatics_ 14, 7 (2013).
#' \doi{10.1186/1471-2105-14-7}
#'
#' Foroutan, M., Bhuva, D.D., Lyu, R. et al. Single sample scoring of molecular
#' phenotypes. _BMC Bioinformatics_ 19, 404 (2018).
#' \doi{10.1186/s12859-018-2435-4}
#'
#' @importFrom dplyr add_count arrange distinct filter inner_join select
#' @importFrom dplyr starts_with
#' @export
#'
#' @examples
#' # my_enrichment <- clustermole_enrichment(
#' #   expr_mat = my_expr_mat, species = "hs"
#' # )
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
  markers_tbl <- filter(markers_tbl, .data$gene %in% rownames(expr_mat))
  markers_tbl <- select(markers_tbl, !"gene_original")
  markers_tbl <- distinct(markers_tbl)
  markers_tbl <- add_count(
    markers_tbl,
    .data$celltype_full,
    name = "n_genes"
  )
  markers_tbl <- filter(markers_tbl, .data$n_genes >= 5)

  # convert markers to a list
  markers_list <- split(x = markers_tbl$gene, f = markers_tbl$celltype_full)

  # create a table of cell types (without genes)
  celltypes_tbl <- select(markers_tbl, !starts_with("gene"))
  celltypes_tbl <- distinct(celltypes_tbl)

  # run the actual enrichment analysis
  scores_tbl <- get_scores(
    expr_mat = expr_mat,
    markers_list = markers_list,
    method = method
  )

  scores_tbl <- filter(scores_tbl, .data$score_rank <= 100)
  scores_tbl <- inner_join(scores_tbl, celltypes_tbl, by = "celltype_full")
  scores_tbl <- arrange(scores_tbl, .data$cluster, .data$score_rank)
  scores_tbl
}

#' @importFrom dplyr full_join select starts_with
#' @importFrom GSEABase GeneSet GeneSetCollection
#' @importFrom GSVA gsva gsvaParam ssgseaParam
#' @importFrom singscore multiScore rankGenes
#' @importFrom stats median
get_scores <- function(
  expr_mat,
  markers_list,
  method = c("gsva", "ssgsea", "singscore", "all")
) {
  method <- match.arg(method)

  if (method == "gsva" || method == "all") {
    gsva_param <- GSVA::gsvaParam(
      exprData = expr_mat,
      geneSets = markers_list,
      kcdf = "Gaussian"
    )
    scores_mat <- GSVA::gsva(gsva_param, verbose = FALSE)
    scores_tbl <- lengthen_scores(scores_mat)
    scores_gsva_tbl <- select(
      scores_tbl,
      "cluster",
      "celltype_full",
      score_rank_gsva = "score_rank"
    )
  }

  if (method == "ssgsea" || method == "all") {
    ssgsea_param <- GSVA::ssgseaParam(
      exprData = expr_mat,
      geneSets = markers_list
    )
    scores_mat <- GSVA::gsva(ssgsea_param, verbose = FALSE)
    scores_tbl <- lengthen_scores(scores_mat)
    scores_ssgsea_tbl <- select(
      scores_tbl,
      "cluster",
      "celltype_full",
      score_rank_ssgsea = "score_rank"
    )
  }

  if (method == "singscore" || method == "all") {
    markers_gsc <- Map(
      function(x, y) GSEABase::GeneSet(x, setName = y),
      markers_list,
      names(markers_list)
    )
    markers_gsc <- GSEABase::GeneSetCollection(markers_gsc)
    scores_mat <- singscore::multiScore(
      rankData = rankGenes(expr_mat),
      upSetColc = markers_gsc
    )
    scores_mat <- scores_mat$Scores
    scores_tbl <- lengthen_scores(scores_mat)
    scores_singscore_tbl <- select(
      scores_tbl,
      "cluster",
      "celltype_full",
      score_rank_singscore = "score_rank"
    )
  }

  if (method == "all") {
    # combine all scores into a single table
    scores_tbl <- scores_gsva_tbl
    scores_tbl <- full_join(
      scores_tbl,
      scores_ssgsea_tbl,
      by = c("cluster", "celltype_full")
    )
    scores_tbl <- full_join(
      scores_tbl,
      scores_singscore_tbl,
      by = c("cluster", "celltype_full")
    )
    # create a score matrix for easier stats
    scores_mat <- select(scores_tbl, starts_with("score_rank_"))
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

#' @importFrom dplyr desc group_by mutate select ungroup
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
lengthen_scores <- function(scores_mat) {
  scores_mat |>
    round(10) |>
    as_tibble(rownames = "celltype_full") |>
    gather(key = "cluster", value = "score", -"celltype_full") |>
    select("cluster", "celltype_full", "score") |>
    group_by(.data$cluster) |>
    mutate(
      score_rank = rank(desc(.data$score), ties.method = "first")
    ) |>
    ungroup()
}
