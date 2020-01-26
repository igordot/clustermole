
#' Retrieve the available cell type markers
#'
#' @param species Species for the appropriate gene symbol format: "hs" for human or "mm" for mouse.
#'
#' @return A data frame of markers with one gene per row.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' markers <- clustermole_markers()
#' head(markers)
clustermole_markers <- function(species = "hs") {
  m_tbl <- clustermole_markers_tbl
  if (species == "hs") {
    m_tbl %>%
      dplyr::select(-.data$gene_mm) %>%
      dplyr::rename(gene = .data$gene_hs)
  } else if (species == "mm") {
    m_tbl %>%
      dplyr::select(-.data$gene_hs) %>%
      dplyr::rename(gene = .data$gene_mm)
  } else {
    stop("unknown species")
  }
}

#' Perform cell type overrepresentation analysis for a set of genes
#'
#' @param genes A vector of genes.
#' @param species Species: "hs" for human or "mm" for mouse.
#'
#' @return A data frame of enrichment results with hypergeometric test p-values.
#'
#' @import methods
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
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

#' Perform cell type enrichment for a given gene expression matrix
#'
#' @param expr_mat Expression matrix (logCPMs, logFPKMs, or logTPMs) with genes as rows.
#' @param species Species: "hs" for human or "mm" for mouse.
#'
#' @return A data frame of enrichment results.
#'
#' @import methods
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom GSVA gsva
#' @export
#'
#' @examples
#' # my_enrichment <- clustermole_enrichment(expr_mat = my_expr_mat, species = "hs")
clustermole_enrichment <- function(expr_mat, species) {

  # check that the expression matrix seems reasonable
  if (!is(expr_mat, "matrix")) {
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
  markers_tbl <- clustermole_markers(species = species)
  markers_list <- split(x = markers_tbl$gene, f = markers_tbl$celltype_full)
  celltypes_tbl <-
    markers_tbl %>%
    dplyr::select(-dplyr::starts_with("gene")) %>%
    dplyr::distinct()

  # run the actual enrichment analysis
  gsva_mat <- GSVA::gsva(
    expr = expr_mat, gset.idx.list = markers_list,
    method = "gsva", kcdf = "Gaussian", parallel.sz = 1, verbose = FALSE
  )

  # clean up the enrichment table
  gsva_tbl <-
    gsva_mat %>%
    round(8) %>%
    tibble::as_tibble(rownames = "celltype_full") %>%
    tidyr::gather(key = "cluster", value = "score", -.data$celltype_full) %>%
    dplyr::select(.data$cluster, .data$celltype_full, .data$score) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::top_n(n = 50) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(celltypes_tbl, by = "celltype_full") %>%
    dplyr::arrange(.data$cluster, desc(.data$score))
  gsva_tbl
}

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
