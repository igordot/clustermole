#' Cell types based on overlap of marker genes
#'
#' Perform overrepresentation analysis for a set of genes compared to all cell
#' type signatures.
#'
#' @param genes A character vector of gene symbols.
#' @inheritParams clustermole_markers species
#'
#' @return A data frame with one row per returned signature:
#'
#'   - `overlap`: Unique gene count shared by the input and signature.
#'   - `p_value`: Hypergeometric test p-value.
#'   - `fdr`: Benjamini-Hochberg adjusted p-value across all tested signatures.
#'   - `n_genes`: Unique gene count in the signature for the requested species.
#'   - Signature metadata (see [clustermole_markers()]).
#'
#' @importFrom dplyr arrange distinct filter inner_join select starts_with
#' @importFrom stats p.adjust phyper
#' @importFrom tibble as_tibble
#' @export
#'
#' @examples
#' my_genes <- c("CD2", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "LTB")
#' my_overlaps <- clustermole_overlaps(genes = my_genes, species = "hs")
#' head(my_overlaps)
clustermole_overlaps <- function(genes, species) {
  # check that the genes vector seems reasonable
  if (!is(genes, "character")) {
    stop("`genes` is not a character vector")
  }
  genes <- unique(sort(genes))
  if (length(genes) < 5) {
    stop("input should be at least 5 genes")
  }
  if (length(genes) > 5000) {
    stop("input should be less than 5,000 genes")
  }

  # retrieve markers
  markers_tbl <- clustermole_markers(species = species)
  markers_tbl <- select(markers_tbl, !"gene_original")
  markers_tbl <- distinct(markers_tbl)
  markers_list <- split(x = markers_tbl$gene, f = markers_tbl$celltype_full)
  celltypes_tbl <- select(markers_tbl, !starts_with("gene"))
  celltypes_tbl <- distinct(celltypes_tbl)

  # check that input genes overlap marker genes for a given species
  all_genes <- unique(markers_tbl$gene)
  input_genes <- genes
  genes <- intersect(genes, all_genes)
  if (length(genes) < max(3, length(input_genes) * 0.2)) {
    problematic_genes <- setdiff(input_genes, genes)
    stop(
      "large fraction of input genes are not known (possibly wrong species): ",
      toString(problematic_genes)
    )
  }

  # run the enrichment analysis
  overlaps_mat <-
    sapply(markers_list, function(celltype_genes) {
      n_overlap <- length(intersect(genes, celltype_genes))
      n_query <- length(genes)
      n_celltype <- length(celltype_genes)
      n_all <- length(all_genes)
      # phyper(success-in-sample, success-in-bg, fail-in-bg, sample-size)
      p_val <- phyper(
        n_overlap - 1,
        n_celltype,
        n_all - n_celltype,
        n_query,
        lower.tail = FALSE
      )
      c("overlap" = n_overlap, "p_value" = p_val, "fdr" = 1)
    })
  overlaps_mat <- t(overlaps_mat)
  overlaps_mat[, "fdr"] <- p.adjust(overlaps_mat[, "p_value"], method = "fdr")

  # clean up the enrichment table
  overlaps_tbl <- as_tibble(overlaps_mat, rownames = "celltype_full")
  overlaps_tbl <- filter(overlaps_tbl, .data$p_value < 0.05)
  overlaps_tbl <- inner_join(celltypes_tbl, overlaps_tbl, by = "celltype_full")
  overlaps_tbl <- arrange(
    overlaps_tbl,
    .data$fdr,
    .data$p_value,
    .data$celltype_full
  )
  overlaps_tbl
}
