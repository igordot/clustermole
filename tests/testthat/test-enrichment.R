
# markers table
markers_hs_tbl <- clustermole_markers(species = "hs")
markers_mm_tbl <- clustermole_markers(species = "mm")

# expression matrix
n_genes <- 10000
expr_mat <- matrix(rnbinom(n_genes * 5, size = 1, mu = 10), nrow = n_genes, ncol = 5)
colnames(expr_mat) <- c("C1", "C2", "C3", "C4", "C5")
cpm_mat <- t(t(expr_mat) / colSums(expr_mat)) * 1e6
log_cpm_mat <- log2(cpm_mat + 0.1)

# add human gene names to the expression matrix
gene_names_hs <- unique(markers_hs_tbl$gene)
gene_names_hs <- sort(sample(gene_names_hs, n_genes))
rownames(expr_mat) <- gene_names_hs
rownames(cpm_mat) <- gene_names_hs
rownames(log_cpm_mat) <- gene_names_hs

test_that("clustermole_enrichment() wrong input", {
  expect_error(clustermole_enrichment(as.data.frame(log_cpm_mat), species = "hs"))
  expect_error(clustermole_enrichment(log_cpm_mat[1:100, ], species = "hs"))
  expect_error(clustermole_enrichment(log_cpm_mat[, 1:3], species = "hs"))
  expect_error(clustermole_enrichment(cpm_mat, species = "hs"))
  expect_error(clustermole_enrichment(log_cpm_mat, species = "hs", method = "x"))
})

# default (gsva)
test_that("clustermole_enrichment() human input default method", {
  enrich_hs_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "hs")
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
})

# gsva
test_that("clustermole_enrichment() human input gsva method", {
  enrich_hs_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "hs", method = "gsva")
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
})

# ssgsea
test_that("clustermole_enrichment() human input ssgsea method", {
  enrich_hs_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "hs", method = "ssgsea")
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
})

# singscore
test_that("clustermole_enrichment() human input singscore method", {
  enrich_hs_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "hs", method = "singscore")
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
})

# combined
test_that("clustermole_enrichment() human input all combined method", {
  enrich_hs_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "hs", method = "all")
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
})

# add mouse gene names to the expression matrix
gene_names_mm <- unique(markers_mm_tbl$gene)
gene_names_mm <- sort(sample(gene_names_mm, n_genes))
rownames(expr_mat) <- gene_names_mm
rownames(cpm_mat) <- gene_names_mm
rownames(log_cpm_mat) <- gene_names_mm

test_that("clustermole_enrichment() mouse input", {
  enrich_mm_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "mm")
  expect_s3_class(enrich_mm_tbl, "tbl_df")
  expect_gt(nrow(enrich_mm_tbl), 100)
  expect_equal(length(unique(enrich_mm_tbl$cluster)), 5)
})
