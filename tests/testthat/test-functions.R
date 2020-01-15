
# clustermole_markers -----

markers_tbl <- clustermole_markers()
test_that("clustermole_markers() default params", {
  expect_s3_class(markers_tbl, "tbl_df")
  expect_gt(nrow(markers_tbl), 100000)
})

test_that("clustermole_markers() wrong input", {
  expect_error(clustermole_markers(species = ""))
  expect_error(clustermole_markers(species = "x"))
  expect_error(clustermole_markers(species = "?"))
  expect_error(clustermole_markers(species = "*"))
})

markers_hs_tbl <- clustermole_markers(species = "hs")
test_that("clustermole_markers() human input", {
  expect_s3_class(markers_hs_tbl, "tbl_df")
  expect_gt(nrow(markers_hs_tbl), 100000)
})

markers_mm_tbl <- clustermole_markers(species = "mm")
test_that("clustermole_markers() mouse input", {
  expect_s3_class(markers_mm_tbl, "tbl_df")
  expect_gt(nrow(markers_mm_tbl), 100000)
})

# clustermole_overlaps -----

# gene list for human overrepresentation tests
gene_names <- unique(markers_hs_tbl$gene)
gene_names <- sample(gene_names)

test_that("clustermole_overlaps() wrong input", {
  expect_error(clustermole_overlaps(gene_names[1:3], species = "hs"))
  expect_error(clustermole_overlaps(gene_names[1:10000], species = "hs"))
  expect_error(clustermole_overlaps(c(gene_names[1:3], "X", "Y", "Z"), species = "hs"))
  expect_error(clustermole_overlaps(as.list(gene_names[1:10]), species = "hs"))
})

test_that("clustermole_overlaps() human input", {
  overlap_tbl <- clustermole_overlaps(genes = gene_names[1:50], species = "hs")
  expect_s3_class(overlap_tbl, "tbl_df")
  expect_gt(nrow(overlap_tbl), 100)
})

# gene list for mouse overrepresentation tests
gene_names <- unique(markers_mm_tbl$gene)
gene_names <- sample(gene_names)

test_that("clustermole_overlaps() mouse input", {
  overlap_tbl <- clustermole_overlaps(genes = gene_names[1:50], species = "mm")
  expect_s3_class(overlap_tbl, "tbl_df")
  expect_gt(nrow(overlap_tbl), 100)
})

# clustermole_enrichment -----

# expression matrix for enrichment tests
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
  expect_error(clustermole_enrichment(cpm_mat, species = "hs"))
})

enrich_hs_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "hs")
test_that("clustermole_enrichment() human input", {
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

enrich_mm_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "mm")
test_that("clustermole_enrichment() mouse input", {
  expect_s3_class(enrich_mm_tbl, "tbl_df")
  expect_gt(nrow(enrich_mm_tbl), 100)
  expect_equal(length(unique(enrich_mm_tbl$cluster)), 5)
})

# test read_gmt()
test_that("read_gmt() output", {
  gmt_tbl <- read_gmt(file = "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt")
  expect_s3_class(gmt_tbl, "tbl_df")
  expect_gt(nrow(gmt_tbl), 1000)
})
