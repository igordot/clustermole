# markers table
markers_hs_tbl <- clustermole_markers(species = "hs")
markers_mm_tbl <- clustermole_markers(species = "mm")

# gene list
gene_names <- unique(markers_hs_tbl$gene)
gene_names <- sample(gene_names)

# generate strings that look like real genes
fake_gene_names <- paste0(sample(LETTERS), sample(LETTERS), sample(LETTERS), sample(1:9, 1000, replace = TRUE))
fake_gene_names <- sample(fake_gene_names)

test_that("clustermole_overlaps() wrong input", {
  expect_error(clustermole_overlaps(gene_names[1:3], species = "hs"))
  expect_error(clustermole_overlaps(gene_names[1:10000], species = "hs"))
  expect_error(clustermole_overlaps(fake_gene_names[1:100], species = "hs"))
  expect_error(clustermole_overlaps(c(gene_names[1:2], fake_gene_names[1:5]), species = "hs"))
  expect_error(clustermole_overlaps(c(gene_names[1:5], fake_gene_names[1:25]), species = "hs"))
  expect_error(clustermole_overlaps(as.list(gene_names[1:10]), species = "hs"))
})

test_that("clustermole_overlaps() human input", {
  overlap_tbl <- clustermole_overlaps(genes = gene_names[1:50], species = "hs")
  expect_s3_class(overlap_tbl, "tbl_df")
  expect_gt(nrow(overlap_tbl), 1)
})

test_that("clustermole_overlaps() human input with fake genes", {
  overlap_tbl <- clustermole_overlaps(genes = c(gene_names[1:10], fake_gene_names[1:20]), species = "hs")
  expect_s3_class(overlap_tbl, "tbl_df")
  expect_gt(nrow(overlap_tbl), 1)
})

# gene list for mouse overrepresentation tests
gene_names <- unique(markers_mm_tbl$gene)
gene_names <- sample(gene_names)

test_that("clustermole_overlaps() mouse input", {
  overlap_tbl <- clustermole_overlaps(genes = gene_names[1:50], species = "mm")
  expect_s3_class(overlap_tbl, "tbl_df")
  expect_gt(nrow(overlap_tbl), 1)
})
