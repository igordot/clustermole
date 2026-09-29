# markers table
markers_hs_tbl <- clustermole_markers(species = "hs")
markers_mm_tbl <- clustermole_markers(species = "mm")

# gene list
gene_names <- unique(markers_hs_tbl$gene)
gene_names <- sample(gene_names)

# generate strings that look like real genes
fake_gene_names <- paste0(
  sample(LETTERS, 1000, replace = TRUE),
  sample(LETTERS, 1000, replace = TRUE),
  sample(LETTERS, 1000, replace = TRUE),
  sample(1:9, 1000, replace = TRUE)
)
fake_gene_names <- setdiff(fake_gene_names, markers_hs_tbl$gene)
fake_gene_names <- sample(fake_gene_names)

test_that("invalid gene list input errors", {
  expect_gt(length(fake_gene_names), 970)
  expect_error(clustermole_overlaps(gene_names[1:3], species = "hs"))
  expect_error(clustermole_overlaps(gene_names[1:10000], species = "hs"))
  expect_error(clustermole_overlaps(fake_gene_names[1:100], species = "hs"))
  expect_error(clustermole_overlaps(
    c(gene_names[1:2], fake_gene_names[1:5]),
    species = "hs"
  ))
  expect_error(clustermole_overlaps(
    c(gene_names[1:5], fake_gene_names[1:25]),
    species = "hs"
  ))
  expect_error(clustermole_overlaps(as.list(gene_names[1:10]), species = "hs"))
})

test_that("human gene list returns overlap results", {
  overlap_tbl <- clustermole_overlaps(genes = gene_names[1:50], species = "hs")
  expect_s3_class(overlap_tbl, "tbl_df")
  expect_gt(nrow(overlap_tbl), 1)
})

test_that("human gene list with fake genes returns overlap results", {
  overlap_tbl <- clustermole_overlaps(
    genes = c(gene_names[1:10], fake_gene_names[1:20]),
    species = "hs"
  )
  expect_s3_class(overlap_tbl, "tbl_df")
  expect_gt(nrow(overlap_tbl), 1)
})

# gene list for mouse overrepresentation tests
gene_names <- unique(markers_mm_tbl$gene)
gene_names <- sample(gene_names)

test_that("mouse gene list returns overlap results", {
  overlap_tbl <- clustermole_overlaps(genes = gene_names[1:50], species = "mm")
  expect_s3_class(overlap_tbl, "tbl_df")
  expect_gt(nrow(overlap_tbl), 1)
})
