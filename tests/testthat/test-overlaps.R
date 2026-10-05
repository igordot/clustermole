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

test_that("non-numeric significance cutoffs error", {
  genes <- c("CD2", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "LTB")
  for (cutoff in list(NULL, "0.05", TRUE, 1i)) {
    expect_error(
      clustermole_overlaps(genes, "hs", max_p = cutoff),
      "`max_p` is not numeric",
      fixed = TRUE
    )
    expect_error(
      clustermole_overlaps(genes, "hs", max_fdr = cutoff),
      "`max_fdr` is not numeric",
      fixed = TRUE
    )
  }
})

test_that("invalid numeric significance cutoffs error", {
  genes <- c("CD2", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "LTB")
  for (cutoff in list(NA_real_, -0.1, 1.1)) {
    expect_error(
      clustermole_overlaps(genes, "hs", max_p = cutoff),
      "`max_p` must be a number between 0 and 1",
      fixed = TRUE
    )
  }
  for (cutoff in list(NA_real_, -0.1, 1.1)) {
    expect_error(
      clustermole_overlaps(genes, "hs", max_fdr = cutoff),
      "`max_fdr` must be a number between 0 and 1",
      fixed = TRUE
    )
  }
})

test_that("significance cutoffs preserve p-values and FDR", {
  genes <- c("CD2", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "LTB")
  all_overlaps <- clustermole_overlaps(genes, "hs", max_p = 1)
  default_overlaps <- clustermole_overlaps(genes, "hs")
  expect_equal(
    default_overlaps,
    all_overlaps[all_overlaps$p_value <= 0.05, ]
  )
  n_signatures <- dplyr::n_distinct(clustermole_markers("hs")$celltype_full)
  expect_equal(
    all_overlaps$fdr,
    p.adjust(all_overlaps$p_value, method = "BH", n = n_signatures)
  )

  fdr_cutoff <- default_overlaps$fdr[ceiling(nrow(default_overlaps) / 2)]
  filtered <- clustermole_overlaps(genes, "hs", max_fdr = fdr_cutoff)
  expect_equal(filtered, default_overlaps[default_overlaps$fdr <= fdr_cutoff, ])
  expect_gt(nrow(filtered), 0)
  expect_lt(nrow(filtered), nrow(default_overlaps))

  p_cutoff <- filtered$p_value[ceiling(nrow(filtered) / 2)]
  expect_equal(
    clustermole_overlaps(
      genes,
      "hs",
      max_p = p_cutoff,
      max_fdr = fdr_cutoff
    ),
    filtered[filtered$p_value <= p_cutoff, ]
  )
  expect_true(all(clustermole_overlaps(genes, "hs", max_p = 0)$p_value == 0))
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
