# markers table
markers_hs_tbl <- clustermole_markers(species = "hs")
markers_mm_tbl <- clustermole_markers(species = "mm")

# expression matrix
n_genes <- 10000
expr_mat <- matrix(
  rnbinom(n_genes * 5, size = 1, mu = 10),
  nrow = n_genes,
  ncol = 5
)
colnames(expr_mat) <- c("C1", "C2", "C3", "C4", "C5")
cpm_mat <- t(t(expr_mat) / colSums(expr_mat)) * 1e6
log_cpm_mat <- log2(cpm_mat + 0.1)

# add human gene names to the expression matrix
gene_names_hs <- unique(markers_hs_tbl$gene)
gene_names_hs <- sort(sample(gene_names_hs, n_genes))
rownames(expr_mat) <- gene_names_hs
rownames(cpm_mat) <- gene_names_hs
rownames(log_cpm_mat) <- gene_names_hs

test_that("invalid expression matrix input errors", {
  expect_error(clustermole_enrichment(
    transform(as.data.frame(log_cpm_mat), C1 = "invalid"),
    species = "hs"
  ))
  expect_error(clustermole_enrichment(as.list(log_cpm_mat), species = "hs"))
  expect_error(clustermole_enrichment(log_cpm_mat[1:100, ], species = "hs"))
  expect_error(clustermole_enrichment(log_cpm_mat[, 1:3], species = "hs"))
  expect_error(clustermole_enrichment(cpm_mat, species = "hs"))
  expect_error(clustermole_enrichment(
    log_cpm_mat,
    species = "hs",
    method = "x"
  ))
})

test_that("duplicate expression row names error for every method", {
  expr_mat <- matrix(seq_len(25000) %% 10, nrow = 5000, ncol = 5)
  rownames(expr_mat) <- c("CD3D", "CD3D", paste0("gene", seq_len(4998)))

  for (method in c("gsva", "ssgsea", "singscore", "all")) {
    expect_error(
      clustermole_enrichment(expr_mat, species = "hs", method = method),
      "^expression matrix row names must be unique$"
    )
  }
})

test_that("non-numeric rank cutoffs error", {
  for (max_rank in list(NULL, "100", TRUE, 1i)) {
    expect_error(
      clustermole_enrichment(log_cpm_mat, species = "hs", max_rank = max_rank),
      "`max_rank` is not numeric",
      fixed = TRUE
    )
  }
})

test_that("invalid numeric rank cutoffs error", {
  for (max_rank in list(NA_real_, 0)) {
    expect_error(
      clustermole_enrichment(log_cpm_mat, species = "hs", max_rank = max_rank),
      "`max_rank` must be a number greater than or equal to 1",
      fixed = TRUE
    )
  }
})

# default (gsva)
test_that("matrix and data frame inputs give the same enrichment", {
  skip_if_not_installed("GSVA")
  # skip this slow test on CRAN to keep the check under 10 minutes
  skip_on_cran()
  enrich_hs_tbl <- clustermole_enrichment(
    expr_mat = log_cpm_mat,
    species = "hs"
  )
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
  expect_equal(
    clustermole_enrichment(as.data.frame(log_cpm_mat), species = "hs"),
    enrich_hs_tbl
  )
})

# gsva
test_that("rank cutoffs filter each cluster without changing scores", {
  skip_if_not_installed("GSVA")
  # skip this slow test on CRAN to keep the check under 10 minutes
  skip_on_cran()
  enrich_hs_tbl <- clustermole_enrichment(
    expr_mat = log_cpm_mat,
    species = "hs",
    method = "gsva"
  )
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)

  all_scores <- clustermole_enrichment(log_cpm_mat, "hs", max_rank = Inf)
  expect_gt(max(all_scores$score_rank), 100)
  expect_equal(enrich_hs_tbl, all_scores[all_scores$score_rank <= 100, ])

  top_scores <- clustermole_enrichment(log_cpm_mat, "hs", max_rank = 2)
  expect_equal(top_scores, all_scores[all_scores$score_rank <= 2, ])
  expect_equal(as.integer(table(top_scores$cluster)), rep(2L, 5))
})

# ssgsea
test_that("ssgsea method returns human enrichment results", {
  skip_if_not_installed("GSVA")
  enrich_hs_tbl <- clustermole_enrichment(
    expr_mat = log_cpm_mat,
    species = "hs",
    method = "ssgsea"
  )
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
})

# singscore
test_that("singscore method returns human enrichment results", {
  skip_if_not_installed("GSEABase")
  skip_if_not_installed("singscore")
  enrich_hs_tbl <- clustermole_enrichment(
    expr_mat = log_cpm_mat,
    species = "hs",
    method = "singscore"
  )
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
})

# combined
test_that("combined methods return human enrichment results", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("GSEABase")
  skip_if_not_installed("singscore")
  # skip this slow test on CRAN to keep the check under 10 minutes
  skip_on_cran()
  enrich_hs_tbl <- clustermole_enrichment(
    expr_mat = log_cpm_mat,
    species = "hs",
    method = "all"
  )
  expect_s3_class(enrich_hs_tbl, "tbl_df")
  expect_gt(nrow(enrich_hs_tbl), 100)
  expect_equal(length(unique(enrich_hs_tbl$cluster)), 5)
  expect_equal(enrich_hs_tbl$score_rank, enrich_hs_tbl$score_ranks_mean)

  all_scores <- clustermole_enrichment(
    log_cpm_mat,
    "hs",
    method = "all",
    max_rank = Inf
  )
  expect_gt(max(all_scores$score_rank), 100)
  expect_equal(enrich_hs_tbl, all_scores[all_scores$score_ranks_mean <= 100, ])
})

# add mouse gene names to the expression matrix
gene_names_mm <- unique(markers_mm_tbl$gene)
gene_names_mm <- sort(sample(gene_names_mm, n_genes))
rownames(expr_mat) <- gene_names_mm
rownames(cpm_mat) <- gene_names_mm
rownames(log_cpm_mat) <- gene_names_mm

test_that("mouse expression matrix returns enrichment results", {
  skip_if_not_installed("GSVA")
  # skip this slow test on CRAN to keep the check under 10 minutes
  skip_on_cran()
  enrich_mm_tbl <- clustermole_enrichment(
    expr_mat = log_cpm_mat,
    species = "mm"
  )
  expect_s3_class(enrich_mm_tbl, "tbl_df")
  expect_gt(nrow(enrich_mm_tbl), 100)
  expect_equal(length(unique(enrich_mm_tbl$cluster)), 5)
})
