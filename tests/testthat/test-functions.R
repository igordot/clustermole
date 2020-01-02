
# markers table tests
markers_tbl <- clustermole_markers()
expect_s3_class(markers_tbl, "tbl_df")
expect_gt(nrow(markers_tbl), 100000)

# expression table for enrichment tests
n_genes <- 10000
gene_names <- markers_tbl
gene_names <- dplyr::filter(gene_names, species == "Human")
gene_names <- unique(gene_names$gene)
gene_names <- sort(sample(gene_names, n_genes))
expr_mat <- matrix(rnbinom(n_genes * 5, size = 1, mu = 10), nrow = n_genes, ncol = 5)
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- c("C1", "C2", "C3", "C4", "C5")
cpm_mat <- t(t(expr_mat) / (colSums(expr_mat) / 1e6))
log_cpm_mat <- log2(cpm_mat + 0.1)

# enrichment tests
enrich_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "hs")
expect_s3_class(enrich_tbl, "tbl_df")
expect_gt(nrow(enrich_tbl), 100)
expect_equal(length(unique(enrich_tbl$cluster)), 5)
expect_error(clustermole_enrichment(expr_mat = cpm_mat, species = "hs"))
expect_error(clustermole_enrichment(expr_mat = log_cpm_mat[1:100, ], species = "hs"))
