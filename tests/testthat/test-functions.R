
# markers table tests
markers_tbl <- clustermole_markers()
expect_s3_class(markers_tbl, "tbl_df")
expect_gt(nrow(markers_tbl), 100000)

# expression matrix for enrichment tests
n_genes <- 10000
expr_mat <- matrix(rnbinom(n_genes * 5, size = 1, mu = 10), nrow = n_genes, ncol = 5)
colnames(expr_mat) <- c("C1", "C2", "C3", "C4", "C5")
cpm_mat <- t(t(expr_mat) / (colSums(expr_mat) / 1e6))
log_cpm_mat <- log2(cpm_mat + 0.1)

# add human gene names to the expression matrix
gene_names <- dplyr::filter(markers_tbl, species == "Human")
gene_names <- unique(gene_names$gene)
gene_names <- sort(sample(gene_names, n_genes))
rownames(expr_mat) <- gene_names
rownames(cpm_mat) <- gene_names
rownames(log_cpm_mat) <- gene_names

# enrichment tests - human
enrich_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "hs")
expect_s3_class(enrich_tbl, "tbl_df")
expect_gt(nrow(enrich_tbl), 100)
expect_equal(length(unique(enrich_tbl$cluster)), 5)
expect_error(clustermole_enrichment(expr_mat = as.data.frame(log_cpm_mat), species = "hs"))
expect_error(clustermole_enrichment(expr_mat = log_cpm_mat[1:100, ], species = "hs"))
expect_error(clustermole_enrichment(expr_mat = cpm_mat, species = "hs"))

# add mouse gene names to the expression matrix
gene_names <- dplyr::filter(markers_tbl, species == "Mouse")
gene_names <- unique(gene_names$gene)
gene_names <- sort(sample(gene_names, n_genes))
rownames(expr_mat) <- gene_names
rownames(cpm_mat) <- gene_names
rownames(log_cpm_mat) <- gene_names

# enrichment tests - mouse
enrich_tbl <- clustermole_enrichment(expr_mat = log_cpm_mat, species = "mm")
expect_s3_class(enrich_tbl, "tbl_df")
expect_gt(nrow(enrich_tbl), 100)
expect_equal(length(unique(enrich_tbl$cluster)), 5)
