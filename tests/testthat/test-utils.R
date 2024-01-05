test_that("read_gmt() output", {
  skip_if_offline(host = "software.broadinstitute.org")
  gmt_tbl <- read_gmt(file = "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt")
  expect_s3_class(gmt_tbl, "tbl_df")
  expect_gt(nrow(gmt_tbl), 1000)
})
