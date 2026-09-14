test_that("read_gmt() output", {
  skip_on_cran()
  skip_if_offline(host = "software.broadinstitute.org")
  gmt_url <- paste0(
    "http://software.broadinstitute.org/gsea/msigdb/",
    "supplemental/scsig.all.v1.0.symbols.gmt"
  )
  gmt_tbl <- read_gmt(file = gmt_url)
  expect_s3_class(gmt_tbl, "tbl_df")
  expect_gt(nrow(gmt_tbl), 1000)
})
