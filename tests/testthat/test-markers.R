
test_that("clustermole_markers() default params", {
  markers_tbl <- clustermole_markers()
  expect_s3_class(markers_tbl, "tbl_df")
  expect_gt(nrow(markers_tbl), 100000)
})

test_that("clustermole_markers() wrong input", {
  expect_error(clustermole_markers(species = ""))
  expect_error(clustermole_markers(species = "x"))
  expect_error(clustermole_markers(species = "?"))
  expect_error(clustermole_markers(species = "*"))
})

test_that("clustermole_markers() human input", {
  markers_hs_tbl <- clustermole_markers(species = "hs")
  expect_s3_class(markers_hs_tbl, "tbl_df")
  expect_gt(nrow(markers_hs_tbl), 100000)
})

test_that("clustermole_markers() mouse input", {
  markers_mm_tbl <- clustermole_markers(species = "mm")
  expect_s3_class(markers_mm_tbl, "tbl_df")
  expect_gt(nrow(markers_mm_tbl), 100000)
})
