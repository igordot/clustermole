test_that("clustermole_markers() default params", {
  markers_tbl <- clustermole_markers()
  expect_s3_class(markers_tbl, "tbl_df")
  expect_equal(ncol(markers_tbl), 8)
  expect_gt(nrow(markers_tbl), 100000)
})

test_that("clustermole_markers() wrong input", {
  expect_error(clustermole_markers(species = ""))
  expect_error(clustermole_markers(species = "x"))
  expect_error(clustermole_markers(species = "?"))
  expect_error(clustermole_markers(species = "*"))
  expect_error(clustermole_markers(species = NA))
  expect_error(clustermole_markers(species = NA_character_))
  expect_error(clustermole_markers(species = 1))
  expect_error(clustermole_markers(species = TRUE))
})

test_that("clustermole_markers() species edge cases fall back to hs", {
  expect_equal(
    clustermole_markers(species = c("hs", "mm")),
    clustermole_markers(species = "hs")
  )
  expect_equal(
    clustermole_markers(species = NULL),
    clustermole_markers(species = "hs")
  )
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

test_that("clustermole_markers() partial species match", {
  expect_equal(
    clustermole_markers(species = "h"),
    clustermole_markers(species = "hs")
  )
  expect_equal(
    clustermole_markers(species = "m"),
    clustermole_markers(species = "mm")
  )
})
