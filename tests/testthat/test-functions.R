
clustermole_tbl <- clustermole_markers()
expect_s3_class(clustermole_tbl, "tbl_df")
expect_gt(nrow(clustermole_tbl), 100000)
