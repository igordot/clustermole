test_that("default call returns expected table shape", {
  markers_tbl <- clustermole_markers()
  expect_s3_class(markers_tbl, "tbl_df")
  expect_equal(ncol(markers_tbl), 8)
  expect_equal(dplyr::n_distinct(markers_tbl$db), 11)
  expect_gt(nrow(markers_tbl), 500000)
  expect_gt(dplyr::n_distinct(markers_tbl$celltype_full), 13000)
  expect_gt(dplyr::n_distinct(markers_tbl$gene), 20000)
  expect_gt(dplyr::n_distinct(markers_tbl$organ), 300)
  expect_equal(sum(is.na(markers_tbl$species)), 0)
})

test_that("invalid species input errors", {
  expect_error(clustermole_markers(species = ""))
  expect_error(clustermole_markers(species = "x"))
  expect_error(clustermole_markers(species = "?"))
  expect_error(clustermole_markers(species = "*"))
  expect_error(clustermole_markers(species = NA))
  expect_error(clustermole_markers(species = NA_character_))
  expect_error(clustermole_markers(species = 1))
  expect_error(clustermole_markers(species = TRUE))
})

test_that("multiple or missing species default to human", {
  expect_equal(
    clustermole_markers(species = c("hs", "mm")),
    clustermole_markers(species = "hs")
  )
  expect_equal(
    clustermole_markers(species = NULL),
    clustermole_markers(species = "hs")
  )
})

test_that("requesting human species returns human markers", {
  markers_hs_tbl <- clustermole_markers(species = "hs")
  expect_s3_class(markers_hs_tbl, "tbl_df")
  expect_gt(nrow(markers_hs_tbl), 500000)
  expect_gt(dplyr::n_distinct(markers_hs_tbl$gene), 20000)
})

test_that("requesting mouse species returns mouse markers", {
  markers_mm_tbl <- clustermole_markers(species = "mm")
  expect_s3_class(markers_mm_tbl, "tbl_df")
  expect_gt(nrow(markers_mm_tbl), 500000)
  expect_gt(dplyr::n_distinct(markers_mm_tbl$gene), 17000)
})

test_that("species-specific marker table excludes unmapped genes", {
  expect_gt(sum(is.na(clustermole_markers_tbl$gene_hs)), 0)
  expect_gt(sum(is.na(clustermole_markers_tbl$gene_mm)), 0)
  expect_equal(sum(is.na(clustermole_markers("hs")$gene)), 0)
  expect_equal(sum(is.na(clustermole_markers("mm")$gene)), 0)
})

test_that("known human markers preserve native symbols without orthologs", {
  genes <- c("CEACAM6", "FCGR2C", "SIGLEC7")
  markers <- dplyr::filter(
    clustermole_markers_tbl,
    species == "HS",
    gene_original %in% genes
  )
  expect_setequal(markers$gene_original, genes)
  expect_identical(markers$gene_hs, markers$gene_original)
  expect_identical(markers$gene_mm, rep(NA_character_, nrow(markers)))
})

test_that("known mouse markers preserve native symbols and mapped orthologs", {
  genes <- c("Retnlg", "Chil3", "Klrb1c")
  markers <- dplyr::filter(
    clustermole_markers_tbl,
    species == "MM",
    gene_original %in% genes
  )
  expect_setequal(markers$gene_original, genes)
  expect_identical(markers$gene_mm, markers$gene_original)
  expect_identical(markers$gene_hs, rep(NA_character_, nrow(markers)))

  markers <- dplyr::filter(
    clustermole_markers_tbl,
    species == "MM",
    gene_original == "Dusp7"
  )
  expect_setequal(markers$gene_original, "Dusp7")
  expect_setequal(markers$gene_hs, "DUSP7")
  expect_setequal(markers$gene_mm, "Dusp7")
})

test_that("known human markers map to multi-ortholog and aliased mouse genes", {
  markers <- dplyr::filter(
    clustermole_markers_tbl,
    species == "HS"
  )
  expect_setequal(markers$gene_hs[markers$gene_original == "DYNLT1"], "DYNLT1")
  expect_setequal(
    markers$gene_mm[markers$gene_original == "DYNLT1"],
    c("Dynlt1a", "Dynlt1b", "Dynlt1c", "Dynlt1f")
  )
  expect_setequal(markers$gene_hs[markers$gene_original == "LHFPL6"], "LHFPL6")
  expect_setequal(markers$gene_mm[markers$gene_original == "LHFPL6"], "Lhfpl6")
  expect_setequal(markers$gene_hs[markers$gene_original == "VPREB1"], "VPREB1")
  expect_setequal(markers$gene_mm[markers$gene_original == "VPREB1"], "Vpreb1a")
})

test_that("known mouse markers map to multi-ortholog and aliased human genes", {
  markers <- dplyr::filter(
    clustermole_markers_tbl,
    species == "MM"
  )
  expect_setequal(markers$gene_mm[markers$gene_original == "Try4"], "Try4")
  expect_setequal(
    markers$gene_hs[
      markers$gene_original %in% c("Dynlt1a", "Dynlt1b", "Dynlt1c", "Dynlt1f")
    ],
    "DYNLT1"
  )
  expect_identical(
    markers$gene_mm[
      markers$gene_original %in% c("Dynlt1a", "Dynlt1b", "Dynlt1c", "Dynlt1f")
    ],
    markers$gene_original[
      markers$gene_original %in% c("Dynlt1a", "Dynlt1b", "Dynlt1c", "Dynlt1f")
    ]
  )
  expect_setequal(markers$gene_hs[markers$gene_original == "Lhfpl6"], "LHFPL6")
  expect_setequal(markers$gene_mm[markers$gene_original == "Lhfpl6"], "Lhfpl6")
  expect_setequal(markers$gene_hs[markers$gene_original == "Vpreb1a"], "VPREB1")
  expect_setequal(
    markers$gene_mm[markers$gene_original == "Vpreb1a"],
    "Vpreb1a"
  )
})

test_that("blank-species markers resolve aliased gene symbols across species", {
  markers <- dplyr::filter(
    clustermole_markers_tbl,
    species == ""
  )
  expect_setequal(markers$gene_hs[markers$gene_original == "CD123"], "IL3RA")
  expect_setequal(markers$gene_mm[markers$gene_original == "CD123"], "Il3ra")
  expect_setequal(markers$gene_hs[markers$gene_original == "CD127"], "IL7R")
  expect_setequal(markers$gene_mm[markers$gene_original == "CD127"], "Il7r")
  expect_setequal(markers$gene_hs[markers$gene_original == "DYNLT1"], "DYNLT1")
  expect_setequal(
    markers$gene_mm[markers$gene_original == "DYNLT1"],
    c("Dynlt1a", "Dynlt1b", "Dynlt1c", "Dynlt1f")
  )
  expect_setequal(markers$gene_hs[markers$gene_original == "LHFP"], "LHFPL6")
  expect_setequal(markers$gene_mm[markers$gene_original == "LHFP"], "Lhfpl6")
  expect_setequal(markers$gene_hs[markers$gene_original == "VPREB1"], "VPREB1")
  expect_setequal(markers$gene_mm[markers$gene_original == "VPREB1"], "Vpreb1a")
})

test_that("blank-species CD16 maps to mismatched human and mouse symbols", {
  markers <- dplyr::filter(
    clustermole_markers_tbl,
    species == "",
    gene_original == "CD16"
  )
  expect_setequal(markers$gene_hs, "FCGR3B")
  expect_setequal(markers$gene_mm, "Fcgr3")
})

test_that("partial species strings match by prefix", {
  expect_equal(
    clustermole_markers(species = "h"),
    clustermole_markers(species = "hs")
  )
  expect_equal(
    clustermole_markers(species = "m"),
    clustermole_markers(species = "mm")
  )
})
