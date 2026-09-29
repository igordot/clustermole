# Download, clean up, and harmonize cell type markers

# Settings and helper functions -----

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(tibble)
library(stringr)
library(janitor)
library(jsonlite)
library(usethis)

options(timeout = 300)
markers <- list()
common_cols <- c(
  "db",
  "species_original",
  "species",
  "organ",
  "celltype",
  "gene"
)

# Extract gene synonyms from the NCBI gene info table
get_alt_symbols <- function(gene_info) {
  symbols <- gene_info$symbol
  gene_info <- distinct(gene_info, symbol, synonyms)
  gene_info <- separate_longer_delim(gene_info, synonyms, delim = "|")
  # Exclude short aliases to reduce nonspecific matches
  gene_info <- filter(gene_info, nchar(synonyms) > 2)
  # Exclude aliases with whitespaces
  gene_info <- filter(gene_info, str_detect(synonyms, "\\s", negate = TRUE))
  gene_info <- distinct(gene_info, symbol, synonyms)
  gene_info <- mutate(gene_info, synonyms_upper = toupper(synonyms))
  # Exclude aliases that match canonical symbols regardless of case
  gene_info <- filter(gene_info, !(synonyms_upper %in% toupper(symbols)))
  # Keep aliases for one canonical gene, allowing case variants
  gene_info <- filter(gene_info, n_distinct(symbol) == 1, .by = synonyms_upper)
  select(gene_info, symbol, alt_symbol = synonyms)
}

# Standardize labels without changing the vector length
standardize_species <- function(x) {
  x <- replace_na(x, "")
  x <- str_squish(x)
  x <- str_to_lower(x)
  x[x %in% c("human", "homo sapiens", "hs")] <- "HS"
  x[x %in% c("mouse", "mus musculus", "mm")] <- "MM"
  # "4" is a data entry error for some rows in PanglaoDB
  x[x %in% c("none", "4")] <- ""
  stopifnot(all(x %in% c("HS", "MM", "")))
  x
}

# Clear species labels that disagree with the genes within each signature
confirm_species <- function(x, genes_hs, genes_mm) {
  # Confirm cell type exists for all entries
  x <- drop_na(x, celltype)
  # Keep only valid genes
  x <- mutate(x, gene = str_trim(gene))
  x <- filter(x, gene %in% c(genes_hs, genes_mm))
  # Standardize species labels
  x <- mutate(x, species_original = replace_na(species, ""))
  x <- mutate(x, species = standardize_species(species_original))
  # Remove cell types with too few markers
  x <- distinct(x)
  x <- mutate(x, soc = paste(species, organ, celltype))
  x <- add_count(x, soc)
  x <- filter(x, n >= 5)
  x <- select(x, !n)
  # Retain the original label only when most markers support it
  x <-
    x |>
    group_by(soc) |>
    mutate(
      species = case_when(
        species == "HS" & mean(gene %in% genes_hs) > 0.6 ~ "HS",
        species == "MM" & mean(gene %in% genes_mm) > 0.6 ~ "MM",
        TRUE ~ ""
      )
    ) |>
    ungroup()
  select(x, !soc)
}

# Human/mouse gene symbols -----

ncbi_hs_source <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
ncbi_mm_source <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"

gene_info_hs <- read_tsv(ncbi_hs_source, show_col_types = FALSE)
gene_info_hs <- clean_names(gene_info_hs)
gene_info_hs <- arrange(gene_info_hs, symbol)
gene_info_hs <- filter(gene_info_hs, type_of_gene != "biological-region")
gene_info_hs <- filter(gene_info_hs, chromosome != "-")
gene_info_hs <- filter(gene_info_hs, map_location != "-")
nrow(gene_info_hs)
stopifnot(nrow(gene_info_hs) > 60000)

gene_info_mm <- read_tsv(ncbi_mm_source, show_col_types = FALSE)
gene_info_mm <- clean_names(gene_info_mm)
gene_info_mm <- arrange(gene_info_mm, symbol)
gene_info_mm <- filter(gene_info_mm, type_of_gene != "biological-region")
gene_info_mm <- filter(gene_info_mm, chromosome != "-")
gene_info_mm <- filter(gene_info_mm, map_location != "-")
gene_info_mm <- filter(gene_info_mm, nomenclature_status != "-")
gene_info_mm <- filter(
  gene_info_mm,
  str_detect(symbol, fixed("("), negate = TRUE)
)
nrow(gene_info_mm)
stopifnot(nrow(gene_info_mm) > 50000)
stopifnot(anyDuplicated(toupper(gene_info_mm$symbol)) == 0)

# Canonical symbols should be distinct across species, with a few exceptions
stopifnot(length(intersect(gene_info_hs$symbol, gene_info_mm$symbol)) < 20)

synonyms_hs <- get_alt_symbols(gene_info_hs)
synonyms_hs <- filter(synonyms_hs, !(alt_symbol %in% gene_info_mm$symbol))
nrow(synonyms_hs)
stopifnot(nrow(synonyms_hs) > 60000)

synonyms_mm <- get_alt_symbols(gene_info_mm)
synonyms_mm <- filter(synonyms_mm, !(alt_symbol %in% gene_info_hs$symbol))
nrow(synonyms_mm)
stopifnot(nrow(synonyms_mm) > 60000)

# Human/mouse orthologs -----

agr_source <- "https://www.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_TSV_COMBINED.tsv.gz"
agr <- read_tsv(agr_source, comment = "#", show_col_types = FALSE)
stopifnot(nrow(agr) > 950000)

# Keep human-mouse ortholog pairs
orthologs <-
  filter(
    agr,
    Gene1SpeciesName %in% c("Homo sapiens", "Mus musculus"),
    Gene2SpeciesName %in% c("Homo sapiens", "Mus musculus"),
    Gene1SpeciesName != Gene2SpeciesName
  )
stopifnot(nrow(orthologs) > 40000)

# Extract human and mouse symbols before filtering ortholog pairs
genes_hs <- c(
  orthologs$Gene1Symbol[orthologs$Gene1SpeciesName == "Homo sapiens"],
  gene_info_hs$symbol,
  synonyms_hs$alt_symbol
)
genes_hs <- sort(unique(genes_hs))
length(genes_hs)
genes_mm <- c(
  orthologs$Gene1Symbol[orthologs$Gene1SpeciesName == "Mus musculus"],
  gene_info_mm$symbol,
  synonyms_mm$alt_symbol
)
genes_mm <- sort(unique(genes_mm))
length(genes_mm)
valid_genes <- sort(unique(c(genes_hs, genes_mm)))
length(valid_genes)
stopifnot(length(valid_genes) > 200000)
stopifnot(length(intersect(genes_hs, genes_mm)) < 8000)

# Keep only the reciprocal best orthologs (ties are possible)
orthologs <- filter(orthologs, IsBestScore == "Yes", IsBestRevScore == "Yes")
stopifnot(nrow(orthologs) > 40000)

orthologs <-
  mutate(
    orthologs,
    gene_hs = if_else(
      Gene1SpeciesName == "Homo sapiens",
      Gene1Symbol,
      Gene2Symbol
    ),
    gene_mm = if_else(
      Gene1SpeciesName == "Mus musculus",
      Gene1Symbol,
      Gene2Symbol
    )
  )
orthologs <- arrange(orthologs, gene_hs, gene_mm)
orthologs <- distinct(orthologs, gene_hs, gene_mm)
nrow(orthologs)
stopifnot(nrow(orthologs) > 20000)

# PanglaoDB -----

# Ref: Franzen et al. Database (2019) https://doi.org/10.1093/database/baz046
# Source: https://panglaodb.se/

panglao_source <- "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
panglao_all <- read_tsv(panglao_source, show_col_types = FALSE)
panglao_all <- clean_names(panglao_all)
nrow(panglao_all)
# 8286 rows, 178 cell types

# Clean up
panglao_clean <-
  panglao_all |>
  mutate(
    db = "PanglaoDB",
    organ = replace_na(organ, ""),
    celltype = cell_type,
    gene = str_trim(official_gene_symbol)
  ) |>
  separate_longer_delim(species, delim = " ") |>
  mutate(
    species_original = species,
    species = standardize_species(species)
  )
panglao_clean <- panglao_clean[, common_cols]

# Split to fix mouse genes (all provided as uppercase)
panglao_hs <- filter(panglao_clean, species == "HS")
stopifnot(nrow(panglao_hs) > nrow(panglao_all) * 0.9)
panglao_mm <- filter(panglao_clean, species == "MM")
stopifnot(nrow(panglao_mm) > nrow(panglao_all) * 0.9)

# Restore mouse symbols when the uppercase version matches a known gene
# Problematic genes: DYNLT1, LHFP, VPREB1
panglao_mm <- mutate(panglao_mm, gene_upper = toupper(gene), .keep = "unused")
panglao_mm <-
  inner_join(
    panglao_mm,
    tibble(
      gene = genes_mm,
      gene_upper = toupper(genes_mm)
    ),
    by = "gene_upper"
  )
panglao_mm <- select(panglao_mm, -gene_upper)
stopifnot(nrow(panglao_mm) > nrow(panglao_all) * 0.95)

# Check that human and mouse genes are valid
stopifnot(mean(unique(panglao_hs$gene) %in% genes_hs) > 0.99)
stopifnot(mean(unique(panglao_mm$gene) %in% genes_mm) > 0.99)

markers$panglao <- bind_rows(panglao_hs, panglao_mm)
nrow(markers$panglao)
n_distinct(markers$panglao$celltype)
# 15749 rows, 178 cell types
stopifnot(nrow(markers$panglao) > 15000)
stopifnot(n_distinct(markers$panglao$celltype) > 150)
stopifnot(n_distinct(markers$panglao$gene) > 8000)

# CellMarker -----

# Ref: Hu et al. Nucleic Acids Research (2023) https://doi.org/10.1093/nar/gkac947
# Source: http://bio-bigdata.hrbmu.edu.cn/CellMarker/

# cellmarker_source <- "http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt"
# cellmarker_source <- "https://bio-bigdata.hrbmu.edu.cn/CellMarker1.0/download/all_cell_markers.txt"
# cellmarker_source <- "https://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/CellMarker_download_files/file/Cell_marker_All.xlsx"
cellmarker_source <- "https://bio-bigdata.hrbmu.edu.cn/CellMarker/file/all_cell_marker.zip"
# cellmarker_source <- "https://zenodo.org/records/22808257/files/all_cell_marker.zip?download=1"
cellmarker_zip <- tempfile(fileext = ".zip")
download.file(
  url = cellmarker_source,
  destfile = cellmarker_zip,
  quiet = TRUE,
  mode = "wb"
)
cellmarker_dir <- tempfile()
unzip(cellmarker_zip, exdir = cellmarker_dir)
cellmarker_all <- read_tsv(
  list.files(cellmarker_dir, full.names = TRUE)[1],
  show_col_types = FALSE
)
cellmarker_all <- clean_names(cellmarker_all)
unlink(c(cellmarker_zip, cellmarker_dir), recursive = TRUE)

nrow(cellmarker_all)
n_distinct(cellmarker_all$cell_name)
# 2537570 rows, 4455 cell types

# Most CellMarker 3.0 rows are algorithm output (marker_source "Method")
cellmarker_clean <- filter(cellmarker_all, marker_source != "Method")
nrow(cellmarker_clean)
n_distinct(cellmarker_clean$cell_name)
# 218643 rows, 4066 cell types

cellmarker_clean <- mutate(cellmarker_clean, symbol = str_trim(symbol))
cellmarker_clean <- filter(cellmarker_clean, symbol %in% valid_genes)
cellmarker_clean <- distinct(cellmarker_clean)
cellmarker_clean <-
  mutate(
    cellmarker_clean,
    db = "CellMarker",
    pmid = replace_na(pmid, 0),
    disease = replace_na(disease, ""),
    disease = str_remove(disease, "Undefined"),
    tissue_type = replace_na(tissue_type, ""),
    organ = str_remove(tissue_type, "Undefined"),
    celltype = str_c(cell_name, " (", disease, ")"),
    celltype = str_remove(celltype, fixed(" (Normal)")),
    celltype = str_remove(celltype, fixed(" ()")),
    celltype = str_c(celltype, " [PMID:", pmid, "]"),
    celltype = str_remove(celltype, fixed(" [PMID:0]")),
    gene = symbol
  )
count(cellmarker_clean, disease, tissue_type, organ, sort = TRUE)
count(cellmarker_clean, disease, organ, celltype, sort = TRUE)
nrow(cellmarker_clean)

# Fix species labels that disagree with the genes
# Problematic PMIDs: 32066997, 34301296, 32286228, 41481707
cellmarker_clean <-
  confirm_species(
    cellmarker_clean,
    genes_hs,
    genes_mm
  )
nrow(cellmarker_clean)
distinct(cellmarker_clean, species_original, species, organ, celltype) |>
  count(species_original, species)

markers$cellmarker <- cellmarker_clean[, common_cols]
nrow(markers$cellmarker)
n_distinct(cellmarker_clean$celltype)
# 126970 rows, 7505 cell types
stopifnot(nrow(markers$cellmarker) > 50000)
stopifnot(n_distinct(markers$cellmarker$celltype) > 5000)
stopifnot(n_distinct(markers$cellmarker$gene) > 25000)

# SaVanT -----

# Ref: Lopez et al. BMC Genomics (2017) https://doi.org/10.1186/s12864-017-4167-7
# Source: http://newpathways.mcdb.ucla.edu/savant-dev/

# The original host is offline
# savant_source <- "http://newpathways.mcdb.ucla.edu/savant-dev/SaVanT_Signatures_Release01.zip"
savant_source <- "https://web.archive.org/web/20210518161907id_/https://newpathways.mcdb.ucla.edu/savant-dev/SaVanT_Signatures_Release01.zip"
savant_txt <- "SaVanT_Signatures_Release01.tab.txt"
savant_tmp <- tempfile(fileext = ".zip")
download.file(
  url = savant_source,
  destfile = savant_tmp,
  quiet = TRUE,
  mode = "wb"
)
savant_list <- strsplit(readLines(unz(savant_tmp, savant_txt)), "\t")
unlink(savant_tmp)
savant_all <- lapply(savant_list, tail, -1)
names(savant_all) <- sapply(savant_list, head, 1)
savant_all <- enframe(savant_all, name = "celltype", value = "gene")
savant_all <- unnest(savant_all, gene)
nrow(savant_all)
# 596248 rows, 619 cell types

savant_clean <-
  savant_all |>
  mutate(
    db = "SaVanT",
    species = case_when(
      str_detect(celltype, "^MBA_") ~ "MM",
      str_detect(celltype, "^IMGN_") ~ "MM",
      str_detect(celltype, "^HBA_") ~ "HS",
      str_detect(celltype, "^HPCA_") ~ "HS",
      str_detect(celltype, "^MA_") ~ "HS",
      TRUE ~ ""
    ),
    organ = ""
  )

savant_clean <- confirm_species(savant_clean, genes_hs, genes_mm)
nrow(savant_clean)
distinct(savant_clean, species_original, species, organ, celltype) |>
  count(species_original, species)

# Select the top 50 genes (default in SaVanT)
savant_clean <-
  savant_clean |>
  group_by(celltype) |>
  slice_head(n = 50) |>
  ungroup()

markers$savant <- savant_clean[, common_cols]
nrow(markers$savant)
n_distinct(markers$savant$celltype)
# 30944 rows, 619 cell types
stopifnot(nrow(markers$savant) > 30000)
stopifnot(n_distinct(markers$savant$celltype) > 600)
stopifnot(n_distinct(markers$savant$gene) > 5000)

# MSigDB C8/M8 (formerly SCSig) -----

# Ref (human): Liberzon et al. Bioinformatics (2011) https://doi.org/10.1093/bioinformatics/btr260
# Ref (mouse): Castanza et al. Nature Methods (2023) https://doi.org/10.1038/s41592-023-02014-7
# Source: http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C8

# msigbd_source <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml"
msigdb_hs_source <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c8.all.v2026.1.Hs.symbols.gmt"
msigdb_mm_source <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Mm/m8.all.v2026.1.Mm.symbols.gmt"

msigdb_hs_all <- clustermole::read_gmt(msigdb_hs_source)
msigdb_hs_all <- mutate(msigdb_hs_all, gene = str_trim(gene))
msigdb_hs_all <- mutate(msigdb_hs_all, species = "HS")
nrow(msigdb_hs_all)
# 157462 rows, 866 cell types
stopifnot(nrow(msigdb_hs_all) > 50000)
stopifnot(n_distinct(msigdb_hs_all$celltype) > 400)
stopifnot(n_distinct(msigdb_hs_all$gene) > 10000)
stopifnot(mean(unique(msigdb_hs_all$gene) %in% genes_hs) > 0.95)

msigdb_mm_all <- clustermole::read_gmt(msigdb_mm_source)
msigdb_mm_all <- mutate(msigdb_mm_all, gene = str_trim(gene))
msigdb_mm_all <- mutate(msigdb_mm_all, species = "MM")
nrow(msigdb_mm_all)
n_distinct(msigdb_mm_all$celltype)
# 47976 rows, 233 cell types
stopifnot(nrow(msigdb_mm_all) > 15000)
stopifnot(n_distinct(msigdb_mm_all$celltype) > 100)
stopifnot(n_distinct(msigdb_mm_all$gene) > 10000)
stopifnot(mean(unique(msigdb_mm_all$gene) %in% genes_mm) > 0.95)

markers$msigdb <- bind_rows(msigdb_hs_all, msigdb_mm_all)
markers$msigdb$db <- "MSigDB"
markers$msigdb$organ <- ""
markers$msigdb$species_original <- markers$msigdb$species
markers$msigdb <- markers$msigdb[, common_cols]
nrow(markers$msigdb)
n_distinct(markers$msigdb$celltype)
# 205438 rows, 1099 cell types

# xCell -----

# Ref: Aran et al. Genome Biology (2017) https://doi.org/10.1186/s13059-017-1349-1
# Source: http://xcell.ucsf.edu/

# NCBI PMC gates file downloads behind a JS puzzle
# xcell_source <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5688663/bin/13059_2017_1349_MOESM3_ESM.xlsx"
xcell_source <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-017-1349-1/MediaObjects/13059_2017_1349_MOESM3_ESM.xlsx"
xcell_tmp <- tempfile(fileext = ".xlsx")
download.file(
  url = xcell_source,
  destfile = xcell_tmp,
  quiet = TRUE,
  mode = "wb"
)
xcell_all <- read_xlsx(xcell_tmp)
xcell_all <- clean_names(xcell_all)
unlink(xcell_tmp)
nrow(xcell_all)
# 489 rows (one row per cell type)

markers$xcell <-
  xcell_all |>
  select(!number_of_genes) |>
  gather(key = "k", value = "gene", !celltype_source_id) |>
  mutate(
    db = "xCell",
    species_original = "HS",
    species = species_original,
    organ = "",
    celltype = celltype_source_id
  ) |>
  drop_na(celltype, gene)
markers$xcell <- markers$xcell[, common_cols]
nrow(markers$xcell)
n_distinct(markers$xcell$celltype)
# 20803 rows, 489 cell types
stopifnot(nrow(markers$xcell) > 20000)
stopifnot(n_distinct(markers$xcell$celltype) > 400)
stopifnot(n_distinct(markers$xcell$gene) > 5000)

# TISSUES -----

# Ref: Palasca et al. Database (2018) https://doi.org/10.1093/database/bay003
# Source: https://tissues.jensenlab.org/

tissues_hs_source <- "https://download.jensenlab.org/human_tissue_knowledge_full.tsv"
tissues_mm_source <- "https://download.jensenlab.org/mouse_tissue_knowledge_full.tsv"

tissues_hs_all <-
  read_tsv(
    tissues_hs_source,
    col_names = FALSE,
    show_col_types = FALSE
  )
tissues_hs_all <- mutate(tissues_hs_all, species_original = "HS")
tissues_hs_all <- mutate(tissues_hs_all, X2 = str_trim(X2))
tissues_hs_all <- filter(tissues_hs_all, X2 %in% genes_hs)
stopifnot(nrow(tissues_hs_all) > 300000)
nrow(tissues_hs_all)
# 338585 rows, 593 cell types

tissues_mm_all <-
  read_tsv(
    tissues_mm_source,
    col_names = FALSE,
    show_col_types = FALSE
  )
tissues_mm_all <- mutate(tissues_mm_all, species_original = "MM")
tissues_mm_all <- mutate(tissues_mm_all, X2 = str_trim(X2))
tissues_mm_all <- filter(tissues_mm_all, X2 %in% genes_mm)
stopifnot(nrow(tissues_mm_all) > 150000)
nrow(tissues_mm_all)
# 176262 rows, 316 cell types

markers$tissues <-
  bind_rows(tissues_hs_all, tissues_mm_all) |>
  filter(X3 != "BTO:0000000") |>
  filter(str_detect(X4, "BTO:", negate = TRUE)) |>
  mutate(
    db = "TISSUES",
    species = species_original,
    organ = "",
    celltype = X4,
    gene = X2
  )
markers$tissues <- markers$tissues[, common_cols]
nrow(markers$tissues)
n_distinct(markers$tissues$celltype)
# 452589 rows, 590 cell types
stopifnot(nrow(markers$tissues) > 100000)
stopifnot(n_distinct(markers$tissues$celltype) > 500)
stopifnot(n_distinct(markers$tissues$gene) > 25000)

# ScType -----

# Ref: Ianevski et al. Nature Communications (2022) https://doi.org/10.1038/s41467-022-28803-w
# Source: https://github.com/IanevskiAleksandr/sc-type

sctype_source <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
sctype_tmp <- tempfile(fileext = ".xlsx")
download.file(
  url = sctype_source,
  destfile = sctype_tmp,
  quiet = TRUE,
  mode = "wb"
)
sctype_all <- read_xlsx(sctype_tmp)
unlink(sctype_tmp)
nrow(sctype_all)
# 268 rows, 193 cell types

sctype_clean <-
  mutate(
    sctype_all,
    db = "ScType",
    species_original = "",
    species = species_original,
    organ = tissueType,
    celltype = cellName,
    gene = geneSymbolmore1
  )
sctype_clean <- separate_longer_delim(sctype_clean, gene, delim = ",")
sctype_clean <- drop_na(sctype_clean, celltype, gene)

markers$sctype <- sctype_clean[, common_cols]
nrow(markers$sctype)
n_distinct(markers$sctype$celltype)
# 4697 rows, 193 cell types
stopifnot(nrow(markers$sctype) > 4000)
stopifnot(n_distinct(markers$sctype$celltype) > 150)
stopifnot(n_distinct(markers$sctype$gene) > 2000)

# HPA -----

# Ref: Karlsson et al. Science Advances (2021) https://doi.org/10.1126/sciadv.abh2169
# Source: https://www.proteinatlas.org/about/download

# "Cell type enriched" >=4x any other single cell type
# "Cell type enhanced" (vs. the mean only)
# "Group enriched" ambiguous within its 2-10 cell types
hpa_source <- "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
hpa_zip <- tempfile(fileext = ".zip")
download.file(
  url = hpa_source,
  destfile = hpa_zip,
  quiet = TRUE,
  mode = "wb"
)
hpa_dir <- tempfile()
unzip(hpa_zip, exdir = hpa_dir)
hpa_all <- read_tsv(
  list.files(hpa_dir, full.names = TRUE)[1],
  show_col_types = FALSE
)
hpa_all <- clean_names(hpa_all)
unlink(c(hpa_zip, hpa_dir), recursive = TRUE)
nrow(hpa_all)
# 20162 rows

hpa_clean <-
  filter(
    hpa_all,
    rna_single_cell_type_specificity == "Cell type enriched"
  )
nrow(hpa_clean)

hpa_clean <-
  hpa_clean |>
  mutate(gene = str_trim(gene)) |>
  filter(gene %in% genes_hs) |>
  mutate(
    db = "HPA",
    species_original = "HS",
    species = species_original,
    organ = "",
    celltype = sub(":.*", "", rna_single_cell_type_specific_n_cpm),
    gene = gene
  )

markers$hpa <- hpa_clean[, common_cols]
nrow(markers$hpa)
n_distinct(markers$hpa$celltype)
# 1882 rows, 127 cell types
stopifnot(nrow(markers$hpa) > 1000)
stopifnot(n_distinct(markers$hpa$celltype) > 100)
stopifnot(n_distinct(markers$hpa$gene) > 1000)

# Cell Taxonomy -----

# Ref: Jiang et al. Nucleic Acids Research (2023) https://doi.org/10.1093/nar/gkac816
# Source: https://ngdc.cncb.ac.cn/celltaxonomy/

celltaxonomy_source <- "https://download.cncb.ac.cn/celltaxonomy/Cell_Taxonomy_resource.txt"
celltaxonomy_all <-
  read_tsv(
    celltaxonomy_source,
    guess_max = 50000,
    show_col_types = FALSE
  )
celltaxonomy_all <- clean_names(celltaxonomy_all)
nrow(celltaxonomy_all)
n_distinct(celltaxonomy_all$cell_standard)
# 226222 rows, 3142 cell types

# Keep only cell types curated by Cell Taxonomy team itself
celltaxonomy_clean <-
  filter(
    celltaxonomy_all,
    source == "Cell Taxonomy",
    species %in% c("Homo sapiens", "Mus musculus")
  )
celltaxonomy_clean <-
  mutate(
    celltaxonomy_clean,
    db = "CellTaxonomy",
    organ = replace_na(tissue_standard, ""),
    condition = replace_na(condition, ""),
    condition = str_remove(condition, "Physiology"),
    celltype = str_c(cell_standard, " (", condition, ")"),
    celltype = str_remove(celltype, fixed(" ()")),
    pmid = replace_na(pmid, 0),
    celltype = str_c(celltype, " [PMID:", pmid, "]"),
    celltype = str_remove(celltype, fixed(" [PMID:0]")),
    gene = cell_marker
  )
nrow(celltaxonomy_clean)
n_distinct(celltaxonomy_clean$celltype)
# 14997 rows, 4911 cell types

# Fix species labels that disagree with the genes
celltaxonomy_clean <-
  confirm_species(
    celltaxonomy_clean,
    genes_hs,
    genes_mm
  )
nrow(celltaxonomy_clean)
distinct(celltaxonomy_clean, species_original, species, organ, celltype) |>
  count(species_original, species)

markers$celltaxonomy <- celltaxonomy_clean[, common_cols]
nrow(markers$celltaxonomy)
n_distinct(markers$celltaxonomy$celltype)
# 4620 rows, 573 cell types
stopifnot(nrow(markers$celltaxonomy) > 4000)
stopifnot(n_distinct(markers$celltaxonomy$celltype) > 500)
stopifnot(n_distinct(markers$celltaxonomy$gene) > 1000)

# CellMatch -----

# Ref: Shao et al. iScience (2020) https://doi.org/10.1016/j.isci.2020.100882
# scCATCH includes CellMatch as its reference database for cell type annotation
# Source: https://github.com/ZJUFanLab/scCATCH

cellmatch_source <- "https://raw.githubusercontent.com/ZJUFanLab/scCATCH/master/data/cellmatch.rda"
cellmatch_tmp <- tempfile(fileext = ".rda")
download.file(
  url = cellmatch_source,
  destfile = cellmatch_tmp,
  quiet = TRUE,
  mode = "wb"
)
cellmatch_env <- new.env()
load(cellmatch_tmp, envir = cellmatch_env)
cellmatch_all <- cellmatch_env$cellmatch
unlink(cellmatch_tmp)
nrow(cellmatch_all)
n_distinct(cellmatch_all$celltype)
# 49560 rows, 353 cell types

# A cell type can cover several biological subtypes or disease contexts.
cellmatch_clean <-
  cellmatch_all |>
  drop_na(celltype) |>
  unite(
    celltype,
    subtype1,
    subtype2,
    subtype3,
    celltype,
    sep = " ",
    na.rm = TRUE
  ) |>
  mutate(
    db = "CellMatch",
    species_original = species,
    organ = replace_na(tissue, ""),
    cancer = replace_na(cancer, ""),
    celltype = str_c(celltype, " (", cancer, ")"),
    celltype = str_remove_all(celltype, " \\(Normal\\)"),
    celltype = str_remove(celltype, fixed(" ()"))
  )

# Fix species labels that disagree with the genes
cellmatch_clean <-
  confirm_species(
    cellmatch_clean,
    genes_hs,
    genes_mm
  )
nrow(cellmatch_clean)
distinct(cellmatch_clean, species_original, species, organ, celltype) |>
  count(species_original, species)

markers$cellmatch <- cellmatch_clean[, common_cols]
nrow(markers$cellmatch)
n_distinct(markers$cellmatch$celltype)
# 47171 rows, 451 cell types
stopifnot(nrow(markers$cellmatch) > 40000)
stopifnot(n_distinct(markers$cellmatch$celltype) > 400)
stopifnot(n_distinct(markers$cellmatch$gene) > 20000)

# DISCO -----

# Ref: Li et al. Nucleic Acids Research (2022) https://doi.org/10.1093/nar/gkab1020
# Source: https://www.immunesinglecell.org/

disco_source <- "https://immunesinglecell.com/disco_v3_api/toolkit/getCellTypeDEGMarkers"
disco_all <- fromJSON(disco_source)
nrow(disco_all)
# 390 rows (one row per cell type)

disco_clean <-
  disco_all |>
  mutate(
    db = "DISCO",
    species_original = "HS",
    species = species_original,
    organ = "",
    celltype = cell_type,
    gene = genes
  ) |>
  drop_na(celltype, genes) |>
  separate_longer_delim(gene, delim = ";") |>
  mutate(gene = str_trim(gene)) |>
  filter(gene %in% valid_genes)

markers$disco <- disco_clean[, common_cols]
nrow(markers$disco)
n_distinct(markers$disco$celltype)
# 48660 rows, 390 cell types
stopifnot(nrow(markers$disco) > 40000)
stopifnot(n_distinct(markers$disco$celltype) > 300)
stopifnot(n_distinct(markers$disco$gene) > 5000)

# Excluded sources -----

# ARCHS4: most signatures are over 2000 genes
# CellTypist Pan-Immune atlas: less than 5 curated markers for most signatures
# Human Cell Landscape (HCL): no bulk marker download
# Human Cell Atlas (HCA): no bulk marker download
# Invitrogen
# scTyper.db
# SHOGoiN
# tinyatlas

# Merge markers data frames -----

# Save the list of data frames
markers_list <- markers

# Check all the sources
names(markers_list)
stopifnot(length(markers_list) == 11)

# Combine a list of data frames into a single data frame
markers <- bind_rows(markers_list)
nrow(markers)
n_distinct(markers$gene)
n_distinct(markers$celltype)
n_celltypes_unfiltered <- n_distinct(markers$celltype)

# Confirm that there are no missing values
stopifnot(!any(is.na(markers)))

# Confirm that the species are correctly set (Human, Mouse, or blank)
stopifnot(setequal(unique(markers$species), c("HS", "MM", "")))

# Keep only valid genes
markers <-
  markers |>
  mutate(gene = str_trim(gene)) |>
  filter(gene %in% valid_genes) |>
  distinct()
nrow(markers)
n_distinct(markers$gene)
n_distinct(markers$celltype)
distinct(markers, db, species, organ, celltype) |> count(db)

stopifnot(n_distinct(markers$celltype) == n_celltypes_unfiltered)

# Clean up cell type signatures -----

# Check the number of cell types per source
markers |>
  distinct(db, species, organ, celltype) |>
  count(db, species)

# Clean up cell type names and create a unique cell type identifier
markers <-
  markers |>
  unite(
    celltype_full,
    celltype,
    organ,
    species,
    db,
    sep = " | ",
    remove = FALSE,
    na.rm = TRUE
  ) |>
  mutate(
    celltype_full = str_replace_all(celltype_full, "\\|  \\|", "\\|"),
    celltype_full = str_replace_all(celltype_full, "\\|  \\|", "\\|")
  ) |>
  add_count(celltype_full, name = "n_genes") |>
  relocate(celltype_full)

# Check the number of signatures per source
markers |>
  distinct(db, species, celltype_full) |>
  count(db, species)

# Check the size of signatures
markers |>
  distinct(celltype_full, n_genes) |>
  pull(n_genes) |>
  quantile(seq(0, 1, 0.1))

# Check large signatures
markers |>
  filter(n_genes > 1000) |>
  distinct(db, celltype_full, n_genes) |>
  arrange(-n_genes)

# Remove very small and large signatures
markers <-
  markers |>
  filter(n_genes >= 5, n_genes < 1000) |>
  relocate(gene, .after = last_col()) |>
  arrange(celltype_full, gene)

stopifnot(n_distinct(markers$celltype) > n_celltypes_unfiltered * 0.95)

# Check the number of signatures per source
markers |>
  distinct(db, species, celltype_full) |>
  count(db, species)

# Check the size of signatures
markers |>
  distinct(celltype_full, n_genes) |>
  pull(n_genes) |>
  quantile(seq(0, 1, 0.1))

n_celltypes_full <- n_distinct(markers$celltype_full)
stopifnot(n_celltypes_full < n_distinct(markers$celltype) * 1.2)

# Resolve human/mouse gene symbols -----

# Use the species label only when a gene belongs to both species
markers_hs <-
  filter(
    markers,
    (gene %in% setdiff(genes_hs, genes_mm)) |
      (gene %in% intersect(genes_hs, genes_mm) & species == "HS")
  )
nrow(markers_hs)
markers_mm <-
  filter(
    markers,
    (gene %in% setdiff(genes_mm, genes_hs)) |
      (gene %in% intersect(genes_hs, genes_mm) & species == "MM")
  )
nrow(markers_mm)
markers_unk <-
  filter(markers, gene %in% intersect(genes_hs, genes_mm) & species == "")
nrow(markers_unk)
stopifnot(nrow(markers_unk) < 200)

# Marker subsets should add up to the initial markers table
stopifnot(
  nrow(markers_hs) + nrow(markers_mm) + nrow(markers_unk) == nrow(markers)
)

# Resolve canonical symbols
markers_hs <-
  markers_hs |>
  left_join(synonyms_hs, by = c("gene" = "alt_symbol")) |>
  rename(gene_hs = symbol) |>
  mutate(gene_hs = coalesce(gene_hs, gene))
markers_mm <-
  markers_mm |>
  left_join(synonyms_mm, by = c("gene" = "alt_symbol")) |>
  rename(gene_mm = symbol) |>
  mutate(gene_mm = coalesce(gene_mm, gene))

# Add orthologs
markers_hs <-
  left_join(
    markers_hs,
    orthologs,
    by = "gene_hs",
    relationship = "many-to-many"
  )
markers_mm <-
  left_join(
    markers_mm,
    orthologs,
    by = "gene_mm",
    relationship = "many-to-many"
  )

# Resolve ambiguous markers independently for each species
markers_unk <-
  markers_unk |>
  left_join(synonyms_hs, by = c("gene" = "alt_symbol")) |>
  rename(gene_hs = symbol) |>
  mutate(gene_hs = coalesce(gene_hs, gene)) |>
  left_join(synonyms_mm, by = c("gene" = "alt_symbol")) |>
  rename(gene_mm = symbol) |>
  mutate(gene_mm = coalesce(gene_mm, gene))

nrow(markers)
markers <- bind_rows(markers_hs, markers_mm, markers_unk)
markers <- rename(markers, gene_original = gene)
nrow(markers)
# 547688

# Confirm this step kept every cell type
stopifnot(n_distinct(markers$celltype_full) == n_celltypes_full)

# Check stats
distinct(markers, db, species, celltype_full) |> count(db, species)
markers |>
  distinct(celltype_full, n_genes) |>
  pull(n_genes) |>
  quantile(seq(0, 1, 0.1))
count(markers, celltype_full, n_genes, sort = TRUE)

# Prepare package -----

# Create package data
clustermole_markers_tbl <- markers
use_data(
  clustermole_markers_tbl,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz",
  version = 3
)
