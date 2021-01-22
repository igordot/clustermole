
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(tibble)
library(stringr)
library(glue)
library(janitor)
library(stringdist)
library(xml2)
library(usethis)


# HCOP gene orthologs -----

# import orthologs table
hcop_source <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz"
hcop_all <- read_tsv(hcop_source, progress = FALSE)

# extract all gene symbols (all valid genes)
valid_genes <- c(hcop_all$human_symbol, hcop_all$mouse_symbol) %>%
  unique() %>%
  sort()

# add the number of sources and string distances between gene names
hcop_clean <-
  hcop_all %>%
  select(gene_hs = human_symbol, gene_mm = mouse_symbol, sources = support) %>%
  mutate(num_sources = str_count(sources, ",") + 1) %>%
  mutate(distance = stringdist(toupper(gene_hs), toupper(gene_mm), method = "jw"))

# check if the sources were properly stored and parsed
nrow(hcop_clean)
hcop_clean %>%
  pull(num_sources) %>%
  quantile()
hcop_clean %>%
  pull(distance) %>%
  quantile()
barplot(table(hcop_clean$num_sources))
boxplot(distance ~ num_sources, data = hcop_clean, outline = FALSE)

# keep only the best orthologs
hcop_clean <-
  hcop_clean %>%
  group_by(gene_hs) %>%
  top_n(1, num_sources) %>%
  top_n(-1, distance) %>%
  ungroup()

# check the filtered results
nrow(hcop_clean)
hcop_clean %>%
  pull(num_sources) %>%
  quantile()
hcop_clean %>%
  pull(distance) %>%
  quantile()
barplot(table(hcop_clean$num_sources))
boxplot(distance ~ num_sources, data = hcop_clean, outline = FALSE)

# keep only the genes
hcop_clean <- hcop_clean %>% distinct(gene_hs, gene_mm)


# PanglaoDB signatures -----

# source: https://panglaodb.se/

# download gene signatures
panglao_source <- "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
panglao_all <- read_tsv(panglao_source, guess_max = 10000, progress = FALSE)

panglao_clean <-
  panglao_all %>%
  clean_names() %>%
  drop_na(gene_type, sensitivity_human, sensitivity_mouse) %>%
  mutate(
    db = "PanglaoDB",
    celltype = cell_type,
    gene = official_gene_symbol
  ) %>%
  separate_rows(species, sep = " ") %>%
  select(db, species, organ, celltype, gene)


# CellMarker signatures -----

# source: http://bio-bigdata.hrbmu.edu.cn/CellMarker/

# download gene signatures
cellmarker_source <- "http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt"
cellmarker_all <- read_tsv(cellmarker_source, guess_max = 10000, progress = FALSE)

cellmarker_clean <-
  cellmarker_all %>%
  clean_names() %>%
  mutate(
    db = "CellMarker",
    species = species_type,
    organ = str_remove(tissue_type, "Undefined"),
    celltype = str_c(cell_name, " (", cancer_type, ")"),
    celltype = str_remove_all(celltype, " \\(Normal\\)"),
    gene = str_remove_all(gene_symbol, "\\[|\\]| ")
  ) %>%
  drop_na(celltype, gene) %>%
  select(db, species, organ, celltype, gene) %>%
  separate_rows(gene, sep = ",")


# SaVanT signatures -----

# source: http://newpathways.mcdb.ucla.edu/savant-dev/

# download gene signatures
savant_source <- "http://newpathways.mcdb.ucla.edu/savant-dev/SaVanT_Signatures_Release01.zip"
savant_txt <- "SaVanT_Signatures_Release01.tab.txt"
savant_tmp <- tempfile(fileext = ".zip")
download.file(url = savant_source, destfile = savant_tmp, quiet = TRUE, mode = "wb")
savant_list <- strsplit(readLines(unz(savant_tmp, savant_txt)), "\t")
unlink(savant_tmp)
savant_all <- lapply(savant_list, tail, -1)
names(savant_all) <- sapply(savant_list, head, 1)
savant_all <- savant_all %>%
  enframe(name = "celltype", value = "gene") %>%
  unnest(gene)

# select the top 50 genes (default in SaVanT)
savant_all <-
  savant_all %>%
  filter(gene %in% valid_genes) %>%
  group_by(celltype) %>%
  slice(1:50) %>%
  ungroup()

savant_clean <-
  savant_all %>%
  mutate(
    db = "SaVanT",
    species = "",
    species =
      case_when(
        str_detect(celltype, "^MBA_") ~ "Mouse",
        str_detect(celltype, "^IMGN_") ~ "Mouse",
        str_detect(celltype, "^HBA_") ~ "Human",
        str_detect(celltype, "^HPCA_") ~ "Human",
        str_detect(celltype, "^MA_") ~ "Human",
        TRUE ~ species
      ),
    organ = ""
  ) %>%
  drop_na(celltype, gene) %>%
  select(db, species, organ, celltype, gene)


# MSigDB C8 (formerly SCSig) signatures -----

# link: http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C8

# download gene signatures
msigbd_source <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/msigdb_v7.2.xml"
msigbd_tmp <- tempfile(fileext = ".xml")
download.file(url = msigbd_source, destfile = msigbd_tmp, mode = "wb")
msigdb_doc <- read_xml(msigbd_tmp)
unlink(msigbd_tmp)

# convert XML document to a data frame
msigdb_records <- xml_find_all(msigdb_doc, xpath = ".//GENESET")
msigdb_all <-
  tibble(
    cat = xml_attr(msigdb_records, attr = "CATEGORY_CODE"),
    celltype = xml_attr(msigdb_records, attr = "STANDARD_NAME"),
    species = xml_attr(msigdb_records, attr = "ORGANISM"),
    members = xml_attr(msigdb_records, attr = "MEMBERS_MAPPING")
  ) %>%
  filter(cat == "C8")

# convert to one gene per row
msigdb_all <- mutate(msigdb_all, members_split = strsplit(members, "|", fixed = TRUE))
msigdb_all <- unnest(msigdb_all, cols = members_split, names_repair = "minimal")
msigdb_all <-
  msigdb_all %>%
  separate(
    col = members_split,
    into = c("source_gene", "gene", "entrez_gene"),
    sep = ","
  ) %>%
  mutate(entrez_gene = as.integer(entrez_gene)) %>%
  filter(entrez_gene > 0)

msigdb_clean <-
  msigdb_all %>%
  mutate(
    db = "MSigDB",
    organ = ""
  ) %>%
  select(db, species, organ, celltype, gene)


# xCell signatures -----

# source: http://xcell.ucsf.edu/

# download gene signatures
xcell_source <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5688663/bin/13059_2017_1349_MOESM3_ESM.xlsx"
xcell_tmp <- tempfile(fileext = ".xlsx")
download.file(url = xcell_source, destfile = xcell_tmp, quiet = TRUE, mode = "wb")
xcell_all <- read_xlsx(xcell_tmp, progress = FALSE)
unlink(xcell_tmp)

xcell_clean <-
  xcell_all %>%
  clean_names() %>%
  select(-number_of_genes) %>%
  gather(key = "k", value = "gene", -celltype_source_id) %>%
  mutate(celltype = celltype_source_id) %>%
  mutate(
    db = "xCell",
    species = "Human",
    organ = "",
    celltype = celltype_source_id
  ) %>%
  drop_na(celltype, gene) %>%
  select(db, species, organ, celltype, gene)


# Enrichr ARCHS4 tissues signatures -----

# source: http://amp.pharm.mssm.edu/archs4

# download gene signatures
archs_source <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=ARCHS4_Tissues"
archs_all <- clustermole::read_gmt(archs_source)

archs_clean <-
  archs_all %>%
  mutate(
    db = "ARCHS4",
    species = "",
    organ = ""
  ) %>%
  select(db, species, organ, celltype, gene)


# TISSUES human signatures -----

# source: https://tissues.jensenlab.org/

# download gene signatures
tissues_h_source <- "https://download.jensenlab.org/human_tissue_knowledge_full.tsv"
tissues_h_all <- read_tsv(tissues_h_source, col_names = FALSE, progress = FALSE)

tissues_h_clean <-
  tissues_h_all %>%
  filter(X3 != "BTO:0000000") %>%
  filter(str_detect(X4, "BTO:", negate = TRUE)) %>%
  mutate(
    db = "TISSUES",
    species = "Human",
    organ = "",
    celltype = X4,
    gene = X2
  ) %>%
  select(db, species, organ, celltype, gene)


# TISSUES mouse signatures -----

# source: https://tissues.jensenlab.org/

# download gene signatures
tissues_m_source <- "https://download.jensenlab.org/mouse_tissue_knowledge_full.tsv"
tissues_m_all <- read_tsv(tissues_m_source, col_names = FALSE, progress = FALSE)

tissues_m_clean <-
  tissues_m_all %>%
  filter(X3 != "BTO:0000000") %>%
  filter(str_detect(X4, "BTO:", negate = TRUE)) %>%
  mutate(
    db = "TISSUES",
    species = "Mouse",
    organ = "",
    celltype = X4,
    gene = X2
  ) %>%
  select(db, species, organ, celltype, gene)


# combine signatures -----

# combine
markers <-
  bind_rows(
    panglao_clean,
    cellmarker_clean,
    savant_clean,
    msigdb_clean,
    xcell_clean,
    archs_clean,
    tissues_h_clean,
    tissues_m_clean
  ) %>%
  filter(gene %in% valid_genes) %>%
  drop_na() %>%
  distinct()

# check the number of signatures per source
markers %>%
  distinct(db, species, celltype) %>%
  count(db, species)

# clean up
markers <-
  markers %>%
  mutate(
    species =
      case_when(
        species == "Human" ~ "Human",
        species == "Homo sapiens" ~ "Human",
        species == "Hs" ~ "Human",
        species == "Mouse" ~ "Mouse",
        species == "Mm" ~ "Mouse",
        species == "None" ~ "",
        TRUE ~ ""
      ),
    celltype_full = str_c(celltype, organ, species, db, sep = " | "),
    celltype_full = str_replace_all(celltype_full, "\\|  \\|", "\\|"),
    celltype_full = str_replace_all(celltype_full, "\\|  \\|", "\\|")
  ) %>%
  add_count(celltype_full, name = "n_genes") %>%
  relocate(celltype_full)

# check the number of signatures per source
markers %>%
  distinct(db, species, celltype) %>%
  count(db, species)

# check the size of signatures
markers %>%
  distinct(celltype_full, n_genes) %>%
  pull(n_genes) %>%
  quantile(seq(0, 1, 0.1))

# remove very small and large signatures
markers <-
  markers %>%
  filter(n_genes >= 5, n_genes <= 2500) %>%
  relocate(gene, .after = last_col()) %>%
  arrange(celltype_full, gene)

# check the number of signatures per source
markers %>%
  distinct(db, species, celltype) %>%
  count(db, species)

# check the size of signatures
markers %>%
  distinct(celltype_full, n_genes) %>%
  pull(n_genes) %>%
  quantile(seq(0, 1, 0.1))

# add human/mouse gene symbols (listed as either in the original table)
markers <- left_join(markers, hcop_clean, by = c("gene" = "gene_mm"))
markers <- left_join(markers, hcop_clean, by = c("gene" = "gene_hs"))
markers <- markers %>% mutate(gene_hs = if_else(is.na(gene_hs), gene, gene_hs))
markers <- markers %>% mutate(gene_mm = if_else(is.na(gene_mm), gene, gene_mm))
markers <- markers %>% rename(gene_original = gene)

# check stats
n_distinct(markers$celltype_full)
markers %>%
  pull(species) %>%
  table()
markers %>%
  distinct(celltype_full, n_genes) %>%
  pull(n_genes) %>%
  hist(breaks = 50, col = "gray20")
markers %>%
  distinct(celltype_full, n_genes) %>%
  arrange(-n_genes)


# prepare package -----

# create package data
clustermole_markers_tbl <- markers
use_data(
  clustermole_markers_tbl,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz",
  version = 3
)
