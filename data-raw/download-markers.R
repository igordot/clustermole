
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(tibble)
library(stringr)
library(glue)
library(stringdist)
library(usethis)


# HCOP gene orthologs -----

# import orthologs table
hcop_source <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz"
hcop_all <- read_tsv(hcop_source, progress = FALSE)

# extract all gene symbols (all valid genes)
valid_genes <- c(hcop_all$human_symbol, hcop_all$mouse_symbol) %>%
  unique() %>%
  sort()

# format orthologs table
hcop_clean <-
  hcop_all %>%
  select(gene_hs = human_symbol, gene_mm = mouse_symbol, sources = support) %>%
  mutate(num_sources = str_count(sources, ",") + 1) %>%
  mutate(distance = stringdist(toupper(gene_hs), toupper(gene_mm), method = "jw")) %>%
  group_by(gene_hs) %>%
  top_n(1, num_sources) %>%
  top_n(-1, distance) %>%
  ungroup() %>%
  distinct() %>%
  select(gene_hs, gene_mm)


# SCSig signatures -----

# link: http://www.gsea-msigdb.org/gsea/msigdb/supplementary_genesets.jsp

# download gene signatures
scsig_source <- "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"
scsig_all <- clustermole::read_gmt(scsig_source)

# download metadata (contains organism and organ)
scsig_metadata_source <- "http://software.broadinstitute.org/gsea/msigdb/supplemental/scsig.v1.0.metadata.xls"
scsig_tmp <- tempfile(fileext = ".xls")
download.file(scsig_metadata_source, scsig_tmp, quiet = TRUE, mode = "wb")
scsig_metadata <- read_xls(scsig_tmp, progress = FALSE)
scsig_metadata <- scsig_metadata %>% rename(celltype = `Gene Set Standard Name`)
unlink(scsig_tmp)

# combine gene signatures with the metadata
scsig_clean <-
  full_join(scsig_all, scsig_metadata, by = "celltype") %>%
  mutate(
    db = "SCSig",
    species = `Source Organism`,
    organ = `Source Organ System`
  ) %>%
  select(db, species, organ, celltype, gene)


# PanglaoDB signatures -----

# link: https://panglaodb.se/

# download gene signatures
panglao_all <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_30_Oct_2019.tsv.gz", progress = FALSE)

panglao_clean <-
  panglao_all %>%
  mutate(
    db = "PanglaoDB",
    celltype = `cell type`,
    gene = `official gene symbol`
  ) %>%
  separate_rows(species, sep = " ") %>%
  select(db, species, organ, celltype, gene)


# CellMarker signatures -----

# link: http://biocc.hrbmu.edu.cn/CellMarker/

# download gene signatures
cellmarker_source <- "http://biocc.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt"
cellmarker_all <- read_tsv(cellmarker_source, guess_max = 10000, progress = FALSE)

cellmarker_clean <-
  cellmarker_all %>%
  mutate(
    db = "CellMarker",
    species = speciesType,
    organ = str_remove(tissueType, "Undefined"),
    celltype = str_c(cellName, " (", cancerType, ")"),
    celltype = str_remove_all(celltype, " \\(Normal\\)"),
    gene = str_remove_all(geneSymbol, "\\[|\\]| ")
  ) %>%
  drop_na(celltype, gene) %>%
  select(db, species, organ, celltype, gene) %>%
  separate_rows(gene, sep = ",")


# SaVanT signatures -----

# link: http://newpathways.mcdb.ucla.edu/savant-dev/

# download gene signatures
savant_source <- "http://newpathways.mcdb.ucla.edu/savant-dev/SaVanT_Signatures_Release01.zip"
savant_txt <- "SaVanT_Signatures_Release01.tab.txt"
savant_tmp <- tempfile(fileext = ".zip")
download.file(savant_source, savant_tmp, quiet = TRUE, mode = "wb")
savant_list <- strsplit(readLines(unz(savant_tmp, savant_txt)), "\t")
unlink(savant_tmp)
savant_all <- lapply(savant_list, tail, -1)
names(savant_all) <- sapply(savant_list, head, 1)
savant_all <- savant_all %>%
  enframe(name = "celltype", value = "gene") %>%
  unnest(gene)
savant_all <- savant_all %>% filter(gene %in% valid_genes)

# select the top 50 genes (default in SaVanT)
savant_all <- savant_all %>%
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


# xCell signatures -----

# link: http://xcell.ucsf.edu/

# download gene signatures
xcell_source <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5688663/bin/13059_2017_1349_MOESM3_ESM.xlsx"
xcell_tmp <- tempfile(fileext = ".xlsx")
download.file(xcell_source, xcell_tmp, quiet = TRUE, mode = "wb")
xcell_all <- read_xlsx(xcell_tmp, progress = FALSE)
unlink(xcell_tmp)

xcell_clean <-
  xcell_all %>%
  select(-`# of genes`) %>%
  gather(key = "k", value = "gene", -Celltype_Source_ID) %>%
  mutate(celltype = Celltype_Source_ID) %>%
  mutate(
    db = "xCell",
    species = "Human",
    organ = "",
    celltype = Celltype_Source_ID
  ) %>%
  drop_na(celltype, gene) %>%
  select(db, species, organ, celltype, gene)


# combine signatures -----

# combine
markers <-
  bind_rows(
    scsig_clean,
    panglao_clean,
    cellmarker_clean,
    savant_clean,
    xcell_clean
  ) %>%
  filter(gene %in% valid_genes) %>%
  drop_na() %>%
  distinct()

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
  add_count(celltype_full, name = "n_genes")

markers <-
  markers %>%
  filter(n_genes >= 5, n_genes <= 2000) %>%
  select(-gene, gene) %>%
  arrange(celltype_full, gene)

# add human/mouse gene symbols (listed as either in the original table)
markers <- left_join(markers, hcop_clean, by = c("gene" = "gene_mm"))
markers <- left_join(markers, hcop_clean, by = c("gene" = "gene_hs"))
markers <- markers %>% mutate(gene_hs = if_else(is.na(gene_hs), gene, gene_hs))
markers <- markers %>% mutate(gene_mm = if_else(is.na(gene_mm), gene, gene_mm))
markers <- markers %>% rename(gene_original = gene)

# markers %>% pull(species) %>% table()
# hist(markers$n_genes, breaks = 50, col = "gray20")
# markers %>% distinct(celltype_full, n_genes) %>% arrange(-n_genes)
# n_distinct(markers$celltype_full)


# prepare package -----

# create package data
clustermole_markers_tbl = markers
use_data(
  clustermole_markers_tbl,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)