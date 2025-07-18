library(tidyverse)

# Extend timeout to avoid failures with large downloads
options(timeout = max(600, getOption("timeout")))

##### SINGLE-CELL DATA
# From original CRC dataset from Schurch et al (2020) - https://data.mendeley.com/datasets/mpjzbtfgfr/1

# Define data directory and file paths
data_dir <- "crc_analysis/data/"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

master_file <- file.path(data_dir, "CRC_master.csv")
alt_file    <- file.path(data_dir, "CRC_clusters_neighborhoods_markers.csv")
cleaned_file <- file.path(data_dir, "CRC_cleaned.csv")

# Download URL (verify this link stays current)
download_url <- "https://data.mendeley.com/public-files/datasets/mpjzbtfgfr/files/c24351b3-76d7-444f-9edf-0246356b0c78/file_downloaded"

# Check if master file exists; if not, download and rename alternate file
if (!file.exists(master_file)) {
  message("CRC_master.csv not found. Downloading dataset...")
  download.file(download_url, destfile = alt_file, mode = "wb")
  file.rename(alt_file, master_file)
}

# Read and process the data
df_raw <- read_csv(master_file, show_col_types = FALSE) %>%
  select(spots,ClusterName,`X:X`,`Y:Y`) %>%
  rename(type=ClusterName,X=`X:X`,Y=`Y:Y`,Spot=spots) %>%
  mutate(type = as.factor(type)) %>%
  filter(!(type %in% c("dirt", "undefined"))) %>%
  mutate(type = fct_collapse(
    type,
    stroma             = c("stroma", "nerves", "lymphatics"),
    "CD163+ macros"    = c("CD68+CD163+ macrophages", "CD163+ macrophages"),
    "CD68+ macros"     = c("CD11b+ monocytes", "CD11b+CD68+ macrophages", "CD68+ macrophages",
                           "CD68+ macrophages GzmB+", "CD11c+ DCs"),
    "generic immune"   = c("tumor cells / immune cells", "immune cells", "immune cells / vasculature"),
    "memory CD4+ T"    = "CD4+ T cells CD45RO+",
    "CD4+ T cells"     = c("CD4+ T cells", "CD4+ T cells GATA3+", "CD3+ T cells")
  )) %>%
  droplevels()

# Save the cleaned data
write_csv(df_raw, cleaned_file)

###### PATIENT METADATA
# From: https://www.cancerimagingarchive.net/collection/crc_ffpe-codex_cellneighs/
patients_xlsx_url  <- "https://www.cancerimagingarchive.net/wp-content/uploads/CRC_TMAs_patient_annotations.xlsx"
patients_xlsx_file <- file.path(data_dir, "CRC_patients.xlsx")
pt_metadata_file   <- file.path(data_dir, "CRC_pt_metadata.csv")

# Download patient annotations XLSX if not present
if (!file.exists(patients_xlsx_file)) {
  download.file(patients_xlsx_url, destfile = patients_xlsx_file, mode = "wb")
}

# Load and process patient metadata
pt_data <- readxl::read_excel(patients_xlsx_file) %>%
  slice(1:35) %>% # remove extra non-patient metadata at bottom
  select(Patient, Group, `TMA spot / region`,) %>%
  rename(Spot = `TMA spot / region`) %>%
  mutate(Group = factor(ifelse(Group == 1, "CRC", "DII"))) %>%
  separate(Spot, into = c("r1","r2")) %>%
  mutate(sp11=paste0(r1,"_A"),sp12=paste0(r1,"_B"),
         sp21=paste0(r2,"_A"),sp22=paste0(r2,"_B")) %>%
  select(-r1,-r2) %>%
  pivot_longer(-c(Patient,Group),values_to = "Spot") %>%
  select(-name)

# Save cleaned metadata
write_csv(pt_data, pt_metadata_file)
