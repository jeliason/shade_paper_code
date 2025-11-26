# ============================================================================
# SETUP: Load Data and Prepare Spatial Point Patterns
# ============================================================================
# Module 0 of method comparison analysis
# Loads CRC data and creates spatial point pattern objects

library(tidyverse)
library(spatstat)
library(Matrix)
library(cmdstanr)
library(SHADE)
library(posterior)
library(mxfda)

source("crc_analysis/utils.R")
source("utils.R")

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./manuscript/images/CRC_analysis_paper/"

# ============================================================================
# LOAD DATA
# ============================================================================

df_raw <- read_csv("crc_analysis/data/CRC_cleaned.csv")
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

df_raw <- df_raw %>%
  dplyr::mutate(type = as.factor(type)) %>%
  mutate(type = fct_recode(type,"CAFs"="smooth muscle","hybrid E/M"="stroma","TAMs"="CD163+ macros","CTLs"="CD8+ T cells"))

targets <- c(
  "CTLs",
  "memory CD4+ T",
  "granulocytes"
)

sources <- c(
  "TAMs",
  "CAFs",
  "vasculature",
  "hybrid E/M",
  "tumor cells"
)

path <- "./crc_analysis/data/"
types_keep <- c(targets, sources)
num_types <- length(types_keep)

# ============================================================================
# PREPARE SPATIAL POINT PATTERNS
# ============================================================================

spots <- pt_data %>% pull(Spot)

dats <- df_raw %>%
  dplyr::filter(Spot %in% spots) %>%
  group_by(Spot) %>%
  group_map(~{
    .x <- .x %>%
      filter(type %in% types_keep) %>%
      droplevels() %>%
      mutate(Spot = .y$Spot)
  })

pats <- lapply(1:length(dats), \(i) {
  df <- dats[[i]]
  tryCatch({
    pat <- make_pat(df$X, df$Y, factor(df$type, levels=types_keep))
    sq_W <- owin(xrange = c(min(df$X), max(df$X)), yrange=c(min(df$Y), max(df$Y)))
    Window(pat) <- sq_W
    pat
  }, error=\(e) {
    print(paste("Error in pattern", i))
    stop(e)
  })
})

# Add metadata to patterns
for(i in seq_along(pats)) {
  spot <- unique(dats[[i]]$Spot)
  pt_info <- pt_data %>% filter(Spot == spot)
  attr(pats[[i]], "Spot") <- spot
  attr(pats[[i]], "Patient") <- pt_info$Patient
  attr(pats[[i]], "Group") <- pt_info$Group
}

cat("✓ Module 0: Loaded", length(pats), "spatial patterns\n")
