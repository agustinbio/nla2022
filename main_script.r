rm(list = ls())
gc()
setwd("/mnt/storage/work/tt1/src")

NA_LIMIT = 100

library(dplyr)
library(tidyr)
library(sf)

#cargamos la base de datos y la especificación de los nombres de las variables
waterchem_raw <- read.csv(gzfile("nla22_waterchem_wide.csv.gz"))
waterchem_specs <-read.delim("nla22_waterchem_wide.txt")

#El objetivo de esta sección de código es dejar los nombres de las variables en forma prolija
#y generar un archivo: variables_description.csv que contiene el nombre de cada variable numérica
#y su correspondiente descripción

#nos quedamos con las columnas que tienen datos numéricos

result_cols <- grep("RESULT", names(waterchem_raw), value = TRUE)
result_cols <- result_cols[!grepl("UNITS", result_cols)]
id_columns <- c("UID", "SITE_ID", "UNIQUE_ID")

waterchem <- waterchem_raw[, c(id_columns, result_cols)]

#le sacamos la palabra RESULT a esas columnas
gen_waterchem_vardesc <- function(waterchem, waterchem_specs){
var_names_orig <- colnames(waterchem)
variables_names <- data.frame(PARAMETER = var_names_orig)
variables_description <- merge(variables_names,waterchem_specs)
variables_description <- variables_description[, colnames(variables_description) %in% c("PARAMETER", "LABEL")]
variables_description$PARAMETER <- gsub("_RESULT", "", variables_description$PARAMETER)
colnames(variables_description)[colnames(variables_description) == "PARAMETER"] <- "COLUMN_NAME"
return(variables_description)
}

waterchem_vardesc <- gen_waterchem_vardesc(waterchem, waterchem_specs)
write.csv(waterchem_vardesc, "waterchem_vardesc.csv", row.names = FALSE)
colnames(waterchem) <- gsub("_RESULT", "", colnames(waterchem))

#sacamos una columna con solo seis valores no nulos: la de nitrato-nitrito
sum(!is.na(waterchem$NITRATE_NITRITE_N))
waterchem <- waterchem[, !colnames(waterchem) %in% c("NITRATE_NITRITE_N")]

secchi_raw <- read.csv(gzfile("nla22_secchi.csv.gz"))
secchi_specs <-read.delim("nla22_secchi.txt")

secchi <- secchi_raw %>%
  select(UID, DISAPPEARS, REAPPEARS) %>%
  rename(
    SECCHI_DISAPPEARS = DISAPPEARS,
    SECCHI_REAPPEARS = REAPPEARS
  )

phabmets_raw <- read.csv(gzfile("nla2022_phabmets_wide_0.csv.gz"))
phabmets_specs <- read.delim("nla2022_phabmets_wide.txt")

phabmets <- phabmets_raw %>%
  # Sacamos la columna VISIT_NO
  select(-VISIT_NO) %>%
  # conservamos UID y solo las columnas que son numéricas y no constantes
  select(where(is.numeric)) %>%
  select(where(~ n_distinct(.) > 1))

generate_vardesc <- function(data, specs) {
var_names_orig <- colnames(data)
variables_names <- data.frame(COLUMN_NAME = var_names_orig)
variables_description <- merge(variables_names,specs)
variables_description <- variables_description[, colnames(variables_description) %in% c("COLUMN_NAME", "LABEL")]
return(variables_description)
}

phabmets_vardesc <- generate_vardesc(phabmets, phabmets_specs)
write.csv(phabmets_vardesc, "phabmets_vardesc.csv", row.names = FALSE)

atrazine_raw <- read.csv(gzfile("nla22_atrazine.csv.gz"))
atrazine_specs <- read.delim("nla22_atrazine.txt")

atrazine <- atrazine_raw[, c("UID", "RESULT")]
names(atrazine)[names(atrazine) == "RESULT"] <- "ATRAZINE"

algaltoxins_raw <- read.csv(gzfile("nla22_algaltoxins.csv.gz"))
algaltoxins_specs <- read.delim("nla22_algaltoxins.txt")

algaltoxins <- algaltoxins_raw %>%
  pivot_wider(
    id_cols = UID,                # Mantenemos UID como identificadora
    names_from = ANALYTE,         # Creamos nuevas columans de los valores de ANALYTE
    values_from = RESULT          # Llenamos las nuevas columnas con los valores de RESULT
  )

benthic_mmi_raw <- read.csv(gzfile("nla22_benthic_mmi_0.csv.gz"))
benthic_mmi_specs <- read.delim("nla22_benthic_mmi.txt")

benthic_mmi <- benthic_mmi_raw %>%
  select(-c(
    PUBLICATION_DATE, 
    UNIQUE_ID, 
    SITE_ID, 
    DATE_COL, 
    VISIT_NO, 
    PSTL_CODE, 
    BENT_MMI_COND_2017
  ))

benthic_mmi_vardesc <- generate_vardesc(benthic_mmi, benthic_mmi_specs)
write.csv(benthic_mmi_vardesc, "benthic_mmi_vardesc.csv", row.names = FALSE)

landscape_raw <- read.csv(gzfile("nla2022_landscape_wide_0.csv.gz"))
landscape_specs <- read.delim("nla2022_landscape_wide_0.txt")

landscape <- landscape_raw %>%
  select(-c(
    PUBLICATION_DATE,
    PSTL_CODE,
    GEOL_HUNT_DOM_DESC_WS,
    GEOL_HUNT_DOM_WS,
    GEOL_REEDBUSH_DOM_WS
  ))

landscape_vardesc <- generate_vardesc(landscape, landscape_specs)
write.csv(landscape_vardesc, "landscape_vardesc.csv", row.names = FALSE)

# Crear variables dummy (sección de código comentada, porque no están en uso en la última versión)
#landscape_dummies <- model.matrix(~ GEOL_REEDBUSH_DOM_WS - 1, data = landscape) # "-1" removes intercept
#landscape_dummies <- as.data.frame(landscape_dummies)

#landscape <- cbind(landscape, landscape_dummies)
#landscape <- landscape[, !(names(landscape) %in% "GEOL_REEDBUSH_DOM_WS")]

waterchem_uid <- waterchem %>%
  select(c(
    UID,
    SITE_ID,
    UNIQUE_ID
  ))

landscape <- waterchem_uid  %>%
  left_join(landscape, by = c("SITE_ID", "UNIQUE_ID"))

landscape <- landscape %>%
  select(-c(
    SITE_ID,
    UNIQUE_ID
  ))

siteinfo_raw <- read.csv(gzfile("nla22_siteinfo.csv.gz"))
siteinfo_specs <- read.delim("nla22_siteinfo.txt")

siteinfo <- siteinfo_raw %>%
  select(c(
    UID,
    SITE_ID,
    UNIQUE_ID,
    XCOORD,
    YCOORD
  ))

siteinfo$UID[siteinfo$UID == 0] <- NA

siteinfo_vardesc <- generate_vardesc(siteinfo, siteinfo_specs)
write.csv(siteinfo_vardesc, "siteinfo_vardesc.csv", row.names = FALSE)

siteinfo_bkp <- siteinfo
rm(siteinfo)

# Crear un archivo temporario
temp_gpkg <- tempfile(fileext = ".gpkg")

# Leer el archivo comprimido .gpkg.gz y descomprimir
compressed_con <- gzfile("siteinfo.gpkg.gz", open = "rb")
raw_data <- readBin(compressed_con, what = "raw", n = 1e8)  # Read all bytes (adjust n if needed)
close(compressed_con)

# Escribir los datos descomprimidos en un GeoPackage temporario
writeBin(raw_data, temp_gpkg)

# Cargar la capa GeoPackage
siteinfo_sf <- st_read(temp_gpkg)
unlink(temp_gpkg)
st_crs(siteinfo_sf)

usa_ec_crs <- "+proj=eqdc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +datum=NAD83 +units=m +no_defs"
site_usaec <- st_transform(siteinfo_sf, crs = usa_ec_crs)

names(site_usaec)[names(site_usaec) == "XCOORD"] <- "XCOORD_OLD"
names(site_usaec)[names(site_usaec) == "YCOORD"] <- "YCOORD_OLD"

coords <- st_coordinates(siteinfo_sf)

siteinfo_equidistant_df <- cbind(
  st_drop_geometry(site_usaec)
)

siteinfo_equidistant_df$XCOORD <- coords[, "X"]
siteinfo_equidistant_df$YCOORD <- coords[, "Y"]

siteinfo <- siteinfo_equidistant_df %>%
  select(c(
    UID,
    SITE_ID,
    UNIQUE_ID,
    XCOORD,
    YCOORD,
    NA_L2CODE
  ))

siteinfo$UID[siteinfo$UID == 0] <- NA

enterococci_raw <- read.csv(gzfile("nla22_enterococci.csv.gz"))
enterococci_specs <- read.delim("nla22_enterococci.txt")

enterococci <- enterococci_raw %>%
  select(c(
    UID,
    SITE_ID,
    UNIQUE_ID,
    ENTEROCOCCI_RESULT
  ))

enterococci_vardesc <- generate_vardesc(enterococci, enterococci_specs)
write.csv(enterococci_vardesc, "enterococci_vardesc.csv", row.names = FALSE)

merged_environment_full <- waterchem %>%
  left_join(secchi, by = "UID") %>%
  left_join(atrazine, by = "UID") %>%
  left_join(algaltoxins, by = "UID") %>%
  left_join(benthic_mmi, by = "UID") %>%
  left_join(landscape, by = "UID") %>%
  left_join(phabmets, by = "UID")

merged_environment <- waterchem %>%
  left_join(secchi, by = "UID") %>%
  left_join(atrazine, by = "UID") %>%
  left_join(algaltoxins, by = "UID") %>%
  left_join(benthic_mmi, by = "UID")

phytoplankton_wide <- read.csv(gzfile("nla2022_phytoplanktoncount_wide.csv.gz"))
phytoplankton_groups_abundance <- phytoplankton_wide %>%
  mutate(ALGAL_GROUP = trimws(ALGAL_GROUP)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = ALGAL_GROUP, 
    values_from = ABUNDANCE,
    values_fill = 0,
    values_fn = sum
  )
phytoplankton_groups_density <- phytoplankton_wide %>%
  mutate(ALGAL_GROUP = trimws(ALGAL_GROUP)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = ALGAL_GROUP, 
    values_from = DENSITY,
    values_fill = 0,
    values_fn = sum
  )
phytoplankton_groups_biovolume <- phytoplankton_wide %>%
  mutate(ALGAL_GROUP = trimws(ALGAL_GROUP)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = ALGAL_GROUP, 
    values_from = BIOVOLUME,
    values_fill = 0,
    values_fn = sum
  )

polish_names_algal_groups <- function(phytoplankton_groups){
names(phytoplankton_groups) <- gsub(" ", "_", names(phytoplankton_groups))
names(phytoplankton_groups) <- gsub("-", "_", names(phytoplankton_groups))
phytoplankton_groups$SITE_ID <- NULL
phytoplankton_groups$UNIQUE_ID <- NULL
return(phytoplankton_groups)
}

phytoplankton_groups_abundance <- polish_names_algal_groups(phytoplankton_groups_abundance)
phytoplankton_groups_density <-polish_names_algal_groups(phytoplankton_groups_density)
phytoplankton_groups_biovolume <-polish_names_algal_groups(phytoplankton_groups_biovolume)

merged_environment_full_w_algal <- merged_environment_full %>%
  left_join(phytoplankton_groups_abundance, by = "UID")

waterchem$NITRATE_N <- NULL
waterchem$NITRITE_N <- NULL

waterchem_w_algal <- waterchem %>%
  left_join(phytoplankton_groups_abundance, by = "UID")

vardesc_remaining <- data.frame(
  COLUMN_NAME = c("SECCHI_DISAPPEARS", "SECCHI_REAPPEARS", 
                  "ATRAZINE", "CYLSPER", "MICX"),
  LABEL = c("Depth at which secchi disk disappears", 
            "Depth at which secchi disk reappears",
            "Atrazine", 
            "Cylindrospermopsin", 
            "Microcystin"),
  stringsAsFactors = FALSE
)

vardesc_merged <- rbind(waterchem_vardesc, vardesc_remaining, phabmets_vardesc, benthic_mmi_vardesc, landscape_vardesc)
write.csv(vardesc_merged, "merged_vardesc.csv", row.names = FALSE)

summarize_na <- function(df, vardesc) {
  
  cols <- names(df)
  summary_list <- lapply(cols, function(col_name) {
    col <- df[[col_name]]
    na_count <- sum(is.na(col))
    
    description <- vardesc$LABEL[vardesc$COLUMN_NAME == col_name]
    
    data.frame(
      Columna = col_name,
      NAs = na_count,
      Descripcion = ifelse(length(description) > 0, description, NA)
    )
  })
  
  do.call(rbind, summary_list)
}

na_count <- summarize_na(merged_environment, vardesc_merged)
write.csv(na_count, "na_count.csv", row.names = FALSE)

na_count_full <- summarize_na(merged_environment_full, vardesc_merged)
write.csv(na_count_full, "na_count-full.csv", row.names = FALSE)

merged_environment_complete_fields <- merged_environment[, colSums(is.na(merged_environment)) < NA_LIMIT]
merged_environment_full_complete_fields <- merged_environment_full[, colSums(is.na(merged_environment_full)) < NA_LIMIT]
merged_environment_full_w_algal_complete_fields <- merged_environment_full_w_algal[, colSums(is.na(merged_environment_full_w_algal)) < NA_LIMIT]

zooplankton_wide <- read.csv(gzfile("nla22_zooplanktoncount_wide.csv.gz"))
benthic_wide <- read.csv(gzfile("nla22_benthic_counts.csv.gz"))


phytoplankton <- phytoplankton_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = ABUNDANCE,
    values_fill = 0,
    values_fn = sum
  )

zooplankton <- zooplankton_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = COUNT,
    values_fill = 0,
    values_fn = sum
  )

benthic <- benthic_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = TOTAL,
    values_fill = 0,
    values_fn = sum
  )


phytoplankton_biovolume <- phytoplankton_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = BIOVOLUME,
    values_fill = 0,
    values_fn = sum
  )

zooplankton_biomass <- zooplankton_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = BIOMASS,
    values_fill = 0,
    values_fn = sum
  )

phytoplankton_density <- phytoplankton_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = DENSITY,
    values_fill = 0,
    values_fn = sum
  )

zooplankton_300 <- zooplankton_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = COUNT_300,
    values_fill = 0,
    values_fn = sum
  )

zooplankton_biomass_300 <- zooplankton_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = BIOMASS_300,
    values_fill = 0,
    values_fn = sum
  )


benthic_300 <- benthic_wide %>%
  mutate(TARGET_TAXON = trimws(TARGET_TAXON)) %>% 
  pivot_wider(
    id_cols = id_columns,
    names_from = TARGET_TAXON, 
    values_from = TOTAL_300,
    values_fill = 0,
    values_fn = sum
  )

add_prefix <- function(df, prefix) {
  df %>% rename_with(~ paste0(prefix, .x), -id_columns)
}

phyto_prefixed <- add_prefix(phytoplankton, "phyto_")
zoo_prefixed <- add_prefix(zooplankton, "zoo_")
benthic_prefixed <- add_prefix(benthic, "benthic_")

merged_species <- phyto_prefixed %>% 
  inner_join(zoo_prefixed, by = id_columns) %>% 
  inner_join(benthic_prefixed, by = id_columns)


source("perform_analysis_new.R")


sink("salida-spatial.txt")

results_all_species <- perform_analysis(merged_species, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "all_species", "full-environment")
results_zooplankton_w_algal <- perform_analysis(zooplankton, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-w-algal")
results_zooplankton <- perform_analysis(zooplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment")
results_phytoplankton <- perform_analysis(phytoplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton", "full-environment")
results_benthic <- perform_analysis(benthic, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "benthic", "full-environment")

results_phytoplankton_biovolume  <- perform_analysis(phytoplankton_biovolume, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton_biovolume", "full-environment")
results_zooplankton_biomass_w_algal  <- perform_analysis(zooplankton_biomass, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass", "full-environment-w-algal")
results_zooplankton_biomass  <- perform_analysis(zooplankton_biomass, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass", "full-environment")
results_phytoplankton_density  <- perform_analysis(phytoplankton_density, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton_density", "full-environment")
results_zooplankton_300_w_algal  <- perform_analysis(zooplankton_300, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_300", "full-environment-w-algal")
results_zooplankton_300  <- perform_analysis(zooplankton_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_300", "full-environment")
results_benthic_300  <- perform_analysis(benthic_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "benthic_300", "full-environment")


results_zooplankton_waterchem <- perform_analysis(zooplankton, waterchem, siteinfo, vardesc_merged, "zooplankton", "waterchem")
results_zooplankton_biomass_waterchem <- perform_analysis(zooplankton_biomass, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem")

results_all_species_waterchem <- perform_analysis(merged_species, waterchem, siteinfo, vardesc_merged, "all_species", "waterchem")
results_zooplankton_waterchem_w_algal <- perform_analysis(zooplankton, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton", "waterchem-w-algal")
results_phytoplankton_waterchem <- perform_analysis(phytoplankton, waterchem, siteinfo, vardesc_merged, "phytoplankton", "waterchem")
results_benthic_waterchem <- perform_analysis(benthic, waterchem, siteinfo, vardesc_merged, "benthic", "waterchem")

results_phytoplankton_biovolume_waterchem <- perform_analysis(phytoplankton_biovolume, waterchem, siteinfo, vardesc_merged, "phytoplankton_biovolume", "waterchem")
results_zooplankton_biomass_waterchem_w_algal <- perform_analysis(zooplankton_biomass, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-w-algal")
results_phytoplankton_density_waterchem <- perform_analysis(phytoplankton_density, waterchem, siteinfo, vardesc_merged, "phytoplankton_density", "waterchem")
results_zooplankton_300_waterchem <- perform_analysis(zooplankton_300, waterchem, siteinfo, vardesc_merged, "zooplankton_300", "waterchem")
results_benthic_300_waterchem <- perform_analysis(benthic_300, waterchem, siteinfo, vardesc_merged, "benthic_300", "waterchem")

results_zooplankton_300_waterchem_w_algal  <- perform_analysis(zooplankton_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-w-algal")
results_zooplankton_biomass_300  <- perform_analysis(zooplankton_biomass_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass_300", "full-environment")
results_zooplankton_biomass_300_w_algal  <- perform_analysis(zooplankton_biomass_300, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass_300", "full-environment-w-algal")

results_zooplankton_biomass_300_waterchem  <- perform_analysis(zooplankton_biomass_300, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem")
results_zooplankton_biomass_300_waterchem_w_algal  <- perform_analysis(zooplankton_biomass_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-w-algal")

results_zooplankton_waterchem_reduced <- perform_analysis(zooplankton, waterchem, siteinfo, vardesc_merged, "zooplankton", "waterchem-reduced")
results_zooplankton_reduced <- perform_analysis(zooplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-reduced")

results_all_species_reduced <- perform_analysis(merged_species, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "all_species", "full-environment-reduced")
results_zooplankton_w_algal_reduced <- perform_analysis(zooplankton, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-w-algal-reduced")
results_zooplankton_reduced <- perform_analysis(zooplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-reduced")
results_phytoplankton_reduced <- perform_analysis(phytoplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton", "full-environment-reduced")
results_benthic_reduced <- perform_analysis(benthic, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "benthic", "full-environment-reduced")

results_phytoplankton_biovolume_reduced  <- perform_analysis(phytoplankton_biovolume, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton_biovolume", "full-environment-reduced")
results_zooplankton_biomass_w_algal_reduced  <- perform_analysis(zooplankton_biomass, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass", "full-environment-w-algal-reduced")
results_zooplankton_biomass_reduced  <- perform_analysis(zooplankton_biomass, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass", "full-environment-reduced")
results_phytoplankton_density_reduced  <- perform_analysis(phytoplankton_density, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton_density", "full-environment-reduced")
results_zooplankton_300_w_algal_reduced  <- perform_analysis(zooplankton_300, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_300", "full-environment-w-algal-reduced")
results_zooplankton_300_reduced  <- perform_analysis(zooplankton_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_300", "full-environment-reduced")
results_benthic_300_reduced  <- perform_analysis(benthic_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "benthic_300", "full-environment-reduced")

results_zooplankton_waterchem_reduced <- perform_analysis(zooplankton, waterchem, siteinfo, vardesc_merged, "zooplankton", "waterchem-reduced")
results_zooplankton_biomass_waterchem_reduced <- perform_analysis(zooplankton_biomass, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-reduced")

results_all_species_waterchem_reduced <- perform_analysis(merged_species, waterchem, siteinfo, vardesc_merged, "all_species", "waterchem-reduced")
results_zooplankton_waterchem_w_algal_reduced <- perform_analysis(zooplankton, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton", "waterchem-w-algal-reduced")
results_phytoplankton_waterchem_reduced <- perform_analysis(phytoplankton, waterchem, siteinfo, vardesc_merged, "phytoplankton", "waterchem-reduced")
results_benthic_waterchem_reduced <- perform_analysis(benthic, waterchem, siteinfo, vardesc_merged, "benthic", "waterchem-reduced")

results_phytoplankton_biovolume_waterchem_reduced <- perform_analysis(phytoplankton_biovolume, waterchem, siteinfo, vardesc_merged, "phytoplankton_biovolume", "waterchem-reduced")
results_zooplankton_biomass_waterchem_w_algal_reduced <- perform_analysis(zooplankton_biomass, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-w-algal-reduced")
results_phytoplankton_density_waterchem_reduced <- perform_analysis(phytoplankton_density, waterchem, siteinfo, vardesc_merged, "phytoplankton_density", "waterchem-reduced")
results_zooplankton_300_waterchem_reduced <- perform_analysis(zooplankton_300, waterchem, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-reduced")
results_benthic_300_waterchem_reduced <- perform_analysis(benthic_300, waterchem, siteinfo, vardesc_merged, "benthic_300", "waterchem-reduced")

results_zooplankton_300_waterchem_w_algal_reduced  <- perform_analysis(zooplankton_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-w-algal-reduced")
results_zooplankton_biomass_300_reduced  <- perform_analysis(zooplankton_biomass_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass_300", "full-environment-reduced")
results_zooplankton_biomass_300_w_algal_reduced  <- perform_analysis(zooplankton_biomass_300, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass_300", "full-environment-w-algal-reduced")

results_zooplankton_biomass_300_waterchem_reduced  <- perform_analysis(zooplankton_biomass_300, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-reduced")
results_zooplankton_biomass_300_waterchem_w_algal_reduced  <- perform_analysis(zooplankton_biomass_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-w-algal-reduced")

results_all_species_no_outliers <- perform_analysis(merged_species, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "all_species", "full-environment-no-outliers")
results_zooplankton_w_algal_no_outliers <- perform_analysis(zooplankton, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-w-algal-no-outliers")
results_zooplankton_no_outliers <- perform_analysis(zooplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-no-outliers")
results_phytoplankton_no_outliers <- perform_analysis(phytoplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton", "full-environment-no-outliers")
results_benthic_no_outliers <- perform_analysis(benthic, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "benthic", "full-environment-no-outliers")

results_phytoplankton_biovolume_no_outliers  <- perform_analysis(phytoplankton_biovolume, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton_biovolume", "full-environment-no-outliers")
results_zooplankton_biomass_w_algal_no_outliers  <- perform_analysis(zooplankton_biomass, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass", "full-environment-w-algal-no-outliers")
results_zooplankton_biomass_no_outliers  <- perform_analysis(zooplankton_biomass, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass", "full-environment-no-outliers")
results_phytoplankton_density_no_outliers  <- perform_analysis(phytoplankton_density, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton_density", "full-environment-no-outliers")
results_zooplankton_300_w_algal_no_outliers  <- perform_analysis(zooplankton_300, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_300", "full-environment-w-algal-no-outliers")
results_zooplankton_300_no_outliers  <- perform_analysis(zooplankton_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_300", "full-environment-no-outliers")
results_benthic_300_no_outliers  <- perform_analysis(benthic_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "benthic_300", "full-environment-no-outliers")

results_zooplankton_waterchem_no_outliers <- perform_analysis(zooplankton, waterchem, siteinfo, vardesc_merged, "zooplankton", "waterchem-no-outliers")
results_zooplankton_biomass_waterchem_no_outliers <- perform_analysis(zooplankton_biomass, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-no-outliers")

results_all_species_waterchem_no_outliers <- perform_analysis(merged_species, waterchem, siteinfo, vardesc_merged, "all_species", "waterchem-no-outliers")
results_zooplankton_waterchem_w_algal_no_outliers <- perform_analysis(zooplankton, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton", "waterchem-w-algal-no-outliers")
results_phytoplankton_waterchem_no_outliers <- perform_analysis(phytoplankton, waterchem, siteinfo, vardesc_merged, "phytoplankton", "waterchem-no-outliers")
results_benthic_waterchem_no_outliers <- perform_analysis(benthic, waterchem, siteinfo, vardesc_merged, "benthic", "waterchem-no-outliers")

results_phytoplankton_biovolume_waterchem_no_outliers <- perform_analysis(phytoplankton_biovolume, waterchem, siteinfo, vardesc_merged, "phytoplankton_biovolume", "waterchem-no-outliers")
results_zooplankton_biomass_waterchem_w_algal_no_outliers <- perform_analysis(zooplankton_biomass, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-w-algal-no-outliers")
results_phytoplankton_density_waterchem_no_outliers <- perform_analysis(phytoplankton_density, waterchem, siteinfo, vardesc_merged, "phytoplankton_density", "waterchem-no-outliers")
results_zooplankton_300_waterchem_no_outliers <- perform_analysis(zooplankton_300, waterchem, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-no-outliers")
results_benthic_300_waterchem_no_outliers <- perform_analysis(benthic_300, waterchem, siteinfo, vardesc_merged, "benthic_300", "waterchem-no-outliers")

results_zooplankton_300_waterchem_w_algal_no_outliers  <- perform_analysis(zooplankton_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-w-algal-no-outliers")
results_zooplankton_biomass_300_no_outliers  <- perform_analysis(zooplankton_biomass_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass_300", "full-environment-no-outliers")
results_zooplankton_biomass_300_w_algal_no_outliers  <- perform_analysis(zooplankton_biomass_300, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass_300", "full-environment-w-algal-no-outliers")

results_zooplankton_biomass_300_waterchem_no_outliers  <- perform_analysis(zooplankton_biomass_300, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-no-outliers")
results_zooplankton_biomass_300_waterchem_w_algal_no_outliers  <- perform_analysis(zooplankton_biomass_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-w-algal-no-outliers")

results_zooplankton_waterchem_reduced_no_outliers <- perform_analysis(zooplankton, waterchem, siteinfo, vardesc_merged, "zooplankton", "waterchem-reduced-no-outliers")
results_zooplankton_reduced_no_outliers <- perform_analysis(zooplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-reduced-no-outliers")

# results_all_species_reduced_no_outliers <- perform_analysis(merged_species, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "all_species", "full-environment-reduced-no-outliers")
results_zooplankton_w_algal_reduced_no_outliers <- perform_analysis(zooplankton, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-w-algal-reduced-no-outliers")
results_zooplankton_reduced_no_outliers <- perform_analysis(zooplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton", "full-environment-reduced-no-outliers")
results_phytoplankton_reduced_no_outliers <- perform_analysis(phytoplankton, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton", "full-environment-reduced-no-outliers")
results_benthic_reduced_no_outliers <- perform_analysis(benthic, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "benthic", "full-environment-reduced-no-outliers")

results_phytoplankton_biovolume_reduced_no_outliers  <- perform_analysis(phytoplankton_biovolume, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton_biovolume", "full-environment-reduced-no-outliers")
results_zooplankton_biomass_w_algal_reduced_no_outliers  <- perform_analysis(zooplankton_biomass, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass", "full-environment-w-algal-reduced-no-outliers")
results_zooplankton_biomass_reduced_no_outliers  <- perform_analysis(zooplankton_biomass, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass", "full-environment-reduced-no-outliers")
results_phytoplankton_density_reduced_no_outliers  <- perform_analysis(phytoplankton_density, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "phytoplankton_density", "full-environment-reduced-no-outliers")
results_zooplankton_300_w_algal_reduced_no_outliers  <- perform_analysis(zooplankton_300, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_300", "full-environment-w-algal-reduced-no-outliers")
results_zooplankton_300_reduced_no_outliers  <- perform_analysis(zooplankton_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_300", "full-environment-reduced-no-outliers")
results_benthic_300_reduced_no_outliers  <- perform_analysis(benthic_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "benthic_300", "full-environment-reduced-no-outliers")

results_zooplankton_waterchem_reduced_no_outliers <- perform_analysis(zooplankton, waterchem, siteinfo, vardesc_merged, "zooplankton", "waterchem-reduced-no-outliers")
results_zooplankton_biomass_waterchem_reduced_no_outliers <- perform_analysis(zooplankton_biomass, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-reduced-no-outliers")

results_all_species_waterchem_reduced_no_outliers <- perform_analysis(merged_species, waterchem, siteinfo, vardesc_merged, "all_species", "waterchem-reduced-no-outliers")
results_zooplankton_waterchem_w_algal_reduced_no_outliers <- perform_analysis(zooplankton, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton", "waterchem-w-algal-reduced-no-outliers")
results_phytoplankton_waterchem_reduced_no_outliers <- perform_analysis(phytoplankton, waterchem, siteinfo, vardesc_merged, "phytoplankton", "waterchem-reduced-no-outliers")
results_benthic_waterchem_reduced_no_outliers <- perform_analysis(benthic, waterchem, siteinfo, vardesc_merged, "benthic", "waterchem-reduced-no-outliers")

results_phytoplankton_biovolume_waterchem_reduced_no_outliers <- perform_analysis(phytoplankton_biovolume, waterchem, siteinfo, vardesc_merged, "phytoplankton_biovolume", "waterchem-reduced-no-outliers")
results_zooplankton_biomass_waterchem_w_algal_reduced_no_outliers <- perform_analysis(zooplankton_biomass, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-w-algal-reduced-no-outliers")
results_phytoplankton_density_waterchem_reduced_no_outliers <- perform_analysis(phytoplankton_density, waterchem, siteinfo, vardesc_merged, "phytoplankton_density", "waterchem-reduced-no-outliers")
results_zooplankton_300_waterchem_reduced_no_outliers <- perform_analysis(zooplankton_300, waterchem, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-reduced-no-outliers")
results_benthic_300_waterchem_reduced_no_outliers <- perform_analysis(benthic_300, waterchem, siteinfo, vardesc_merged, "benthic_300", "waterchem-reduced-no-outliers")

results_zooplankton_300_waterchem_w_algal_reduced_no_outliers  <- perform_analysis(zooplankton_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-w-algal-reduced-no-outliers")
results_zooplankton_biomass_300_reduced_no_outliers  <- perform_analysis(zooplankton_biomass_300, merged_environment_full_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass_300", "full-environment-reduced-no-outliers")
results_zooplankton_biomass_300_w_algal_reduced_no_outliers  <- perform_analysis(zooplankton_biomass_300, merged_environment_full_w_algal_complete_fields, siteinfo, vardesc_merged, "zooplankton_biomass_300", "full-environment-w-algal-reduced-no-outliers")

results_zooplankton_biomass_300_waterchem_reduced_no_outliers  <- perform_analysis(zooplankton_biomass_300, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-reduced-no-outliers")
results_zooplankton_biomass_300_waterchem_w_algal_reduced_no_outliers  <- perform_analysis(zooplankton_biomass_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-w-algal-reduced-no-outliers")

results_zooplankton_waterchem_minimal <- perform_analysis(zooplankton, waterchem, siteinfo, vardesc_merged, "zooplankton", "waterchem-minimal")
results_zooplankton_biomass_waterchem_minimal <- perform_analysis(zooplankton_biomass, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-minimal")

results_all_species_waterchem_minimal <- perform_analysis(merged_species, waterchem, siteinfo, vardesc_merged, "all_species", "waterchem-minimal")
results_zooplankton_waterchem_w_algal_minimal <- perform_analysis(zooplankton, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton", "waterchem-w-algal-minimal")
results_phytoplankton_waterchem_minimal <- perform_analysis(phytoplankton, waterchem, siteinfo, vardesc_merged, "phytoplankton", "waterchem-minimal")
results_benthic_waterchem_minimal <- perform_analysis(benthic, waterchem, siteinfo, vardesc_merged, "benthic", "waterchem-minimal")

results_phytoplankton_biovolume_waterchem_minimal <- perform_analysis(phytoplankton_biovolume, waterchem, siteinfo, vardesc_merged, "phytoplankton_biovolume", "waterchem-minimal")
results_zooplankton_biomass_waterchem_w_algal_minimal <- perform_analysis(zooplankton_biomass, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-w-algal-minimal")
results_phytoplankton_density_waterchem_minimal <- perform_analysis(phytoplankton_density, waterchem, siteinfo, vardesc_merged, "phytoplankton_density", "waterchem-minimal")
results_zooplankton_300_waterchem_minimal <- perform_analysis(zooplankton_300, waterchem, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-minimal")
results_benthic_300_waterchem_minimal <- perform_analysis(benthic_300, waterchem, siteinfo, vardesc_merged, "benthic_300", "waterchem-minimal")

results_zooplankton_300_waterchem_w_algal_minimal  <- perform_analysis(zooplankton_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-w-algal-minimal")

results_zooplankton_biomass_300_waterchem_minimal  <- perform_analysis(zooplankton_biomass_300, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-minimal")
results_zooplankton_biomass_300_waterchem_w_algal_minimal  <- perform_analysis(zooplankton_biomass_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-w-algal-minimal")

results_zooplankton_waterchem_minimal_no_outliers <- perform_analysis(zooplankton, waterchem, siteinfo, vardesc_merged, "zooplankton", "waterchem-minimal-no-outliers")
results_zooplankton_biomass_waterchem_minimal_no_outliers <- perform_analysis(zooplankton_biomass, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-minimal-no-outliers")

results_all_species_waterchem_minimal_no_outliers <- perform_analysis(merged_species, waterchem, siteinfo, vardesc_merged, "all_species", "waterchem-minimal-no-outliers")
results_zooplankton_waterchem_w_algal_minimal_no_outliers <- perform_analysis(zooplankton, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton", "waterchem-w-algal-minimal-no-outliers")
results_phytoplankton_waterchem_minimal_no_outliers <- perform_analysis(phytoplankton, waterchem, siteinfo, vardesc_merged, "phytoplankton", "waterchem-minimal-no-outliers")
results_benthic_waterchem_minimal_no_outliers <- perform_analysis(benthic, waterchem, siteinfo, vardesc_merged, "benthic", "waterchem-minimal-no-outliers")

results_phytoplankton_biovolume_waterchem_minimal_no_outliers <- perform_analysis(phytoplankton_biovolume, waterchem, siteinfo, vardesc_merged, "phytoplankton_biovolume", "waterchem-minimal-no-outliers")
results_zooplankton_biomass_waterchem_w_algal_minimal_no_outliers <- perform_analysis(zooplankton_biomass, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass", "waterchem-w-algal-minimal-no-outliers")
results_phytoplankton_density_waterchem_minimal_no_outliers <- perform_analysis(phytoplankton_density, waterchem, siteinfo, vardesc_merged, "phytoplankton_density", "waterchem-minimal-no-outliers")
results_zooplankton_300_waterchem_minimal_no_outliers <- perform_analysis(zooplankton_300, waterchem, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-minimal-no-outliers")
results_benthic_300_waterchem_minimal_no_outliers <- perform_analysis(benthic_300, waterchem, siteinfo, vardesc_merged, "benthic_300", "waterchem-minimal-no-outliers")

results_zooplankton_300_waterchem_w_algal_minimal_no_outliers  <- perform_analysis(zooplankton_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_300", "waterchem-w-algal-minimal-no-outliers")

results_zooplankton_biomass_300_waterchem_minimal_no_outliers  <- perform_analysis(zooplankton_biomass_300, waterchem, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-minimal-no-outliers")
results_zooplankton_biomass_300_waterchem_w_algal_minimal_no_outliers  <- perform_analysis(zooplankton_biomass_300, waterchem_w_algal, siteinfo, vardesc_merged, "zooplankton_biomass_300", "waterchem-w-algal-minimal-no-outliers")

all_results <- list(
  all_species = results_all_species,
  zooplankton_w_algal = results_zooplankton_w_algal,
  zooplankton = results_zooplankton,
  phytoplankton = results_phytoplankton,
  benthic = results_benthic,
  phytoplankton_biovolume = results_phytoplankton_biovolume,
  zooplankton_biomass_w_algal = results_zooplankton_biomass_w_algal,
  zooplankton_biomass = results_zooplankton_biomass,
  phytoplankton_density = results_phytoplankton_density,
  zooplankton_300_w_algal = results_zooplankton_300_w_algal,
  zooplankton_300 = results_zooplankton_300,
  benthic_300 = results_benthic_300,
  zooplankton_waterchem = results_zooplankton_waterchem,
  zooplankton_biomass_waterchem = results_zooplankton_biomass_waterchem,
  all_species_waterchem = results_all_species_waterchem,
  zooplankton_waterchem_w_algal = results_zooplankton_waterchem_w_algal,
  phytoplankton_waterchem = results_phytoplankton_waterchem,
  benthic_waterchem = results_benthic_waterchem,
  phytoplankton_biovolume_waterchem = results_phytoplankton_biovolume_waterchem,
  zooplankton_biomass_waterchem_w_algal = results_zooplankton_biomass_waterchem_w_algal,
  phytoplankton_density_waterchem = results_phytoplankton_density_waterchem,
  zooplankton_300_waterchem = results_zooplankton_300_waterchem,
  benthic_300_waterchem = results_benthic_300_waterchem,
  zooplankton_300_waterchem_w_algal = results_zooplankton_300_waterchem_w_algal,
  zooplankton_biomass_300 = results_zooplankton_biomass_300,
  zooplankton_biomass_300_w_algal = results_zooplankton_biomass_300_w_algal,
  zooplankton_biomass_300_waterchem = results_zooplankton_biomass_300_waterchem,
  zooplankton_biomass_300_waterchem_w_algal = results_zooplankton_biomass_300_waterchem_w_algal,
  all_species_reduced = results_all_species_reduced,
  zooplankton_w_algal_reduced = results_zooplankton_w_algal_reduced,
  zooplankton_reduced = results_zooplankton_reduced,
  phytoplankton_reduced = results_phytoplankton_reduced,
  benthic_reduced = results_benthic_reduced,
  phytoplankton_biovolume_reduced = results_phytoplankton_biovolume_reduced,
  zooplankton_biomass_w_algal_reduced = results_zooplankton_biomass_w_algal_reduced,
  zooplankton_biomass_reduced = results_zooplankton_biomass_reduced,
  phytoplankton_density_reduced = results_phytoplankton_density_reduced,
  zooplankton_300_w_algal_reduced = results_zooplankton_300_w_algal_reduced,
  zooplankton_300_reduced = results_zooplankton_300_reduced,
  benthic_300_reduced = results_benthic_300_reduced,
  zooplankton_waterchem_reduced = results_zooplankton_waterchem_reduced,
  zooplankton_biomass_waterchem_reduced = results_zooplankton_biomass_waterchem_reduced,
  all_species_waterchem_reduced = results_all_species_waterchem_reduced,
  zooplankton_waterchem_w_algal_reduced = results_zooplankton_waterchem_w_algal_reduced,
  phytoplankton_waterchem_reduced = results_phytoplankton_waterchem_reduced,
  benthic_waterchem_reduced = results_benthic_waterchem_reduced,
  phytoplankton_biovolume_waterchem_reduced = results_phytoplankton_biovolume_waterchem_reduced,
  zooplankton_biomass_waterchem_w_algal_reduced = results_zooplankton_biomass_waterchem_w_algal_reduced,
  phytoplankton_density_waterchem_reduced = results_phytoplankton_density_waterchem_reduced,
  zooplankton_300_waterchem_reduced = results_zooplankton_300_waterchem_reduced,
  benthic_300_waterchem_reduced = results_benthic_300_waterchem_reduced,
  zooplankton_300_waterchem_w_algal_reduced = results_zooplankton_300_waterchem_w_algal_reduced,
  zooplankton_biomass_300_reduced = results_zooplankton_biomass_300_reduced,
  zooplankton_biomass_300_w_algal_reduced = results_zooplankton_biomass_300_w_algal_reduced,
  zooplankton_biomass_300_waterchem_reduced = results_zooplankton_biomass_300_waterchem_reduced,
  zooplankton_biomass_300_waterchem_w_algal_reduced = results_zooplankton_biomass_300_waterchem_w_algal_reduced,
  all_species_no_outliers = results_all_species_no_outliers,
  zooplankton_w_algal_no_outliers = results_zooplankton_w_algal_no_outliers,
  zooplankton_no_outliers = results_zooplankton_no_outliers,
  phytoplankton_no_outliers = results_phytoplankton_no_outliers,
  benthic_no_outliers = results_benthic_no_outliers,
  phytoplankton_biovolume_no_outliers = results_phytoplankton_biovolume_no_outliers,
  zooplankton_biomass_w_algal_no_outliers = results_zooplankton_biomass_w_algal_no_outliers,
  zooplankton_biomass_no_outliers = results_zooplankton_biomass_no_outliers,
  phytoplankton_density_no_outliers = results_phytoplankton_density_no_outliers,
  zooplankton_300_w_algal_no_outliers = results_zooplankton_300_w_algal_no_outliers,
  zooplankton_300_no_outliers = results_zooplankton_300_no_outliers,
  benthic_300_no_outliers = results_benthic_300_no_outliers,
  zooplankton_waterchem_no_outliers = results_zooplankton_waterchem_no_outliers,
  zooplankton_biomass_waterchem_no_outliers = results_zooplankton_biomass_waterchem_no_outliers,
  all_species_waterchem_no_outliers = results_all_species_waterchem_no_outliers,
  zooplankton_waterchem_w_algal_no_outliers = results_zooplankton_waterchem_w_algal_no_outliers,
  phytoplankton_waterchem_no_outliers = results_phytoplankton_waterchem_no_outliers,
  benthic_waterchem_no_outliers = results_benthic_waterchem_no_outliers,
  phytoplankton_biovolume_waterchem_no_outliers = results_phytoplankton_biovolume_waterchem_no_outliers,
  zooplankton_biomass_waterchem_w_algal_no_outliers = results_zooplankton_biomass_waterchem_w_algal_no_outliers,
  phytoplankton_density_waterchem_no_outliers = results_phytoplankton_density_waterchem_no_outliers,
  zooplankton_300_waterchem_no_outliers = results_zooplankton_300_waterchem_no_outliers,
  benthic_300_waterchem_no_outliers = results_benthic_300_waterchem_no_outliers,
  zooplankton_300_waterchem_w_algal_no_outliers = results_zooplankton_300_waterchem_w_algal_no_outliers,
  zooplankton_biomass_300_no_outliers = results_zooplankton_biomass_300_no_outliers,
  zooplankton_biomass_300_w_algal_no_outliers = results_zooplankton_biomass_300_w_algal_no_outliers,
  zooplankton_biomass_300_waterchem_no_outliers = results_zooplankton_biomass_300_waterchem_no_outliers,
  zooplankton_biomass_300_waterchem_w_algal_no_outliers = results_zooplankton_biomass_300_waterchem_w_algal_no_outliers,
  zooplankton_waterchem_reduced_no_outliers = results_zooplankton_waterchem_reduced_no_outliers,
  zooplankton_reduced_no_outliers = results_zooplankton_reduced_no_outliers,
  #all_species_reduced_no_outliers = results_all_species_reduced_no_outliers,
  zooplankton_w_algal_reduced_no_outliers = results_zooplankton_w_algal_reduced_no_outliers,
  zooplankton_reduced_no_outliers = results_zooplankton_reduced_no_outliers,
  phytoplankton_reduced_no_outliers = results_phytoplankton_reduced_no_outliers,
  benthic_reduced_no_outliers = results_benthic_reduced_no_outliers,
  phytoplankton_biovolume_reduced_no_outliers = results_phytoplankton_biovolume_reduced_no_outliers,
  zooplankton_biomass_w_algal_reduced_no_outliers = results_zooplankton_biomass_w_algal_reduced_no_outliers,
  zooplankton_biomass_reduced_no_outliers = results_zooplankton_biomass_reduced_no_outliers,
  phytoplankton_density_reduced_no_outliers = results_phytoplankton_density_reduced_no_outliers,
  zooplankton_300_w_algal_reduced_no_outliers = results_zooplankton_300_w_algal_reduced_no_outliers,
  zooplankton_300_reduced_no_outliers = results_zooplankton_300_reduced_no_outliers,
  benthic_300_reduced_no_outliers = results_benthic_300_reduced_no_outliers,
  zooplankton_waterchem_reduced_no_outliers = results_zooplankton_waterchem_reduced_no_outliers,
  zooplankton_biomass_waterchem_reduced_no_outliers = results_zooplankton_biomass_waterchem_reduced_no_outliers,
  all_species_waterchem_reduced_no_outliers = results_all_species_waterchem_reduced_no_outliers,
  zooplankton_waterchem_w_algal_reduced_no_outliers = results_zooplankton_waterchem_w_algal_reduced_no_outliers,
  phytoplankton_waterchem_reduced_no_outliers = results_phytoplankton_waterchem_reduced_no_outliers,
  benthic_waterchem_reduced_no_outliers = results_benthic_waterchem_reduced_no_outliers,
  phytoplankton_biovolume_waterchem_reduced_no_outliers = results_phytoplankton_biovolume_waterchem_reduced_no_outliers,
  zooplankton_biomass_waterchem_w_algal_reduced_no_outliers = results_zooplankton_biomass_waterchem_w_algal_reduced_no_outliers,
  phytoplankton_density_waterchem_reduced_no_outliers = results_phytoplankton_density_waterchem_reduced_no_outliers,
  zooplankton_300_waterchem_reduced_no_outliers = results_zooplankton_300_waterchem_reduced_no_outliers,
  benthic_300_waterchem_reduced_no_outliers = results_benthic_300_waterchem_reduced_no_outliers,
  zooplankton_300_waterchem_w_algal_reduced_no_outliers = results_zooplankton_300_waterchem_w_algal_reduced_no_outliers,
  zooplankton_biomass_300_reduced_no_outliers = results_zooplankton_biomass_300_reduced_no_outliers,
  zooplankton_biomass_300_w_algal_reduced_no_outliers = results_zooplankton_biomass_300_w_algal_reduced_no_outliers,
  zooplankton_biomass_300_waterchem_reduced_no_outliers = results_zooplankton_biomass_300_waterchem_reduced_no_outliers,
  zooplankton_biomass_300_waterchem_w_algal_reduced_no_outliers = results_zooplankton_biomass_300_waterchem_w_algal_reduced_no_outliers,
  zooplankton_waterchem_minimal = results_zooplankton_waterchem_minimal,
  zooplankton_biomass_waterchem_minimal = results_zooplankton_biomass_waterchem_minimal,
  all_species_waterchem_minimal = results_all_species_waterchem_minimal,
  zooplankton_waterchem_w_algal_minimal = results_zooplankton_waterchem_w_algal_minimal,
  phytoplankton_waterchem_minimal = results_phytoplankton_waterchem_minimal,
  benthic_waterchem_minimal = results_benthic_waterchem_minimal,
  phytoplankton_biovolume_waterchem_minimal = results_phytoplankton_biovolume_waterchem_minimal,
  zooplankton_biomass_waterchem_w_algal_minimal = results_zooplankton_biomass_waterchem_w_algal_minimal,
  phytoplankton_density_waterchem_minimal = results_phytoplankton_density_waterchem_minimal,
  zooplankton_300_waterchem_minimal = results_zooplankton_300_waterchem_minimal,
  benthic_300_waterchem_minimal = results_benthic_300_waterchem_minimal,
  zooplankton_300_waterchem_w_algal_minimal = results_zooplankton_300_waterchem_w_algal_minimal,
  zooplankton_biomass_300_waterchem_minimal = results_zooplankton_biomass_300_waterchem_minimal,
  zooplankton_biomass_300_waterchem_w_algal_minimal = results_zooplankton_biomass_300_waterchem_w_algal_minimal,
  zooplankton_waterchem_minimal_no_outliers = results_zooplankton_waterchem_minimal_no_outliers,
  zooplankton_biomass_waterchem_minimal_no_outliers = results_zooplankton_biomass_waterchem_minimal_no_outliers,
  all_species_waterchem_minimal_no_outliers = results_all_species_waterchem_minimal_no_outliers,
  zooplankton_waterchem_w_algal_minimal_no_outliers = results_zooplankton_waterchem_w_algal_minimal_no_outliers,
  phytoplankton_waterchem_minimal_no_outliers = results_phytoplankton_waterchem_minimal_no_outliers,
  benthic_waterchem_minimal_no_outliers = results_benthic_waterchem_minimal_no_outliers,
  phytoplankton_biovolume_waterchem_minimal_no_outliers = results_phytoplankton_biovolume_waterchem_minimal_no_outliers,
  zooplankton_biomass_waterchem_w_algal_minimal_no_outliers = results_zooplankton_biomass_waterchem_w_algal_minimal_no_outliers,
  phytoplankton_density_waterchem_minimal_no_outliers = results_phytoplankton_density_waterchem_minimal_no_outliers,
  zooplankton_300_waterchem_minimal_no_outliers = results_zooplankton_300_waterchem_minimal_no_outliers,
  benthic_300_waterchem_minimal_no_outliers = results_benthic_300_waterchem_minimal_no_outliers,
  zooplankton_300_waterchem_w_algal_minimal_no_outliers = results_zooplankton_300_waterchem_w_algal_minimal_no_outliers,
  zooplankton_biomass_300_waterchem_minimal_no_outliers = results_zooplankton_biomass_300_waterchem_minimal_no_outliers,
  zooplankton_biomass_300_waterchem_w_algal_minimal_no_outliers = results_zooplankton_biomass_300_waterchem_w_algal_minimal_no_outliers
)

save(all_results, file = "all_results.RData")

sink()

