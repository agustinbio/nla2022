library(dplyr)
library(tidyr)

# Cargar los resultados guardados
load("all_results.RData")

# Traducción de los nombres
generar_nombre <- function(spp, env) {
  env_translations <- list(
    "full-environment" = "COMPL",
    "full-environment-w-algal" = "COMPL+ALG",
    "waterchem" = "QUIM",
    "waterchem-w-algal" = "QUIM+ALG",
    "full-environment-reduced" = "COMPL (VIF)",
    "full-environment-w-algal-reduced" = "COMPL+ALG (VIF)",
    "waterchem-reduced" = "QUIM (VIF)",
    "waterchem-w-algal-reduced" = "QUIM+ALG (VIF)",
    "full-environment-no-outliers" = "COMPL (-OA)",
    "full-environment-w-algal-no-outliers" = "COMPL+ALG (-OA)",
    "waterchem-no-outliers" = "QUIM (-OA)",
    "waterchem-w-algal-no-outliers" = "QUIM+ALG (-OA)",
    "full-environment-reduced-no-outliers" = "COMPL (VIF-OA)",
    "full-environment-w-algal-reduced-no-outliers" = "COMPL+ALG (VIF-OA)",
    "waterchem-reduced-no-outliers" = "QUIM (VIF-OA)",
    "waterchem-w-algal-reduced-no-outliers" = "QUIM+ALG (VIF-OA)"
  )
  
  spp_translations <- list(
    "all_species" = "TODAS ABUND",
    "phytoplankton" = "FITO ABUND",
    "zooplankton" = "ZOO ABUND",
    "benthic" = "BENTOS ABUND",
    "phytoplankton_density" = "FITO DENS",
    "phytoplankton_biovolume" = "FITO BIOVOL",
    "zooplankton_300" = "ZOO 300",
    "benthic_300" = "BENTOS 300",
    "zooplankton_biomass" = "ZOO BIOM",
    "zooplankton_biomass_300" = "ZOO BIOM 300"
  )
  
  translated_env <- ifelse(tolower(env) %in% names(env_translations),
                           env_translations[[tolower(env)]],
                           env)
  
  translated_spp <- ifelse(tolower(spp) %in% names(spp_translations),
                           spp_translations[[tolower(spp)]],
                           spp)
  
  paste0(translated_spp, " (", translated_env, ")")
}

# Inicializar las variables
all_outliers <- list()
overlap_summary <- list()

# procesar cada resultado en all_results
for (i in seq_along(all_results)) {
  res <- all_results[[i]]
  
  # obetner el nombre traducido
  analysis_name <- generar_nombre(res$spp, res$env)
  
  # obtener los outliers de isolation forest
  iso_species <- res$outliers$isolation$species$UID
  iso_env <- res$outliers$isolation$env$UID
  
  # Obtener los outliers de CCA
  cca_species <- res$outliers$cca$species$UID
  cca_env <- res$outliers$cca$env$UID
  
  # Combinar los UIDs de los outliers para este análisis
  all_uids <- unique(c(iso_species, iso_env, cca_species, cca_env))
  
  # crear un data frame con el resumen
  analysis_df <- data.frame(
    UID = all_uids,
    Analysis = analysis_name,
    ISO_Species = as.integer(all_uids %in% iso_species),
    ISO_Env = as.integer(all_uids %in% iso_env),
    CCA_Species = as.integer(all_uids %in% cca_species),
    CCA_Env = as.integer(all_uids %in% cca_env)
  )
  
  # Calcular las superposiciones
  analysis_df$ISO_Any <- as.integer(analysis_df$ISO_Species | analysis_df$ISO_Env)
  analysis_df$CCA_Any <- as.integer(analysis_df$CCA_Species | analysis_df$CCA_Env)
  analysis_df$Any_Overlap <- as.integer(analysis_df$ISO_Any & analysis_df$CCA_Any)
  
  # Guarda los resultados
  all_outliers[[analysis_name]] <- analysis_df
  
  # Calcula la estadística de resumen (conteos)
  overlap_summary[[analysis_name]] <- data.frame(
    Analysis = analysis_name,
    ISO_Count = length(unique(c(iso_species, iso_env))),
    CCA_Count = length(unique(c(cca_species, cca_env))),
    Overlap_Count = sum(analysis_df$Any_Overlap),
    Overlap_UIDs = paste(all_uids[analysis_df$Any_Overlap == 1], collapse = ", ")
  )
}

# combina todos los análisis
combined_outliers <- bind_rows(all_outliers)

# Crea las tablas de frecuencias
outlier_frequency <- combined_outliers %>%
  group_by(UID) %>%
  summarise(
    Total_Analyses = n(),
    ISO_Count = sum(ISO_Any),
    CCA_Count = sum(CCA_Any),
    Analyses = paste(unique(Analysis), collapse = ", ")
  ) %>%
  arrange(desc(Total_Analyses))

# Crea un resumen global
global_summary <- data.frame(
  Metric = c("Total de outliers únicos (sitios)", 
             "Sitios detectados por Isolation Forest",
             "Sitios detectados por CCA",
             "Sitios detectados por ambos métodos"),
  Count = c(
    n_distinct(combined_outliers$UID),
    sum(combined_outliers$ISO_Any > 0),
    sum(combined_outliers$CCA_Any > 0),
    sum(combined_outliers$ISO_Any & combined_outliers$CCA_Any)
  )
)

# Guarda las salidas
write.csv(outlier_frequency, "outlier_frequency_summary.csv", row.names = FALSE)
write.csv(bind_rows(overlap_summary), "analysis_overlap_summary.csv", row.names = FALSE)
write.csv(global_summary, "global_outlier_summary.csv", row.names = FALSE)

