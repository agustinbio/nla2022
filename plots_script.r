library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Cargar los datos de inercia
df <- read.csv("enhanced_variance_decomposition.csv")
  
# función de traducción
translate_spp <- function(spp) {
  spp_translations <- list(
    "all_species" = "TODAS ABUND",
    "phytoplankton" = "FITO ABUND",
    "zooplankton" = "ZOO ABUND",
    "benthic" = "BENTOS ABUND",
    "phytoplankton_density" = "FITO DENS",
    "phytoplankton_biovolume" = "FITO BIOVOL",
    "zooplankton_300" = "ZOO ABUND 300",
    "benthic_300" = "BENTOS ABUND 300",
    "zooplankton_biomass" = "ZOO BIOM",
    "zooplankton_biomass_300" = "ZOO BIOM 300"
  )
  
  sapply(spp, function(x) {
    if (tolower(x) %in% names(spp_translations)) {
      spp_translations[[tolower(x)]]
    } else {
      x
    }
  })
}

translate_env_base <- function(env) {
  env_translations <- c(
    "full-environment" = "COMPL",
    "full-environment-w-algal" = "COMPL+ALG",
    "waterchem" = "QUIM",
    "waterchem-w-algal" = "QUIM+ALG"
  )
  
  # Vectorized translation
  result <- env_translations[tolower(env)]
  result[is.na(result)] <- env[is.na(result)]
  as.character(result)
}

translate_env_type <- function(env_type) {
  case_when(
    env_type == "standard" ~ "estándar",
    env_type == "reduced" ~ "VIF",
    env_type == "minimal" ~ "STEP",
    env_type == "standard-no-outliers" ~ "estándar-OA",
    env_type == "reduced-no-outliers" ~ "VIF-OA",
    env_type == "minimal-no-outliers" ~ "STEP-OA",
    TRUE ~ env_type
  )
}

translate_variance <- function(var) {
  case_when(
    var == "CONDITIONED" ~ "Condicionada",
    var == "FULL_CONSTRAINED" ~ "Restringida (Todas)",
    var == "SEL_CONSTRAINED" ~ "Restringida (Selección)",
    var == "SEL_CONSTRAINED_GLM" ~ "Restringida (Selección-GLM)",
    TRUE ~ var
  )
}

# Funciones de preparación de datos
get_dataset_type <- function(env_label) {
  if (grepl("-minimal-no-outliers$", env_label)) {
    "minimal-no-outliers"
  } else if (grepl("-minimal$", env_label)) {
    "minimal"
  } else if (grepl("-reduced-no-outliers$", env_label)) {
    "reduced-no-outliers"
  } else if (grepl("-no-outliers$", env_label)) {
    "standard-no-outliers"
  } else if (grepl("-reduced$", env_label)) {
    "reduced"
  } else {
    "standard"
  }
}

extract_base_env <- function(env_label) {
  env_label <- gsub("-reduced-no-outliers$", "", env_label)
  env_label <- gsub("-no-outliers$", "", env_label)
  env_label <- gsub("-reduced$", "", env_label)
  env_label <- gsub("-minimal-no-outliers$", "", env_label)  
  env_label <- gsub("-minimal$", "", env_label)
  env_label
}

# Traduce el dataframe
df <- df %>%
  mutate(
    SPP_TRANS = translate_spp(SPP),
    ENV_TYPE = sapply(ENV, get_dataset_type) %>% translate_env_type(),
    BASE_ENV = sapply(ENV, extract_base_env) %>% translate_env_base()
  )

# Convierte al formato largo
df_long <- df %>%
  pivot_longer(
    cols = c(CONDITIONED, FULL_CONSTRAINED, SEL_CONSTRAINED),
    names_to = "Variance_Component",
    values_to = "Porcentaje"
  ) %>%
  mutate(
    Variance_Component = translate_variance(Variance_Component)
  )


# Crea el formato largo para inercia restringida
df_long_constrained <- df %>%
  select(ID, SPP, ENV, BASE_ENV, ENV_TYPE, NUM_SEL, FULL_CONSTRAINED, SEL_CONSTRAINED) %>%
  pivot_longer(
    cols = c(FULL_CONSTRAINED, SEL_CONSTRAINED),
    names_to = "Variance_Type",
    values_to = "Constrained_Variance"
  ) %>%
  mutate(
    Variance_Type = translate_variance(Variance_Type)
  )

plot_df <- df %>%
  mutate(
    ENV_TYPE = factor(
      ENV_TYPE,
      levels = c("estándar", "estándar-OA", "VIF", "VIF-OA", "STEP", "STEP-OA")
    )
  )


# ==============================================================================
# BLOQUE 1: Gráficos de Porcentajes Absolutos (Escala Unificada y Combinados)
# ==============================================================================

# 1. Calcular el rango global para las 4 métricas
vals_const_std <- plot_df$SEL_CONSTRAINED
vals_const_glm <- plot_df$SEL_CONSTRAINED_GLM
vals_expl_std  <- plot_df$SEL_CONSTRAINED + plot_df$CONDITIONED
vals_expl_glm  <- plot_df$SEL_CONSTRAINED_GLM + plot_df$CONDITIONED

all_values_abs <- c(vals_const_std, vals_const_glm, vals_expl_std, vals_expl_glm)
rng_global_abs <- range(all_values_abs, na.rm = TRUE)
expand_amt_abs <- 0.05 * diff(rng_global_abs)
ylim_global_abs <- c(rng_global_abs[1] - expand_amt_abs, rng_global_abs[2] + expand_amt_abs)

# 2. Generar los gráficos (Con títulos centrados y en negrita)

# A. Inercia Restringida
p2_species <- ggplot(plot_df, aes(x = ENV_TYPE, y = SEL_CONSTRAINED, color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(df$SPP_TRANS))) +
  labs(title = "Flujo de trabajo", x = NULL, y = "Inercia restringida (%)",
       color = "Variables ambientales", shape = "Especies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) + # Estilo de título agregado
  coord_cartesian(ylim = ylim_global_abs)

p2_species_glm <- ggplot(plot_df, aes(x = ENV_TYPE, y = SEL_CONSTRAINED_GLM, color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(df$SPP_TRANS))) +
  labs(title = "GLM", x = NULL, y = NULL,
       color = "Variables ambientales", shape = "Especies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) + # Estilo de título agregado
  coord_cartesian(ylim = ylim_global_abs)

# B. Inercia Explicada (Restringida + Condicionada)
p2_species_w_cond <- ggplot(plot_df, aes(x = ENV_TYPE, y = (SEL_CONSTRAINED + CONDITIONED), color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(df$SPP_TRANS))) +
  labs(title = "Flujo de trabajo", x = "Tipo de procesamiento", y = "Inercia restringida + condicionada (%)",
       color = "Variables ambientales", shape = "Especies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) + # Estilo de título agregado
  coord_cartesian(ylim = ylim_global_abs)

p2_species_glm_w_cond <- ggplot(plot_df, aes(x = ENV_TYPE, y = (SEL_CONSTRAINED_GLM + CONDITIONED), color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(df$SPP_TRANS))) +
  labs(title = "GLM", x = "Tipo de procesamiento", y = NULL,
       color = "Variables ambientales", shape = "Especies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) + # Estilo de título agregado
  coord_cartesian(ylim = ylim_global_abs)

# 3. Combinar y Guardar

# Restringida: SIN leyenda
combined_restr <- p2_species + p2_species_glm + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "none") # Quita la leyenda

ggsave("restringida-combinado.svg", combined_restr, width = 16, height = 7)

# Explicada: Leyenda ABAJO
combined_expl <- p2_species_w_cond + p2_species_glm_w_cond + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") # Leyenda abajo

ggsave("explicada-combinado.svg", combined_expl, width = 16, height = 8) # Height ajustado un poco para la leyenda


# ==============================================================================
# BLOQUE 2: Gráficos Normalizados por N (Escala Unificada y Combinados)
# ==============================================================================

# 1. Calcular el rango global para las versiones normalizadas
norm_const_std <- plot_df$SEL_CONSTRAINED / plot_df$NUM_SEL
norm_const_glm <- plot_df$SEL_CONSTRAINED_GLM / plot_df$NUM_SEL_GLM
norm_expl_std  <- (plot_df$SEL_CONSTRAINED + plot_df$CONDITIONED) / (plot_df$NUM_SEL + plot_df$NUM_MEM)
norm_expl_glm  <- (plot_df$SEL_CONSTRAINED_GLM + plot_df$CONDITIONED) / (plot_df$NUM_SEL_GLM + plot_df$NUM_MEM)

all_values_norm <- c(norm_const_std, norm_const_glm, norm_expl_std, norm_expl_glm)
rng_global_norm <- range(all_values_norm, na.rm = TRUE)
expand_amt_norm <- 0.05 * diff(rng_global_norm)
ylim_global_norm <- c(rng_global_norm[1] - expand_amt_norm, rng_global_norm[2] + expand_amt_norm)

# 2. Generar gráficos normalizados con títulos centrados y en negrita

p2_species_n <- ggplot(plot_df, aes(x = ENV_TYPE, y = (SEL_CONSTRAINED/NUM_SEL), color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(df$SPP_TRANS))) +
  labs(title = "Flujo de trabajo", x = NULL, y = "Inercia restringida (%) / N",
       color = "Variables ambientales", shape = "Especies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(ylim = ylim_global_norm)

p2_species_glm_n <- ggplot(plot_df, aes(x = ENV_TYPE, y = (SEL_CONSTRAINED_GLM/NUM_SEL_GLM), color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(df$SPP_TRANS))) +
  labs(title = "GLM", x = NULL, y = NULL,
       color = "Variables ambientales", shape = "Especies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(ylim = ylim_global_norm)

p2_species_w_cond_n <- ggplot(plot_df, aes(x = ENV_TYPE, y = ((SEL_CONSTRAINED + CONDITIONED)/(NUM_SEL+NUM_MEM)), color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(df$SPP_TRANS))) +
  labs(title = "Flujo de trabajo", x = "Tipo de procesamiento", y = "(Inercia restringida + condicionada) (%) / N",
       color = "Variables ambientales", shape = "Especies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(ylim = ylim_global_norm)

p2_species_glm_w_cond_n <- ggplot(plot_df, aes(x = ENV_TYPE, y = ((SEL_CONSTRAINED_GLM + CONDITIONED)/(NUM_SEL_GLM+NUM_MEM)), color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(df$SPP_TRANS))) +
  labs(title = "GLM", x = "Tipo de procesamiento", y = NULL,
       color = "Variables ambientales", shape = "Especies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(ylim = ylim_global_norm)

# 3. Combinar y Guardar (Normalizados)

# Restringida / N: SIN leyenda
combined_restr_n <- p2_species_n + p2_species_glm_n + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "none")

ggsave("n-restringida-combinado.svg", combined_restr_n, width = 16, height = 7)

# Explicada / N: Leyenda ABAJO
combined_expl_n <- p2_species_w_cond_n + p2_species_glm_w_cond_n + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("n-explicada-combinado.svg", combined_expl_n, width = 16, height = 8)



# Carga las métricas de los modelos
metrics_df <- read.csv("comprehensive_model_metrics.csv")

# Aplica las traducciones
metrics_df <- metrics_df %>%
  mutate(
    SPP_TRANS = translate_spp(SPP),
    ENV_TYPE = sapply(ENV, get_dataset_type) %>% translate_env_type(),
    BASE_ENV = sapply(ENV, extract_base_env) %>% translate_env_base()
  )

create_rmse_plots <- function() {
  # calcula el rango intercuartílico para las proporciones de RMSE del GAM
  gam_rmse <- metrics_df$gam_test_rmse_ratio
  q1 <- quantile(gam_rmse, 0.25, na.rm = TRUE)
  q3 <- quantile(gam_rmse, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  upper_bound <- q3 + 1.5 * iqr
  
  # gráfico de GLM
  p_glmnet <- ggplot(metrics_df, aes(x = ENV_TYPE, y = glmnet_test_rmse_ratio, 
                                     color = BASE_ENV, shape = SPP_TRANS)) +
    geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
    scale_shape_manual(values = 1:length(unique(metrics_df$SPP_TRANS))) +
    labs(title = "GLM",
         x = "Tipo de procesamiento",
         y = "RMSE/Abundancia total",
         color = "Variables ambientales",
         shape = "Especies") +
    ylim(0, max(metrics_df$glmnet_test_rmse_ratio, na.rm = TRUE)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # gráfico de GAM
  p_gam_log <- ggplot(metrics_df, aes(x = ENV_TYPE, y = gam_test_rmse_ratio, 
                                      color = BASE_ENV, shape = SPP_TRANS)) +
    geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8, na.rm = TRUE) +
    scale_shape_manual(values = 1:length(unique(metrics_df$SPP_TRANS))) +
    scale_y_continuous(trans='log10') +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.1) +
    labs(title = "GAM",
         x = "Tipo de procesamiento",
         y = "RMSE/Abundancia total (escala logarítmica)",
         color = "Variables ambientales",
         shape = "Especies") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Combina los gráficos
  combined <- p_glmnet + p_gam_log +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave("rmse-ratio.svg", combined, width = 16, height = 8)
}

# Crea los gráficos de R²
create_r2_plots <- function() {
  
  all_r2_values <- c(metrics_df$glmnet_test_r2_weighted, 
                     metrics_df$gam_test_r2_weighted)
  rng <- range(all_r2_values, na.rm = TRUE)
  expand_amount <- 0.05 * diff(rng)
  ylim_unified <- c(rng[1] - expand_amount, rng[2] + expand_amount)
  
  p_glmnet <- ggplot(metrics_df, aes(x = ENV_TYPE, y = glmnet_test_r2_weighted, 
                                     color = BASE_ENV, shape = SPP_TRANS)) +
    geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
    scale_shape_manual(values = 1:length(unique(metrics_df$SPP_TRANS))) +
    labs(title = "GLM",
         x = NULL,
         y = "R² Ponderado",
         color = "Variables ambientales",
         shape = "Especies") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    coord_cartesian(ylim = ylim_unified)  # Apply unified scale
  
  p_gam <- ggplot(metrics_df, aes(x = ENV_TYPE, y = gam_test_r2_weighted, 
                                  color = BASE_ENV, shape = SPP_TRANS)) +
    geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
    scale_shape_manual(values = 1:length(unique(metrics_df$SPP_TRANS))) +
    labs(title = "GAM",
         x = NULL,
         y = "R² Ponderado",
         color = "Variables ambientales",
         shape = "Especies") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    coord_cartesian(ylim = ylim_unified)  # Apply same unified scale
  
  combined <- p_glmnet + p_gam +
    plot_layout(guides = "collect") &
    theme(legend.position = "none")
  
  ggsave("r2-ponderado.svg", combined, width = 16, height = 8)

}

# Genera los gráficos
rmse_plots <- create_rmse_plots()
r2_plots <- create_r2_plots()

p_rmse_comp <- ggplot(metrics_df, aes(x = ENV_TYPE, y = test_rmse_comp, 
                                      color = BASE_ENV, shape = SPP_TRANS)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +
  scale_shape_manual(values = 1:length(unique(metrics_df$SPP_TRANS))) +
  scale_y_continuous(trans='log10') +
 # annotation_logticks(sides = "l") +  # Agrega los ticks de la escala logarítmica
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.1) +
  labs(#title = "Comparación de RMSE entre GAM y GLM",
       #subtitle = "Y = RMSE (GAM) / RMSE (GLM) (Escala logarítmica)",
       x = "Tipo de procesamiento",
       y = "Proporción de RMSE (GAM/GLM)",
       color = "Variables ambientales",
       shape = "Especies") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 0, size = 8),  # fuente más chica para las etiquetas del eje Y
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()  # Saca las grillas menos importantes para mayor claridad
  ) +
  guides(shape = guide_legend(ncol = 2))

#Guarda el gráfico con más altura
ggsave("rmse-comparacion.svg", p_rmse_comp, width = 10, height = 9)

