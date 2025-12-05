library(ggplot2)
library(patchwork)
library(dplyr)
library(purrr)

# Cargar los resultados
load("all_results.RData")
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
    "full-environment-minimal" = "COMPL (STEP)",
    "full-environment-w-algal-minimal" = "COMPL+ALG (STEP)",
    "waterchem-minimal" = "QUIM (STEP)",
    "waterchem-w-algal-minimal" = "QUIM+ALG (STEP)",
    "full-environment-no-outliers" = "COMPL (-OA)",
    "full-environment-w-algal-no-outliers" = "COMPL+ALG (-OA)",
    "waterchem-no-outliers" = "QUIM (-OA)",
    "waterchem-w-algal-no-outliers" = "QUIM+ALG (-OA)",
    "full-environment-reduced-no-outliers" = "COMPL (VIF-OA)",
    "full-environment-w-algal-reduced-no-outliers" = "COMPL+ALG (VIF-OA)",
    "waterchem-reduced-no-outliers" = "QUIM (VIF-OA)",
    "waterchem-w-algal-reduced-no-outliers" = "QUIM+ALG (VIF-OA)",
    "full-environment-minimal-no-outliers" = "COMPL (STEP-OA)",
    "full-environment-w-algal-minimal-no-outliers" = "COMPL+ALG (STEP-OA)",
    "waterchem-minimal-no-outliers" = "QUIM (STEP-OA)",
    "waterchem-w-algal-minimal-no-outliers" = "QUIM+ALG (STEP-OA)"
  )
  
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
  
  translated_env <- ifelse(tolower(env) %in% names(env_translations),
                           env_translations[[tolower(env)]],
                           env)
  
  translated_spp <- ifelse(tolower(spp) %in% names(spp_translations),
                           spp_translations[[tolower(spp)]],
                           spp)
  
  paste0(translated_spp, " (", translated_env, ")")
}


# Función auxiliar para categorizar los datasets
get_dataset_type <- function(env_label) {
  if (grepl("-reduced-no-outliers$", env_label)) {
    "reduced-no-outliers"
  } else if (grepl("-no-outliers$", env_label)) {
    "standard-no-outliers"
  } else if (grepl("-reduced$", env_label)) {
    "reduced"
  } else {
    "standard"
  }
}

# Categorizar todos los resultados en una lista de categorías
categorized_results <- list(
  standard = list(),
  reduced = list(),
  `standard-no-outliers` = list(),
  `reduced-no-outliers` = list()
)

# Poblar las categorías
for (res in all_results) {
  type <- get_dataset_type(res$env)
  categorized_results[[type]] <- c(categorized_results[[type]], list(res))
}

# Tabla con estadísticas de R² (train)
r2_stats_df <- map_dfr(all_results, function(res) {
  r2_vals <- res$glmnet$train_metrics$r2
  data.frame(
    ID = generar_nombre(res$spp, res$env),
    SPP = res$spp,
    ENV = res$env,
    R2_MEAN = mean(r2_vals, na.rm = TRUE),
    R2_MEDIAN = median(r2_vals, na.rm = TRUE),
    R2_SD = sd(r2_vals, na.rm = TRUE)
  )
})

write.csv(r2_stats_df, "train-r2_statistics.csv", row.names = FALSE)

# Tabla con estadísticas de R² (test)
r2_stats_df <- map_dfr(all_results, function(res) {
  r2_vals <- res$glmnet$test_metrics$r2
  data.frame(
    ID = generar_nombre(res$spp, res$env),
    SPP = res$spp,
    ENV = res$env,
    R2_MEAN = mean(r2_vals, na.rm = TRUE),
    R2_MEDIAN = median(r2_vals, na.rm = TRUE),
    R2_SD = sd(r2_vals, na.rm = TRUE)
  )
})

write.csv(r2_stats_df, "test-r2_statistics.csv", row.names = FALSE)



# Tabla con la descomposición de la varianza
var_decomp_df <- map_dfr(all_results, function(res) {
  # Calcular porcentajes
  total_full <- res$var_full$total
  cond_full <- ifelse(is.null(res$var_full$cond), 0, res$var_full$cond)
  constr_full <- res$var_full$constr
  
  total_sel <- res$var_selected$total
  constr_sel <- res$var_selected$constr
  
  total_sel_glmnet <- res$var_selected_glmnet$total
  constr_sel_glmnet <- res$var_selected_glmnet$constr
  
  data.frame(
    ID = generar_nombre(res$spp, res$env),
    SPP = res$spp,
    ENV = res$env,
    CONDITIONED = (cond_full / total_full) * 100,
    FULL_CONSTRAINED = (constr_full / total_full) * 100,
    SEL_CONSTRAINED = (constr_sel / total_sel) * 100,
    SEL_CONSTRAINED_GLMNET = (constr_sel_glmnet / total_sel_glmnet) * 100
  )
})

write.csv(var_decomp_df, "variance_decomposition.csv", row.names = FALSE)

# ========================================================================
# función SAFE_EXTRACT para obtener las métricas y las proporciones
# ========================================================================

safe_extract <- function(metrics, abund = NULL) {
  result <- list(
    r2_mean = NA,
    r2_weighted = NA,
    total_rmse = NA,
    total_mae = NA,
    total_abund = NA,
    rmse_ratio = NA,
    mae_ratio = NA
  )
  
  # Calcular la abundancia total si está disponible
  if (!is.null(abund) && all(!is.na(abund))) {
    result$total_abund <- sum(abund, na.rm = TRUE)
  }
  
  # Procesar las métricas si las hay
  if (!is.null(metrics)) {
    # cálculo de R²
    if (!is.null(metrics$r2) && length(metrics$r2) > 0) {
      result$r2_mean <- mean(metrics$r2, na.rm = TRUE)
      if (!is.null(abund) && length(abund) == length(metrics$r2)) {
        result$r2_weighted <- weighted.mean(metrics$r2, w = abund, na.rm = TRUE)
      }
    }
    
    # cálculo de RMSE
    if (!is.null(metrics$rmse) && length(metrics$rmse) > 0) {
      result$total_rmse <- sum(metrics$rmse, na.rm = TRUE)
      if (!is.na(result$total_abund) && result$total_abund > 0) {
        result$rmse_ratio <- result$total_rmse / result$total_abund
      }
    }
    
    # cálculo de MAE
    if (!is.null(metrics$mae) && length(metrics$mae) > 0) {
      result$total_mae <- sum(metrics$mae, na.rm = TRUE)
      if (!is.na(result$total_abund) && result$total_abund > 0) {
        result$mae_ratio <- result$total_mae / result$total_abund
      }
    }
  }
  
  return(result)
}

# ========================================================================
# Tablas de métricas con proporciones y comparaciones
# ========================================================================

metrics_df <- map_dfr(all_results, function(res) {
  # Extraer las métricas
  glmnet_train <- safe_extract(res$glmnet$train_metrics, res$spp_abund$train)
  glmnet_test <- safe_extract(res$glmnet$test_metrics, res$spp_abund$test)
  
  gam_train <- if (!is.null(res$gam)) {
    safe_extract(res$gam$train_metrics, res$spp_abund$train)
  } else {
    list(r2_mean = NA, r2_weighted = NA, total_rmse = NA, total_mae = NA,
         total_abund = NA, rmse_ratio = NA, mae_ratio = NA)
  }
  
  gam_test <- if (!is.null(res$gam)) {
    safe_extract(res$gam$test_metrics, res$spp_abund$test)
  } else {
    list(r2_mean = NA, r2_weighted = NA, total_rmse = NA, total_mae = NA,
         total_abund = NA, rmse_ratio = NA, mae_ratio = NA)
  }
  
  # Calcular las comparaciones (GAM vs. GLM)
  train_rmse_comp <- ifelse(
    !is.na(gam_train$rmse_ratio) & !is.na(glmnet_train$rmse_ratio) & glmnet_train$rmse_ratio > 0,
    gam_train$rmse_ratio / glmnet_train$rmse_ratio,
    NA
  )
  
  train_mae_comp <- ifelse(
    !is.na(gam_train$mae_ratio) & !is.na(glmnet_train$mae_ratio) & glmnet_train$mae_ratio > 0,
    gam_train$mae_ratio / glmnet_train$mae_ratio,
    NA
  )
  
  test_rmse_comp <- ifelse(
    !is.na(gam_test$rmse_ratio) & !is.na(glmnet_test$rmse_ratio) & glmnet_test$rmse_ratio > 0,
    gam_test$rmse_ratio / glmnet_test$rmse_ratio,
    NA
  )
  
  test_mae_comp <- ifelse(
    !is.na(gam_test$mae_ratio) & !is.na(glmnet_test$mae_ratio) & glmnet_test$mae_ratio > 0,
    gam_test$mae_ratio / glmnet_test$mae_ratio,
    NA
  )
  
  # Construir la fila-resultado
  data.frame(
    ID = generar_nombre(res$spp, res$env),
    SPP = res$spp,
    ENV = res$env,
    
    # métricas glmnet train
    glmnet_train_r2_mean = glmnet_train$r2_mean,
    glmnet_train_r2_weighted = glmnet_train$r2_weighted,
    glmnet_train_total_rmse = glmnet_train$total_rmse,
    glmnet_train_total_mae = glmnet_train$total_mae,
    glmnet_train_rmse_ratio = glmnet_train$rmse_ratio,
    glmnet_train_mae_ratio = glmnet_train$mae_ratio,
    
    # métricas glmnet test
    glmnet_test_r2_mean = glmnet_test$r2_mean,
    glmnet_test_r2_weighted = glmnet_test$r2_weighted,
    glmnet_test_total_rmse = glmnet_test$total_rmse,
    glmnet_test_total_mae = glmnet_test$total_mae,
    glmnet_test_rmse_ratio = glmnet_test$rmse_ratio,
    glmnet_test_mae_ratio = glmnet_test$mae_ratio,
    
    # métricas gam train 
    gam_train_r2_mean = gam_train$r2_mean,
    gam_train_r2_weighted = gam_train$r2_weighted,
    gam_train_total_rmse = gam_train$total_rmse,
    gam_train_total_mae = gam_train$total_mae,
    gam_train_rmse_ratio = gam_train$rmse_ratio,
    gam_train_mae_ratio = gam_train$mae_ratio,
    
    # métricas gam test 
    gam_test_r2_mean = gam_test$r2_mean,
    gam_test_r2_weighted = gam_test$r2_weighted,
    gam_test_total_rmse = gam_test$total_rmse,
    gam_test_total_mae = gam_test$total_mae,
    gam_test_rmse_ratio = gam_test$rmse_ratio,
    gam_test_mae_ratio = gam_test$mae_ratio,
    
    # Métricas de comparación
    train_rmse_comp = train_rmse_comp,
    train_mae_comp = train_mae_comp,
    test_rmse_comp = test_rmse_comp,
    test_mae_comp = test_mae_comp
  )
})

write.csv(metrics_df, "comprehensive_model_metrics.csv", row.names = FALSE)

# ========================================================================
# Tabla de descomposición de la varianza ampliada con metadata
# ========================================================================

var_decomp_df <- map_dfr(all_results, function(res) {
  # Inicializar con valores NA por defecto
  result <- list(
    CONDITIONED = NA,
    FULL_CONSTRAINED = NA,
    SEL_CONSTRAINED = NA,
    SEL_CONSTRAINED_GLM = NA,
    NUM_SEL = NA,
    NUM_SEL_GLM = NA,
    NUM_SPP = NA,
    NUM_ENV = NA,
    NUM_MEM = NA,
    TRAIN_OBS = NA,
    TEST_OBS = NA,
    TOTAL_OBS = NA
  )
  
  # Calcular los porcentajes de varianza si la información está disponible
  if (!is.null(res$var_full) && !is.null(res$var_selected)) {
    total_full <- res$var_full$total
    cond_full <- if (!is.null(res$var_full$cond)) res$var_full$cond else 0
    constr_full <- res$var_full$constr
    
    total_sel <- res$var_selected$total
    constr_sel <- res$var_selected$constr
    
    total_sel_glmnet <- res$var_selected_glmnet$total
    constr_sel_glmnet <- res$var_selected_glmnet$constr
    
    if (!is.na(total_full) && total_full > 0) {
      result$CONDITIONED <- (cond_full / total_full) * 100
      result$FULL_CONSTRAINED <- (constr_full / total_full) * 100
    }
    
    if (!is.na(total_sel) && total_sel > 0) {
      result$SEL_CONSTRAINED <- (constr_sel / total_sel) * 100
    }
    
    if (!is.na(total_sel_glmnet) && total_sel_glmnet > 0) {
      result$SEL_CONSTRAINED_GLM <- (constr_sel_glmnet / total_sel_glmnet) * 100
    }
  }
  
  # Extraer la metadata del componente n de resultados
  tryCatch({
    # Extraer la cantidad de especies
    result$NUM_SPP <- res$n$vars$spp
    
    # Extraer la cantidad de variables ambientales o MEMs
    result$NUM_ENV <- res$n$vars$env
    result$NUM_MEM <- res$n$vars$mem
    
    result$NUM_SEL <- nrow(res$loadings_selected)
    result$NUM_SEL_GLM <- nrow(res$loadings_selected_glmnet)
    
    # Extraer la cantidad de observaciones
    result$TRAIN_OBS <- res$n$obs$train
    result$TEST_OBS <- res$n$obs$test
    result$TOTAL_OBS <- res$n$obs$total
  }, error = function(e) {
    # La metadata tendrá el valor NA si la extracción fallase
    message("Error extracting metadata for: ", res$spp, " - ", res$env)
    message("Error message: ", e$message)
  })
  
  # devolver como una fila de dataframe
  data.frame(
    ID = generar_nombre(res$spp, res$env),
    SPP = res$spp,
    ENV = res$env,
    result
  )
})

write.csv(var_decomp_df, "enhanced_variance_decomposition.csv", row.names = FALSE)
