# Cargar librerías necesarias
library(dplyr)
library(missMDA)
library(vegan)
library(psych)
library(ggvegan)
library(ggplot2)
library(ggbiplot)
library(corrplot)
library(purrr)
library(patchwork)
library(GGally)
library(ggrepel)
library(reshape2)
library(adespatial)
library(spdep)
library(mcga)
library(digest)
library(VIM)
library(glmnet)
library(mgcv)
library(isotree)


OUTLIERS_UID <- c(2020495,
                 2020502) #, 
#si se quitan los dos outliers que figuran abajo disminuye el rendimiento de los análisis.
                 #2021400,
                #2021692)

#transformación de las variables porcentuales
PCT_TRANSFORM <- "logit" #Opciones: "asinsqrt", "logit", "none"

EXCLUDE_OUTLIERS_KEYWORD <- "no-outliers"

VIF_ANALYSIS_KEYWORD <- "reduced"
VIF_THRESHOLD <- 10

# Cantidad máxima de pasos para stepwise selection de los datasets "minimal"
MAX_STEPS <- 200

CCA_OUTLIER_QUANTILE <- 0.999

KNN_IMPUTATION_N <- 5

DETECT_ENV_OUTLIERS <- TRUE 
DETECT_SPECIES_OUTLIERS <- TRUE
OUTLIER_QUANTILE <- 0.99
ISOLATION_FOREST_NTREES_ENV <- 500
ISOLATION_FOREST_NTREES_SPECIES <- 500

TEST_PROP <- 0.1

NO_LOG_TRANSFORM_GLMNET <- FALSE
GLMNET_FAMILY <- "mgaussian"
GLMNET_LINK <- NULL #"log"

NO_LOG_TRANSFORM_GAM <- FALSE

PERFORM_GAM <- TRUE
GAM_FAMILY <- "gaussian"  # Opciones: "gaussian", "poisson", "binomial", "nb" (negative binomial)
NOT_BAM <- FALSE
GAM_LINK <- NULL # Opciones: "log", "sqrt", NULL
GAM_K <- 3  # Número de nudos por spline
GAM_XY_K <- GAM_K
GAM_ENV_K <- GAM_K
GAM_EVAL_METRICS <- TRUE  # Si el valor es FALSE nos se calculan las métricas en test
GAM_CONTROL = list(
  maxit = 50,          # reduce las iteraciones máximas
  mgcv.tol = 1e-4,     # relaja la tolerancia
  mgcv.half = 10       # Limita el step halving
)
GAM_CACHE <- TRUE  # si FALSE se deshabilita la cache del análisis GAM

PERFORM_RDA <- FALSE
if (PERFORM_RDA) {
TRANSFORM <- "log" # opciones: "sqrt", "log"
} else {TRANSFORM <- "none"}
SCALING <- 2
MIN_SITES <- 2
MEMS_NUMBER <- 20
NO_TOP_VARS <- FALSE
VAR_CUTOFF <- 50 # máximo número de variables admitido para no pasar por el proceso de selección
TOP_VARS <- 100 # 64 es un valor recomendable, no debería valer mucho más que 100
NBITS <- TOP_VARS
SELECTED_VARS <- 20
PROP_SEEDED <- 0.25
N_CORES <- FALSE # permite el procesamiento en paralelo, funciona mejor en FALSE
SEED <- 1
POPSIZE <- 250
MAXITER <- 100
ELITISM <- 10
PCROSSOVER <- 0.8
PMUTATION <- 0.1


strip_gam_for_prediction <- function(fit) {
  fit_stripped <- list(
    coefficients = fit$coefficients,
    smooth = fit$smooth,
    pterms = fit$pterms,
    terms = fit$terms,
    var.summary = fit$var.summary,
    family = fit$family,
    formula = fit$formula,
    nsdf = fit$nsdf,
    Vp = fit$Vp,               # necesario para el error estándar
    cmX = fit$cmX,             # necesario para el centrado
    n.para = fit$n.para,
    n.smooth = fit$n.smooth,
    n.fixed = fit$n.fixed
  )
  class(fit_stripped) <- c("bam", "gam", "glm", "lm")
  return(fit_stripped)
}

#Función encargada de realizar el modelado de los GAM
fit_gam_models_xy <- function(species_train, env_train, coord_train, 
                              species_test = NULL, env_test = NULL, coord_test = NULL,
                              family = GAM_FAMILY, 
                              xy_k = GAM_XY_K, 
                              env_k = GAM_ENV_K,
                              eval_metrics = GAM_EVAL_METRICS,
                              use_cache = GAM_CACHE,
                              control = GAM_CONTROL,
                              species_label = "", 
                              environment_label = "") {
  
  train_obs_orig <- species_train
  test_obs_orig <- species_test
  
  if(!NO_LOG_TRANSFORM_GAM){
    species_train = log(species_train+1)
    species_test = log(species_test+1)
  }
  
  # Crea un hash único para la cache
  cache_params <- list(
    species_train = digest(species_train),
    env_train = digest(env_train),
    coord_train = digest(coord_train),
    family = family,
    control = control,
    xy_k = xy_k,
    env_k = env_k
  )
  
  if(eval_metrics && !is.null(species_test)) {
    cache_params$species_test <- digest(species_test)
    cache_params$env_test <- digest(env_test)
    cache_params$coord_test <- digest(coord_test)
  }
  
  cache_hash <- digest(cache_params)
  cache_file <- paste0("gam_cache_", species_label, "-", environment_label, "-", cache_hash, ".rds")
  
  
  # Prepara el almacenamiento de los resultados
  gam_results <- list()
  train_metrics <- list(
    r2 = numeric(ncol(species_train)),
    rmse = numeric(ncol(species_train)),
    mae = numeric(ncol(species_train))
  )
  
  if(eval_metrics && !is.null(species_test)) {
    test_metrics <- list(
      r2 = numeric(ncol(species_train)),
      rmse = numeric(ncol(species_train)),
      mae = numeric(ncol(species_train))
    )
  } else {
    test_metrics <- NULL
  }

  # Crea la fórmula de los modelos GAM con splines para las coordenadas X e Y
  env_terms <- if(ncol(env_train) > 0) {
    paste0("s(", colnames(env_train), ", k = ", env_k, ")", collapse = " + ")
  } else { "" }
  
  xy_terms <- paste0("s(", colnames(coord_train), ", k = ", xy_k, ")", collapse = " + ")
  
  gam_formula <- if(ncol(env_train) > 0) {
    as.formula(paste("abundance ~", env_terms, "+", xy_terms))
  } else {
    as.formula(paste("abundance ~", xy_terms))
  }

  # Recupera la cache si está guardada
  if(use_cache && file.exists(cache_file)) {
    cached_data <- readRDS(cache_file)
    cat("Loaded GAM fits from cache\n")
    gam_results <- cached_data$models
  } else {
  # ajusta el GAM para cada especie
  for (i in seq_len(ncol(species_train))) {
    species_name <- colnames(species_train)[i]

    # crea un dataframe del modelo con nombres de columna explícitos
    model_data <- cbind(
      abundance = species_train[, i],  # Variable respuesta
      env_train,
      coord_train
    )
    
    cat("Ajustando especie:", species_name, "\n")
    
    # Realiza el ajuste
    tryCatch({
      # Permite llamar a la familia de funciones a partir de una cadena de caracteres
      function_family <- match.fun(family)
      if(NOT_BAM){
        gam_fit <- gam(gam_formula, 
                       data = model_data,
                       family = if(!is.null(GAM_LINK)){function_family(link=GAM_LINK)}else{function_family()},
                       discrete = TRUE,
                       method = "fREML",
                       control = control,
                       strip = TRUE)}
      else {
        gam_fit <- bam(gam_formula, 
                       data = model_data,
                       family = if(!is.null(GAM_LINK)){function_family(link=GAM_LINK)}else{function_family()},
                       discrete = TRUE,
                       method = "fREML",
                       control = control,
                       strip = TRUE)}
      
      fit_stripped <- strip_gam_for_prediction(gam_fit)
      # Guarda el modelo
      gam_results[[species_name]] <- fit_stripped
    }, error = function(e) {
      message("Error en el ajuste del GAM para la especie ", species_name, ": ", e$message)
    })}
    
    results_cache <- list(
      models = gam_results,
      cache_hash = cache_hash
    )
    
    # Guarda en la cache el ajuste (si la cache está habilitada)
    if(use_cache && !file.exists(cache_file)) {
      saveRDS(results_cache, file = cache_file, compress = "xz")
      cat("Se guardaron los ajustes GAM a la caché\n")
    }
    
    }
  
if(!exists("cached_data") || is.null(cached_data$train_metrics)){
  
  cat("Calculando las métricas de los modelos GAM...\n")
  
  # Precalular los dataframes con las predictoras =============================================
  train_predictor_df <- as.data.frame(cbind(env_train, coord_train))
  if(eval_metrics && !is.null(species_test)) {
    test_predictor_df <- as.data.frame(cbind(env_test, coord_test))
  }
  
  # Define la inversa de la trasformación logarítmica =======================================
  inv_transform <- if(NO_LOG_TRANSFORM_GAM) {
    function(x) x
  } else {
    function(x) exp(x) - 1
  }
  
  # Precálculo de las matrices de respuesta originales ========================================
  train_obs_orig_mat <- train_obs_orig
  if(eval_metrics && !is.null(species_test)) {
    test_obs_orig_mat <- test_obs_orig
  }
  
  # Inicializa el guardado de las métricas ==================================================
  n_species <- ncol(species_train)
  train_metrics <- list(
    r2 = rep(NA_real_, n_species),
    rmse = rep(NA_real_, n_species),
    mae = rep(NA_real_, n_species)
  )
  
  if(eval_metrics && !is.null(species_test)) {
    test_metrics <- list(
      r2 = rep(NA_real_, n_species),
      rmse = rep(NA_real_, n_species),
      mae = rep(NA_real_, n_species)
    )
  } else {
    test_metrics <- NULL
  }
  
  # Calcula las métricas =============================================================
  for (i in seq_len(n_species)) {
    species_name <- colnames(species_train)[i]
    gam_fit <- gam_results[[species_name]]
    
    print(paste("Calculando las métricas para la especie:", species_name))
    if(is.null(gam_fit)){
             train_metrics$r2[i] <- NA
             train_metrics$rmse[i] <- NA
             train_metrics$mae[i] <- NA
             test_metrics$r2[i] <- NA
             test_metrics$rmse[i] <- NA
             test_metrics$mae[i] <- NA    
             print(paste("No hubo ajuste para la especie:", species_name))
             next
    }
    # MÉTRICAS DE ENTRENAMIENTO ----------------------------------------------------------
    train_preds <- predict(gam_fit, newdata = train_predictor_df, type = "response")
    train_obs <- species_train[, i]
    
    # Manejar los datos faltantes
    complete_train <- !is.na(train_preds) & !is.na(train_obs)
    n_complete_train <- sum(complete_train)
    
    if(n_complete_train >= 2) {
      # convierte los predichos a la escala original
      train_preds_orig <- inv_transform(train_preds[complete_train])
      train_obs_orig <- train_obs_orig_mat[complete_train, i]
      
      # calcula las distintas métricas
      train_metrics$r2[i] <- cor(train_preds[complete_train], train_obs[complete_train])^2
      train_metrics$rmse[i] <- sqrt(mean((train_obs_orig - train_preds_orig)^2))
      train_metrics$mae[i] <- mean(abs(train_obs_orig - train_preds_orig))
    }
    
    # MÉTRICAS EN TEST (si se solicita su cálculo) ----------------------------------------------
    if(eval_metrics && !is.null(species_test)) {
      test_preds <- predict(gam_fit, newdata = test_predictor_df, type = "response")
      test_obs <- species_test[, i]
      
      # Manejo de valores faltantes
      complete_test <- !is.na(test_preds) & !is.na(test_obs)
      n_complete_test <- sum(complete_test)
      
      if(n_complete_test >= 2) {
        # Convierte las predicciones a la escala original
        test_preds_orig <- inv_transform(test_preds[complete_test])
        test_obs_orig <- test_obs_orig_mat[complete_test, i]
        
        # Calcula las métricas
        test_metrics$r2[i] <- cor(test_preds[complete_test], test_obs[complete_test])^2
        test_metrics$rmse[i] <- sqrt(mean((test_obs_orig - test_preds_orig)^2))
        test_metrics$mae[i] <- mean(abs(test_obs_orig - test_preds_orig))
      }
    }
  }
  
  # Prepara los resultados
  results <- list(
    models = gam_results,
    train_metrics = train_metrics,
    test_metrics = test_metrics,
    cache_hash = cache_hash
  )
  # Los guarda en la caché si está habilitada
  if(use_cache) {
    saveRDS(results, file = cache_file, compress = "xz")
    cat("Saved GAM fits to cache\n")
  }
  
  return(results)
} else{
  return(cached_data)
}
}


#Esta es la función principal que realiza los análisis y las salidas gráficas

perform_analysis <- function(species_df, environment_df, coord_df, vardesc, species_label, environment_label) {

print(match.call())  

test_prop <- TEST_PROP
  
# Setea la semilla para el split en train y test
set.seed(SEED)

# Corrobora que esté la columna para realizar el muestreo estratificado
if (!"NA_L2CODE" %in% colnames(coord_df)) {
  stop("coord_df debe contener la columna 'NA_L2CODE' para el muestreo estratificado")
}


if(grepl(EXCLUDE_OUTLIERS_KEYWORD, environment_label, fixed = TRUE)){
  print(paste("Número de filas en el dataset de variables ambientales, antes de sacar los outliers:",nrow(environment_df)))
  environment_df <- environment_df %>%
    filter(!(UID %in% OUTLIERS_UID))
  print(paste("Número de filas en el dataset de variables ambientales, después de sacar los outliers:",nrow(environment_df)))
  print(paste("Número de filas en el dataset de especies, antes de sacar los outliers:",nrow(species_df)))
  species_df <- species_df %>%
    filter(!(UID %in% OUTLIERS_UID))
  print(paste("Número de filas en el dataset de especies, después de sacar los outliers:",nrow(species_df)))
}

if(PCT_TRANSFORM != "none"){
  if(PCT_TRANSFORM == "logit")
{#Transformación logit de las variables porcentuales en el dataset de variables ambientales
 environment_df <- environment_df %>%
 mutate(across(starts_with("PCT"), ~ qlogis((. + 0.01) / 100.02)))}
  else if(PCT_TRANSFORM == "asinsqrt"){
    #Transformación arcoseno de la raíz cuadrada de las variables porcentuales en el dataset de variables ambientales
    environment_df <- environment_df %>%
      mutate(across(starts_with("PCT"),
                    ~ asin(sqrt((. + 0.01) / 100.02))))}
  else{stop("Transformación de variables porcentuales no reconocida")}
  }

# convierte el UID a entero en todos los dataframes fuente
coord_df_converted <- coord_df %>% mutate(UID = as.integer(UID))
species_df_converted <- species_df %>% mutate(UID = as.integer(UID))
environment_df_converted <- environment_df %>% mutate(UID = as.integer(UID))

# Crea un dataframe unificado para realizar el muestreo estratificado
all_data <- coord_df_converted %>%
  select(UID, stratum = NA_L2CODE) %>%
  inner_join(environment_df_converted %>% select(UID), by = "UID") %>%
  inner_join(species_df_converted %>% select(UID), by = "UID")

# Realiza el muestreo estratificado
test_uids <- all_data %>%
  group_by(stratum) %>%
  mutate(stratum_size = n()) %>%
  ungroup() %>%
  mutate(sample_prob = ifelse(
    stratum_size > 1, 
    test_prop, 
    min(1, test_prop)  # Manejo de los estratos de una única muestra
  )) %>%
  group_by(stratum) %>%
  sample_frac(size = unique(sample_prob)) %>%
  pull(UID)


test_n <- length(test_uids)
print(paste("Total test UIDs:", test_n))

#función que procesa los dataframes
process_dfs <- function(species_df, environment_df, coord_df)  
{
id_columns <- c("UID", "SITE_ID", "UNIQUE_ID")

# Asegura el alineamiento de UIDs
coord_df$NA_L2CODE <- NULL
env_reduced <- environment_df[environment_df$UID %in% species_df$UID, ]
coord_reduced <- coord_df[coord_df$UID %in% env_reduced$UID,]

env_data <- env_reduced[order(env_reduced$UID), ]
coord_data <- coord_reduced[order(coord_reduced$UID), ]
species_data <- species_df[order(species_df$UID), ]


# Toma un subconjunto de los datos para la posterior imputación
env_numeric <- env_data[, !names(env_data) %in% id_columns]

# Cache para los outliers del dataset de variables ambientales
env_outliers <- NULL
if (DETECT_ENV_OUTLIERS) {
  env_cache_params <- list(
    data = digest(env_numeric),
    params = list(
      ntrees = ISOLATION_FOREST_NTREES_ENV,
      quantile = OUTLIER_QUANTILE
    )
  )
  env_cache_hash <- digest(env_cache_params)
  env_cache_file <- paste0("iso_env_cache_", env_cache_hash, ".rds")
  
  if (file.exists(env_cache_file)) {
    cached_env <- readRDS(env_cache_file)
    outlier_scores_env <- cached_env$scores
    env_outliers <- cached_env$outliers
    cat("Se cargaron los outliers del dataset de variables ambientales de la caché\n")
  } else {
    iso_model <- isotree::isolation.forest(
      env_numeric, 
      ntrees = ISOLATION_FOREST_NTREES_ENV
    )
    outlier_scores_env <- predict(iso_model, env_numeric)
    env_outliers <- outlier_scores_env > quantile(outlier_scores_env, OUTLIER_QUANTILE)
    saveRDS(list(scores = outlier_scores_env, outliers = env_outliers), env_cache_file)
    cat("Se guardaron los outliers del dataset de variables ambientales a la caché\n")
  }
}

# Realiza la imputación en el dataset de variables ambientales

env_imputed <- imputePCA(env_numeric, ncp = 2, scale = TRUE, method = "Regularized",
                         coeff.ridge = 1, threshold = 1e-06, seed = NULL, nb.init = 1,
                         maxiter = 1000)

# Accede al dataset completado (con las imputaciones)
env_completed <- as.data.frame(env_imputed$completeObs)

# Si es necesario vuelve a introducir las columnas de ID en el dataframe
env_final <- cbind(env_data[id_columns], env_completed)

species_values <- species_data[, !(names(species_data) %in% id_columns)]


# Caché para los outliers en el dataset de especies
species_outliers <- NULL
if (DETECT_SPECIES_OUTLIERS) {
  species_cache_params <- list(
    data = digest(species_values),
    params = list(
      ntrees = ISOLATION_FOREST_NTREES_SPECIES,
      quantile = OUTLIER_QUANTILE
    )
  )
  species_cache_hash <- digest(species_cache_params)
  species_cache_file <- paste0("iso_species_cache_", species_cache_hash, ".rds")
  
  if (file.exists(species_cache_file)) {
    cached_species <- readRDS(species_cache_file)
    outlier_scores_species <- cached_species$scores
    species_outliers <- cached_species$outliers
    cat("Se cargaron los outliers del dataset de especies de la caché\n")
  } else {
    iso_model <- isotree::isolation.forest(
      species_values, 
      ntrees = ISOLATION_FOREST_NTREES_SPECIES
    )
    outlier_scores_species <- predict(iso_model, species_values)
    species_outliers <- outlier_scores_species > quantile(outlier_scores_species, OUTLIER_QUANTILE)
    saveRDS(list(scores = outlier_scores_species, outliers = species_outliers), species_cache_file)
    cat("Se guardaron los outliers del dataset de especies de la caché\n")
  }
}

filename <- paste0(environment_label, "-", sprintf("%0.2f", OUTLIER_QUANTILE), "-outliers.svg")
if(!file.exists(filename)){
  # Visualizar los outliers en el dataset de variables ambientales
  if (DETECT_ENV_OUTLIERS && !is.null(env_outliers)) {
    # PCA para la visualización
    pca_env <- prcomp(env_completed, scale. = TRUE)
    pca_scores <- as.data.frame(pca_env$x[, 1:2])

    # Agrega el UID y si es outlier o no
    pca_scores$UID <- env_data$UID
    pca_scores$observaciones <- factor(ifelse(env_outliers, "Atípicas", "Normales"))

    p_env <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = observaciones)) +
      geom_point(alpha = 0.7) +
      # Agrega las etiquetas a los outliers exclusivamente
      geom_text_repel(
        data = subset(pca_scores, observaciones == "Atípicas"),
        aes(label = UID),
        color = "red",
        size = 3,
        max.overlaps = 20
      ) +
      scale_color_manual(values = c("Normales" = "blue", "Atípicas" = "red")) +
      theme_bw()

    ggsave(filename, plot = p_env, width = 8, height = 6)
  }}



# Función para obtener un hash de un dataframe
get_data_hash <- function(param_list) {
  digest::digest(param_list, algo = "md5")
}

# Función que chequea si la versión cacheada existe
get_cached_imputation <- function(data_hash, cols) {
  cache_file <- paste0("knn_cache_", data_hash, "_", paste(species_label, environment_label, collapse="_"), ".rds")
  if(file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  return(NULL)
}

# Función para guardar los datos imputados a la caché
save_to_cache <- function(data, data_hash, cols) {
  cache_file <- paste0("knn_cache_", data_hash, "_", paste(species_label, environment_label, collapse="_"), ".rds")
  saveRDS(data, cache_file)
}

# Código de imputación (con manejo de caché)
if(sum(is.na(species_values)) > 0) {
  # Obtener el hash único
  data_hash <- get_data_hash(list(species = species_values, knn_n = KNN_IMPUTATION_N))
  cols_to_impute <- colnames(species_values)
  
  cached_data <- get_cached_imputation(data_hash, cols_to_impute)
  
  if(!is.null(cached_data)) {
    species_values <- cached_data
  } else {
    # Calcular el kNN si no está en la caché
    species_values <- VIM::kNN(species_values, variable = cols_to_impute, k = KNN_IMPUTATION_N, imp_var = FALSE)
    
    # Guardar en la caché
    save_to_cache(species_values, data_hash, cols_to_impute)
  }
}

species_data_imputed <- as.data.frame(cbind(species_data[,id_columns], species_values))


filename <- paste0(species_label, "-", sprintf("%0.2f", OUTLIER_QUANTILE), "-outliers-pca.svg")

if(!file.exists(filename)){
  # Visualizar los outliers del dataset de especies
  if (DETECT_SPECIES_OUTLIERS && !is.null(species_outliers)) {
    # Se aplica la transformación de Hellinger para poder correr el PCA
    species_hel <- decostand(species_values, "hellinger")
    pca_species <- prcomp(species_hel)
    pca_scores <- as.data.frame(pca_species$x[, 1:2])

    # Se agrega el UID y si es outlier o no
    pca_scores$UID <- species_data_imputed$UID
    pca_scores$observaciones <- factor(ifelse(species_outliers, "Atípicas", "Normales"))

    p_species <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = observaciones)) +
      geom_point(alpha = 0.7) +
      # Se agregan las etiquetas para los outliers solamente
      geom_text_repel(
        data = subset(pca_scores, observaciones == "Atípicas"),
        aes(label = UID),
        color = "red",
        size = 3,
        max.overlaps = 20
      ) +
      scale_color_manual(values = c("Normales" = "darkgreen", "Atípicas" = "red")) +
      theme_bw()

    ggsave(filename, plot = p_species, width = 8, height = 6)
  }}


zero_rows <- rowSums(species_values != 0) == 0

species_filtered <- species_data_imputed[!zero_rows, "UID", drop = FALSE ]

env_filtered <- env_final[env_final$UID %in% species_filtered$UID, ]
coord_filtered <- coord_data[coord_data$UID %in% species_filtered$UID, ]

# Verificar que los dataframes estén alineados
stopifnot(all(env_filtered$UID == species_filtered$UID))

# Crear el split entre train y test
species_train <- species_data_imputed[!zero_rows,] %>% filter(!UID %in% test_uids)
species_test <- species_data_imputed[!zero_rows,] %>% filter(UID %in% test_uids)

env_train <- env_filtered %>% filter(!UID %in% test_uids)
env_test <- env_filtered %>% filter(UID %in% test_uids)

coord_train <- coord_filtered %>% filter(!UID %in% test_uids)
coord_test <- coord_filtered %>% filter(UID %in% test_uids)

env_train <- env_train[, !colnames(env_train) %in% id_columns]
species_train <- species_train[, !colnames(species_train) %in% id_columns]
coord_train <- coord_train[, !colnames(coord_train) %in% id_columns]

env_test <- env_test[, !colnames(env_test) %in% id_columns]
species_test <- species_test[, !colnames(species_test) %in% id_columns]
coord_test <- coord_test[, !colnames(coord_test) %in% id_columns]

# quitar la columna UID (no necesaria para el análisis)
env <- env_filtered[, !colnames(env_filtered) %in% id_columns]
species <- species_values[!zero_rows, ]
coord <- coord_filtered[, !colnames(coord_filtered) %in% id_columns]

print("Total de especies")
print(ncol(species))

#Sacar las especies presentes en menos de MIN_SITES
print(paste("Especies presentes en al menos", MIN_SITES, "sitios"))
species <- species[, colSums(species > 0) >= MIN_SITES]

species_train <- species_train[, c(colnames(species))]
species_test <- species_test[, c(colnames(species))]


processed_dfs <- list(full_dataset = list(species = species, env = env, coord = coord, uid = coord_filtered$UID),
                      train_dataset = list(species = species_train, env = env_train, coord = coord_train),
                      test_dataset= list(species = species_test, env = env_test, coord = coord_test),
                      outliers_isof = list(species = species_data[species_outliers,] , env = env_data[env_outliers,]))

return(processed_dfs)
}


datasets <- process_dfs(species_df, environment_df, coord_df)

attach(datasets)
attach(full_dataset)

spp_nvars <- ncol(species)
print(spp_nvars)

train_n = nrow(train_dataset$species)
total_n = nrow(species)


env_std <- as.data.frame(scale(env))


if(PERFORM_RDA & TRANSFORM != "none"){if(TRANSFORM == "log"){
species <- log(species+1)
}else if (TRANSFORM == "sqrt")
species <- sqrt(species)
}

filename <- paste0(environment_label, "-all-vars-corrplot.svg")
if(!file.exists(filename))
{
  cov_matrix_std <- cov(env_std)
  if(ncol(env_std)>100){svg(filename, width=40, height=40, onefile = FALSE)}
  else {svg(filename, width=20, height=20, onefile = FALSE)}
  
  corrplot(cov_matrix_std, type = 'lower', order = 'hclust', tl.col = 'black',
           cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))
  
  dev.off()
}

filename <- paste0(environment_label, "-all-vars-pcorrplot.svg")
if(!file.exists(filename))
{
  partial_cor_matrix <- partial.r(data=env_std, method = "pearson")
  
  if(ncol(env_std)>100){svg(filename, width=40, height=40, onefile = FALSE)}
  else {svg(filename, width=20, height=20, onefile = FALSE)}
  
  corrplot(partial_cor_matrix, type = 'lower', order = 'hclust', tl.col = 'black',
           cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))
  
  dev.off()
}

# --------------------------------------------------
# Cachear la generación de MEMs y el test I de Moran
# --------------------------------------------------

# Crear un hash único para los parámetros de la generación de MEMs
mem_params <- list(
  species = species,
  coord = coord,          # coordenadas espaciales
  thresh = 0.1,           # umbral dbmem
  k = 5,                  # Valor de K de Kneardigh
  alpha = 0.05,           # Nivel de significancia del test I de Moran
  mems_number = MEMS_NUMBER
)
mem_hash <- digest(mem_params)
mem_cache_file <- paste0("mem_cache_", mem_hash, ".rds")

if (file.exists(mem_cache_file)) {
  # Cargar los MEMs y los resultados del test I de Moran
  cached_data <- readRDS(mem_cache_file)
  mem.signif <- cached_data$mem.signif
  mem.selected <- cached_data$mem.selected
  cat("Loaded MEMs from cache\n")
} else {
  # Calcularlo de cero
  distmat <- dist(coord)
  mem <- dbmem(distmat, thresh = 0.1)
  mem.signif <- mem[, which(attr(mem, "values") > 0)]
  
  # Test I de Moran
  nb <- knn2nb(knearneigh(coord, k = 5))
  lw <- nb2listw(nb)
  moran_results <- apply(mem.signif, 2, function(x) moran.test(x, lw)$p.value)
  mem.signif <- mem.signif[, moran_results < 0.05]

  # Usar tryCatch para manejar los errores de forward.sel
  mod.sel <- tryCatch(
    expr = {
      forward.sel(
        Y = mem_params$species,
        X = mem.signif,
        K = mem_params$mems_number,
        R2thresh = 0.99,
        nperm = 99,
        verbose = FALSE
      )
    },
    error = function(e) {
      # Verifica si el error se debe a que no se eligió ningún MEM
      if (grepl("No variables selected", e$message)) {
        message("No MEMs selected. Continuing with empty set.")
        # Devolver objeto dummy con 'order' seteado en 0 y con el elemento variables
        list(variables = character(0), order = integer(0))
      } else {
        # Arrojar otros errores no manejados
        stop(e)
      }
    }
  )

  # Usar los nombres de las variables para la selección (más seguro que usar índices)
  selected_vars <- mod.sel$variables
  
  mem.selected <- as.data.frame(mem.signif[, selected_vars, drop = FALSE])

  # Guardar a la caché
  saveRDS(list(
    mem.signif = mem.signif,
    mem.selected = mem.selected
  ), file = mem_cache_file)
  cat("Saved MEMs to cache\n")
}

# --------------------------------------------------
# Realizar una copia en la caché de la matriz de distancias
# --------------------------------------------------
dist_hash <- digest(list(coord = coord, method = "euclidean"))
dist_cache_file <- paste0("dist_cache_", dist_hash, ".rds")

if (file.exists(dist_cache_file)) {
  distmat <- readRDS(dist_cache_file)
  cat("Se cargó la matriz de distancias de la caché\n")
} else {
  distmat <- dist(coord)
  saveRDS(distmat, dist_cache_file)
  cat("Se guardó la matriz de distancias a la caché\n")
}

# --------------------------------------------------
# Verificar la integridad de la caché
# --------------------------------------------------
# Después de cargar de la caché, verificar que las dimensiones son compatibles con la de los datos presentes
stopifnot(
  nrow(mem.signif) == nrow(coord),
  ncol(mem.selected) <= ncol(mem.signif)
)

env_std_backup <- env_std

# Parámetros de la caché del proceso de reducción por VIFs
vif_cache_params <- list(
  env_std = env_std,
  mem.selected = mem.selected,
  environment_label = environment_label,
  species_label = species_label,
  vif_threshold = VIF_THRESHOLD,
  perform_rda = PERFORM_RDA
)
vif_cache_hash <- digest(vif_cache_params)
vif_cache_file <- paste0("vif_cache_", vif_cache_hash, ".rds")


# Verificar si la reducción por VIFs es solicitada y corroborar si existe la caché
if (grepl(VIF_ANALYSIS_KEYWORD, environment_label, fixed = TRUE)) {
  if (file.exists(vif_cache_file)) {
    cat("Cargando las variables ambientales reducidas por VIF de la caché\n")
    env_mems_vif_reduced <- readRDS(vif_cache_file)
  } else {
    cat("Realizando la reducción por VIF y guardando una caché de los resultados\n")
    repeat {
      env_std_with_mem_reduced <- cbind(env_std, mem.selected)      
      if(ncol(mem.selected) > 0) {
        formula_str <- paste(
          "species ~", 
          paste(names(env_std), collapse = " + "), 
          "+ Condition(", 
          paste(names(mem.selected), collapse = " + "), 
          ")"
        )} else {
          formula_str <- paste(
            "species ~", 
            paste(names(env_std), collapse = " + "))
        }
      
      if(PERFORM_RDA){
        canon_model <- rda(as.formula(formula_str), data = env_std_with_mem_reduced)
      } else {
        canon_model <- cca(as.formula(formula_str), data = env_std_with_mem_reduced)
      }
      
      # Calcular los VIFs
      vif_vals <- vif.cca(canon_model)

      # Verificar si todos los VIFs se encuentran por
      if (all(vif_vals < VIF_THRESHOLD)) break

      # Sacar la variable con VIF más alto
      var_to_remove <- names(which.max(vif_vals))
      message("Quitando variable por VIF elevado: ", var_to_remove)
      env_std[[var_to_remove]] <- NULL
      mem.selected[[var_to_remove]] <- NULL
    }
    env_mems_vif_reduced = list(env_std = env_std, mem.selected = mem.selected)
    saveRDS(env_mems_vif_reduced, vif_cache_file)
  }
  rm(env_std, mem.selected)
  env_std <- env_mems_vif_reduced$env_std
  mem.selected <- env_mems_vif_reduced$mem.selected
}

env_nvars <- ncol(env_std)
mem_nvars <- ncol(mem.selected)

# CCA parcial: efectos del ambiente luego de tener en consideración el espacio

env_std_with_mem <- cbind(env_std, mem.selected)
if(ncol(mem.selected) > 0) {
formula_str <- paste(
  "species ~", 
  paste(names(env_std), collapse = " + "), 
  "+ Condition(", 
  paste(names(mem.selected), collapse = " + "), 
  ")"
)}else{formula_str <- paste(
  "species ~", 
  paste(names(env_std), collapse = " + "))}

if(PERFORM_RDA){
  rda.res <- rda(as.formula(formula_str), data = env_std_with_mem)
  print(summary(rda.res))
}
else{
cca.res <- cca(as.formula(formula_str), data = env_std_with_mem)
print(summary(cca.res))

}

# --------------------------------------------------
# Calcular la contribución de las variables a la inercia explicada
# --------------------------------------------------
if(PERFORM_RDA){loadings <- scores(rda.res, display = "bp", choices = 1:rda.res$CCA$rank, correlation = TRUE, scaling = 0)
eigenvalues <- rda.res$CCA$eig}
else{loadings <- scores(cca.res, display = "bp", choices = 1:cca.res$CCA$rank, correlation = TRUE, scaling = 0)
eigenvalues <- cca.res$CCA$eig}

var_contrib <- apply(loadings^2, 1, function(x) sum(x * eigenvalues)) %>%
  sort(decreasing = TRUE)

print(paste("Suma de la contribución de todas las variables:",sum(var_contrib)))

selected_vars <- NULL

# Caso 1: Selección stepwise si 'minimal' está en environment_label
if (grepl("minimal", environment_label, fixed = TRUE)) {
  # Preparar fórmulas para ordistep
  if (ncol(mem.selected) > 0) {
    base_formula <- as.formula(paste("species ~ 1 + Condition(", paste(names(mem.selected), collapse = " + "), ")"))
    full_formula <- as.formula(paste("species ~", paste(names(env_std), collapse = " + "), "+ Condition(", paste(names(mem.selected), collapse = " + "), ")"))
  } else {
    base_formula <- as.formula("species ~ 1")
    full_formula <- as.formula(paste("species ~", paste(names(env_std), collapse = " + ")))
  }
  
  # Cache para modelos stepwise
  step_cache_params <- list(
    species = species,
    env_std = env_std,
    mem.selected = mem.selected,
    seed = SEED,
    Pin = 0.5,
    Pout = 0.5
  )
  step_hash <- digest(step_cache_params)
  step_cache_file <- paste0("minimal_cache_", step_hash, ".rds")
  
  if (file.exists(step_cache_file)) {
    step_model <- readRDS(step_cache_file)
    cat("Se cargó modelo stepwise de caché\n")
  } else {
    base_model <- cca(base_formula, data = cbind(env_std, mem.selected))
    full_model <- cca(full_formula, data = cbind(env_std, mem.selected))
    set.seed(step_cache_params$seed)
    step_model <- ordistep(base_model,
                           scope = formula(full_model),
                           direction = "both",
                           trace = 2,
                           Pin = step_cache_params$Pin,
                           Pout = step_cache_params$Pout)
    saveRDS(step_model, step_cache_file)
    cat("Se guardó modelo stepwise en caché\n")
  }
  
  # Extraer variables seleccionadas
  term_labels <- attr(step_model$terminfo$terms, "term.labels")
  ordered_env_terms <- term_labels[term_labels %in% names(env_std)]
  selected_vars <- head(ordered_env_terms, SELECTED_VARS)
  
} 

#Caso 2: Selección con MCGA si hay muchas variables

else if(ncol(env_std) > VAR_CUTOFF){
print(paste("Suma de la contribución de las variables en el top ranking (ranking de ",TOP_VARS,"variables):", sum(var_contrib[1:TOP_VARS])))

if(NO_TOP_VARS){
  TOP_VARS <- ncol(env_std)
  NBITS <- TOP_VARS}

top_vars <- names(var_contrib)[1:TOP_VARS]


# 1. Definir los parámetros y el hash único
input_params <- list(
  variables = top_vars,
  seed = SEED,
  nBits = NBITS,
  popSize = POPSIZE,
  maxiter = MAXITER,
  elitism = ELITISM,
  pcrossover = PCROSSOVER,
  pmutation = PMUTATION,
  mem.selected = mem.selected,
  species = species,
  env_std = env_std,
  prop_seeded = PROP_SEEDED,
  perform_rda = PERFORM_RDA,
  selected_vars = SELECTED_VARS
)
unique_hash <- digest(input_params)
filename_rds <- paste0("mcga_result_", unique_hash, ".rds")

# 1. Crear un vector mapeando los nombres de las variables con sus posiciones
var_rank <- rank(-var_contrib)  # Rankear por contribución (más alto = 1)
names(var_rank) <- names(var_contrib)

# 2. Crear soluciones semilla para el algoritmo genético según este ranking
n_seeded <- floor(input_params$popSize * input_params$prop_seeded)
seeded_solutions <- matrix(0, nrow = n_seeded, ncol = input_params$nBits)

# Tomar los índices de las mejores SELECTED_VARS en el conjunto completo
top_var_indices <- which(names(var_rank) %in% top_vars)
next_range_end <- min(length(var_contrib), input_params$selected_vars + 10)
next_best_indices <- which(names(var_rank) %in% names(sort(var_contrib, decreasing = TRUE)[(input_params$selected_vars+1):next_range_end]))

for (i in 1:n_seeded) {
    # Empezar con las variables que más contribuyen
    base_solution <- rep(0, input_params$nBits)
    base_solution[top_var_indices] <- 1
    
    # Introducir diversidad controlada
    n_swaps <- sample(1:3, 1)
    if(length(top_var_indices) >= n_swaps && length(next_best_indices) >= n_swaps) {
        swap_out <- sample(top_var_indices, n_swaps)
        swap_in <- sample(next_best_indices, n_swaps)
        base_solution[swap_out] <- 0
        base_solution[swap_in] <- 1
    }
    
    seeded_solutions[i, ] <- base_solution
}


# 3. Función de inicialización modificada que conserva las posiciones de las variables

custom_init <- function(object) {
  
  # Crear población al azar usando los parámetros externos
  n <- input_params$popSize * input_params$nBits
  population <- matrix(
    runif(n),
    nrow = input_params$popSize,
    ncol = input_params$nBits
  )
  
  # Reemplazar los valores NA si existiecen con valores aleatorios
  if (any(is.na(population))) {
    population[is.na(population)] <- runif(sum(is.na(population)))
  }
  
  # Introducir las soluciones semilla
  if (n_seeded > 0) {
    population[1:n_seeded, ] <- seeded_solutions
  }
  
  for (i in (n_seeded + 1):input_params$popSize) {
    ones_pos <- sample(1:input_params$nBits, input_params$selected_vars)
    population[i, ones_pos] <- 1
  }
  
  return(population)
  
}

# Función de fitness

fitness_mcga <- function(x) {
  
  # Verificar si hay NAs en los cromosomas
  if (any(is.na(x))) {
    return(-1e9)  # Penalización fuerte si existen NAs
  }
  
  # Convertir los valores reales a binarios usando un umbral
  subset <- ifelse(x >= 0.5, 1, 0)
  
  # aplicar restricción: solo aceptar SELECTED_VARS número de variables
  if (sum(subset) != input_params$selected_vars) {
    # Penalización proporcional a la desviación
    return(-1e6 * abs(sum(subset) - input_params$selected_vars))
  }
  
  selected <- input_params$variables[as.logical(subset)]
  
  # Corroborar si ocurre una selección inválida
  if (length(selected) == 0 || !all(selected %in% colnames(input_params$env_std))) {
    return(-1e6)  # Fuerte penalización
  }
  
  env_selected <- input_params$env_std[, selected]
  env_selected_with_mem <- cbind(env_selected, input_params$mem.selected)
  
  
  if(ncol(mem.selected) > 0) {
    rerun_formula_str <- paste(
      "species ~", 
      paste(selected, collapse = " + "), 
      "+ Condition(", 
      paste(names(input_params$mem.selected), collapse = " + "), 
      ")"
    )}else{rerun_formula_str <- paste(
      "species ~", 
      paste(selected, collapse = " + "))}
  

  if(input_params$perform_rda){
  rda_temp <- rda(as.formula(rerun_formula_str), data = env_selected_with_mem)
  return(sum(rda_temp$CCA$eig))
  }
  else{
    cca_temp <- cca(as.formula(rerun_formula_str), data = env_selected_with_mem)
    return(cca_temp$CCA$tot.chi)
  }
  
}


# Corroborar si los resultados fueron guardados o correr MCGA
if (file.exists(filename_rds)) {
  mcga_result <- readRDS(filename_rds)
} else {
  
  # correr MCGA2 con parámetros específicos del problema
  mcga_result <- mcga2(
    fitness = fitness_mcga,
    min = rep(0, input_params$nBits),  # Límite inferior (rango de 0 a 1)
    max = rep(1, input_params$nBits),  # Límite superior
    popSize = input_params$popSize,
    maxiter = input_params$maxiter,
    pcrossover = input_params$pcrossover,
    pmutation = input_params$pmutation,
    elitism = input_params$elitism,
    crossover = byte_crossover,
    mutation = byte_mutation,
    parallel = N_CORES,
    # esta función de inicialización personalilzada introduce las soluciones semilla
    population = custom_init,
    seed = input_params$seed
  )
  saveRDS(mcga_result, filename_rds)
}

# Extraer la mejor solución
best_solution <- ifelse(mcga_result@solution[1, ] >= 0.5, 1, 0)
selected_vars <- input_params$variables[as.logical(best_solution)]

# Corroborar si alguna variable fue seleccionada
if (length(selected_vars) == 0) {
  stop("Ninguna variable fue seleccionada. Ajustar parámetros")
}
}
# Caso 3: Selección forward a SELECTED_VARS si hay pocas variables
else if((ncol(env_std) <= VAR_CUTOFF) & (ncol(env_std)>SELECTED_VARS)){
  # Paso 1: Modelo base (solo MEMs)
  if (ncol(mem.selected) > 0) {
    base_formula <- as.formula(paste("species ~ 1 + Condition(", paste(names(mem.selected), collapse = " + "), ")"))
    full_formula <- as.formula(paste("species ~", paste(names(env_std), collapse = " + "), "+ Condition(", paste(names(mem.selected), collapse = " + "), ")"))
  } else {
    base_formula <- as.formula("species ~ 1")
    full_formula <- as.formula(paste("species ~", paste(names(env_std), collapse = " + ")))
  }
  
  # Paso 2: Forward selection con límite de SELECTED_VARS variables
  step_cache_params <- list(
    species = species,
    env_std = env_std,
    mem.selected = mem.selected,
    steps = SELECTED_VARS,
    permutations = 199,
    Pin = 0.5,
    seed = SEED
  )
  step_hash <- digest(step_cache_params)
  step_cache_file <- paste0("reduced_cache_", step_hash, ".rds")
  
  if (file.exists(step_cache_file)) {
    step_model <- readRDS(step_cache_file)
  } else {
    base_model <- cca(base_formula, data = cbind(env_std, mem.selected))
    full_model <- cca(full_formula, data = cbind(env_std, mem.selected))
    set.seed(step_cache_params$seed)
    step_model <- ordistep(
      object    = base_model,
      scope     = list(lower = base_formula, upper = full_formula),
      direction = "forward",
      steps = step_cache_params$steps,
      trace = 2,
      permutations = step_cache_params$permutations,
      Pin = step_cache_params$Pin
    )
    saveRDS(step_model, step_cache_file)
  }
  
  # Paso 3: Extraer variables seleccionadas
  selected_vars <- intersect(names(env_std), attr(step_model$terminfo$terms, "term.labels"))
  }
# Caso 4: Hay menor o igual cantidad de variables que SELECTED_VARS
else {selected_vars <- colnames(env_std)}

# --------------------------------------------------
# Crear la tabla de salida con las contribuciones
# --------------------------------------------------
var_contrib_df <- as.data.frame(var_contrib)
names(var_contrib_df)[1] <- "Total_Contribution"

# Unir las contribuciones con loadings y descripciones

if(PERFORM_RDA){selected_vars_df <- data.frame(
  COLUMN_NAME = selected_vars,
  loadings[selected_vars, ],
  Total_Contribution = var_contrib_df$Total_Contribution[match(selected_vars, var_contrib_df$COLUMN_NAME)]
) %>%
  left_join(vardesc[, c("COLUMN_NAME", "LABEL")], by = "COLUMN_NAME") %>%
  mutate(
    # Calcular la proporción de la varianza total explicada
    Prop_Explained = Total_Contribution / sum(rda.res$CCA$eig[1:4])
  ) %>%  select(COLUMN_NAME, RDA1, RDA2, RDA3, RDA4, Prop_Explained, LABEL)} else{
selected_vars_df <- data.frame(
  COLUMN_NAME = selected_vars,
  loadings[selected_vars, ],
  Total_Contribution = var_contrib_df$Total_Contribution[match(selected_vars, var_contrib_df$COLUMN_NAME)]
) %>%
  left_join(vardesc[, c("COLUMN_NAME", "LABEL")], by = "COLUMN_NAME") %>%
  mutate(
    # Calcular la proporción de la incercia restringida
    Prop_Explained = Total_Contribution / sum(cca.res$CCA$eig[1:4])
  ) %>%  select(COLUMN_NAME, CCA1, CCA2, CCA3, CCA4, Prop_Explained, LABEL)
  }


# Imprimir las variables seleccionadas
print(paste("Variables seleccionadas:", paste(selected_vars, collapse = ", ")))
if(PERFORM_RDA){
  selected_loadings_df <- data.frame(
    COLUMN_NAME = selected_vars,
    RDA1 = loadings[selected_vars, 1],
    RDA2 = loadings[selected_vars, 2],
    RDA3 = loadings[selected_vars, 3],
    RDA4 = loadings[selected_vars, 4]
  )
  print(names(selected_loadings_df))
  
  # Unir con las descripciones de las variables
  selected_loadings_df <- merge(
    selected_loadings_df,
    vardesc[, c("COLUMN_NAME", "LABEL")],
    by = "COLUMN_NAME",
    all.x = TRUE
  )
  # Reordenar las columnas
  selected_loadings_df <- selected_loadings_df[, c("COLUMN_NAME", "RDA1", "RDA2", "RDA3", "RDA4", "LABEL")]
  filename <- paste0(species_label, "-", environment_label, "-", TRANSFORM ,"-spatial-rda-selected_variables_loadings.csv")
}else{
# Crear un dataframe con los loadings y las descripciones
selected_loadings_df <- data.frame(
  COLUMN_NAME = selected_vars,
  CCA1 = loadings[selected_vars, 1],
  CCA2 = loadings[selected_vars, 2],
  CCA3 = loadings[selected_vars, 3],
  CCA4 = loadings[selected_vars, 4]
)
# 2. Unir con las descripciones de las variables
selected_loadings_df <- merge(
  selected_loadings_df,
  vardesc[, c("COLUMN_NAME", "LABEL")],
  by = "COLUMN_NAME",
  all.x = TRUE
)
# Reordenar las columnas
selected_loadings_df <- selected_loadings_df[, c("COLUMN_NAME", "CCA1", "CCA2", "CCA3", "CCA4", "LABEL")]

filename <- paste0(species_label, "-", environment_label, "-", "spatial-cca-selected_variables_loadings.csv")
}
write.csv(selected_loadings_df, filename, row.names = FALSE)



# crear el dataframe de variables ambientales seleccionadas
env_selected <- env_std[, selected_vars, drop = FALSE]
env_selected_with_mem <- cbind(env_selected, mem.selected)

if(ncol(mem.selected) > 0) {
rerun_formula_str <- paste(
  "species ~", 
  paste(names(env_selected), collapse = " + "), 
  "+ Condition(", 
  paste(names(mem.selected), collapse = " + "), 
  ")"
)}else{rerun_formula_str <- paste(
  "species ~", 
  paste(names(env_selected), collapse = " + "))}


if(PERFORM_RDA){
  # Volver a correr el RDA con las variables ambientales seleccionadas (y los MEMs)
  rda_rerun <- rda(as.formula(rerun_formula_str), data = env_selected_with_mem)

  print(summary(rda_rerun))
  
  cat("Varianza explicada del RDA original:", rda.res$CCA$tot.chi / rda.res$tot.chi, "\n")
  cat("Varianza explicada del RDA con selección de variables:", rda_rerun$CCA$tot.chi / rda_rerun$tot.chi, "\n")
  
}else{
# Volver a correr el CCA con las variables ambientales seleccionadas (y los MEMs)
cca_rerun <- cca(as.formula(rerun_formula_str), data = env_selected_with_mem)

print(summary(cca_rerun))

cat("Inercia restringida del CCA original:", cca.res$CCA$tot.chi / cca.res$tot.chi, "\n")
cat("Inercia restringida del CCA con selección de variables:", cca_rerun$CCA$tot.chi / cca_rerun$tot.chi, "\n")
}


if(PERFORM_RDA){
  site_scores <- scores(rda.res, display = "wa", choices = 1:4)
  colnames(site_scores) <- paste0("RDA", 1:ncol(site_scores))
  # env_std tiene las variables ambientales estandarizadas
  cor_loadings <- cor(env_std, site_scores, method = "pearson")
  
  # Convertir al formato largo
  cor_melt <- melt(cor_loadings)
  colnames(cor_melt) <- c("Variable", "Eje_RDA", "Correlacion")
  
  # Graficar
  p <- ggplot(cor_melt, aes(x = RDA_Axis, y = Variable, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    geom_text(aes(label = round(Correlation, 2)), color = "black", size = 4) +
    theme_minimal() +
    labs(title = "Loadings de las variables ambientales en los ejes del RDA")
  
  filename <- paste0(species_label, "-", environment_label, "-", TRANSFORM,"-spatial-rda-loadings_plot.svg")
  # Guardar como SVG
  ggsave(filename, plot = p, width = 8, height = 80, device = "svg", limitsize = FALSE)
  
  arrow_scale <- 50  # Aumentar este valor hace las flechas más largas
  
  env_scores <- scores(rda_rerun, display = "bp", scaling = SCALING, choices=c(1:4)) %>% 
    as.data.frame() %>% 
    mutate(
      RDA1 = RDA1 * arrow_scale,  # Scale arrow length
      RDA2 = RDA2 * arrow_scale,
      RDA3 = RDA3 * arrow_scale,
      RDA4 = RDA4 * arrow_scale
    ) %>% 
    tibble::rownames_to_column("Variable")
  
  site_scores <- scores(rda_rerun, display = "wa", scaling = SCALING, choices=c(1:4)) %>% 
    as.data.frame()
  
  species_scores <- scores(rda_rerun, display = "sp", scaling = SCALING, choices=c(1:4)) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Species")
  
  # Crear gráfico
  p <- ggplot() +
    # Graficar sitios como círculos rojos
    geom_point(data = site_scores, 
               aes(x = RDA1, y = RDA2), 
               shape = 19, color = "red", size = 2, alpha = 0.5) +
    
    # Graficar las especies como triángulos verdes
    geom_point(data = species_scores,
               aes(x = RDA1, y = RDA2, color = "darkgreen"),
               shape = 17, size = 2, alpha = 0.7) +
  
    # Flechas de variables ambientales (escaladas)
    geom_segment(data = env_scores, 
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 arrow = arrow(length = unit(0.03, "npc")),      # Punta de la flecha más grande
                 color = "darkblue", linewidth = 0.8) +          # Línea más gruesa
    
    # Etiquetas de variables ambientales (en la punta de las flechas)
    geom_text_repel(data = env_scores,
                    aes(x = RDA1, y = RDA2, label = Variable),
                    color = "black", size = 4, max.overlaps = 20,
                    point.padding = 0.5) +  # Agregar espacio entre la flecha y la etiqueta
    
    # Líneas de los ejes y etiquetas
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    labs(x = paste("RDA1 (", round(rda_rerun$CCA$eig[1]/sum(rda_rerun$CCA$eig)*100, 1), "%)"),
         y = paste("RDA2 (", round(rda_rerun$CCA$eig[2]/sum(rda_rerun$CCA$eig)*100, 1), "%)")) +
    
    # Relación de aspecto y tema
    coord_fixed() +
    theme_bw()
  
  filename <- paste0(species_label, "-", environment_label, "-", TRANSFORM, "-spatial-rda-1-2_plot.svg")
  # guardar como SVG
  ggsave(filename, plot = p, width = 8, height = 6, device = "svg")
  
  # Create plot
  p <- ggplot() +
    # Graficar sitios como círculos rojos
    geom_point(data = site_scores, 
               aes(x = RDA3, y = RDA4), 
               shape = 19, color = "red", size = 2, alpha = 0.5) +
    
    # Graficar las especies como triángulos verdes
    geom_point(data = species_scores, 
               aes(x = RDA3, y = RDA4), 
               shape = 17, color = "darkgreen", size = 2, alpha = 0.7) +
    
    # Flechas de variables ambientales (escaladas)
    geom_segment(data = env_scores, 
                 aes(x = 0, y = 0, xend = RDA3, yend = RDA4),
                 arrow = arrow(length = unit(0.03, "npc")),  # Larger arrowhead
                 color = "darkblue", linewidth = 0.8) +          # Thicker line
    
    # Etiquetas de variables ambientales (en la punta de las flechas)
    geom_text_repel(data = env_scores,
                    aes(x = RDA3, y = RDA4, label = Variable),
                    color = "black", size = 4, max.overlaps = 20,
                    point.padding = 0.5) +  # Add space between label and arrow
    
    # Líneas de los ejes y etiquetas
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    labs(x = paste("RDA3 (", round(rda_rerun$CCA$eig[3]/sum(rda_rerun$CCA$eig)*100, 1), "%)"),
         y = paste("RDA4 (", round(rda_rerun$CCA$eig[4]/sum(rda_rerun$CCA$eig)*100, 1), "%)")) +
    
    # Relación de aspecto y tema
    coord_fixed() +
    theme_bw()
  
  filename <- paste0(species_label, "-", environment_label, "-", TRANSFORM, "-spatial-rda-3-4_plot.svg")
  # Guardar como SVG
  ggsave(filename, plot = p, width = 8, height = 6, device = "svg")
  
  outliers_cca <- NULL
} else{
site_scores <- scores(cca.res, display = "wa", choices = 1:4)
colnames(site_scores) <- paste0("CCA", 1:ncol(site_scores))     # Rename columns
# env_std tiene las variables ambientales estandarizadas
cor_loadings <- cor(env_std, site_scores, method = "pearson")

# Convertir al formato largo
cor_melt <- melt(cor_loadings)
colnames(cor_melt) <- c("Variable", "CCA_Axis", "Correlation")

# Graficar
p <- ggplot(cor_melt, aes(x = CCA_Axis, y = Variable, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 4) +
  theme_minimal() +
  labs(title = "Environmental Variable Loadings on CCA Axes")

filename <- paste0(species_label, "-", environment_label, "-", "spatial-cca-loadings_plot.svg")
# Guardar como SVG
ggsave(filename, plot = p, width = 8, height = 80, device = "svg", limitsize = FALSE)

# extraer los scores y hacer las flechas de las variables más largas
arrow_scale <- 20  # Aumentar este valor para hacer las flechas más largas

env_scores <- scores(cca_rerun, display = "bp", scaling = SCALING, choices=c(1:4)) %>% 
  as.data.frame() %>% 
  mutate(
    CCA1 = CCA1 * arrow_scale,  # Escalar el largo de la flecha
    CCA2 = CCA2 * arrow_scale,
    CCA3 = CCA3 * arrow_scale,
    CCA4 = CCA4 * arrow_scale
  ) %>% 
  tibble::rownames_to_column("Variable")

site_scores <- scores(cca_rerun, display = "wa", scaling = SCALING, choices=c(1:4)) %>% 
  as.data.frame()

# Agregar los UIDs del dataset completo
site_scores$UID <- full_dataset$uid
outliers_cca <- list(species = cbind(full_dataset$uid,full_dataset$species), env = cbind(full_dataset$uid,full_dataset$env))

species_scores <- scores(cca_rerun, display = "sp", scaling = SCALING, choices=c(1:4)) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Species")

# Calcular las distancias y los umbrales para las especies
species_scores$dist12 <- sqrt(species_scores$CCA1^2 + species_scores$CCA2^2)
species_scores$dist34 <- sqrt(species_scores$CCA3^2 + species_scores$CCA4^2)
threshold_species12 <- quantile(species_scores$dist12, CCA_OUTLIER_QUANTILE)
threshold_species34 <- quantile(species_scores$dist34, CCA_OUTLIER_QUANTILE)

# Calcular la distancia para la detección de outliers
site_scores$dist12 <- sqrt(site_scores$CCA1^2 + site_scores$CCA2^2)
site_scores$dist34 <- sqrt(site_scores$CCA3^2 + site_scores$CCA4^2)

# Determinar los umbrales de los outliers
threshold12 <- quantile(site_scores$dist12, CCA_OUTLIER_QUANTILE)
threshold34 <- quantile(site_scores$dist34, CCA_OUTLIER_QUANTILE)

# después de calcular los umbrales, extraer los UIDs de los outliers
outlier_uids_12 <- site_scores$UID[site_scores$dist12 > threshold12]
outlier_uids_34 <- site_scores$UID[site_scores$dist34 > threshold34]
all_outlier_uids <- unique(c(outlier_uids_12, outlier_uids_34))

# extraer la data original de estos outliers
species_outliers <- species_df[species_df$UID %in% all_outlier_uids, ]
env_outliers <- environment_df[environment_df$UID %in% all_outlier_uids, ]

# crear la lista de outliers del CCA
outliers_cca <- list(
  species = species_outliers,
  env = env_outliers
)

# Graficar
p <- ggplot() +
  # Graficar sitios como círculos rojos
  geom_point(data = site_scores, 
             aes(x = CCA1, y = CCA2), 
             shape = 19, color = "red", size = 2, alpha = 0.5) +
  
  # Graficar las especies como triángulos verdes
  geom_point(data = species_scores, 
             aes(x = CCA1, y = CCA2), 
             shape = 17, color = "darkgreen", size = 2, alpha = 0.7) +
  
  # Flechas de variables ambientales (escaladas)
  geom_segment(data = env_scores, 
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.03, "npc")),      # Punta de la flecha más grande
               color = "darkblue", linewidth = 0.8) +          # Línea más gruesa
  
  # Etiquetas de variables ambientales (en la punta de las flechas)
  geom_text_repel(data = env_scores,
                  aes(x = CCA1, y = CCA2, label = Variable),
                  color = "black", size = 4, max.overlaps = 20,
                  point.padding = 0.5) +  # Agregar espacio entre la flecha y la etiqueta
  geom_text_repel(
    data = subset(site_scores, dist12 > threshold12),
    aes(x = CCA1, y = CCA2, label = UID),
    color = "red", size = 3, max.overlaps = 20
  ) + 
  geom_text_repel(
    data = subset(species_scores, dist12 > threshold_species12),
    aes(x = CCA1, y = CCA2, label = Species),
    color = "darkgreen", size = 3, max.overlaps = 20
  ) +
  # Líneas de los ejes y etiquetas
  geom_hline(yintercept = 0, linetype = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  labs(x = paste("CCA1 (", round(cca_rerun$CCA$eig[1]/sum(cca_rerun$CCA$eig)*100, 1), "%)"),
       y = paste("CCA2 (", round(cca_rerun$CCA$eig[2]/sum(cca_rerun$CCA$eig)*100, 1), "%)")) +
  
  # Relación de aspecto y tema
  coord_fixed() +
  theme_bw()

filename <- paste0(species_label, "-", environment_label, "-", "spatial-cca-1-2_plot.svg")
# guardar como SVG
ggsave(filename, plot = p, width = 8, height = 6, device = "svg")

# Graficar
p <- ggplot() +
  # Graficar sitios como círculos rojos
  geom_point(data = site_scores, 
             aes(x = CCA3, y = CCA4), 
             shape = 19, color = "red", size = 2, alpha = 0.5) +
  
  # Graficar las especies como triángulos verdes
  geom_point(data = species_scores, 
             aes(x = CCA3, y = CCA4), 
             shape = 17, color = "darkgreen", size = 2, alpha = 0.7) +
  
  # Flechas de variables ambientales (escaladas)
  geom_segment(data = env_scores, 
               aes(x = 0, y = 0, xend = CCA3, yend = CCA4),
               arrow = arrow(length = unit(0.03, "npc")),     # Punta de la flecha más grande
               color = "darkblue", linewidth = 0.8) +         # Línea más gruesa
  
  # Etiquetas de variables ambientales (en la punta de las flechas)
  geom_text_repel(data = env_scores,
                  aes(x = CCA3, y = CCA4, label = Variable),
                  color = "black", size = 4, max.overlaps = 20,
                  point.padding = 0.5) +  # Agregar espacio entre la flecha y la etiqueta
  geom_text_repel(
    data = subset(site_scores, dist34 > threshold34),
    aes(x = CCA3, y = CCA4, label = UID),
    color = "red", size = 3, max.overlaps = 20
  ) + 
  geom_text_repel(
    data = subset(species_scores, dist34 > threshold_species34),
    aes(x = CCA3, y = CCA4, label = Species),
    color = "darkgreen", size = 3, max.overlaps = 20,
    fontface = "bold"
  ) +
  # Líneas de los ejes y etiquetas
  geom_hline(yintercept = 0, linetype = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  labs(x = paste("CCA3 (", round(cca_rerun$CCA$eig[3]/sum(cca_rerun$CCA$eig)*100, 1), "%)"),
       y = paste("CCA4 (", round(cca_rerun$CCA$eig[4]/sum(cca_rerun$CCA$eig)*100, 1), "%)")) +
  
  # Relación de aspecto y tema
  coord_fixed() +
  theme_bw()

filename <- paste0(species_label, "-", environment_label, "-", "spatial-cca-3-4_plot.svg")
# Guardar como SVG
ggsave(filename, plot = p, width = 8, height = 6, device = "svg")
}

env_pca <- env_std[,selected_vars]

env.pca <- rda(env_pca)

pca_result <- prcomp(env_pca)

filename <- paste0(environment_label, "-ggbiplot.svg")
if(!file.exists(filename)){
theme <- theme(text = element_text(size=10),
               plot.title = element_text(size=12, face="bold.italic",
                                         hjust = 0.5),
               axis.title.x = element_text(size=10, face="bold", colour='black'),
               axis.title.y = element_text(size=10, face="bold"),
               panel.border = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title = element_text(face="bold"))

#Comparación de biplots
#Biplot del PCA clásico
pca_biplot <- ggbiplot(pca_result, var.factor = 5, alpha=0.25, varname.adjust = 1.5, var.alpha = 0.5, varname.color = "red") + theme
  

# Guardar como SVG
ggsave(filename, plot = pca_biplot, width = 8, height = 6, device = "svg")
}

filename <- paste0(environment_label, "-corrplot.svg")
if(!file.exists(filename)){

cov_matrix_std <- cov(env_std[,selected_vars])

svg(filename, onefile = FALSE)

corrplot(cov_matrix_std, type = 'lower', order = 'hclust', tl.col = 'black',
         cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))

dev.off()
}

filename <- paste0(environment_label, "-pcorrplot.svg")
if(!file.exists(filename)){
partial_cor_matrix <- partial.r(data=env_std[,selected_vars], method = "pearson")

svg(filename, onefile = FALSE)

corrplot(partial_cor_matrix, type = 'lower', order = 'hclust', tl.col = 'black',
         cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))

dev.off()
}


# setear el umbral de correlaciones
CORR_THRESHOLD <- 0.8
cov_matrix_std <- cov(env_std[,selected_vars])
# Identificar las variables altamente correlacionadas
cor_matrix <- cov_matrix_std  # Esta es una matriz de covarianzas
if(ncol(cor_matrix) > 1) {  # Solo proceder si tenemos múltiples variables
  # Convertir la matriz de covarianzas en correlación
  corr_matrix <- cov2cor(cor_matrix)
  
  # Encontrar pares altamente correlacionados (valor absoluto > umbral)
  high_corr_pairs <- which(abs(corr_matrix) > CORR_THRESHOLD & 
                             upper.tri(corr_matrix, diag = FALSE), 
                           arr.ind = TRUE)
  
  # Extraer los nombres de las variables con alta correlación
  high_corr_vars <- unique(c(
    rownames(corr_matrix)[high_corr_pairs[, 1]],
    colnames(corr_matrix)[high_corr_pairs[, 2]]
  ))
} else {
  high_corr_vars <- character(0)
}

# Crear el ggpairs plot para variables correlacionadas, solo si las hallamos
if(length(high_corr_vars) > 1) {
  filename_focused <- paste0(environment_label, "-highcorr-ggpairs.png")
  
  if(!file.exists(filename_focused)) {
    cat("Creando gráfico de ggpairs para", length(high_corr_vars), "variables altamente correlacionadas\n")
    
    # Calcular el tamaõ del gráfico en función de la cantidad de variables
    plot_size <- max(6, 2 * length(high_corr_vars))
    text_scale <- 1.2  # Factor de escalado del texto
    
    focused_plot <- env_std[, high_corr_vars, drop = FALSE] %>%
      ggpairs(
        upper = list(continuous = wrap("cor", 
                                       size = 5 * text_scale,  #Agrandar el texto
                                       hjust = 0.5)),
        diag = list(continuous = wrap("densityDiag", size = 0.8)),
        lower = list(continuous = wrap("points", size = 1.5, alpha = 0.5))
      ) +
      theme_bw(base_size = 14 * text_scale) +  # Tamaño base de la fuente tipográfica
      theme(
        axis.text = element_text(size = 14 * text_scale),
        axis.title = element_text(size = 14 * text_scale, face = "bold"),
        strip.text = element_text(
          size = 14 * text_scale, 
          face = "bold",
          margin = margin(0.2, 0, 0.2, 0, "cm")  # agregar espaciado
        ),
        panel.spacing = unit(1, "lines")  # Más espacio entre los paneles
      )
    
    ggsave(filename_focused, 
           plot = focused_plot, 
           width = plot_size,
           height = plot_size,
           dpi = 300)
  }
} else {
  cat("No se hallaron variables altamente correlacionadas por encima del umbral:", CORR_THRESHOLD, "\n")
}


filename <- paste0(environment_label, "-ggpairs.png")

if(!file.exists(filename)){

grafico_corr <- env %>% 
  dplyr::select(all_of(selected_vars)) %>%
  ggpairs(., upper = list(continuous = wrap("cor", size = 3, hjust=0.5)), legend = 25) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5), legend.position = "bottom")

# Guardar como SVG
ggsave(filename, plot = grafico_corr, width = 20, height = 20, dpi = 600, device = "png")
}

source("cleanplot.pca.R")


filename <- paste0(environment_label, "-pca_scaling1.svg")
if(!file.exists(filename)){
svg(paste0(environment_label, "-pca_scaling1.svg"), onefile = FALSE)
cleanplot.pca(env.pca, scaling = 1, mar.percent = 0.08)

dev.off()

svg(paste0(environment_label, "-pca_scaling2.svg"), onefile = FALSE)
cleanplot.pca(env.pca, scaling = 2, mar.percent = 0.04)

dev.off()
}

if(ncol(env[,selected_vars]) <= 20){
  N_COL <- 4
  N_ROW <- 5
}else{
  N_COL <- 5
  N_ROW <- 6
}
  
filename <-paste0(environment_label, "-boxplots.svg")
if(!file.exists(filename)){
long_data <- environment_df[,selected_vars] %>%
  pivot_longer(cols = names(env_pca), names_to = "Variable", values_to = "Value")

p <- ggplot(long_data, aes(x = "", y = Value)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  facet_wrap(~ Variable, scales = "free_y", ncol = N_COL, nrow = N_ROW) +  # Ajustar el número de columnas para el diseño
  labs(x = NULL, y = NULL) +  # Sacar el eje X
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  # Ocultar el texto del eje X
    axis.ticks.x = element_blank(), # Ocultar los ticks del eje X
    strip.text = element_text(size = 8)  # Ajustar el tamaño de texto de las etiquetas de las facetas
  )
ggsave(filename, plot = p, device = "svg")
}


waterchem_vars <- c("ANC","COND", "MAGNESIUM", "POTASSIUM", "PTL_DISS", "PTL", "SODIUM", "SULFATE", "TURB")
filename <-paste0(environment_label, "-density-plots.svg")
if(!file.exists(filename)){
  
  # Función auxiliar para calcular el sesgo (skewness)
  my_skewness <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 3) return(0)  # Datos insuficientes
    m <- mean(x)
    s <- sd(x)
    if (s == 0) return(0) # Evitar la división por cero
    (sum((x - m)^3) / n) / (s^3)
  }
  
  # Generar una lista de los gráficos
  plot_list <- map(names(env[, selected_vars]), function(var_name) {
    x <- env[[var_name]]
    p <- ggplot(env, aes(x = .data[[var_name]])) +
      geom_density(fill = "steelblue", alpha = 0.7) +
      labs(title = var_name, x = "", y = "") +
      theme_minimal()
    
    # Aplicar transformaciones para las funciones sesgadas a la derecha
    if (min(x, na.rm = TRUE) >= 0) {
      skew_val <- my_skewness(x)
      
      if (skew_val > 1) {  # Sesgo positivo (a la derecha)
        if (min(x, na.rm = TRUE) > 0) {
          p <- p + scale_x_log10() + 
            xlab(paste("log10(", var_name, ")", sep = ""))
        } else {
          p <- p + scale_x_sqrt() + 
            xlab(paste("sqrt(", var_name, ")", sep = ""))
        }
      }
    }else if(var_name %in% waterchem_vars) {p <- p + scale_x_log10() + 
      xlab(paste("log10(", var_name, ")", sep = ""))}  
    return(p)
  })
  
  # Organizar los gráficos en una grilla
  density_plot <- wrap_plots(plot_list, ncol = N_COL)
  
  # Guardar el gráfico
  ggsave(filename, plot = density_plot, height = 14, width = 12, device = "svg")
}
  

original_observed_train <- as.matrix(train_dataset$species)
original_observed_test <- as.matrix(test_dataset$species)

if(PERFORM_RDA || NO_LOG_TRANSFORM_GLMNET){
  observed_train <- as.matrix(train_dataset$species)} else
{ species_log_train <- log((train_dataset$species)+1)
  observed_train <- as.matrix(species_log_train)}

newx_train = as.matrix(train_dataset$env[,colnames(env_selected)])


glmnet_params <- list(
  x = newx_train,
  y = observed_train,
  family = GLMNET_FAMILY,
  nlambda = 10,
  alpha = 0.5,                 # Elastic net (0 = ridge, 1 = lasso)
  nfolds = 4,                  # CV de 4 folds
  maxit = 1e5,                 # por defecto
  thresh = 1e-7                # por defecto
)
glmnet_hash <- digest(glmnet_params)
glmnet_cache_file <- paste0("glmnet_cache_", glmnet_hash, ".rds")

if(!is.null(GLMNET_LINK)){
function_family <- match.fun(glmnet_params$family)
}else{function_family <- NULL}

if (file.exists(glmnet_cache_file)) {
  cached_data <- readRDS(glmnet_cache_file)
  fit <- cached_data$fit
  cat("Se cargó el ajuste GLMNet de la caché\n")
} else {
  fit <- cv.glmnet(
    x = glmnet_params$x,
    y = glmnet_params$y,
    family = if(!is.null(GLMNET_LINK)){function_family(link=GLMNET_LINK)}else{glmnet_params$family},
    nlambda = glmnet_params$nlambda,
    alpha = glmnet_params$alpha,              
    nfolds = glmnet_params$nfolds,             
    maxit = glmnet_params$maxit,         
    thresh = glmnet_params$thresh,    
    trace.it=1
  )
  saveRDS(list(
    fit = fit
  ), file = glmnet_cache_file)
  cat("Se guardó el ajuste GLMNet a la caché\n")
}

# Predecir usando el mejor lambda
preds_train <- predict(fit, newx_train, s = "lambda.min")

if(PERFORM_RDA || NO_LOG_TRANSFORM_GLMNET) {
  # Sin transformación: usar el valor tal cual
  preds_train_orig <- preds_train[, , 1]
} else {
  # Caso de datos transformados, invertir la transformación (exp(x)-1)
  preds_train_orig <- exp(preds_train[, , 1]) - 1
}

compute_r2 <- function(obs, pred) {
  # Si los valores observados son constantes (varianza cero)
  if (sd(obs) < .Machine$double.eps) {
    # Si las predicciones se ajustan al valor observado, R² = 1; caso contrario NaN
    if (sd(pred) < .Machine$double.eps && abs(mean(pred) - mean(obs)) < .Machine$double.eps) {
      return(1)
    } else {
      return(NaN)
    }
  }
  # Caso estándar: usar cor()^2
  cor(obs, pred)^2
}

#  Calcular el R² para cada especie
r2_per_species_train <- sapply(1:ncol(observed_train), function(i) {
  compute_r2(observed_train[, i], preds_train[, i, 1])
})

# Calcular la abundancia total para cada especie
species_abundance_train <- colSums(original_observed_train)

# Calcular el R² promedio ponderado por las abundancias de las especies
weighted_mean_r2_train <- weighted.mean(r2_per_species_train, species_abundance_train, na.rm = TRUE)

print("R² promedio por especie (train)")
print(mean(r2_per_species_train,na.rm = TRUE))  # Mean R² across species
print("R² promedio ponderado por especie (train)")
print(weighted_mean_r2_train)
print("R² mínimo (train)")
print(min(r2_per_species_train,na.rm = TRUE))
print("R² máximo (train)")
print(max(r2_per_species_train,na.rm = TRUE))

# Calcular RMSE/MAE en la escala original
rmse_per_species_train <- sqrt(colMeans((original_observed_train - preds_train_orig)^2))
mae_per_species_train <- colMeans(abs(original_observed_train - preds_train_orig))

print("RMSE total (train)")
print(sum(rmse_per_species_train))
print("MAE total (train)")
print(sum(mae_per_species_train))
print("Abundancia total de las especies (train)")
print(sum(species_abundance_train))

if(PERFORM_RDA){
  observed_test <- as.matrix(test_dataset$species)} else
{ species_log_test <- log(test_dataset$species+1)
  observed_test <- as.matrix(species_log_test)}

newx_test = as.matrix(test_dataset$env[,colnames(env_selected)])

# Predecir usando el mejor lambda
preds_test <- predict(fit, newx_test, s = "lambda.min")
if(PERFORM_RDA || NO_LOG_TRANSFORM_GLMNET) {
  # Sin transformación: usar las predicciones tal cual están
  preds_test_orig <- preds_test[, , 1]
} else {
  # Caso de transformación, calcular la inversa: (exp(x)-1)
  preds_test_orig <- exp(preds_test[, , 1]) - 1
}


# Calcular el R² para cada especie
r2_per_species_test <- sapply(1:ncol(observed_test), function(i) {
  compute_r2(observed_test[, i], preds_test[, i, 1])
})

# Calcular la abundancia total de cada especia
species_abundance_test <- colSums(original_observed_test)

# Calcular el R² ponderado por las abundancias totales de las especies
weighted_mean_r2_test <- weighted.mean(r2_per_species_test, species_abundance_test, na.rm = TRUE)


print("R² promedio (test)")
print(mean(r2_per_species_test,na.rm = TRUE))  # Mean R² across species
print("R² promedio ponderado (test)")
print(weighted_mean_r2_test)
print("R² mínimo (test)")
print(min(r2_per_species_test,na.rm = TRUE))
print("R² máximo (test)")
print(max(r2_per_species_test,na.rm = TRUE))


rmse_per_species_test <- sqrt(colMeans((original_observed_test - preds_test_orig)^2))
mae_per_species_test <- colMeans(abs(original_observed_test - preds_test_orig))

print("RMSE total (test)")
print(sum(rmse_per_species_test))
print("MAE total (test)")
print(sum(mae_per_species_test))
print("Abundancia total de las especies (test)")
print(sum(species_abundance_test))

total_abundances_per_species <- list(train = species_abundance_train, test = species_abundance_test)

glmnet_results <- list()
glmnet_results$train_metrics <- list(r2 = r2_per_species_train, rmse = rmse_per_species_train, mae = mae_per_species_train)
glmnet_results$test_metrics <- list(r2 = r2_per_species_test, rmse = rmse_per_species_test, mae = mae_per_species_test)

if (PERFORM_GAM) {
  env_train <- as.data.frame(train_dataset$env[, colnames(env_selected)])
  env_test <- as.data.frame(test_dataset$env[, colnames(env_selected)])
  
  # Ajustar los modelos GAM con evaluación de métricas en train y test
  
  gam_results <- fit_gam_models_xy(
    species_train = train_dataset$species, 
    env_train = env_train,
    coord_train = train_dataset$coord,
    species_test = test_dataset$species,
    env_test = env_test,
    coord_test = test_dataset$coord,
    family = GAM_FAMILY,
    eval_metrics = GAM_EVAL_METRICS,
    species_label = species_label,
    environment_label = environment_label
  )
  
  species_train <- train_dataset$species
  # Calcular las métricas pesadas
  total_abundance_per_species_train <- colSums(species_train[, !names(species_train) %in% id_columns])
  
  # Métricas de entrenamiento
  weighted_mean_r2_train <- weighted.mean(gam_results$train_metrics$r2, 
                                          total_abundance_per_species_train, 
                                          na.rm = TRUE)
  
  # Imprimir los resultados en train
  print("Resultados de GAM (train):")
  print(paste("R² promedio:", mean(gam_results$train_metrics$r2, na.rm = TRUE)))
  print(paste("R² promedio ponderado:", weighted_mean_r2_train))
  print(paste("RMSE total:", sum(gam_results$train_metrics$rmse, na.rm = TRUE)))
  print(paste("MAE total:", sum(gam_results$train_metrics$mae, na.rm = TRUE)))
  
  species_test <- test_dataset$species
  # Calcular las métricas pesadas
  total_abundance_per_species_test <- colSums(species_test[, !names(species_test) %in% id_columns])
  
  # Mostrar las métricas en test si es necesario
  if(GAM_EVAL_METRICS && !is.null(gam_results$test_metrics)) {
    weighted_mean_r2_test <- weighted.mean(gam_results$test_metrics$r2, 
                                           total_abundance_per_species_test, 
                                           na.rm = TRUE)
    
    print("Resultados del GAM (test):")
    print(paste("R² promedio:", mean(gam_results$test_metrics$r2, na.rm = TRUE)))
    print(paste("R² promedio ponderado:", weighted_mean_r2_test))
    print(paste("RMSE total:", sum(gam_results$test_metrics$rmse, na.rm = TRUE)))
    print(paste("MAE total:", sum(gam_results$test_metrics$mae, na.rm = TRUE)))
  }
  
  gam_results$models <- NULL
  gam_results$cache_hash <- NULL
} else {gam_results = NULL}

# =================================================================
# Selección alternativa con LASSO (siempre se ejecuta)
# =================================================================
# Preparar datos para LASSO

X_lasso <- as.matrix(cbind(env_std, mem.selected))


if(PERFORM_RDA){
  Y_lasso <- as.matrix(species)} else
  { species_log <- log(species+1)
  Y_lasso <- as.matrix(species_log)}

# Configurar factores de penalización (0 para MEMs, 1 para variables ambientales)
penalty_factors <- c(rep(1, ncol(env_std)), rep(0, ncol(mem.selected)))

# Cache para modelo LASSO
lasso_cache_params <- list(
  X = X_lasso,
  Y = Y_lasso,
  penalty_factors = penalty_factors,
  alpha = 1
)
lasso_hash <- digest(lasso_cache_params)
lasso_cache_file <- paste0("lasso_cache_", lasso_hash, ".rds")

if (file.exists(lasso_cache_file)) {
  lasso_fit <- readRDS(lasso_cache_file)
  cat("Se cargó modelo LASSO de caché\n")
} else {
  lasso_fit <- glmnet(
    x = X_lasso, 
    y = Y_lasso, 
    family = "mgaussian", 
    alpha = 1,
    penalty.factor = penalty_factors,
    trace.it=1
  )
  saveRDS(lasso_fit, lasso_cache_file)
  cat("Se guardó modelo LASSO en caché\n")
}

max_n_vars <- length(selected_vars)

# Encontrar lambda que produce ~SELECTED_VARS variables
coef_matrix <- coef(lasso_fit)[[1]]
non_zero_counts <- colSums(coef_matrix[-1, ] != 0)

# Seleccionar lambda con al menos SELECTED_VARS variables
target_idx <- which(non_zero_counts >= max_n_vars)
if (length(target_idx) > 0) {
  chosen_idx <- max(target_idx)  # Más regularizado con suficientes variables
} else {
  chosen_idx <- which.max(non_zero_counts)  # Mejor disponible
}

# Calcular importancia de variables
coef_at_lambda <- coef(lasso_fit, s = lasso_fit$lambda[chosen_idx])
norms <- sapply(1:ncol(env_std), function(j) {
  sqrt(sum(sapply(coef_at_lambda, function(resp_coef) resp_coef[j+1]^2)))
})
  
  # Seleccionar top SELECTED_VARS variables
  selected_vars_lasso_idx <- order(norms, decreasing = TRUE)[1:min(max_n_vars, length(norms))]
  selected_vars_lasso <- colnames(env_std)[selected_vars_lasso_idx]
  
  # Ejecutar CCA con variables seleccionadas por LASSO
  if (ncol(mem.selected) > 0) {
    formula_lasso <- paste(
      "species ~", 
      paste(selected_vars_lasso, collapse = " + "), 
      "+ Condition(", 
      paste(names(mem.selected), collapse = " + "), 
      ")"
    )
  } else {
    formula_lasso <- paste("species ~", paste(selected_vars_lasso, collapse = " + "))
  }

if (PERFORM_RDA){canonical_lasso <- rda(as.formula(formula_lasso), data = cbind(env_std, mem.selected))}
else {canonical_lasso <- cca(as.formula(formula_lasso), data = cbind(env_std, mem.selected))}

print(summary(canonical_lasso))

# Calcular contribuciones de variables para LASSO
loadings_lasso <- scores(canonical_lasso, display = "bp", choices = 1:4, correlation = TRUE, scaling = 0)
eigenvalues_lasso <- canonical_lasso$CCA$eig
var_contrib_lasso <- apply(loadings_lasso^2, 1, function(x) sum(x * eigenvalues_lasso))

# Crear dataframe de resultados para LASSO
selected_loadings_df_lasso <- data.frame(
  COLUMN_NAME = names(var_contrib_lasso),
  loadings_lasso,
  Total_Contribution = var_contrib_lasso,
  Prop_Explained = var_contrib_lasso / sum(canonical_lasso$CCA$eig[1:4])
) %>% 
  left_join(vardesc[, c("COLUMN_NAME", "LABEL")], by = "COLUMN_NAME")

var_selected_glmnet = list(
  total = canonical_lasso$tot.chi,
  constr = canonical_lasso$CCA$tot.chi,
  cond = canonical_lasso$pCCA$tot.chi
)

outliers <- list(isolation = outliers_isof, cca = outliers_cca)
nvars = list(spp = spp_nvars, env = env_nvars, mem = mem_nvars)
nobs = list(train = train_n, test = test_n, total = total_n)
n = list(vars = nvars, obs = nobs)

if(PERFORM_RDA){
  svg(paste0(species_label,"-",environment_label, "-" ,TRANSFORM,"-rda-glmnet.svg"), onefile = FALSE)
  plot(fit)
  
  dev.off()
  
  rda_full <- list(total = rda.res$tot.chi, constr = rda.res$CCA$tot.chi, cond = rda.res$pCCA$tot.chi)
  rda_selected <- list(total = rda_rerun$tot.chi, constr = rda_rerun$CCA$tot.chi, cond = rda_rerun$pCCA$tot.chi)
  
  
  results <- list(is_rda = TRUE, tform = TRANSFORM, spp = species_label, env = environment_label, n = n, outliers = outliers, var_full = rda_full, var_selected = rda_selected, loadings_selected = selected_loadings_df, spp_abund = total_abundances_per_species, glmnet = glmnet_results, gam = gam_results,
                  var_selected_glmnet = var_selected_glmnet,
                  loadings_selected_glmnet = selected_loadings_df_lasso
                  )

  saveRDS(results, paste0(species_label,"-", environment_label, "-",TRANSFORM, "-rda-results.rds"))
} else{
  svg(paste0(species_label,"-",environment_label, "-cca-glmnet.svg"), onefile = FALSE)
  plot(fit)
  
  dev.off()
  
  cca_full <- list(total = cca.res$tot.chi, constr = cca.res$CCA$tot.chi, cond = cca.res$pCCA$tot.chi)
  cca_selected <- list(total = cca_rerun$tot.chi, constr = cca_rerun$CCA$tot.chi, cond = cca_rerun$pCCA$tot.chi)
  
  results <- list(is_rda = FALSE, tform = TRANSFORM, spp = species_label, env = environment_label, n = n, outliers = outliers, var_full = cca_full, var_selected = cca_selected, loadings_selected = selected_loadings_df, spp_abund = total_abundances_per_species, glmnet = glmnet_results, gam = gam_results,
                  var_selected_glmnet = var_selected_glmnet,
                  loadings_selected_glmnet = selected_loadings_df_lasso)

  saveRDS(results, paste0(species_label,"-", environment_label, "-cca-results.rds"))
}


return(results)
}

