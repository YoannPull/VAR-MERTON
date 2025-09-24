# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Préambule : Installation et chargement des packages ---
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("lubridate")) install.packages("lubridate")
library(lubridate)
if (!require("parallel")) install.packages("parallel")
library(parallel)
if (!require("data.table")) install.packages("data.table")
library(data.table)
if (!require("glmnet")) install.packages("glmnet")
library(glmnet)
if (!require("lmtest")) install.packages("lmtest")
library(lmtest)
if (!require("sandwich")) install.packages("sandwich")
library(sandwich)
if (!require("MASS")) install.packages("MASS")
library(MASS)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}

# --- Paramètres de modélisation ---
P_VALUE_THRESHOLD      <- 0.05      # Seuil de significativité pour les variables
NB_LAGS                <- 4         # Lags 0..NB_LAGS
MAX_VARIABLES_IN_MODEL <- 4         # Nombre max de prédicteurs par modèle
constraint_sign        <- FALSE
center_data            <- TRUE

set.seed(12345)

################################################################################
# ÉTAPE 1 : CALCUL DU FACTEUR DE RISQUE SYSTÉMIQUE Z
################################################################################
cat("--- Étape 1: Calcul du Facteur Z pour les USA ---\n")

f_Z_estimation <- function(vect_TD) {
  var_Z_rho <- function(rho) {
    if (rho <= 0 || rho >= 1) return(Inf)
    Z_estim <- (qnorm(mean(vect_TD, na.rm=TRUE)) - qnorm(vect_TD) * sqrt(1 - rho)) / sqrt(rho)
    var(Z_estim, na.rm=TRUE) - 1
  }
  result <- tryCatch(uniroot(var_Z_rho, lower = 1e-6, upper = 1 - 1e-6), error = function(e) NULL)
  if (is.null(result)) return(list(Z = rep(NA, length(vect_TD)), rho = NA))
  rho_opt <- result$root
  Z_estim <- (qnorm(mean(vect_TD, na.rm=TRUE)) - qnorm(vect_TD) * sqrt(1 - rho_opt)) / sqrt(rho_opt)
  list(Z = Z_estim, rho = rho_opt)
}

data_risk_all <- read.csv("data/processed/data_risk_EBA.csv", header = TRUE, sep = ";")
data_risk_usa <- data_risk_all %>% filter(Country == "United States")
Z_result_usa <- f_Z_estimation(data_risk_usa$Corpo_DR_WA)
Z_result_usa$Z <- Z_result_usa$Z - mean(Z_result_usa$Z)
dates_risk <- as.Date(data_risk_usa$Date, format = "%d/%m/%Y")
Y_df <- data.frame(Date = dates_risk, Y = Z_result_usa$Z)

if (all(is.na(Y_df$Y))) stop("Calcul du facteur Z a échoué.")
cat("Facteur Z calculé avec succès.\n\n")

################################################################################
# ÉTAPE 2 : PRÉPARATION DES VARIABLES EXPLICATIVES
################################################################################
cat("--- Étape 2: Préparation des variables explicatives ---\n")

data_macro_usa <- read.csv("data/processed/data_var_for_model.csv", header = TRUE)
data_macro_usa$Date <- as.Date(data_macro_usa$Date)
data_macro_usa$Date <- data_macro_usa$Date %m+% months(2)

if (center_data){
  setDT(data_macro_usa)
  # Centrer toutes les colonnes numériques sauf Date
  num_vars <- setdiff(names(data_macro_usa), "Date")
  data_macro_usa[, (num_vars) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
                 .SDcols = num_vars]
  setDF(data_macro_usa)
}
cat("Correction de l'alignement des dates macroéconomiques effectuée.\n")

# Ensemble des variables macro testées
vars_to_lag <- c("vix","log_sp500_real","log_oil_real","log_hours_pc",
                 "log_gdp_pc")

lag_indices <- 0:NB_LAGS
lagged_vars_list <- lapply(vars_to_lag, function(var_name) {
  sapply(lag_indices, function(k) dplyr::lag(data_macro_usa[[var_name]], n = k))
})

data_lags <- as.data.frame(do.call(cbind, lagged_vars_list))
names(data_lags) <- paste0(rep(vars_to_lag, each = length(lag_indices)), "_lag", lag_indices)
data_lags$Date <- data_macro_usa$Date
cat("Variables décalées créées avec succès.\n\n")

################################################################################
# ÉTAPE 3 : MODÉLISATION SUR L'ENSEMBLE DES DONNÉES (+ ETA temps réel)
################################################################################
cat("--- Étape 3: Recherche du meilleur modèle (parallèle) ---\n")

Y_df$Date      <- floor_date(Y_df$Date, unit = "month")
data_lags$Date <- floor_date(data_lags$Date, unit = "month")
model_df       <- inner_join(Y_df, data_lags, by = "Date") %>% na.omit()

if (nrow(model_df) < 20) stop("Pas assez de données valides (< 20) après fusion.")
cat(sprintf("Données prêtes pour la modélisation: %d observations.\n", nrow(model_df)))

library(glmnet)
X_names <- names(model_df)[!names(model_df) %in% c("Date","Y")]
X = as.matrix(model_df[X_names])
y = as.matrix(Y_df$Y)
lasso_model <- cv.glmnet(X,y, alpha = 1)
# Afficher les résultats
print(lasso_model)
# Coefficients du modèle optimal (selon la validation croisée)
coef(lasso_model, s = "lambda.min")

# Prédictions avec le modèle Lasso optimal (lambda.min)
y_pred <- predict(lasso_model, s = "lambda.min", newx = X)

# Calcul du R²
rss <- sum((y - y_pred)^2)  # Residual Sum of Squares
tss <- sum((y - mean(y))^2)  # Total Sum of Squares
r_squared <- 1 - rss / tss

# Afficher le R²
print(paste("R²: ", round(r_squared, 4)))


# Charger le package ggplot2 pour la visualisation
library(ggplot2)

# Créer un DataFrame avec les valeurs observées, prédites et la date
results_time_df <- data.frame(Date = model_df$Date, 
                              Observed = y, 
                              Predicted = as.vector(y_pred))

# Tracer l'évolution des valeurs observées et prédites au fil du temps
ggplot(results_time_df, aes(x = Date)) +
  geom_line(aes(y = Observed, color = "Observed"), linewidth = 1) +  # valeurs observées
  geom_line(aes(y = Predicted, color = "Predicted"), linewidth = 1, linetype = "dashed") +  # valeurs prédites
  labs(title = "Evolution of Observed vs Predicted Values Over Time", 
       x = "Date", 
       y = "Value") +
  scale_color_manual(name = "Legend", values = c("Observed" = "blue", "Predicted" = "red")) +
  theme_minimal() +
  theme(legend.position = "top")


# Étape 1 : Extraire les coefficients non nuls du modèle Lasso
coefficients_lasso <- coef(lasso_model, s = "lambda.min")

# Convertir les coefficients en un vecteur nommé
coefficients_lasso <- as.vector(coefficients_lasso)
names(coefficients_lasso) <- rownames(coef(lasso_model))

# Retirer l'intercept et les coefficients nuls (représentés par '.')
non_zero_coeffs <- coefficients_lasso[coefficients_lasso != 0]
non_zero_coeffs <- non_zero_coeffs[names(non_zero_coeffs) != "(Intercept)"]

# Classer les coefficients par valeur absolue décroissante
sorted_coeffs <- sort(abs(non_zero_coeffs), decreasing = TRUE)

# Sélectionner les 4 plus significatifs
top_4_vars <- names(sorted_coeffs)[1:4]
cat("Les 4 variables les plus significatives sont :", top_4_vars, "\n")
