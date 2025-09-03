# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Préambule : Installation et chargement des packages ---
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("lubridate")) install.packages("lubridate")
library(lubridate)
if (!require("parallel")) install.packages("parallel")
library(parallel)


unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}

# --- Paramètres de modélisation ---
P_VALUE_THRESHOLD      <- 0.05      # Seuil de significativité pour les variables
NB_LAGS                <- 5         # Lags 0..NB_LAGS
MAX_VARIABLES_IN_MODEL <- 5         # Nombre max de prédicteurs par modèle
constraint_sign <- FALSE

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
cat("Correction de l'alignement des dates macroéconomiques effectuée.\n")

# Ensemble des variables macro testées
vars_to_lag <- c( "log_gdp_pc", "log_hours_pc", "log_sp500_real", "log_oil_real", "vix",
                  "log_payems")

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

if (constraint_sign) {
  expected_signs <- list(
    # Risque & marchés
    vix              = -1,  # ↑ volatilité => ↓ Z
    log_sp500_real   =  1,  # ↑ equity => ↑ Z
    log_oil_real     = -1,  # ↑ pétrole => ↓ Z (coûts/infla/resserrement)
    # Taux & spreads
    gs2              =  1,  # ↑ 2Y => (interprétation pro-cyclique modérée) => ↑ Z
    t10Y2Y           =  1,  # ↑ pente => ↑ Z
    t10Y3M           =  1,  # ↑ pente => ↑ Z
    # Macro réel
    log_inv_pc       =  1,  log_gdp_pc   =  1,  log_gdp      =  1,
    log_hours_pc     =  1,  log_payems   =  1,  unrate       = -1,
    # Cond. financières
    nfci             = -1,  epu          = -1,
    # Inflation
    infl_annualized_pct = 0, infl_yoy_pct = 0
  )
} else {
  expected_signs <- NULL
}

Y_df$Date      <- floor_date(Y_df$Date, unit = "month")
data_lags$Date <- floor_date(data_lags$Date, unit = "month")
model_df       <- inner_join(Y_df, data_lags, by = "Date") %>% na.omit()

if (nrow(model_df) < 20) stop("Pas assez de données valides (< 20) après fusion.")
cat(sprintf("Données prêtes pour la modélisation: %d observations.\n", nrow(model_df)))

# -------------------------------------------------------------
# ÉTAPE 2bis : TEST DE COLINÉARITÉ SUR LE DESIGN (X = toutes var_lag0..NB_LAGS)
# -------------------------------------------------------------
cat("--- Étape 2bis: Test de colinéarité ---\n")

# Fonctions utilitaires
drop_zero_variance <- function(X, tol = 1e-12) {
  sds <- apply(X, 2, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > tol
  list(X = X[, keep, drop = FALSE], dropped = colnames(X)[!keep], kept = colnames(X)[keep])
}
vif_from_matrix <- function(X) {
  p <- ncol(X); out <- numeric(p)
  for (j in 1:p) {
    y <- X[, j]; Z <- X[, -j, drop = FALSE]
    mod <- lm(y ~ Z)
    r2 <- summary(mod)$r.squared
    out[j] <- if (is.finite(r2) && r2 < 1) 1/(1 - r2) else Inf
  }
  setNames(out, colnames(X))
}
top_abs_cor_pairs <- function(C, top_k = 20, thresh = 0.95) {
  cn <- colnames(C); p <- ncol(C)
  res <- data.frame(i = character(), j = character(), corr = numeric(), stringsAsFactors = FALSE)
  for (a in 1:(p-1)) {
    for (b in (a+1):p) {
      res <- rbind(res, data.frame(i = cn[a], j = cn[b], corr = C[a,b], stringsAsFactors = FALSE))
    }
  }
  res$abs_corr <- abs(res$corr)
  res <- res[order(-res$abs_corr), ]
  list(top = head(res, top_k), above_thresh = subset(res, abs_corr >= thresh))
}
condition_index <- function(X) {
  Xs <- scale(X) # centre + scale
  s <- svd(Xs, nu = 0, nv = 0)$d
  list(kappa = max(s) / min(s), singvals = s)
}

# Construction du design et diagnostics
X_names <- names(model_df)[!names(model_df) %in% c("Date","Y")]
X_raw   <- as.matrix(model_df[, X_names, drop = FALSE])
nzv     <- drop_zero_variance(X_raw)
X2      <- nzv$X

if (length(nzv$dropped)) {
  cat("Colonnes à variance nulle supprimées :", paste(nzv$dropped, collapse = ", "), "\n")
}

Xs  <- scale(X2)
C   <- cor(Xs)
tp  <- top_abs_cor_pairs(C, top_k = 40, thresh = 0.8)
vif <- vif_from_matrix(Xs)
ci  <- condition_index(X2)

# Résumés console
cat("Top 10 |corr| :\n"); print(head(tp$top[, c("i","j","corr")], 10))
if (nrow(tp$above_thresh) > 0) {
  cat(sprintf("Paires |corr| >= 0.95 : %d (voir export)\n", nrow(tp$above_thresh)))
} else cat("Aucune paire |corr| >= 0.95\n")
vif_dt <- data.frame(variable = names(vif), VIF = as.numeric(vif), stringsAsFactors = FALSE)
vif_dt <- vif_dt[order(-vif_dt$VIF), ]
cat(sprintf("Indice de condition (kappa, SVD) = %.2f\n", ci$kappa))
if (ci$kappa > 30)  cat("  Colinéarité forte (kappa > 30)\n")
if (ci$kappa > 100) cat(" Colinéarité sévère (kappa > 100)\n")

# Exports
output_dir <- "scripts/outputs_USA_Z_Factor"  # défini aussi plus bas
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
colld <- file.path(output_dir, "collinearity")
dir.create(colld, showWarnings = FALSE, recursive = TRUE)
write.csv(vif_dt,               file.path(colld, "vif.csv"),              row.names = FALSE)
write.csv(tp$top,               file.path(colld, "top_correlations.csv"),  row.names = FALSE)
write.csv(tp$above_thresh,      file.path(colld, "corr_above_0p95.csv"),   row.names = FALSE)
write.csv(data.frame(singval=ci$singvals), file.path(colld, "singular_values.csv"), row.names = FALSE)
cat("Diagnostics de colinéarité écrits dans 'scripts/outputs_USA_Z_Factor/collinearity/'.\n\n")

# -------------------------------------------------------------
# Suite Étape 3 : enumeration des modèles + ETA
# -------------------------------------------------------------
X_names <- names(model_df)[!names(model_df) %in% c("Date","Y")]

# Toutes les combinaisons de 1 à MAX_VARIABLES_IN_MODEL
combs <- unlist(lapply(1:min(MAX_VARIABLES_IN_MODEL, length(X_names)), function(k) {
  combn(X_names, k, simplify = FALSE)
}), recursive = FALSE)
total_models <- length(combs)
cat(sprintf("Test de %d combinaisons de modèles...\n", total_models))

# -------- ESTIMATEUR de temps (warm-up sur 500 modèles max) -------- #
num_cores <- max(1, floor(detectCores() * 0.75))  # Mac OK
sample_size <- min(500, total_models)
cat(sprintf("Estimation du temps basée sur %d modèles (mono-thread)...\n", sample_size))
t0 <- Sys.time()
for (vars in combs[1:sample_size]) {
  formula <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
  invisible(lm(formula, data = model_df))
}
t1 <- Sys.time()
time_per_model <- as.numeric(difftime(t1, t0, units = "secs")) / sample_size
est_time_sec   <- (time_per_model * total_models) / num_cores
fmt_dur <- function(sec) {
  sec <- max(0, as.numeric(sec))
  h <- floor(sec/3600); m <- floor((sec - 3600*h)/60); s <- round(sec - 3600*h - 60*m)
  sprintf("%02d:%02d:%02d", h, m, s)
}
cat(sprintf("Temps moyen par modèle: %.4f sec | Estimation totale: ~%s (sur %d cœurs)\n",
            time_per_model, fmt_dur(est_time_sec), num_cores))

# -------- PROGRESS/ETA en temps réel (par batch) -------- #
evaluate_combo <- function(vars) {
  f <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
  model <- try(lm(f, data = model_df), silent = TRUE)
  if (inherits(model, "try-error")) return(NULL)
  
  p_values <- try(summary(model)$coefficients[-1, "Pr(>|t|)"], silent = TRUE)
  if (inherits(p_values, "try-error") || any(is.na(p_values))) return(NULL)
  if (any(p_values > P_VALUE_THRESHOLD)) return(NULL)
  
  # Vérif signes (si liste fournie)
  if (!is.null(expected_signs)) {
    coefs <- coef(model)[-1]
    for (i in seq_along(coefs)) {
      base_var_name <- sub("_lag[0-9]+$", "", names(coefs)[i])
      want <- expected_signs[[base_var_name]]
      if (!is.null(want) && want != 0) {
        if (sign(coefs[i]) != want) return(NULL)
      }
    }
  }
  
  list(model = model, vars = vars, aic = AIC(model), bic = BIC(model))
}

# Création du cluster
cl <- makeCluster(num_cores)
on.exit({
  try(stopCluster(cl), silent = TRUE)
  unregister_dopar()  # fonction pour nettoyer l'environnement parallèle
}, add = TRUE)
# Exporter les paquets nécessaires dans les cœurs parallèles
clusterExport(cl, varlist = c("P_VALUE_THRESHOLD", "model_df", "expected_signs"))
clusterEvalQ(cl, {
  library(dplyr)
  library(lubridate)
  library(stats)
})

# Découpage en batches pour afficher l’ETA régulièrement (faible overhead)
BATCH_SIZE <- max(2000, min(25000, ceiling(total_models / (num_cores * 6))))
idx_all    <- seq_len(total_models)
batches    <- split(idx_all, ceiling(idx_all / BATCH_SIZE))

cat(sprintf("Traitement en %d batch(es) de ~%d modèles...\n", length(batches), 
            BATCH_SIZE))
flush.console()

start_time <- Sys.time()
count_done <- 0L
results    <- vector("list", total_models)

for (b in seq_along(batches)) {
  ids <- batches[[b]]
  combs_batch <- combs[ids]
  
  t_b0 <- Sys.time()
  ans  <- parLapply(cl, combs_batch, evaluate_combo)
  t_b1 <- Sys.time()
  
  results[ids] <- ans
  count_done   <- count_done + length(ids)
  
  # Progress + ETA basés sur le débit réel
  elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  speed     <- if (elapsed > 0) count_done / elapsed else NA_real_
  remaining <- if (is.finite(speed) && speed > 0) (total_models - count_done) / speed else NA_real_
  pct       <- 100 * count_done / total_models
  
  cat(sprintf("[%.1f%%] %d/%d faits | écoulé: %s | ETA: %s | batch: %.2fs\n",
              pct, count_done, total_models, fmt_dur(elapsed), fmt_dur(remaining),
              as.numeric(difftime(t_b1, t_b0, units = "secs"))))
  flush.console()
}

cat("Calculs en parallèle terminés.\n")
stopCluster(cl)
unregister_dopar() 
gc()  # forcer libération mémoire

# Sélection du meilleur modèle par AIC
valid_models <- Filter(Negate(is.null), results)
if (length(valid_models) == 0) stop("Aucun modèle trouvé respectant toutes les contraintes.")

best_model_info <- valid_models[[ which.min(sapply(valid_models, `[[`, "aic")) ]]
best_model <- best_model_info$model
cat(sprintf("Meilleur modèle trouvé pour le Facteur Z avec un AIC de: %.5f\n\n", best_model_info$aic))

# --- Sauvegarde du meilleur modèle ---
dir.create("scripts/model", showWarnings = FALSE, recursive = TRUE)
saveRDS(best_model, file = "scripts/model/best_model.rds")
cat("Meilleur modèle sauvegardé dans 'best_model.rds'\n")

################################################################################
# ÉTAPE 4 : GÉNÉRATION DES RÉSULTATS
################################################################################
cat("--- Étape 4: Écriture des résultats et du graphique ---\n")
output_dir <- "scripts/outputs_USA_Z_Factor"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Écriture du résumé texte
summary_path <- file.path(output_dir, "model_summary.txt")
sink(summary_path)
cat("Résumé du Meilleur Modèle pour le Facteur Z (USA)\n\n")
cat("Param lag =", NB_LAGS, " | MAX VARIABLES =", MAX_VARIABLES_IN_MODEL, "\n\n")
cat("Variables testées :", paste(vars_to_lag, collapse = ", "), "\n\n")
cat("Contrainte de signe :", constraint_sign, "\n\n")
cat("NOTE: Le modèle est estimé sur l'ensemble de l'échantillon disponible.\n\n")
cat("Critère de sélection : AIC le plus faible parmi les modèles valides.\n\n")
cat("AIC du modèle sélectionné :", best_model_info$aic, "\n\n")
print(summary(best_model))
sink()
cat(sprintf("Résumé du modèle sauvegardé dans '%s'.\n", summary_path))

# Graphique d'ajustement
plot_path <- file.path(output_dir, "model_fit_diagnostic.png")
png(plot_path, width = 1000, height = 600)
fitted_values <- predict(best_model, newdata = model_df)
plot(model_df$Date, model_df$Y, type = "l", col = "grey50", lwd = 2,
     ylab = "Valeur du Facteur Z", xlab = "Date",
     main = "Ajustement du Modèle pour le Facteur Z (USA)",
     sub  = paste("Modèle sélectionné par AIC =", round(best_model_info$aic, 3)))
lines(model_df$Date, fitted_values, col = "darkgreen", lwd = 2)
legend("topleft", bty = "n",
       legend = c("Z Observé", "Ajustement du Modèle"),
       col = c("grey50", "darkgreen"), lwd = 2, lty = 1)
dev.off()
cat(sprintf("Graphique d'ajustement sauvegardé dans '%s'.\n\n", plot_path))

cat("--- Analyse terminée ! ---\n")
