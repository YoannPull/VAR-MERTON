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
NB_LAGS                <- 4         # Lags 0..NB_LAGS
MAX_VARIABLES_IN_MODEL <- 4         # Nombre max de prédicteurs par modèle
constraint_sign        <- FALSE
center_data            <- FALSE

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





################################################################################
# ÉTAPE 3bis : Régressions régularisées (Ridge & Lasso) pour robustesse
################################################################################
cat("--- Étape 3bis: Régressions Ridge & Lasso (robustesse) ---\n")

# --- Packages ---
if (!require("glmnet")) install.packages("glmnet")
library(glmnet)

# --- Switches (activez/désactivez) ---
RUN_RIDGE <- TRUE
RUN_LASSO <- TRUE
USE_LAMBDA <- "lambda.1se"   # "lambda.min" ou "lambda.1se"

set.seed(12345)

# --- Données d'entrée (mêmes que le design précédemment construit) ---
# On repart de model_df déjà créé plus haut
stopifnot(all(c("Date","Y") %in% names(model_df)))
y_vec <- model_df$Y
X_all <- as.matrix(model_df[, setdiff(names(model_df), c("Date","Y")), drop = FALSE])

# Suppression variance nulle (réutilise la fonction définie plus haut)
nzv2      <- drop_zero_variance(X_all)
X_mat     <- nzv2$X
kept_cols <- colnames(X_mat)
if (length(nzv2$dropped)) {
  cat("Ridge/Lasso - colonnes à variance nulle supprimées :",
      paste(nzv2$dropped, collapse = ", "), "\n")
}

# --- Contraintes de signes (optionnel) ---
# On mappe les contraintes de 'expected_signs' (définies plus haut) aux colonnes (avec lags)
build_limits <- function(col_names, expected_signs, use_constraints = FALSE) {
  p <- length(col_names)
  ll <- rep(-Inf, p); ul <- rep(Inf, p)
  if (!use_constraints || is.null(expected_signs)) return(list(ll = ll, ul = ul))
  base_of <- function(nm) sub("_lag[0-9]+$", "", nm)
  for (j in seq_len(p)) {
    base <- base_of(col_names[j])
    want <- expected_signs[[base]]
    if (!is.null(want) && want != 0) {
      if (want ==  1) ll[j] <- 0   # coefficient >= 0
      if (want == -1) ul[j] <- 0   # coefficient <= 0
    }
  }
  list(ll = ll, ul = ul)
}
lims <- build_limits(kept_cols, expected_signs, constraint_sign)

# --- Cross-validation par blocs temporels (respect de l'ordre) ---
make_blocked_foldid <- function(n, K = 5) {
  K <- max(3, min(K, n))  # sécurité
  rep(rep(1:K, each = ceiling(n / K)), length.out = n)
}
K_FOLDS <- 5
foldid  <- make_blocked_foldid(nrow(X_mat), K = K_FOLDS)

# --- Répertoire de sortie ---
output_dir <- "scripts/outputs_USA_Z_Factor"
reg_dir    <- file.path(output_dir, "regularization")
dir.create(reg_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helpers de sortie ---
write_nonzero_coefs <- function(fit_cv, which_lambda, xnames, path_csv) {
  cf <- as.matrix(coef(fit_cv, s = which_lambda))
  out <- data.frame(variable = rownames(cf), coef = as.numeric(cf[,1]), row.names = NULL)
  out$nonzero <- out$coef != 0
  out <- out[order(-abs(out$coef)), ]
  write.csv(out, path_csv, row.names = FALSE)
  invisible(out)
}
compute_fit_stats <- function(y, yhat, df, name) {
  n   <- length(y)
  rss <- sum((y - yhat)^2)
  mse <- rss / n
  rmse <- sqrt(mse)
  aic <- n * log(mse) + 2 * df
  bic <- n * log(mse) + log(n) * df
  data.frame(model = name, n = n, df = df, RMSE = rmse, MSE = mse, RSS = rss, AIC = aic, BIC = bic)
}

metrics_list <- list()

# =========================
# RIDGE (alpha = 0)
# =========================
if (RUN_RIDGE) {
  cat("• Estimation RIDGE (alpha=0)...\n")
  cv_ridge <- cv.glmnet(
    x = X_mat, y = y_vec, family = "gaussian", alpha = 0,
    standardize = TRUE, foldid = foldid,
    lower.limits = lims$ll, upper.limits = lims$ul,
    type.measure = "mse", nfolds = length(unique(foldid))
  )
  saveRDS(cv_ridge, file = file.path(reg_dir, "ridge_cvglmnet.rds"))
  
  # Plot CV
  png(file.path(reg_dir, "ridge_cv_curve.png"), width = 1000, height = 600)
  plot(cv_ridge); title("Ridge CV (MSE)")
  dev.off()
  
  # Coefficients
  ridge_coefs <- write_nonzero_coefs(cv_ridge, USE_LAMBDA, kept_cols,
                                     file.path(reg_dir, paste0("ridge_coefs_", USE_LAMBDA, ".csv")))
  
  # Ajustement in-sample
  yhat_ridge <- as.numeric(predict(cv_ridge, newx = X_mat, s = USE_LAMBDA))
  df_ridge   <- cv_ridge$glmnet.fit$df[which(cv_ridge$glmnet.fit$lambda == get(USE_LAMBDA, cv_ridge))]
  if (length(df_ridge) == 0) df_ridge <- cv_ridge$glmnet.fit$df[which.min(abs(cv_ridge$glmnet.fit$lambda - get(USE_LAMBDA, cv_ridge)))]
  metrics_list[["ridge"]] <- compute_fit_stats(y_vec, yhat_ridge, df = df_ridge, name = paste0("RIDGE_", USE_LAMBDA))
  
  # Plot fit
  png(file.path(reg_dir, paste0("ridge_fit_", USE_LAMBDA, ".png")), width = 1000, height = 600)
  plot(model_df$Date, y_vec, type = "l", lwd = 2, col = "grey50",
       ylab = "Valeur du Facteur Z", xlab = "Date",
       main = paste0("Ajustement Ridge (", USE_LAMBDA, ")"))
  lines(model_df$Date, yhat_ridge, lwd = 2)
  legend("topleft", bty = "n", legend = c("Z Observé", "Ridge (fitted)"),
         lwd = 2, col = c("grey50", "black"))
  dev.off()
}

# =========================
# LASSO (alpha = 1)
# =========================
if (RUN_LASSO) {
  cat("• Estimation LASSO (alpha=1)...\n")
  cv_lasso <- cv.glmnet(
    x = X_mat, y = y_vec, family = "gaussian", alpha = 1,
    standardize = TRUE, foldid = foldid,
    lower.limits = lims$ll, upper.limits = lims$ul,
    type.measure = "mse", nfolds = length(unique(foldid))
  )
  saveRDS(cv_lasso, file = file.path(reg_dir, "lasso_cvglmnet.rds"))
  
  # Plot CV
  png(file.path(reg_dir, "lasso_cv_curve.png"), width = 1000, height = 600)
  plot(cv_lasso); title("Lasso CV (MSE)")
  dev.off()
  
  # Coefficients
  lasso_coefs <- write_nonzero_coefs(cv_lasso, USE_LAMBDA, kept_cols,
                                     file.path(reg_dir, paste0("lasso_coefs_", USE_LAMBDA, ".csv")))
  
  # Ajustement in-sample
  yhat_lasso <- as.numeric(predict(cv_lasso, newx = X_mat, s = USE_LAMBDA))
  df_lasso   <- cv_lasso$glmnet.fit$df[which(cv_lasso$glmnet.fit$lambda == get(USE_LAMBDA, cv_lasso))]
  if (length(df_lasso) == 0) df_lasso <- cv_lasso$glmnet.fit$df[which.min(abs(cv_lasso$glmnet.fit$lambda - get(USE_LAMBDA, cv_lasso)))]
  metrics_list[["lasso"]] <- compute_fit_stats(y_vec, yhat_lasso, df = df_lasso, name = paste0("LASSO_", USE_LAMBDA))
  
  # Plot fit
  png(file.path(reg_dir, paste0("lasso_fit_", USE_LAMBDA, ".png")), width = 1000, height = 600)
  plot(model_df$Date, y_vec, type = "l", lwd = 2, col = "grey50",
       ylab = "Valeur du Facteur Z", xlab = "Date",
       main = paste0("Ajustement Lasso (", USE_LAMBDA, ")"))
  lines(model_df$Date, yhat_lasso, lwd = 2)
  legend("topleft", bty = "n", legend = c("Z Observé", "Lasso (fitted)"),
         lwd = 2, col = c("grey50", "black"))
  dev.off()
}

# --- Export des métriques comparatives ---
if (length(metrics_list)) {
  metrics <- do.call(rbind, metrics_list)
  write.csv(metrics, file.path(reg_dir, "regularization_metrics.csv"), row.names = FALSE)
  print(metrics)
}

cat("Sorties ridge/lasso écrites dans 'scripts/outputs_USA_Z_Factor/regularization/'.\n\n")



################################################################################
# ÉTAPE 3ter : Résumé texte des résultats Ridge/Lasso
################################################################################
cat("--- Étape 3ter: Rédaction d'un résumé texte (regularization_summary.txt) ---\n")

# Sécurité : s'il n'existe pas, reg_dir a été défini plus haut
if (!exists("reg_dir")) {
  output_dir <- "scripts/outputs_USA_Z_Factor"
  reg_dir    <- file.path(output_dir, "regularization")
  dir.create(reg_dir, showWarnings = FALSE, recursive = TRUE)
}

# Helpers
idx_for_lambda <- function(cvfit, lam) which.min(abs(cvfit$lambda - lam))
take_lambda <- function(cvfit, use) if (use == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
nonzero_df <- function(cvfit, s_label){
  cf <- as.matrix(coef(cvfit, s = s_label))
  df <- data.frame(variable = rownames(cf), coef = as.numeric(cf[,1]), row.names = NULL)
  df <- subset(df, variable != "(Intercept)")
  df$nonzero <- df$coef != 0
  df[df$nonzero, c("variable","coef")]
}

# Prépare le fichier
txt_path <- file.path(reg_dir, "regularization_summary.txt")
con <- file(txt_path, open = "wt"); on.exit(close(con), add = TRUE)

wln <- function(...) cat(paste0(..., collapse = ""), "\n", file = con, append = TRUE)

wln("Résumé des régressions régularisées — Facteur Z (USA)")
wln("Généré le : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "")
wln("Paramètres :")
if (exists("NB_LAGS"))                wln("  NB_LAGS = ", NB_LAGS)
if (exists("MAX_VARIABLES_IN_MODEL")) wln("  MAX_VARIABLES_IN_MODEL = ", MAX_VARIABLES_IN_MODEL)
wln("  K_FOLDS = ", if (exists("K_FOLDS")) K_FOLDS else "5",
    " | USE_LAMBDA = ", if (exists("USE_LAMBDA")) USE_LAMBDA else "lambda.1se")
wln("  constraint_sign = ", if (exists("constraint_sign")) constraint_sign else FALSE, "")
wln("Données : n = ", nrow(model_df), " obs ; p (après filtrage variance nulle) = ",
    if (exists("X_mat")) ncol(X_mat) else NA, "")
if (exists("nzv2") && length(nzv2$dropped)) {
  wln("Colonnes supprimées (variance nulle) : ", paste(nzv2$dropped, collapse = ", "), "")
}

# Comparaison avec OLS si dispo
if (exists("best_model")) {
  yhat_ols <- as.numeric(predict(best_model, newdata = model_df))
  rmse_ols <- sqrt(mean((model_df$Y - yhat_ols)^2))
  wln("=== OLS (meilleur modèle par AIC, déjà estimé) ===")
  wln("AIC = ", sprintf("%.3f", AIC(best_model)),
      " | BIC = ", sprintf("%.3f", BIC(best_model)),
      " | RMSE (in-sample) = ", sprintf("%.5f", rmse_ols))
  ols_vars <- setdiff(names(coef(best_model)), "(Intercept)")
  wln("Variables OLS : ", paste(ols_vars, collapse = ", "), "")
}

# Résumé RIDGE
if (exists("RUN_RIDGE") && RUN_RIDGE && exists("cv_ridge")) {
  lam_ridge <- take_lambda(cv_ridge, USE_LAMBDA)
  idx_r     <- idx_for_lambda(cv_ridge, lam_ridge)
  yhat_r    <- as.numeric(predict(cv_ridge, newx = X_mat, s = USE_LAMBDA))
  rmse_r    <- sqrt(mean((y_vec - yhat_r)^2))
  cor_r     <- suppressWarnings(cor(y_vec, yhat_r))
  df_r      <- cv_ridge$glmnet.fit$df[idx_r]  # nb de coef non nuls (glmnet)
  cv_mse_r  <- cv_ridge$cvm[idx_r]
  cv_sd_r   <- cv_ridge$cvsd[idx_r]
  
  nz_r      <- nonzero_df(cv_ridge, USE_LAMBDA)
  nz_r      <- nz_r[order(-abs(nz_r$coef)), , drop = FALSE]
  
  wln("=== RIDGE (alpha = 0) — ", USE_LAMBDA, " ===")
  wln("lambda = ", sprintf("%.6f", lam_ridge),
      " | df (glmnet) = ", df_r,
      " | CV-MSE = ", sprintf("%.6f", cv_mse_r), " ± ", sprintf("%.6f", cv_sd_r),
      " | RMSE (in-sample) = ", sprintf("%.5f", rmse_r),
      " | corr(Y, Ŷ) = ", sprintf("%.4f", cor_r))
  wln("Top coefficients (|coef| décroissant) :")
  if (nrow(nz_r)) {
    topk <- head(nz_r, 20)
    apply(topk, 1, function(row) wln("  ", sprintf("%-28s", row[["variable"]]), " = ", sprintf("% .6f", row[["coef"]])))
  } else {
    wln("  (Ridge donne souvent tous les coefficients non nuls ; cf. CSV détaillé.)")
  }
  wln("")
  
  # Overlap OLS si dispo
  if (exists("best_model")) {
    ols_vars <- setdiff(names(coef(best_model)), "(Intercept)")
    inter_r  <- intersect(ols_vars, nz_r$variable)
    wln("Chevauchement avec OLS (variables) : ", if (length(inter_r)) paste(inter_r, collapse = ", ") else "(aucune)", "")
  }
}

# Résumé LASSO
if (exists("RUN_LASSO") && RUN_LASSO && exists("cv_lasso")) {
  lam_lasso <- take_lambda(cv_lasso, USE_LAMBDA)
  idx_l     <- idx_for_lambda(cv_lasso, lam_lasso)
  yhat_l    <- as.numeric(predict(cv_lasso, newx = X_mat, s = USE_LAMBDA))
  rmse_l    <- sqrt(mean((y_vec - yhat_l)^2))
  cor_l     <- suppressWarnings(cor(y_vec, yhat_l))
  df_l      <- cv_lasso$glmnet.fit$df[idx_l]  # nb de coef non nuls (lasso)
  cv_mse_l  <- cv_lasso$cvm[idx_l]
  cv_sd_l   <- cv_lasso$cvsd[idx_l]
  
  nz_l      <- nonzero_df(cv_lasso, USE_LAMBDA)
  nz_l      <- nz_l[order(-abs(nz_l$coef)), , drop = FALSE]
  
  wln("=== LASSO (alpha = 1) — ", USE_LAMBDA, " ===")
  wln("lambda = ", sprintf("%.6f", lam_lasso),
      " | df (glmnet) = ", df_l,
      " | CV-MSE = ", sprintf("%.6f", cv_mse_l), " ± ", sprintf("%.6f", cv_sd_l),
      " | RMSE (in-sample) = ", sprintf("%.5f", rmse_l),
      " | corr(Y, Ŷ) = ", sprintf("%.4f", cor_l))
  wln("Coefficients non nuls (|coef| décroissant) :")
  if (nrow(nz_l)) {
    topk <- head(nz_l, 30)
    apply(topk, 1, function(row) wln("  ", sprintf("%-28s", row[["variable"]]), " = ", sprintf("% .6f", row[["coef"]])))
  } else {
    wln("  (Aucun coefficient non nul à ce lambda.)")
  }
  wln("")
  
  # Overlap OLS si dispo
  if (exists("best_model")) {
    ols_vars <- setdiff(names(coef(best_model)), "(Intercept)")
    inter_l  <- intersect(ols_vars, nz_l$variable)
    wln("Chevauchement avec OLS (variables) : ", if (length(inter_l)) paste(inter_l, collapse = ", ") else "(aucune)", "")
  }
}

wln("Notes d'interprétation :")
wln("- Les métriques CV (MSE ± écart-type) proviennent d'une validation croisée par BLOCS temporels,")
wln("  ce qui est adapté aux séries temporelles mais ne remplace pas une vraie évaluation out-of-sample (rolling).")
wln("- Les coefficients glmnet sont reportés sur l'échelle originale des variables (pas standardisés dans le fichier).")
wln("- Pour une décision finale, compare aussi les performances hors-échantillon (ex. rolling window/expanding).")

cat(sprintf("Résumé texte sauvegardé dans '%s'.\n\n", txt_path))
