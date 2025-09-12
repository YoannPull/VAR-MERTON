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
center_data            <- FALSE

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
output_dir <- "scripts/outputs_USA_Z_Factor"
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
# ÉTAPE 3bis : Régressions régularisées (Ridge & Lasso) pour robustesse
################################################################################
cat("--- Étape 3bis: Régressions Ridge & Lasso (robustesse) ---\n")

# Données pour glmnet
stopifnot(all(c("Date","Y") %in% names(model_df)))
y_vec <- model_df$Y
X_all <- as.matrix(model_df[, setdiff(names(model_df), c("Date","Y")), drop = FALSE])

nzv2      <- drop_zero_variance(X_all)
X_mat     <- nzv2$X
kept_cols <- colnames(X_mat)
if (length(nzv2$dropped)) {
  cat("Ridge/Lasso - colonnes à variance nulle supprimées :",
      paste(nzv2$dropped, collapse = ", "), "\n")
}

# Contraintes de signes (optionnel)
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

# CV par blocs temporels
make_blocked_foldid <- function(n, K = 5) {
  K <- max(3, min(K, n))
  rep(rep(1:K, each = ceiling(n / K)), length.out = n)
}
K_FOLDS <- 5
foldid  <- make_blocked_foldid(nrow(X_mat), K = K_FOLDS)

# Répertoire
reg_dir    <- file.path(output_dir, "regularization")
dir.create(reg_dir, showWarnings = FALSE, recursive = TRUE)

# Helpers
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
take_lambda <- function(cvfit, use) if (use == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
idx_for_lambda <- function(cvfit, lam) which.min(abs(cvfit$lambda - lam))

metrics_list <- list()
RUN_RIDGE <- TRUE
RUN_LASSO <- TRUE
USE_LAMBDA <- "lambda.1se"   # "lambda.min" ou "lambda.1se"

# ===== RIDGE =====
if (RUN_RIDGE) {
  cat("• Estimation RIDGE (alpha=0)...\n")
  cv_ridge <- cv.glmnet(
    x = X_mat, y = y_vec, family = "gaussian", alpha = 0,
    standardize = TRUE, foldid = foldid,
    lower.limits = lims$ll, upper.limits = lims$ul,
    type.measure = "mse", nfolds = length(unique(foldid))
  )
  saveRDS(cv_ridge, file = file.path(reg_dir, "ridge_cvglmnet.rds"))
  
  png(file.path(reg_dir, "ridge_cv_curve.png"), width = 1000, height = 600)
  plot(cv_ridge); title("Ridge CV (MSE)")
  dev.off()
  
  ridge_coefs <- write_nonzero_coefs(cv_ridge, USE_LAMBDA, kept_cols,
                                     file.path(reg_dir, paste0("ridge_coefs_", USE_LAMBDA, ".csv")))
  
  yhat_ridge <- as.numeric(predict(cv_ridge, newx = X_mat, s = USE_LAMBDA))
  lam_r <- take_lambda(cv_ridge, USE_LAMBDA)
  idx_r <- idx_for_lambda(cv_ridge, lam_r)
  df_ridge <- cv_ridge$glmnet.fit$df[idx_r]
  metrics_list[["ridge"]] <- compute_fit_stats(y_vec, yhat_ridge, df = df_ridge, name = paste0("RIDGE_", USE_LAMBDA))
  
  png(file.path(reg_dir, paste0("ridge_fit_", USE_LAMBDA, ".png")), width = 1000, height = 600)
  plot(model_df$Date, y_vec, type = "l", lwd = 2, col = "grey50",
       ylab = "Valeur du Facteur Z", xlab = "Date",
       main = paste0("Ajustement Ridge (", USE_LAMBDA, ")"))
  lines(model_df$Date, yhat_ridge, lwd = 2)
  legend("topleft", bty = "n", legend = c("Z Observé", "Ridge (fitted)"),
         lwd = 2, col = c("grey50", "black"))
  dev.off()
}

# ===== LASSO =====
if (RUN_LASSO) {
  cat("• Estimation LASSO (alpha=1)...\n")
  cv_lasso <- cv.glmnet(
    x = X_mat, y = y_vec, family = "gaussian", alpha = 1,
    standardize = TRUE, foldid = foldid,
    lower.limits = lims$ll, upper.limits = lims$ul,
    type.measure = "mse", nfolds = length(unique(foldid))
  )
  saveRDS(cv_lasso, file = file.path(reg_dir, "lasso_cvglmnet.rds"))
  
  png(file.path(reg_dir, "lasso_cv_curve.png"), width = 1000, height = 600)
  plot(cv_lasso); title("Lasso CV (MSE)")
  dev.off()
  
  lasso_coefs <- write_nonzero_coefs(cv_lasso, USE_LAMBDA, kept_cols,
                                     file.path(reg_dir, paste0("lasso_coefs_", USE_LAMBDA, ".csv")))
  
  yhat_lasso <- as.numeric(predict(cv_lasso, newx = X_mat, s = USE_LAMBDA))
  lam_l <- take_lambda(cv_lasso, USE_LAMBDA)
  idx_l <- idx_for_lambda(cv_lasso, lam_l)
  df_lasso <- cv_lasso$glmnet.fit$df[idx_l]
  metrics_list[["lasso"]] <- compute_fit_stats(y_vec, yhat_lasso, df = df_lasso, name = paste0("LASSO_", USE_LAMBDA))
  
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
print_top <- function(df, top_n = 20) {
  if (is.null(df) || !nrow(df)) { wln("  (aucun)"); return(invisible(NULL)) }
  df2 <- df[order(-abs(df$coef)), , drop = FALSE]
  df2 <- head(df2, top_n)
  for (i in seq_len(nrow(df2))) {
    v <- as.character(df2$variable[i])
    cval <- suppressWarnings(as.numeric(df2$coef[i]))
    if (is.na(cval)) {
      wln("  ", sprintf("%-28s", v), " = NA")
    } else {
      wln("  ", sprintf("%-28s", v), " = ", sprintf("% .6f", cval))
    }
  }
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

# Résumé RIDGE @ lambda.1se
if (exists("RUN_RIDGE") && RUN_RIDGE && exists("cv_ridge")) {
  lam_ridge <- take_lambda(cv_ridge, "lambda.1se")
  idx_r     <- idx_for_lambda(cv_ridge, lam_ridge)
  yhat_r    <- as.numeric(predict(cv_ridge, newx = X_mat, s = "lambda.1se"))
  rmse_r    <- sqrt(mean((y_vec - yhat_r)^2))
  cor_r     <- if (sd(yhat_r) > 0) suppressWarnings(cor(y_vec, yhat_r)) else NA
  df_r      <- cv_ridge$glmnet.fit$df[idx_r]
  cv_mse_r  <- cv_ridge$cvm[idx_r]
  cv_sd_r   <- cv_ridge$cvsd[idx_r]
  nz_r      <- nonzero_df(cv_ridge, "lambda.1se")
  
  wln("=== RIDGE (alpha = 0) — lambda.1se ===")
  wln("lambda = ", sprintf("%.6f", lam_ridge),
      " | df (glmnet) = ", df_r,
      " | CV-MSE = ", sprintf("%.6f", cv_mse_r), " ± ", sprintf("%.6f", cv_sd_r),
      " | RMSE (in-sample) = ", sprintf("%.5f", rmse_r),
      " | corr(Y, Ŷ) = ", ifelse(is.na(cor_r), "NA", sprintf("%.4f", cor_r)))
  wln("Top coefficients (|coef| décroissant) :"); print_top(nz_r, top_n = 20); wln("")
  if (exists("best_model")) {
    ols_vars <- setdiff(names(coef(best_model)), "(Intercept)")
    inter_r  <- intersect(ols_vars, nz_r$variable)
    wln("Chevauchement avec OLS (variables) : ", if (length(inter_r)) paste(inter_r, collapse = ", ") else "(aucune)", "")
  }
}

# Résumé LASSO @ lambda.1se
if (exists("RUN_LASSO") && RUN_LASSO && exists("cv_lasso")) {
  lam_lasso <- take_lambda(cv_lasso, "lambda.1se")
  idx_l     <- idx_for_lambda(cv_lasso, lam_lasso)
  yhat_l    <- as.numeric(predict(cv_lasso, newx = X_mat, s = "lambda.1se"))
  rmse_l    <- sqrt(mean((y_vec - yhat_l)^2))
  cor_l     <- if (sd(yhat_l) > 0) suppressWarnings(cor(y_vec, yhat_l)) else NA
  df_l      <- cv_lasso$glmnet.fit$df[idx_l]
  cv_mse_l  <- cv_lasso$cvm[idx_l]
  cv_sd_l   <- cv_lasso$cvsd[idx_l]
  nz_l      <- nonzero_df(cv_lasso, "lambda.1se")
  
  wln("=== LASSO (alpha = 1) — lambda.1se ===")
  wln("lambda = ", sprintf("%.6f", lam_lasso),
      " | df (glmnet) = ", df_l,
      " | CV-MSE = ", sprintf("%.6f", cv_mse_l), " ± ", sprintf("%.6f", cv_sd_l),
      " | RMSE (in-sample) = ", sprintf("%.5f", rmse_l),
      " | corr(Y, Ŷ) = ", ifelse(is.na(cor_l), "NA", sprintf("%.4f", cor_l)))
  wln("Coefficients non nuls (|coef| décroissant) :"); print_top(nz_l, top_n = 30); wln("")
  if (exists("best_model")) {
    ols_vars <- setdiff(names(coef(best_model)), "(Intercept)")
    inter_l  <- intersect(ols_vars, nz_l$variable)
    wln("Chevauchement avec OLS (variables) : ", if (length(inter_l)) paste(inter_l, collapse = ", ") else "(aucune)", "")
  }
}

# --- RIDGE @ lambda.min ---
if (exists("cv_ridge")) {
  lam_r_min <- cv_ridge$lambda.min
  yhat_r_min <- as.numeric(predict(cv_ridge, newx = X_mat, s = "lambda.min"))
  rmse_r_min <- sqrt(mean((y_vec - yhat_r_min)^2))
  cor_r_min  <- if (sd(yhat_r_min) > 0) cor(y_vec, yhat_r_min) else NA
  idx_r_min  <- which.min(abs(cv_ridge$lambda - lam_r_min))
  df_r_min   <- cv_ridge$glmnet.fit$df[idx_r_min]
  cv_mse_r_min <- cv_ridge$cvm[idx_r_min]; cv_sd_r_min <- cv_ridge$cvsd[idx_r_min]
  
  wln("=== RIDGE (alpha = 0) — lambda.min ===")
  wln("lambda = ", sprintf("%.6f", lam_r_min),
      " | df = ", df_r_min,
      " | CV-MSE = ", sprintf("%.6f", cv_mse_r_min), " ± ", sprintf("%.6f", cv_sd_r_min),
      " | RMSE (in-sample) = ", sprintf("%.5f", rmse_r_min),
      " | corr(Y, Ŷ) = ", ifelse(is.na(cor_r_min), "NA", sprintf("%.4f", cor_r_min)))
  wln("")
}

# --- LASSO @ lambda.min ---
if (exists("cv_lasso")) {
  lam_l_min <- cv_lasso$lambda.min
  yhat_l_min <- as.numeric(predict(cv_lasso, newx = X_mat, s = "lambda.min"))
  rmse_l_min <- sqrt(mean((y_vec - yhat_l_min)^2))
  cor_l_min  <- if (sd(yhat_l_min) > 0) cor(y_vec, yhat_l_min) else NA
  idx_l_min  <- which.min(abs(cv_lasso$lambda - lam_l_min))
  df_l_min   <- cv_lasso$glmnet.fit$df[idx_l_min]
  cv_mse_l_min <- cv_lasso$cvm[idx_l_min]; cv_sd_l_min <- cv_lasso$cvsd[idx_l_min]
  
  wln("=== LASSO (alpha = 1) — lambda.min ===")
  wln("lambda = ", sprintf("%.6f", lam_l_min),
      " | df = ", df_l_min,
      " | CV-MSE = ", sprintf("%.6f", cv_mse_l_min), " ± ", sprintf("%.6f", cv_sd_l_min),
      " | RMSE (in-sample) = ", sprintf("%.5f", rmse_l_min),
      " | corr(Y, Ŷ) = ", ifelse(is.na(cor_l_min), "NA", sprintf("%.4f", cor_l_min)))
  wln("")
}

# --- OLS en CV sur les mêmes folds ---
if (exists("best_model")) {
  # Si foldid/X_mat n'existent pas, on les reconstruit
  if (!exists("foldid")) foldid <- make_blocked_foldid(nrow(model_df), 5)
  uniq_folds <- sort(unique(foldid))
  ols_form <- formula(best_model)
  cv_mses <- numeric(length(uniq_folds))
  for (ii in seq_along(uniq_folds)) {
    k <- uniq_folds[ii]
    test_idx  <- which(foldid == k)
    train_idx <- setdiff(seq_len(nrow(model_df)), test_idx)
    fit_k <- try(lm(ols_form, data = model_df[train_idx, ]), silent = TRUE)
    if (inherits(fit_k, "try-error")) { cv_mses[ii] <- NA; next }
    pred_k <- predict(fit_k, newdata = model_df[test_idx, ])
    cv_mses[ii] <- mean((model_df$Y[test_idx] - pred_k)^2)
  }
  ols_cv_mse <- mean(cv_mses, na.rm = TRUE)
  ols_cv_sd  <- sd(cv_mses, na.rm = TRUE)
  
  wln("=== OLS — CV (mêmes folds que glmnet) ===")
  wln("CV-MSE = ", sprintf("%.6f", ols_cv_mse), " ± ", sprintf("%.6f", ols_cv_sd))
  wln("Note: cette CV peut inclure du 'futur' dans l'entraînement ; préférer un rolling-origin pour du time-safe.")
  wln("")
}

wln("Notes d'interprétation :")
wln("- Les métriques CV (MSE ± écart-type) proviennent d'une validation croisée par BLOCS temporels,")
wln("  ce qui est adapté aux séries temporelles mais ne remplace pas une vraie évaluation out-of-sample (rolling).")
wln("- Les coefficients glmnet sont reportés sur l'échelle originale des variables (pas standardisés dans le fichier).")
wln("- Pour une décision finale, compare aussi les performances hors-échantillon (ex. rolling window/expanding).")

cat(sprintf("Résumé texte sauvegardé dans '%s'.\n\n", txt_path))

################################################################################
# ÉTAPE 4 : GÉNÉRATION DES RÉSULTATS OLS (texte et graphique)
################################################################################
cat("--- Étape 4: Écriture des résultats et du graphique OLS ---\n")
output_dir <- "scripts/outputs_USA_Z_Factor"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Écriture du résumé texte OLS
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

# Graphique d'ajustement OLS
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

################################################################################
# ÉTAPE 5 : ROBUSTESSES ADDITIONNELLES (HAC/HC3, Rolling OOS, Model averaging, Influence)
################################################################################
cat("--- Étape 5: Robustesses additionnelles ---\n")

# (A) SE robustes (HAC + HC3)
hac <- coeftest(best_model, vcov = NeweyWest(best_model, lag = NB_LAGS, prewhite = TRUE, adjust = TRUE))
hc3 <- coeftest(best_model, vcov = vcovHC(best_model, type = "HC3"))
write.table(capture.output(hac),  file.path(output_dir, "hac_newey_west.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(capture.output(hc3),  file.path(output_dir, "hc3_robust_se.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("SE robustes HAC/HC3 exportées.\n")

# (B) Rolling-origin OOS (expanding)
rolling_oos <- function(df, form, initial, h = 1, step = 1) {
  n <- nrow(df); preds <- rep(NA_real_, n)
  for (t in seq(initial, n - h, by = step)) {
    fit <- lm(form, data = df[1:t, ])
    preds[t + h] <- predict(fit, newdata = df[t + h, , drop = FALSE])
  }
  ok <- !is.na(preds)
  rmse <- sqrt(mean((df$Y[ok] - preds[ok])^2))
  mae  <- mean(abs(df$Y[ok] - preds[ok]))
  list(rmse = rmse, mae = mae, preds = preds)
}
oos <- rolling_oos(model_df, formula(best_model), initial = max(24, ceiling(nrow(model_df)/2)), h = 1)
cat(sprintf("Rolling OOS (h=1) — RMSE=%.4f | MAE=%.4f\n", oos$rmse, oos$mae))
write.csv(data.frame(Date=model_df$Date, Y=model_df$Y, Yhat=oos$preds),
          file.path(output_dir, "rolling_oos_predictions.csv"), row.names = FALSE)

# (C) Model averaging & fréquences d’inclusion
vm <- Filter(Negate(is.null), results)
if (length(vm) > 0) {
  tbl <- do.call(rbind, lapply(vm, function(x) data.frame(aic = x$aic, model = paste(x$vars, collapse = " + "))))
  tbl$delta <- tbl$aic - min(tbl$aic)
  tbl$wAIC  <- exp(-0.5 * tbl$delta); tbl$wAIC <- tbl$wAIC / sum(tbl$wAIC)
  vars <- sort(unique(unlist(lapply(vm, `[[`, "vars"))))
  incl <- sapply(vars, function(v) sum(sapply(vm, function(x) v %in% x$vars)))
  incl_w <- sapply(vars, function(v) sum(sapply(seq_along(vm), function(i) (v %in% vm[[i]]$vars) * tbl$wAIC[i])))
  incl_df <- data.frame(variable = vars, freq = incl/length(vm), weight = incl_w)
  write.csv(tbl[order(-tbl$wAIC), ], file.path(output_dir, "akaike_weights.csv"), row.names = FALSE)
  write.csv(incl_df[order(-incl_df$weight), ], file.path(output_dir, "variable_inclusion.csv"), row.names = FALSE)
  cat("Model averaging (poids AIC) exporté.\n")
} else {
  cat("Aucun modèle valide pour le model averaging.\n")
}

# (D) Robustesse à l’influence / outliers
inf <- influence.measures(best_model)
infl_idx <- which(apply(inf$is.inf, 1, any))
write.csv(data.frame(obs = infl_idx), file.path(output_dir, "influential_points.csv"), row.names = FALSE)

fit_rlm <- rlm(formula(best_model), data = model_df, psi = psi.huber)
write.table(capture.output(summary(fit_rlm)),
            file.path(output_dir, "robust_rlm_summary.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Diagnostics d'influence et régression robuste exportés.\n")

cat("--- Analyse terminée ! ---\n")