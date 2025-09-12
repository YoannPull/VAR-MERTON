# ------------------------------------------------------------
# Script complet : Facteur Z (USA) + sélection AIC & sélection OOS
# ------------------------------------------------------------

# --- Préambule : Installation et chargement des packages ---
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("lubridate")) install.packages("lubridate"); library(lubridate)
if (!require("parallel")) install.packages("parallel"); library(parallel)
if (!require("data.table")) install.packages("data.table"); library(data.table)
if (!require("glmnet")) install.packages("glmnet"); library(glmnet)
if (!require("lmtest")) install.packages("lmtest"); library(lmtest)
if (!require("sandwich")) install.packages("sandwich"); library(sandwich)
if (!require("MASS")) install.packages("MASS"); library(MASS)
if (!require("zoo")) install.packages("zoo"); library(zoo)

# Nettoyage sécurisé d'un éventuel environnement foreach (si présent)
unregister_dopar <- function() {
  if (requireNamespace("foreach", quietly = TRUE)) {
    env <- get(".foreachGlobals", envir = asNamespace("foreach"))
    rm(list = ls(name = env), pos = env)
  }
}

# --- Paramètres de modélisation ---
P_VALUE_THRESHOLD      <- 0.05      # Seuil de significativité pour les variables
NB_LAGS                <- 4         # Lags 0..NB_LAGS
MAX_VARIABLES_IN_MODEL <- 6         # Nombre max de prédicteurs par modèle
constraint_sign        <- TRUE
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

load_csv_dt <- function(filename) {
  dt <- fread(file.path(data_dir, filename))
  dt[, DATE := as.Date(observation_date)]
  dt[, observation_date := NULL]
  return(dt)
}

# Conversion de yearqtr en décimal
convert_qtr_to_decimal <- function(qtr) {
  year <- floor(as.numeric(format(qtr, "%Y")))
  qnum <- as.numeric(cycle(qtr))  # 1..4
  return(year + (qnum - 1) * 0.25)
}

data_dir <- "data/raw"

data_risk_usa <- load_csv_dt("default/DRALACBN.csv")
# data_risk_usa <- load_csv_dt("default/CORALACBN.csv")
data_risk_usa[, quarter := as.yearqtr(DATE)]
data_risk_usa <- data_risk_usa[, .(default = DRALACBN), by = quarter]
data_risk_usa[, Date_quarter := convert_qtr_to_decimal(quarter)]
data_risk_usa[, default := default/100]
data_risk_usa[, Date := as.Date(paste0(
  floor(Date_quarter), "-",
  sprintf("%02d", 3 * (Date_quarter %% 1 / 0.25) + 1),
  "-01"
))]

Z_result_usa <- f_Z_estimation(data_risk_usa$default)
Y_df <- data.frame(Date = data_risk_usa$Date,
                   Y = Z_result_usa$Z)

if (all(is.na(Y_df$Y))) stop("Calcul du facteur Z a échoué.")
cat("Facteur Z calculé avec succès.\n\n")

################################################################################
# ÉTAPE 2 : PRÉPARATION DES VARIABLES EXPLICATIVES
################################################################################
cat("--- Étape 2: Préparation des variables explicatives ---\n")

data_macro_usa <- read.csv("data/processed/data_var_for_model.csv", header = TRUE)
data_macro_usa$Date <- as.Date(data_macro_usa$Date)
# data_macro_usa$Date <- data_macro_usa$Date %m+% months(2)

if (center_data) {
  setDT(data_macro_usa)
  num_vars <- setdiff(names(data_macro_usa), "Date")
  data_macro_usa[, (num_vars) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), .SDcols = num_vars]
  setDF(data_macro_usa)
}
cat("Correction de l'alignement des dates macroéconomiques effectuée.\n")

# Ensemble des variables macro testées
vars_to_lag <- c("t10Y2Y","gs2","log_gdp_pc","log_hours_pc","log_sp500_real")

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
    vix              = -1,
    log_sp500_real   =  1,
    log_oil_real     = -1,
    # Taux & spreads
    gs2              =  1,
    t10Y2Y           =  1,
    t10Y3M           =  1,
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
# ÉTAPE 2bis : TEST DE COLINÉARITÉ SUR LE DESIGN
# -------------------------------------------------------------
cat("--- Étape 2bis: Test de colinéarité ---\n")

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
  for (a in 1:(p-1)) for (b in (a+1):p) {
    res <- rbind(res, data.frame(i = cn[a], j = cn[b], corr = C[a,b], stringsAsFactors = FALSE))
  }
  res$abs_corr <- abs(res$corr)
  res <- res[order(-res$abs_corr), ]
  list(top = head(res, top_k), above_thresh = subset(res, abs_corr >= thresh))
}
condition_index <- function(X) {
  Xs <- scale(X)
  s <- svd(Xs, nu = 0, nv = 0)$d
  list(kappa = max(s) / min(s), singvals = s)
}

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

cat("Top 10 |corr| :\n"); print(head(tp$top[, c("i","j","corr")], 10))
if (nrow(tp$above_thresh) > 0) {
  cat(sprintf("Paires |corr| >= 0.95 : %d (voir export)\n", nrow(tp$above_thresh)))
} else cat("Aucune paire |corr| >= 0.95\n")
vif_dt <- data.frame(variable = names(vif), VIF = as.numeric(vif), stringsAsFactors = FALSE)
vif_dt <- vif_dt[order(-vif_dt$VIF), ]
cat(sprintf("Indice de condition (kappa, SVD) = %.2f\n", ci$kappa))
if (ci$kappa > 30)  cat("  Colinéarité forte (kappa > 30)\n")
if (ci$kappa > 100) cat(" Colinéarité sévère (kappa > 100)\n")

output_dir <- "scripts/scripts_robustness/outputs_USA_Z_Factor"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
colld <- file.path(output_dir, "collinearity")
dir.create(colld, showWarnings = FALSE, recursive = TRUE)
write.csv(vif_dt,                        file.path(colld, "vif.csv"),              row.names = FALSE)
write.csv(tp$top,                        file.path(colld, "top_correlations.csv"),  row.names = FALSE)
write.csv(tp$above_thresh,               file.path(colld, "corr_above_0p95.csv"),   row.names = FALSE)
write.csv(data.frame(singval=ci$singvals), file.path(colld, "singular_values.csv"), row.names = FALSE)
cat("Diagnostics de colinéarité écrits dans '.../collinearity/'.\n\n")

# ---------- Utilitaires OOS (rolling expanding, h=1) ----------
rolling_oos_metrics <- function(df, form, initial, h = 1, step = 1, benchmark = c("mean","rw")) {
  benchmark <- match.arg(benchmark)
  n <- nrow(df)
  preds  <- rep(NA_real_, n)
  bench  <- rep(NA_real_, n)
  for (t in seq(initial, n - h, by = step)) {
    fit <- lm(form, data = df[1:t, , drop = FALSE])
    preds[t + h] <- predict(fit, newdata = df[t + h, , drop = FALSE])
    if (benchmark == "mean") bench[t + h] <- mean(df$Y[1:t], na.rm = TRUE) else bench[t + h] <- df$Y[t]
  }
  ok <- !is.na(preds)
  if (!any(ok)) return(list(rmse = Inf, mae = Inf, r2_oos = -Inf, preds = preds, ok = ok))
  err  <- df$Y[ok] - preds[ok]
  rmse <- sqrt(mean(err^2))
  mae  <- mean(abs(err))
  sse_model <- sum(err^2)
  sse_bench <- sum((df$Y[ok] - bench[ok])^2)
  r2_oos <- 1 - sse_model / sse_bench
  list(rmse = rmse, mae = mae, r2_oos = r2_oos, preds = preds, ok = ok)
}

# Fenêtre initiale OOS (>= 50% de l'échantillon ou >= 24 obs)
initial_oos <- max(24, ceiling(nrow(model_df) / 2))

# -------------------------------------------------------------
# Suite Étape 3 : énumération des modèles + ETA
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

# -------- Fonction d'évaluation avec métriques OOS -------- #
evaluate_combo <- function(vars) {
  f <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
  model <- try(lm(f, data = model_df), silent = TRUE)
  if (inherits(model, "try-error")) return(NULL)
  
  p_values <- try(summary(model)$coefficients[-1, "Pr(>|t|)"], silent = TRUE)
  if (inherits(p_values, "try-error") || any(is.na(p_values))) return(NULL)
  if (any(p_values > P_VALUE_THRESHOLD)) return(NULL)
  
  # Vérif signes (si contrainte)
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
  
  # OOS metrics (rolling expanding, h=1)
  oos <- try(rolling_oos_metrics(model_df, f, initial = initial_oos, h = 1, benchmark = "mean"),
             silent = TRUE)
  if (inherits(oos, "try-error")) return(NULL)
  
  list(
    model = model,
    vars  = vars,
    aic   = AIC(model),
    bic   = BIC(model),
    rmse  = oos$rmse,
    mae   = oos$mae,
    r2    = oos$r2_oos
  )
}

# -------- PROGRESS/ETA en temps réel (par batch) -------- #
cl <- makeCluster(num_cores)
on.exit({
  try(stopCluster(cl), silent = TRUE)
  unregister_dopar()
}, add = TRUE)

clusterExport(cl, varlist = c("P_VALUE_THRESHOLD", "model_df", "expected_signs",
                              "rolling_oos_metrics", "initial_oos"))
clusterEvalQ(cl, { library(dplyr); library(lubridate); library(stats) })

BATCH_SIZE <- max(2000, min(25000, ceiling(total_models / (num_cores * 6))))
idx_all    <- seq_len(total_models)
batches    <- split(idx_all, ceiling(idx_all / BATCH_SIZE))

cat(sprintf("Traitement en %d batch(es) de ~%d modèles...\n", length(batches), BATCH_SIZE))
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
gc()

# Sélection du meilleur modèle par AIC
valid_models <- Filter(Negate(is.null), results)
if (length(valid_models) == 0) stop("Aucun modèle trouvé respectant toutes les contraintes.")

best_model_info <- valid_models[[ which.min(sapply(valid_models, `[[`, "aic")) ]]
best_model <- best_model_info$model
cat(sprintf("Meilleur modèle trouvé (AIC) avec un AIC de: %.5f\n\n", best_model_info$aic))

# --- Sauvegarde du meilleur modèle AIC ---
dir.create("scripts/scripts_robustness/model", showWarnings = FALSE, recursive = TRUE)
saveRDS(best_model, file = "scripts/scripts_robustness/model/best_model.rds")
cat("Meilleur modèle (AIC) sauvegardé dans 'best_model.rds'\n")

################################################################################
# ÉTAPE 4 : RÉSULTATS OLS (AIC) — texte et graphique
################################################################################
cat("--- Étape 4: Écriture des résultats et du graphique OLS (AIC) ---\n")
output_dir <- "scripts/scripts_robustness/outputs_USA_Z_Factor"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Résumé texte OLS (AIC)
summary_path <- file.path(output_dir, "model_summary.txt")
sink(summary_path)
cat("Résumé du Meilleur Modèle (sélection AIC) pour le Facteur Z (USA)\n\n")
cat("Param lag =", NB_LAGS, " | MAX VARIABLES =", MAX_VARIABLES_IN_MODEL, "\n\n")
cat("Variables testées :", paste(vars_to_lag, collapse = ", "), "\n\n")
cat("Contrainte de signe :", constraint_sign, "\n\n")
cat("NOTE: Le modèle est estimé sur l'ensemble de l'échantillon disponible.\n\n")
cat("Critère de sélection : AIC le plus faible parmi les modèles valides.\n\n")
cat("AIC du modèle sélectionné :", best_model_info$aic, "\n\n")
print(summary(best_model))
sink()
cat(sprintf("Résumé du modèle (AIC) sauvegardé dans '%s'.\n", summary_path))

# Graphique d'ajustement OLS (AIC)
plot_path <- file.path(output_dir, "model_fit_diagnostic.png")
png(plot_path, width = 1000, height = 600)
fitted_values <- predict(best_model, newdata = model_df)
plot(model_df$Date, model_df$Y, type = "l", col = "grey50", lwd = 2,
     ylab = "Valeur du Facteur Z", xlab = "Date",
     main = "Ajustement du Modèle (sélection AIC) — Facteur Z (USA)",
     sub  = paste("AIC =", round(best_model_info$aic, 3)))
lines(model_df$Date, fitted_values, col = "darkgreen", lwd = 2)
legend("topleft", bty = "n",
       legend = c("Z Observé", "Ajustement du Modèle"),
       col = c("grey50", "darkgreen"), lwd = 2, lty = 1)
dev.off()
cat(sprintf("Graphique d'ajustement (AIC) sauvegardé dans '%s'.\n\n", plot_path))

################################################################################
# ÉTAPE 4bis : SÉLECTION OUT-OF-SAMPLE (RMSE, tie-breaker AIC)
################################################################################
cat("--- Étape 4bis: Sélection OUT-OF-SAMPLE (RMSE) ---\n")

# Tableau des métriques (AIC/BIC + OOS) pour tous les modèles valides
metrics_tbl <- do.call(rbind, lapply(valid_models, function(x) {
  data.frame(
    model = paste(x$vars, collapse = " + "),
    k     = length(x$vars),
    aic   = x$aic,
    bic   = x$bic,
    rmse  = x$rmse,
    mae   = x$mae,
    r2_oos= x$r2,
    stringsAsFactors = FALSE
  )
}))
metrics_tbl <- metrics_tbl[order(metrics_tbl$rmse, metrics_tbl$aic), ]
write.csv(metrics_tbl, file.path(output_dir, "oos_model_metrics.csv"), row.names = FALSE)
cat("Tableau OOS exporté: 'oos_model_metrics.csv'\n")

# Meilleur modèle OOS par RMSE (tie-breaker AIC)
best_oos_idx <- which.min(sapply(valid_models, `[[`, "rmse"))
best_model_oos_info <- valid_models[[best_oos_idx]]
best_model_oos_vars <- best_model_oos_info$vars
best_model_oos_form <- as.formula(paste("Y ~", paste(best_model_oos_vars, collapse = " + ")))

# Refit sur tout l'échantillon pour sauvegarde propre
best_model_oos <- lm(best_model_oos_form, data = model_df)

# Recalcule prédictions OOS détaillées pour export/graphique
oos_best <- rolling_oos_metrics(model_df, best_model_oos_form, initial = initial_oos, h = 1, benchmark = "mean")
cat(sprintf("Meilleur modèle OOS — RMSE=%.4f | MAE=%.4f | R2_OOS=%.3f\n",
            oos_best$rmse, oos_best$mae, oos_best$r2_oos))

# Sauvegardes
dir.create("scripts/scripts_robustness/model", showWarnings = FALSE, recursive = TRUE)
saveRDS(best_model_oos, file = "scripts/scripts_robustness/model/best_model_oos.rds")
cat("Meilleur modèle (sélection OOS) sauvegardé dans 'best_model_oos.rds'\n")

# Résumé du modèle OOS
summary_oos_path <- file.path(output_dir, "model_summary_oos.txt")
sink(summary_oos_path)
cat("Résumé du Meilleur Modèle OUT-OF-SAMPLE pour le Facteur Z (USA)\n\n")
cat("Param lag =", NB_LAGS, " | MAX VARIABLES =", MAX_VARIABLES_IN_MODEL, "\n\n")
cat("Critère de sélection OOS : RMSE (tie-breaker AIC)\n\n")
cat("Variables du modèle OOS :", paste(best_model_oos_vars, collapse = " + "), "\n\n")
cat(sprintf("OOS — RMSE=%.6f | MAE=%.6f | R2_OOS=%.4f\n\n",
            oos_best$rmse, oos_best$mae, oos_best$r2_oos))
print(summary(best_model_oos))
sink()
cat(sprintf("Résumé OOS sauvegardé dans '%s'.\n", summary_oos_path))

# Exports prédictions OOS (rolling) du meilleur modèle OOS
preds_oos_path <- file.path(output_dir, "rolling_oos_predictions_best_model.csv")
write.csv(data.frame(Date = model_df$Date, Y = model_df$Y, Yhat = oos_best$preds),
          preds_oos_path, row.names = FALSE)
cat(sprintf("Prédictions OOS (rolling) exportées dans '%s'.\n", preds_oos_path))

# Graphique OOS (rolling, points prévus)
plot_oos_path <- file.path(output_dir, "oos_best_model_plot.png")
png(plot_oos_path, width = 1000, height = 600)
ok <- oos_best$ok
plot(model_df$Date, model_df$Y, type = "l", col = "grey60", lwd = 2,
     ylab = "Valeur du Facteur Z", xlab = "Date",
     main = "Prédictions OUT-OF-SAMPLE (rolling, h=1) — Meilleur modèle OOS")
if (any(ok)) lines(model_df$Date[ok], oos_best$preds[ok], lwd = 2)
legend("topleft", bty = "n",
       legend = c("Z Observé", "Prédiction OOS (rolling)"),
       col = c("grey60", "black"), lwd = 2, lty = 1)
dev.off()
cat(sprintf("Graphique OOS sauvegardé dans '%s'.\n", plot_oos_path))

# Petite synthèse console sur la dernière prévision OOS disponible
if (any(ok)) {
  last_idx <- max(which(ok))
  cat(sprintf("Dernière prévision OOS : Date=%s | Yhat=%.4f | Y=%.4f | Erreur=%.4f\n",
              format(model_df$Date[last_idx]), oos_best$preds[last_idx],
              model_df$Y[last_idx], model_df$Y[last_idx] - oos_best$preds[last_idx]))
} else {
  cat("Aucune prédiction OOS disponible (vérifier initial_oos vs taille d'échantillon).\n")
}

################################################################################
# ÉTAPE 5 : ROBUSTESSES ADDITIONNELLES (HAC/HC3, OOS pour AIC, Model averaging, Influence)
################################################################################
cat("--- Étape 5: Robustesses additionnelles ---\n")

# (A) SE robustes (HAC + HC3) pour le modèle AIC
hac <- coeftest(best_model, vcov = NeweyWest(best_model, lag = NB_LAGS, prewhite = TRUE, adjust = TRUE))
hc3 <- coeftest(best_model, vcov = vcovHC(best_model, type = "HC3"))
write.table(capture.output(hac),  file.path(output_dir, "hac_newey_west.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(capture.output(hc3),  file.path(output_dir, "hc3_robust_se.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("SE robustes HAC/HC3 exportées.\n")

# (B) Rolling-origin OOS (expanding) pour le modèle AIC
oos_aic <- rolling_oos_metrics(model_df, formula(best_model), initial = initial_oos, h = 1)
cat(sprintf("Rolling OOS AIC (h=1) — RMSE=%.4f | MAE=%.4f | R2_OOS=%.3f\n",
            oos_aic$rmse, oos_aic$mae, oos_aic$r2_oos))
write.csv(data.frame(Date=model_df$Date, Y=model_df$Y, Yhat=oos_aic$preds),
          file.path(output_dir, "rolling_oos_predictions_AIC_model.csv"), row.names = FALSE)

# (C) Model averaging & fréquences d’inclusion (poids AIC)
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

# (D) Robustesse à l’influence / outliers (modèle AIC)
inf <- influence.measures(best_model)
infl_idx <- which(apply(inf$is.inf, 1, any))
write.csv(data.frame(obs = infl_idx), file.path(output_dir, "influential_points.csv"), row.names = FALSE)

fit_rlm <- rlm(formula(best_model), data = model_df, psi = psi.huber)
write.table(capture.output(summary(fit_rlm)),
            file.path(output_dir, "robust_rlm_summary.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Diagnostics d'influence et régression robuste exportés.\n")

cat("--- Analyse terminée ! ---\n")
