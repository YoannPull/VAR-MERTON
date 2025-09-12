# --- NOTE IMPORTANTE ---
# (1) Limite clef : petit échantillon (~30 obs) => faible puissance des tests.
#     "Ne pas rejeter" un test != "validation".
# (2) Model averaging : pas de filtrage par p-value avant l'AIC.
#     On pondère toutes les spécifications estimables, puis on applique les SE robustes
#     sur le modèle choisi (par AIC).

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Préambule : Installation et chargement des packages ---
if (!require("dplyr")) install.packages("dplyr");        library(dplyr)
if (!require("lubridate")) install.packages("lubridate"); library(lubridate)
if (!require("parallel")) install.packages("parallel");   library(parallel)
if (!require("data.table")) install.packages("data.table"); library(data.table)
if (!require("glmnet")) install.packages("glmnet");       library(glmnet)
if (!require("lmtest")) install.packages("lmtest");       library(lmtest)
if (!require("sandwich")) install.packages("sandwich");   library(sandwich)
if (!require("MASS")) install.packages("MASS");           library(MASS)
if (!require("car")) install.packages("car");             library(car)
if (!require("strucchange")) install.packages("strucchange"); library(strucchange)

unregister_dopar <- function() {
  if ("foreach" %in% rownames(installed.packages())) {
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
  }
}

# --- Paramètres de modélisation ---
P_VALUE_THRESHOLD      <- 0.05      # (conservé si tu veux produire un tableau "filtré", mais pas utilisé pour l'AIC)
NB_LAGS                <- 4         # Lags 0..NB_LAGS
MAX_VARIABLES_IN_MODEL <- 4         # Nb max de prédicteurs par modèle
constraint_sign        <- FALSE     # Par défaut, pas de contrainte de signe (on calcule sign_ok mais on ne filtre pas)
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
data_risk_usa <- data_risk_all %>% dplyr::filter(Country == "United States")
Z_result_usa  <- f_Z_estimation(data_risk_usa$Corpo_DR_WA)
dates_risk    <- as.Date(data_risk_usa$Date, format = "%d/%m/%Y")
Y_df          <- data.frame(Date = dates_risk, Y = Z_result_usa$Z)

if (all(is.na(Y_df$Y))) stop("Calcul du facteur Z a échoué.")
cat("Facteur Z calculé avec succès.\n\n")

################################################################################
# ÉTAPE 2 : PRÉPARATION DES VARIABLES EXPLICATIVES
################################################################################
cat("--- Étape 2: Préparation des variables explicatives ---\n")

data_macro_usa <- read.csv("data/processed/data_var_for_model.csv", header = TRUE)
data_macro_usa$Date <- as.Date(data_macro_usa$Date)

# Alignement : mettre à fin de trimestre si besoin ; ici +2 mois comme chez toi
data_macro_usa$Date <- data_macro_usa$Date %m+% months(2)

if (center_data){
  setDT(data_macro_usa)
  num_vars <- setdiff(names(data_macro_usa), "Date")
  data_macro_usa[, (num_vars) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), .SDcols = num_vars]
  setDF(data_macro_usa)
}
cat("Correction de l'alignement des dates macroéconomiques effectuée.\n")

# Ensemble des variables macro testées
vars_to_lag <- c("log_gdp_pc", "log_hours_pc", "log_sp500_real", "log_oil_real", "vix", "log_payems")

lag_indices <- 0:NB_LAGS
lagged_vars_list <- lapply(vars_to_lag, function(var_name) {
  sapply(lag_indices, function(k) dplyr::lag(data_macro_usa[[var_name]], n = k))
})

data_lags <- as.data.frame(do.call(cbind, lagged_vars_list))
names(data_lags) <- paste0(rep(vars_to_lag, each = length(lag_indices)), "_lag", lag_indices)
data_lags$Date   <- data_macro_usa$Date
cat("Variables décalées créées avec succès.\n\n")

################################################################################
# ÉTAPE 3 : MODÉLISATION SUR L'ENSEMBLE DES DONNÉES (+ ETA temps réel)
################################################################################
cat("--- Étape 3: Recherche du meilleur modèle (parallèle) ---\n")

expected_signs <- if (constraint_sign) list(
  # Risque & marchés
  vix              = -1,  # ↑ volatilité => ↓ Z
  log_sp500_real   =  1,  # ↑ equity => ↑ Z
  log_oil_real     = -1,  # ↑ pétrole => ↓ Z
  # Macro réel
  log_inv_pc       =  1,  log_gdp_pc   =  1,  log_gdp      =  1,
  log_hours_pc     =  1,  log_payems   =  1,  unrate       = -1
) else NULL  # sinon NULL (pas de contrainte)

Y_df$Date      <- floor_date(Y_df$Date, unit = "month")
data_lags$Date <- floor_date(data_lags$Date, unit = "month")
model_df       <- inner_join(Y_df, data_lags, by = "Date") %>% na.omit()

if (nrow(model_df) < 20) stop("Pas assez de données valides (< 20) après fusion.")
cat(sprintf("Données prêtes pour la modélisation: %d observations.\n", nrow(model_df)))

# ---- Étape 2bis : DIAGNOSTICS DE COLINEARITÉ SUR LE DESIGN ----
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
  Xs <- scale(X) # centre + scale
  s  <- svd(Xs, nu = 0, nv = 0)$d
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

# Exports colinéarité
output_dir <- "scripts/outputs_USA_Z_Factor"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
colld <- file.path(output_dir, "collinearity")
dir.create(colld, showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(variable = names(vif), VIF = as.numeric(vif))[order(-vif), ],
          file.path(colld, "vif.csv"), row.names = FALSE)
write.csv(tp$top,          file.path(colld, "top_correlations.csv"),  row.names = FALSE)
write.csv(tp$above_thresh, file.path(colld, "corr_above_0p95.csv"),   row.names = FALSE)
write.csv(data.frame(singval=ci$singvals), file.path(colld, "singular_values.csv"), row.names = FALSE)
cat(sprintf("Indice de condition (kappa) = %.2f\n", ci$kappa))
if (ci$kappa > 30)  cat("  Colinéarité forte (kappa > 30)\n")
if (ci$kappa > 100) cat(" Colinéarité sévère (kappa > 100)\n")
cat("Diagnostics de colinéarité écrits dans 'scripts/outputs_USA_Z_Factor/collinearity/'.\n\n")

# ---- Énumération des modèles (sans filtrage par p-value) ----
X_names <- names(model_df)[!names(model_df) %in% c("Date","Y")]
combs <- unlist(lapply(1:min(MAX_VARIABLES_IN_MODEL, length(X_names)), function(k) {
  combn(X_names, k, simplify = FALSE)
}), recursive = FALSE)
total_models <- length(combs)
cat(sprintf("Test de %d combinaisons de modèles...\n", total_models))

# --- Estimation du temps (mono-thread sur un échantillon) ---
num_cores <- max(1, floor(detectCores() * 0.75))  # Mac OK
sample_size <- min(500, total_models)
cat(sprintf("Estimation du temps basée sur %d modèles (mono-thread)...\n", sample_size))
t0 <- Sys.time()
invisible(lapply(combs[1:sample_size], function(vars) {
  formula <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
  try(lm(formula, data = model_df), silent = TRUE)
}))
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

# --- Fonction d'évaluation : AUCUN FILTRE PAR P-VALUE ---
evaluate_combo <- function(vars) {
  f <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
  model <- try(lm(f, data = model_df), silent = TRUE)
  if (inherits(model, "try-error")) return(NULL)
  coefs <- coef(model)[-1]
  # sign check informatif (non bloquant)
  sign_ok <- TRUE
  if (!is.null(expected_signs)) {
    for (nm in names(coefs)) {
      base_var_name <- sub("_lag[0-9]+$", "", nm)
      want <- expected_signs[[base_var_name]]
      if (!is.null(want) && want != 0) {
        if (sign(coefs[[nm]]) != want) { sign_ok <- FALSE; break }
      }
    }
  }
  pvals <- try(summary(model)$coefficients[-1, "Pr(>|t|)"], silent = TRUE)
  if (inherits(pvals, "try-error")) pvals <- rep(NA_real_, length(coefs))
  list(model = model,
       vars  = vars,
       aic   = AIC(model),
       bic   = BIC(model),
       pvals = pvals,
       sign_ok = sign_ok,
       p_ok  = all(!is.na(pvals) & pvals <= P_VALUE_THRESHOLD))
}

# --- Parallélisation par batchs ---
cl <- makeCluster(num_cores)
on.exit({
  try(stopCluster(cl), silent = TRUE)
  unregister_dopar()
}, add = TRUE)

clusterExport(cl, varlist = c("model_df", "expected_signs", "P_VALUE_THRESHOLD"), envir = environment())
clusterEvalQ(cl, { library(stats) })

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
stopCluster(cl); unregister_dopar(); gc()

# --- Agrégation des modèles estimés (non NULL) ---
valid_models <- Filter(Negate(is.null), results)
if (length(valid_models) == 0) stop("Aucun modèle estimable.")

# --- Sélection du meilleur modèle par AIC (SUR TOUS LES MODÈLES) ---
best_model_info <- valid_models[[ which.min(sapply(valid_models, `[[`, "aic")) ]]
best_model <- best_model_info$model
cat(sprintf("Meilleur modèle (AIC global) pour Z (USA) : AIC = %.5f\n", best_model_info$aic))

# --- Sauvegarde du meilleur modèle ---
dir.create("scripts/model", showWarnings = FALSE, recursive = TRUE)
saveRDS(best_model, file = "scripts/model/best_model.rds")
cat("Meilleur modèle sauvegardé dans 'scripts/model/best_model.rds'\n\n")

################################################################################
# ÉTAPE 4 : ROBUSTESSES ADDITIONNELLES
################################################################################
cat("--- Étape 4: Robustesses additionnelles ---\n")

# (A) SE robustes (HAC + HC3)
hac <- coeftest(best_model, vcov = NeweyWest(best_model, lag = NB_LAGS, prewhite = TRUE, adjust = TRUE))
hc3 <- coeftest(best_model, vcov = vcovHC(best_model, type = "HC3"))
write.table(capture.output(hac),  file.path(output_dir, "hac_newey_west.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(capture.output(hc3),  file.path(output_dir, "hc3_robust_se.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("SE robustes HAC/HC3 exportées.\n")

# (A2) Sensibilité HAC aux lags (0..6)
hac_sens <- list()
coef_names <- names(coef(best_model))
for (L in 0:6) {
  tt <- try(coeftest(best_model, vcov = NeweyWest(best_model, lag = L, prewhite = TRUE, adjust = TRUE)), silent = TRUE)
  if (!inherits(tt, "try-error")) {
    df <- as.data.frame(tt)
    df$term <- rownames(df)
    df$lag  <- L
    hac_sens[[length(hac_sens)+1]] <- df[, c("lag","term","t value")]
  }
}
if (length(hac_sens)) {
  hac_sens_df <- do.call(rbind, hac_sens)
  write.csv(hac_sens_df, file.path(output_dir, "hac_sensitivity.csv"), row.names = FALSE)
  cat("Sensibilité HAC aux lags exportée (hac_sensitivity.csv).\n")
}

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

# (C) Model averaging & fréquences d’inclusion (SUR TOUS LES MODÈLES)
tbl <- do.call(rbind, lapply(valid_models, function(x) {
  data.frame(aic = x$aic, bic = x$bic, model = paste(x$vars, collapse = " + "),
             p_ok = x$p_ok, sign_ok = x$sign_ok, stringsAsFactors = FALSE)
}))
tbl$delta <- tbl$aic - min(tbl$aic)
tbl$wAIC  <- exp(-0.5 * tbl$delta); tbl$wAIC <- tbl$wAIC / sum(tbl$wAIC)

vars <- sort(unique(unlist(lapply(valid_models, `[[`, "vars"))))
incl <- sapply(vars, function(v) sum(sapply(valid_models, function(x) v %in% x$vars)))
incl_w <- sapply(vars, function(v) sum(sapply(seq_along(valid_models), function(i) (v %in% valid_models[[i]]$vars) * tbl$wAIC[i])))

incl_df <- data.frame(variable = vars, freq = incl/length(valid_models), weight = incl_w)

# (Option) version "filtrée" si tu veux pour annexe : p_ok &/ou sign_ok
tbl_filt <- subset(tbl, p_ok & sign_ok)
if (nrow(tbl_filt) > 0) {
  vars_f <- sort(unique(unlist(strsplit(tbl_filt$model, " \\+ "))))
  incl_f <- sapply(vars_f, function(v) sum(grepl(paste0("\\b", v, "\\b"), tbl_filt$model)))
  incl_wf <- sapply(vars_f, function(v) sum(tbl_filt$wAIC[grepl(paste0("\\b", v, "\\b"), tbl_filt$model)]))
  incl_df_f <- data.frame(variable = vars_f, freq = incl_f/nrow(tbl_filt), weight = incl_wf)
  write.csv(incl_df_f[order(-incl_df_f$weight), ], file.path(output_dir, "variable_inclusion_filtered.csv"), row.names = FALSE)
}

write.csv(tbl[order(-tbl$wAIC), ], file.path(output_dir, "akaike_weights.csv"), row.names = FALSE)
write.csv(incl_df[order(-incl_df$weight), ], file.path(output_dir, "variable_inclusion.csv"), row.names = FALSE)
cat("Model averaging (poids AIC) exporté.\n")

# (D) Robustesse à l’influence / outliers
inf  <- influence.measures(best_model)
infl_idx <- which(apply(inf$is.inf, 1, any))
write.csv(data.frame(obs = infl_idx), file.path(output_dir, "influential_points.csv"), row.names = FALSE)

fit_rlm <- rlm(formula(best_model), data = model_df, psi = psi.huber)
write.table(capture.output(summary(fit_rlm)),
            file.path(output_dir, "robust_rlm_summary.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Diagnostics d'influence et régression robuste exportés.\n")

# (E) Diagnostics de résidus & spécification
bg1 <- bgtest(best_model, order = 1)
bg2 <- bgtest(best_model, order = 2)
bp  <- bptest(best_model)
rst <- resettest(best_model, power = 2:3)

diag_res <- data.frame(
  test = c("Breusch-Godfrey AR(1)", "Breusch-Godfrey AR(2)", "Breusch-Pagan", "RESET(2:3)"),
  pval = c(bg1$p.value, bg2$p.value, bp$p.value, rst$p.value)
)
write.csv(diag_res, file.path(output_dir, "residual_diagnostics.csv"), row.names = FALSE)
cat("Diagnostics de résidus exportés.\n")

# (F) Stabilité des coefficients : CUSUM
fm  <- formula(best_model)
cus <- efp(fm, data = model_df, type = "Rec-CUSUM")
cus_test <- sctest(cus)
writeLines(c(paste("CUSUM test statistic:", as.numeric(cus_test$statistic)),
             paste("CUSUM p-value:", cus_test$p.value)),
           con = file.path(output_dir, "cusum_test.txt"))
# Plot CUSUM
png(file.path(output_dir, "cusum_plot.png"), width = 900, height = 600)
plot(cus, main = "CUSUM (stabilité des coefficients)")
dev.off()
cat("CUSUM test + plot exportés.\n")

cat("\n--- Analyse terminée ! ---\n")
cat("NOTE : Le principal facteur limitant est la TAILLE DE L'ÉCHANTILLON (~30 obs).\n")
cat("Les tests ont une puissance limitée ; lire les non-rejets avec prudence.\n")
