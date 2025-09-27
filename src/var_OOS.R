suppressPackageStartupMessages({
  library(data.table)
  library(vars)
  library(urca)
})

##======================= Fonctions fournies =======================##
# Ton test Johansen (inchangé, encapsulé ci-dessous)
test_johansen <- function(data, variables) {
  # Sous-ensemble des variables pour cette combinaison
  df <- data[,-1] # do not take the key
  df_subset <- df[, ..variables]
  
  # Test de Johansen (type "trace" et 2 lags, à ajuster si nécessaire)
  johan_test <- tryCatch({
    ca.jo(df_subset, type = "trace", ecdet = "none", K = 2)
  }, error = function(e) {
    NULL  # Retourner NULL en cas d'erreur (par exemple, si trop de variables)
  })
  
  # Si le test a réussi et qu'il y a des relations de cointégration, on retourne TRUE
  if (!is.null(johan_test)) {
    # Résumé du test de cointégration
    summary_johan <- summary(johan_test)
    
    # Vérification de la cointégration
    if (summary_johan@teststat[1] > summary_johan@cval[1,3]) {  # Compare statistique et valeur critique
      return(list(variables = variables, result = summary_johan))
    }
  }
  return(NULL)  # Aucun résultat de cointégration trouvé
}

##======================= Aides robustes VAR =======================##
drop_constant_cols <- function(dt, tol = 1e-12) {
  keep <- vapply(dt, function(x) {
    x <- as.numeric(x)
    if (all(!is.finite(x))) return(FALSE)
    v <- var(x, na.rm = TRUE)
    is.finite(v) && v > tol
  }, logical(1))
  dt[, which(keep), with = FALSE]
}

prep_dt_for_var <- function(dt_sub) {
  # nettoie NA/Inf et colonnes constantes
  dt_clean <- copy(dt_sub)
  M <- as.matrix(dt_clean)
  bad <- !is.finite(M)
  if (any(bad)) {
    # on retire les lignes contenant NA/Inf
    dt_clean <- dt_clean[rowSums(!is.finite(as.matrix(dt_clean))) == 0]
  }
  dt_clean <- na.omit(dt_clean)
  dt_clean <- drop_constant_cols(dt_clean)
  dt_clean
}

calc_ic_mv <- function(fit, p, det_terms_per_eq = 1L) {
  # IC multivariés basés sur ΣU (fallback si logLik indisponible)
  SigmaU <- tryCatch(fit$SigmaU, error = function(e) NULL)
  if (is.null(SigmaU)) return(list(bic = NA_real_, aic = NA_real_))
  detS <- tryCatch(det(SigmaU), error = function(e) NA_real_)
  if (!is.finite(detS) || detS <= 0) return(list(bic = NA_real_, aic = NA_real_))
  T_eff <- fit$obs
  k     <- length(fit$varresult)
  q     <- k * k * p + k * det_terms_per_eq
  aic_mv <- log(detS) + (2 * q) / T_eff
  bic_mv <- log(detS) + (log(T_eff) * q) / T_eff
  list(bic = bic_mv, aic = aic_mv)
}

safe_logLik_BIC_AIC <- function(fit) {
  ll <- tryCatch(logLik(fit), error = function(e) NA)
  if (!is.finite(ll)) return(list(bic = NA_real_, aic = NA_real_))
  npar <- attr(ll, "df")
  nobs <- fit$obs
  bic  <- -2 * as.numeric(ll) + log(nobs) * npar
  aic  <- -2 * as.numeric(ll) + 2 * npar
  list(bic = bic, aic = aic)
}

enough_sample <- function(T_eff, k, p, det_terms_per_eq = 1L) {
  # Règle de pouce pour éviter les fits trop paramétrés
  T_eff > (k * p + det_terms_per_eq) * 4
}

safe_test_johansen_on_levels <- function(data, varset) {
  lvl_vars <- varset[!grepl("^diff", varset)]
  if (length(lvl_vars) >= 2) return(test_johansen(data, lvl_vars))
  NULL
}

##======================= Contrôles de base =======================##
stopifnot(exists("df_train"))
stopifnot("Date_quarter" %in% names(df_train))
# s'assurer que Date_quarter est en 1re position (pour test_johansen)
setcolorder(df_train, c("Date_quarter", setdiff(names(df_train), "Date_quarter")))

##=================== Construction du pool de vars ================##
all_cols <- setdiff(names(df_train), "Date_quarter")

# Contraintes :
keep_mandatory <- "diff1_log_GPRD"  # doit être présent
if (!keep_mandatory %in% all_cols) stop("'diff1_log_GPRD' est absent de df_train.")

ban_exact <- c("cpi")               # ne jamais prendre 'cpi' (mais log_cpi/diff1_log_cpi ok)
gpr_cols  <- grep("GPR", all_cols, value = TRUE) # toutes colonnes GPR*

# Retirer toutes les GPR* sauf 'diff1_log_GPRD', et retirer 'cpi'
allowed_pool <- setdiff(all_cols, union(ban_exact, setdiff(gpr_cols, keep_mandatory)))
# garantir présence du mandatory
allowed_pool <- unique(c(keep_mandatory, setdiff(allowed_pool, keep_mandatory)))

##================ Paramètres de recherche de modèles ==============##
min_k <- 3L
max_k <- 6L
set.seed(42)
max_models <- 400L   # échantillonne si trop de combinaisons

# (Optionnel) interdire de mélanger niveau/log/diff d'une même "famille"
enforce_family_unique <- FALSE
base_of <- function(v) {
  # retire préfixes diff1_/diff2_/log_
  gsub("^log_", "", gsub("^diff[12]_", "", v))
}
family_ok <- function(vars) {
  if (!enforce_family_unique) return(TRUE)
  bases <- sapply(vars, base_of)
  # interdit d'avoir >1 forme pour une même base
  all(tabulate(factor(bases)) <= 1)
}

##===================== Génération des combinaisons ================##
base_rest <- setdiff(allowed_pool, keep_mandatory)

gen_combos <- function(rest, k, mandatory) {
  if ((k - 1L) > length(rest)) return(list())
  cmb <- combn(rest, k - 1L, simplify = FALSE)
  combos <- lapply(cmb, function(x) sort(c(mandatory, x)))
  # filtre familles si activé
  combos[vapply(combos, family_ok, logical(1))]
}

all_combos <- list()
for (k in min_k:max_k) {
  all_combos <- c(all_combos, gen_combos(base_rest, k, keep_mandatory))
}
# déduplication
all_combos <- unique(all_combos)

# échantillonnage si trop de modèles
if (length(all_combos) > max_models) {
  all_combos <- sample(all_combos, max_models)
}

cat(sprintf("Nombre de combinaisons à tester: %d\n", length(all_combos)))

##===================== Boucle d'évaluation VAR ====================##
results <- vector("list", length(all_combos))
names(results) <- sapply(all_combos, function(v) paste(v, collapse = ","))

i <- 0L
for (vars_set in all_combos) {
  i <- i + 1L
  
  # sous-dataset (sans la date)
  dt_sub_raw <- df_train[, ..vars_set]
  dt_sub <- prep_dt_for_var(dt_sub_raw)
  
  if (ncol(dt_sub) < 2L || nrow(dt_sub) < 30L) {
    results[[i]] <- list(vars = vars_set, k = ncol(dt_sub_raw), p = NA_integer_,
                         bic = NA_real_, aic = NA_real_, has_coint = NA,
                         note = "too_few_vars_or_obs_after_clean")
    next
  }
  
  kvars <- ncol(dt_sub)
  # lag max raisonnable selon l'échantillon nettoyé
  lag_max <- min(6L, max(1L, floor((nrow(dt_sub) - 20L) / max(1L, kvars))))
  lag_max <- max(lag_max, 2L)  # au moins 2
  
  # Sélection par BIC
  sel <- tryCatch(VARselect(as.data.frame(dt_sub), lag.max = lag_max, type = "const"),
                  error = function(e) NULL)
  if (is.null(sel)) {
    results[[i]] <- list(vars = vars_set, k = kvars, p = NA_integer_,
                         bic = NA_real_, aic = NA_real_, has_coint = NA,
                         note = "VARselect_failed")
    next
  }
  
  p_sc <- sel$selection["SC(n)"]
  if (is.na(p_sc) || !is.finite(p_sc)) p_sc <- which.min(sel$SC) - 1L
  p_sc <- as.integer(max(1L, p_sc))
  
  # Check capacité (simple règle de pouce)
  if (!enough_sample(nrow(dt_sub), kvars, p_sc, det_terms_per_eq = 1L)) {
    results[[i]] <- list(vars = vars_set, k = kvars, p = p_sc,
                         bic = NA_real_, aic = NA_real_, has_coint = NA,
                         note = "insufficient_sample_for_params")
    next
  }
  
  fit <- tryCatch(VAR(as.data.frame(dt_sub), p = p_sc, type = "const"),
                  error = function(e) NULL)
  if (is.null(fit)) {
    results[[i]] <- list(vars = vars_set, k = kvars, p = p_sc,
                         bic = NA_real_, aic = NA_real_, has_coint = NA,
                         note = "VAR_fit_failed")
    next
  }
  
  # IC via logLik
  ic_std <- safe_logLik_BIC_AIC(fit)
  used_ic <- "logLik"
  
  # Fallback IC multivariés si besoin
  if (!is.finite(ic_std$bic) || !is.finite(ic_std$aic)) {
    ic_mv <- calc_ic_mv(fit, p_sc, det_terms_per_eq = 1L)
    ic_std$bic <- ic_mv$bic
    ic_std$aic <- ic_mv$aic
    used_ic <- "mv"
  }
  
  note_msg <- if (!is.finite(ic_std$bic) || !is.finite(ic_std$aic)) "singular_covariance" else used_ic
  
  # Johansen sur niveaux uniquement (si >=2)
  jres <- safe_test_johansen_on_levels(df_train, vars_set)
  has_coint <- !is.null(jres)
  
  results[[i]] <- list(
    vars = vars_set,
    k = kvars,
    p = p_sc,
    bic = ic_std$bic,
    aic = ic_std$aic,
    has_coint = has_coint,
    note = note_msg,
    johan = jres,
    fit = fit
  )
}

##===================== Résultats & meilleur modèle =================##
res_dt <- rbindlist(lapply(results, function(r) {
  data.table(
    variables = paste(r$vars, collapse = ","),
    k = r$k,
    p = r$p,
    BIC = r$bic,
    AIC = r$aic,
    has_coint = r$has_coint,
    note = r$note
  )
}), fill = TRUE)

res_dt_ok <- res_dt[is.finite(BIC)]
setorder(res_dt_ok, BIC)

if (nrow(res_dt_ok) == 0L) {
  message("Aucun modèle avec BIC fini. Voici les 10 premiers avec diagnostics :")
  print(head(res_dt, 10))
  best_fit <- NULL
  best_vars <- NULL
  best_p <- NA_integer_
  best_coint_summary <- NULL
} else {
  print(head(res_dt_ok, 10))
  best_row <- res_dt_ok[1]
  best_idx <- which(names(results) == best_row$variables[1])
  best_fit <- results[[best_idx]]$fit
  best_vars <- results[[best_idx]]$vars
  best_p    <- results[[best_idx]]$p
  best_coint_summary <- results[[best_idx]]$johan
  
  cat("\nMeilleur modèle (BIC):\n")
  print(best_row)
  cat("\nVariables:", paste(best_vars, collapse = ", "),
      "\nLag (p):", best_p,
      "\nIC_source:", results[[best_idx]]$note, "\n")
  if (!is.null(best_coint_summary)) {
    cat("\nJohansen: cointégration détectée (niveaux) sur ->\n",
        paste(best_coint_summary$variables, collapse = ", "), "\n")
  } else {
    cat("\nJohansen: pas de cointégration détectée (ou non applicable).\n")
  }
}

##===================== (Optionnel) Sauvegardes ====================##
# fwrite(res_dt,    "var_model_search_all.csv")
# fwrite(res_dt_ok, "var_model_search_ok.csv")
# saveRDS(list(results = results, best_fit = best_fit, table = res_dt_ok), "var_model_search.rds")
# Exemple d'usage : summary(best_fit); serial.test(best_fit, lags.pt = 12, type = "PT.asymptotic")






