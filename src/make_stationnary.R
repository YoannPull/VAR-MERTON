## ========================== SETUP ==========================
suppressPackageStartupMessages({
  library(data.table)
  library(caret)
  library(tseries)
  library(urca)
})

# IMPORTANT : split_data() doit effectuer un split chronologique (sans shuffle).
source("src/utils/split.R")
# Doit fournir check_stationarity() et check_I1()
source("src/utils/timeseries_test.R")

transform_level <- TRUE

## ===================== ROBUST HELPERS ======================
to_num <- function(x) {
  # Conversion sûre vers numérique (gère character, integer64, "1,234.56")
  if (is.numeric(x)) return(x)
  if (inherits(x, "integer64")) return(as.numeric(x))
  x <- as.character(x)
  x <- trimws(gsub(",", "", x, fixed = FALSE))
  x[x %in% c("", "NA", "NaN")] <- NA
  suppressWarnings(as.numeric(x))
}

apply_log_transform <- function(dt, vars, eps = .Machine$double.eps) {
  # Ajoute log_* pour chaque variable présente ; eps pour éviter log(0)
  missing <- setdiff(vars, names(dt))
  if (length(missing)) {
    message("Log variables missing (ignored): ", paste(missing, collapse = ", "))
  }
  present <- intersect(vars, names(dt))
  if (!length(present)) return(invisible(dt))
  # Coercition
  dt[, (present) := lapply(.SD, to_num), .SDcols = present]
  # Logs sûrs
  for (v in present) {
    nv <- paste0("log_", v)
    x  <- dt[[v]]
    if (any(!is.na(x) & x <= 0)) {
      nbad <- sum(!is.na(x) & x <= 0)
      warning(sprintf("'%s' contains %d value(s) <= 0: using pmax(x, eps) before log.", v, nbad))
    }
    dt[[nv]] <- log(pmax(x, eps))
  }
  invisible(dt)
}

## ======================= LOAD DATA =========================
data_macro      <- fread("data/processed/data_macro_no_na.csv")
data_macro_full <- fread("data/processed/data_macro_full.csv")

# Colonne temporelle
key_col_candidates <- c("Date_quarter", "Date")
key_col <- key_col_candidates[key_col_candidates %in% names(data_macro)][1]
if (is.na(key_col)) stop("No time column 'Date_quarter' or 'Date' in data_macro.")
if (!key_col %in% names(data_macro_full)) {
  stop(sprintf("Time column '%s' missing in data_macro_full.", key_col))
}

# Coercition numérique globale (hors colonnes temps)
nn_macro      <- setdiff(names(data_macro),      key_col_candidates)
nn_macro_full <- setdiff(names(data_macro_full), key_col_candidates)
data_macro[,      (nn_macro)      := lapply(.SD, to_num), .SDcols = nn_macro]
data_macro_full[, (nn_macro_full) := lapply(.SD, to_num), .SDcols = nn_macro_full]

# Optionnel : log-transform (transfo fixe → OK avant split)
if (transform_level) {
  level_variables <- c("sp500","gdp","vix","cpi","inv","wti","payems","priv_emp",
                       "epu","GPRD","GPRD_ACT","GPRD_THREAT","pop16")
  data_macro      <- apply_log_transform(data_macro,      level_variables)
  data_macro_full <- apply_log_transform(data_macro_full, level_variables)
}

## ==================== TRI CHRONOLOGIQUE ====================
# Sanity checks
if (anyNA(data_macro[[key_col]]))      stop(sprintf("'%s' has NA in data_macro.", key_col))
if (anyNA(data_macro_full[[key_col]])) stop(sprintf("'%s' has NA in data_macro_full.", key_col))

# Tri (programmable)
setorderv(data_macro,      key_col, na.last = TRUE)
setorderv(data_macro_full, key_col, na.last = TRUE)

## ================== SPLIT (chrono requis) ==================
split_res <- split_data(data_macro, train_prop = 0.8, validation_prop = 0.1, test_prop = 0.1)
train_dt  <- copy(split_res$trainData)
valid_dt  <- copy(split_res$validData)
test_dt   <- copy(split_res$testData)

# Frontières temporelles pour re-split après transfo
train_end <- max(train_dt[[key_col]])
valid_end <- max(valid_dt[[key_col]])

## ===== FIT : stationnarité sur TRAIN UNIQUEMENT ===========
# Drop colonnes non informatives si présentes
cols_drop <- c("Date", "GPRD_THREAT", "GPRD_ACT")
cols_drop <- cols_drop[cols_drop %in% names(train_dt)]
if (length(cols_drop)) train_dt[, (cols_drop) := NULL]

# Numériques uniquement
num_cols_train <- names(train_dt)[sapply(train_dt, is.numeric)]
if (!length(num_cols_train)) stop("No numeric columns in TRAIN to test stationarity.")
train_num <- train_dt[, ..num_cols_train]

# Tests
statio_res     <- check_stationarity(train_num)
list_I0_train  <- statio_res$statio_col
not_statio_col <- statio_res$not_statio_col

i1_res         <- check_I1(train_num, not_statio_col)
list_I1_train  <- i1_res$list_I1
list_I2_train  <- i1_res$list_I2

# Spec gelée
stationarity_spec <- list(
  I0 = list_I0_train,
  I1 = list_I1_train,
  I2 = list_I2_train
)

cat("FIT spec — I(0):", paste(stationarity_spec$I0, collapse=", "), "\n")
cat("FIT spec — I(1):", paste(stationarity_spec$I1, collapse=", "), "\n")
cat("FIT spec — I(2):", paste(stationarity_spec$I2, collapse=", "), "\n")

## ===== TRANSFORM : appliquer la spec au FULL (causal) ======
apply_stationarity_spec_to_full <- function(full_dt, spec, key_col) {
  exist_full <- names(full_dt)
  I0 <- intersect(spec$I0, exist_full)
  I1 <- intersect(spec$I1, exist_full)
  I2 <- intersect(spec$I2, exist_full)
  
  dt_full <- copy(full_dt)
  
  # Sortie avec clé normalisée
  dt_out <- dt_full[, ..key_col]
  setnames(dt_out, key_col, "Date_quarter")
  
  # 0) Conserver **toutes** les séries brutes I(0), I(1) et I(2) (pour cointégration)
  keep_raw <- unique(c(I0, I1, I2))
  if (length(keep_raw)) {
    dt_out <- cbind(dt_out, dt_full[, ..keep_raw])
  }
  
  # 1) Features stationnaires : I(0) (déjà stationnaires) => rien à faire
  
  # 2) I(1) -> ajouter diff1_*
  if (length(I1)) {
    for (col in I1) {
      if (!is.numeric(dt_full[[col]])) dt_full[[col]] <- to_num(dt_full[[col]])
      new_name <- paste0("diff1_", col)
      dt_out[[new_name]] <- dt_full[[col]] - shift(dt_full[[col]], 1L)
    }
  }
  
  # 3) I(2) -> ajouter diff2_*
  if (length(I2)) {
    for (col in I2) {
      if (!is.numeric(dt_full[[col]])) dt_full[[col]] <- to_num(dt_full[[col]])
      new_name <- paste0("diff2_", col)
      dt_out[[new_name]] <- dt_full[[col]] - 2*shift(dt_full[[col]], 1L) + shift(dt_full[[col]], 2L)
    }
  }
  
  dt_out[]
}

dt_statio_full <- apply_stationarity_spec_to_full(data_macro_full, stationarity_spec, key_col)

## ===== Dictionnaires pour sélection / traçabilité ==========
# Dictionnaire d'intégration "de base"
integration_dict <- rbindlist(list(
  data.table(variable = sort(stationarity_spec$I0), integration_order = "I(0)"),
  data.table(variable = sort(stationarity_spec$I1), integration_order = "I(1)"),
  data.table(variable = sort(stationarity_spec$I2), integration_order = "I(2)")
), fill = TRUE)

# Carte des features générées (utile pour filtrer les colonnes stationnaires)
feature_map <- rbindlist(list(
  if (length(stationarity_spec$I0)) data.table(
    feature = sort(stationarity_spec$I0),
    based_on = sort(stationarity_spec$I0),
    integration_order = "I(0)",
    type = "raw_stationary"
  ),
  if (length(stationarity_spec$I1)) rbindlist(lapply(stationarity_spec$I1, function(v)
    data.table(feature = paste0("diff1_", v), based_on = v, integration_order = "I(1)", type = "diff1")
  )),
  if (length(stationarity_spec$I2)) rbindlist(lapply(stationarity_spec$I2, function(v)
    data.table(feature = paste0("diff2_", v), based_on = v, integration_order = "I(2)", type = "diff2")
  ))
), fill = TRUE)

# Colonnes stationnaires à utiliser par défaut dans les modèles
stationary_feature_cols <- feature_map$feature

## ===== Post-check SANITY (sur TRAIN transformé) ============
# On ne teste la stationnarité que sur les features stationnaires
split_statio_by_bounds <- function(dt, key_col_norm = "Date_quarter", train_end, valid_end) {
  list(
    train = dt[get(key_col_norm) <= train_end],
    valid = dt[get(key_col_norm) >  train_end & get(key_col_norm) <= valid_end],
    test  = dt[get(key_col_norm) >  valid_end]
  )
}
splits_st <- split_statio_by_bounds(dt_statio_full, "Date_quarter", train_end, valid_end)

num_cols_st <- intersect(stationary_feature_cols,
                         names(splits_st$train)[sapply(splits_st$train, is.numeric)])
if (length(num_cols_st)) {
  train_st_num <- splits_st$train[, ..num_cols_st]
  res_train_st <- check_stationarity(train_st_num)
  if (length(res_train_st$not_statio_col)) {
    message("Still non-stationary on TRAIN after transforms (unexpected): ",
            paste(res_train_st$not_statio_col, collapse = ", "))
  } else {
    message("All transformed TRAIN features look stationary.")
  }
} else {
  message("No stationary features found to sanity-check on TRAIN (check spec).")
}

## ========================= SAVE ============================
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
fwrite(dt_statio_full,                    "data/processed/data_macro_statio.csv")
fwrite(na.omit(splits_st$train),          "data/processed/data_macro_statio_train.csv")
fwrite(na.omit(splits_st$valid),          "data/processed/data_macro_statio_valid.csv")
fwrite(na.omit(splits_st$test),           "data/processed/data_macro_statio_test.csv")
fwrite(integration_dict,                  "data/processed/integration_dict.csv")
fwrite(feature_map,                       "data/processed/feature_map.csv")

# Résumé
message("Summary:")
message("- FIT (train) — I(0): ", ifelse(length(stationarity_spec$I0), paste(stationarity_spec$I0, collapse=", "), "none"))
message("- FIT (train) — I(1): ", ifelse(length(stationarity_spec$I1), paste(stationarity_spec$I1, collapse=", "), "none"))
message("- FIT (train) — I(2): ", ifelse(length(stationarity_spec$I2), paste(stationarity_spec$I2, collapse=", "), "none"))
message("- Stationary features: ", ifelse(length(stationary_feature_cols), paste(stationary_feature_cols, collapse=", "), "none"))
message("Files written: data_macro_statio{,_train,_valid,_test}.csv, integration_dict.csv, feature_map.csv")
