# ==================== CHARGEMENT & STATIONNARITÉ GPR ==================== #

library(data.table)
library(lubridate)
library(tseries)
library(zoo)


# Conversion de yearqtr en décimal
convert_qtr_to_decimal <- function(qtr) {
  year <- floor(as.numeric(format(qtr, "%Y")))
  qnum <- as.numeric(cycle(qtr))  # retourne 1 à 4
  return(year + (qnum - 1) * 0.25)
}

# 1. Chargement du fichier GPR avec séparateur décimal en virgule
df_gpr_daily <- read.csv2("data/raw/data_gpr_daily.csv", dec = ",",
                          stringsAsFactors = FALSE)
df_gpr_daily$date <- NULL

# 2. Conversion des colonnes numériques
cols_to_skip <- c("DAY", "event")
df_gpr_daily[] <- lapply(names(df_gpr_daily), function(col) {
  x <- df_gpr_daily[[col]]
  if (col %in% cols_to_skip) return(x)
  if (is.character(x)) x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
})
names(df_gpr_daily)[names(df_gpr_daily) == "DAY"] <- "Date"
df_gpr_daily$Date <- as.Date(as.character(df_gpr_daily$Date), format = "%Y%m%d")

# 3. Ajout des trimestres
setDT(df_gpr_daily)
df_gpr_daily[, `:=`(
  Year = year(Date),
  Quarter = quarter(Date)
)]

df_gpr_daily[, `:=`(
  Date_quarter = as.Date(paste0(Year, "-", (Quarter - 1) * 3 + 1, "-01"))
)]

df_gpr_daily <- df_gpr_daily[Date >= as.Date("1990-01-01")]

# 4. Calcul des moyennes trimestrielles
cols_num <- setdiff(names(df_gpr_daily), c("Date", "event", "Year", 
                                           "Quarter", "Date_quarter"))
gpr_quarterly <- df_gpr_daily[, lapply(.SD, mean, na.rm = TRUE),
                              by = Date_quarter, .SDcols = cols_num]


# 5. Liste des colonnes GPR
gpr_columns <- grep("(?i)^gpr", names(gpr_quarterly), value = TRUE)

# 6. Fonction de test de stationnarité
test_stationarity <- function(x) {
  ts_clean <- na.omit(x)
  adf <- tryCatch(adf.test(ts_clean, k = 2, alternative = "stationary"), error = function(e) NULL)
  kpss <- tryCatch(kpss.test(ts_clean, null = "Level"), error = function(e) NULL)
  
  list(
    ADF_Stat = if (!is.null(adf)) round(adf$statistic, 3) else NA,
    ADF_p = if (!is.null(adf)) round(adf$p.value, 3) else NA,
    KPSS_Stat = if (!is.null(kpss)) round(kpss$statistic, 3) else NA,
    KPSS_p = if (!is.null(kpss)) round(kpss$p.value, 3) else NA,
    Conclusion = if (!is.null(adf) && !is.null(kpss)) {
      if (adf$p.value < 0.05 && kpss$p.value >= 0.05) "Stationnaire"
      else if (adf$p.value >= 0.05 && kpss$p.value < 0.05) "Non stationnaire"
      else "Ambigu"
    } else {
      "Erreur"
    }
  )
}

# 7. Application des tests à toutes les colonnes GPR
results_list <- lapply(gpr_columns, function(col) {
  res <- test_stationarity(gpr_quarterly[[col]])
  c(Variable = col, res)
})
results_dt <- rbindlist(lapply(results_list, as.data.table), fill = TRUE)

# 8. Affichage
cat("=== Résumé des tests de stationnarité sur les indices GPR ===\n\n")
print(results_dt)


# ==================== STATIONNARITÉ DE GPRD (diff & log-diff) ==================== #


# --- 2. Transformations GPR ---
gpr_quarterly[, log_GPRD := 100 * log(GPRD)]
gpr_quarterly[, logdiff_GPRD := 100 * (log_GPRD - shift(log_GPRD))]
gpr_quarterly[, diff_GPRD := 100 * (GPRD - shift(GPRD))]

gpr_quarterly[, log_GPRD_ACT := 100 * log(GPRD_ACT)]
gpr_quarterly[, logdiff_GPRD_ACT := 100 * (log_GPRD_ACT - shift(log_GPRD_ACT))]
gpr_quarterly[, diff_GPRD_ACT := 100 * (GPRD_ACT - shift(GPRD_ACT))]

gpr_quarterly[, log_GPRD_THREAT := 100 * log(GPRD_THREAT)]
gpr_quarterly[, logdiff_GPRD_THREAT := 100 * (log_GPRD_THREAT - shift(log_GPRD_THREAT))]
gpr_quarterly[, diff_GPRD_THREAT := 100 * (GPRD_THREAT - shift(GPRD_THREAT))]


# --- 3. Fonction de test de stationnarité ---
test_stationarity <- function(ts_vector, name) {
  ts_clean <- na.omit(ts_vector)
  adf <- tryCatch(adf.test(ts_clean, k = 2, alternative = "stationary"), error = function(e) NULL)
  kpss <- tryCatch(kpss.test(ts_clean, null = "Level"), error = function(e) NULL)
  
  adf_stat <- if (!is.null(adf)) round(adf$statistic, 3) else NA
  adf_p <- if (!is.null(adf)) round(adf$p.value, 3) else NA
  kpss_stat <- if (!is.null(kpss)) round(kpss$statistic, 3) else NA
  kpss_p <- if (!is.null(kpss)) round(kpss$p.value, 3) else NA
  
  conclusion <- if (!is.null(adf) && !is.null(kpss)) {
    if (adf$p.value < 0.05 && kpss$p.value >= 0.05) "Stationnaire"
    else if (adf$p.value >= 0.05 && kpss$p.value < 0.05) "Non stationnaire"
    else "Ambigu"
  } else {
    "Erreur"
  }
  
  list(
    Variable = name,
    ADF_Stat = adf_stat,
    ADF_p = adf_p,
    KPSS_Stat = kpss_stat,
    KPSS_p = kpss_p,
    Conclusion = conclusion
  )
}

# --- 4. Application aux séries GPR ---
gpr_vars <- c("log_GPRD","logdiff_GPRD", "diff_GPRD", "GPRD","GPRD_ACT", "log_GPRD_ACT",
              "logdiff_GPRD_ACT","diff_GPRD_ACT","GPRD_THREAT",
              "log_GPRD_THREAT","logdiff_GPRD_THREAT","diff_GPRD_THREAT")
gpr_results <- lapply(gpr_vars, 
                      function(v) test_stationarity(gpr_quarterly[[v]], v))
gpr_results_dt <- rbindlist(gpr_results)

# --- 5. Affichage du résumé ---
cat("\n=== Résumé des tests de stationnarité - Indices GPR ===\n\n")
print(gpr_results_dt)

# --- 6. Ajout des colonnes de date trimestrielle ---
gpr_quarterly$Date <- gpr_quarterly$Date_quarter
gpr_quarterly[, `:=`(
  Annee = year(Date),
  Mois = month(Date)
)]
gpr_quarterly[, `:=`(
  Date_quarter = Annee + (floor((Mois - 1) / 3) + 1) * 0.25
)]
gpr_quarterly[, Date_quarter := as.yearqtr(Date_quarter) - 0.25]
gpr_quarterly[, Date_quarter := convert_qtr_to_decimal(Date_quarter)]
# --- 7. Sélection des colonnes finales pour modélisation ---
gpr_quarterly <- gpr_quarterly[, .(Date, Date_quarter, log_GPRD, logdiff_GPRD,
                                   diff_GPRD, GPRD,
                                   GPRD_ACT, log_GPRD_ACT,logdiff_GPRD_ACT,
                                   diff_GPRD_ACT,GPRD_THREAT,
                                   log_GPRD_THREAT,logdiff_GPRD_THREAT,
                                   diff_GPRD_THREAT
)]



fwrite(gpr_quarterly, file = "data/processed/data_gpr_quarterly.csv")