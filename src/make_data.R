# ─────────────────────────────────────────────────────────────────────────────
# Chargement des packages
# ─────────────────────────────────────────────────────────────────────────────

library(tidyquant)
library(data.table)
library(readxl)
library(lubridate)
library(zoo)  # pour as.yearqtr et as.yearmon

# ─────────────────────────────────────────────────────────────────────────────
# Paramètres globaux
# ─────────────────────────────────────────────────────────────────────────────

# Répertoire contenant les fichiers bruts
data_dir <- "data/raw"

# ─────────────────────────────────────────────────────────────────────────────
# Fonctions utilitaires
# ─────────────────────────────────────────────────────────────────────────────

# Chargement depuis CSV FRED
load_csv_dt <- function(filename) {
  dt <- fread(file.path(data_dir, filename))
  dt[, DATE := as.Date(observation_date)]
  dt[, observation_date := NULL]
  return(dt)
}

# Conversion de yearqtr en décimal
convert_qtr_to_decimal <- function(qtr) {
  year <- floor(as.numeric(format(qtr, "%Y")))
  qnum <- as.numeric(cycle(qtr))  # retourne 1 à 4
  return(year + (qnum - 1) * 0.25)
}

# ─────────────────────────────────────────────────────────────────────────────
# Données mensuelles/quotidiennes/hebdomadaires -> quarterly décimal
# ─────────────────────────────────────────────────────────────────────────────

# S&P 500 (tidyquant → Yahoo Finance)
dt_sp500 <- tq_get("^GSPC", get = "stock.prices", from = "1900-01-01")
dt_sp500 <- as.data.table(dt_sp500)
dt_sp500[, quarter := as.yearqtr(date)]
dt_sp500_quarterly <- dt_sp500[, .(sp500 = mean(adjusted, na.rm = TRUE)), by = quarter]
# dt_sp500_quarterly[, quarter := quarter + 0.25]
dt_sp500_quarterly[, quarter := convert_qtr_to_decimal(quarter)]



# VIX
# Charger et préparer dt_vix
dt_vix <- load_csv_dt("VIXCLS.csv")
dt_vix[, quarter := as.yearqtr(DATE)]
dt_vix <- na.omit(dt_vix)
dt_vix_quarterly <- dt_vix[, .(vix = mean(VIXCLS, na.rm = TRUE)), by = quarter]
dt_vix_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# EMPL
# Charger et préparer dt_vix
dt_payems <- load_csv_dt("PAYEMS.csv")
dt_payems[, quarter := as.yearqtr(DATE)]
dt_payems <- na.omit(dt_payems)
dt_payems_quarterly <- dt_payems[, .(payems = mean(PAYEMS, na.rm = TRUE)), by = quarter]
dt_payems_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Charger et préparer dt_vix_caldara
dt_vix_caldara <- fread("data/raw/vix_caldara.csv")

# Fusionner les deux tables
dt_merged <- merge(
  dt_vix_quarterly,
  dt_vix_caldara,
  by = "quarter",
  all = TRUE
)
all_quarters <- seq(from = min(dt_merged$quarter),
                    to = max(dt_merged$quarter),
                    by = 0.25)

# Compléter les trimestres manquants
dt_complete <- data.table(quarter = all_quarters)
dt_merged <- merge(dt_complete, dt_merged, by = "quarter", all.x = TRUE)

setorder(dt_merged, quarter)
dt_merged[, vix_composite := fifelse(!is.na(vix), vix, SPVXO)]

dt_vix_quarterly <- dt_merged[,.(quarter,vix = vix_composite)]


# # TEST VIX FOR COMPARE_DATA.R
# dt_vix <- fread("data/raw/VIX_History.csv")
# dt_vix[, DATE := as.Date(DATE, format = "%m/%d/%Y")]
# dt_vix[, quarter := as.yearqtr(DATE)]
# dt_vix <- na.omit(dt_vix)
# dt_vix_quarterly <- dt_vix[, .(HIGH = mean(HIGH, na.rm = TRUE),
#                                CLOSE = mean(CLOSE, na.rm = TRUE),
#                                OPEN = mean(OPEN, na.rm = TRUE)), by = quarter]
# # Récupération de la valeur de fin de trimestre (CLOSE du dernier jour de chaque quarter)
# vix_eoq <- dt_vix[, .SD[.N], by = quarter][, .(quarter, EOQ = CLOSE)]
# # Fusion des deux tables
# dt_vix_quarterly <- merge(dt_vix_quarterly, vix_eoq, by = "quarter", all.x = TRUE)
# dt_vix_quarterly[, quarter_decimal := convert_qtr_to_decimal(quarter)]
# 
# setcolorder(dt_vix_quarterly, c("quarter_decimal", "quarter"))
# fwrite(dt_vix_quarterly,"data/processed/all_vix.csv")

# CPI
dt_cpi <- load_csv_dt("CPIAUCSL.csv")
dt_cpi[, quarter := as.yearqtr(DATE)]
dt_cpi_quarterly <- dt_cpi[, .(cpi = mean(CPIAUCSL, na.rm = TRUE)), by = quarter]
# dt_cpi_quarterly[, quarter := quarter + 0.25]
dt_cpi_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Population 16+
dt_pop16 <- load_csv_dt("CNP16OV.csv")
dt_pop16[, quarter := as.yearqtr(DATE)]
dt_pop16_quarterly <- dt_pop16[, .(pop16 = mean(CNP16OV, na.rm = TRUE)), by = quarter]
# dt_pop16_quarterly[, quarter := quarter + 0.25]
dt_pop16_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Emploi privé
dt_priv_emp <- load_csv_dt("USPRIV.csv")
dt_priv_emp[, quarter := as.yearqtr(DATE)]
dt_priv_emp_quarterly <- dt_priv_emp[, .(priv_emp = mean(USPRIV, na.rm = TRUE)), by = quarter]
# dt_priv_emp_quarterly[, quarter := quarter + 0.25]
dt_priv_emp_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Prix du pétrole WTI
dt_wti <- load_csv_dt("WTISPLC.csv")
dt_wti[, quarter := as.yearqtr(DATE)]
dt_wti_quarterly <- dt_wti[, .(wti = mean(WTISPLC, na.rm = TRUE)), by = quarter]
# dt_wti_quarterly[, quarter := quarter + 0.25]
dt_wti_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Taux 2 ans
dt_gs2 <- load_csv_dt("GS2.csv")
dt_gs2[, quarter := as.yearqtr(DATE)]
dt_gs2_quarterly <- dt_gs2[, .(gs2 = mean(GS2, na.rm = TRUE)), by = quarter]
# dt_gs2_quarterly[, quarter := quarter + 0.25]
dt_gs2_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Spread T10 ans T2 ans
dt_T10Y2Y <- load_csv_dt("T10Y2Y.csv")
dt_T10Y2Y[, quarter := as.yearqtr(DATE)]
dt_t10Y2Y_quarterly <- dt_T10Y2Y[, .(t10Y2Y = mean(T10Y2Y, na.rm = TRUE)),
                                 by = quarter]
# dt_gs2_quarterly[, quarter := quarter + 0.25]
dt_t10Y2Y_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Spread T10 ans T2 ans
dt_T10Y3M <- load_csv_dt("T10Y3M.csv")
dt_T10Y3M[, quarter := as.yearqtr(DATE)]
dt_t10Y3M_quarterly <- dt_T10Y3M[, .(t10Y3M = mean(T10Y3M, na.rm = TRUE)),
                                 by = quarter]
# dt_gs2_quarterly[, quarter := quarter + 0.25]
dt_t10Y3M_quarterly[, quarter := convert_qtr_to_decimal(quarter)]


# NFCI Index
dt_nfci <- load_csv_dt("NFCI.csv")
dt_nfci[, quarter := as.yearqtr(DATE)]
dt_nfci_quarterly <- dt_nfci[, .(nfci = mean(NFCI, na.rm = TRUE)), by = quarter]
# dt_nfci_quarterly[, quarter := quarter + 0.25]
dt_nfci_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# EPU Index (Excel)
dt_epu <- as.data.table(read_excel(file.path(data_dir, "EPU_US.xlsx"), sheet = 1))
dt_epu <- dt_epu[!is.na(Month) & !is.na(`News_Based_Policy_Uncert_Index`)]
dt_epu[, Year := as.integer(Year)]
dt_epu[, date := as.Date(as.yearmon(paste(Year, Month), format = "%Y %m"))]
dt_epu[, quarter := as.yearqtr(date)]
dt_epu_quarterly <- dt_epu[, .(epu = mean(`News_Based_Policy_Uncert_Index`, na.rm = TRUE)), by = quarter]
# dt_epu_quarterly[, quarter := quarter + 0.25]
dt_epu_quarterly[, quarter := convert_qtr_to_decimal(quarter)]


# Données GPR
dt_gpr <- fread("data/raw/data_gpr_daily.csv")
dt_gpr$date  <- NULL
dt_gpr$event <- NULL
dt_gpr$N10D  <- NULL
# 2. Conversion des colonnes numériques
cols_to_skip <- c("DAY", "event")
dt_gpr[] <- lapply(names(dt_gpr), function(col) {
  x <- dt_gpr[[col]]
  if (col %in% cols_to_skip) return(x)
  if (is.character(x)) x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
})
names(dt_gpr)[names(dt_gpr) == "DAY"] <- "date"
dt_gpr$date <- as.Date(as.character(dt_gpr$date), format = "%Y%m%d")
dt_gpr[, quarter := as.yearqtr(date)]
dt_gpr_quarterly <- dt_gpr[,.(GPRD        = mean(GPRD, na.rm = TRUE),
                              GPRD_ACT    = mean(GPRD_ACT, na.rm = TRUE),
                              GPRD_THREAT = mean(GPRD_THREAT, na.rm = TRUE)),
                           by = quarter]
dt_gpr_quarterly[, quarter := convert_qtr_to_decimal(quarter)]


# Taux d'Unemployment
dt_unemp <- load_csv_dt("UNRATE.csv")
dt_unemp[, quarter := as.yearqtr(DATE)]
dt_unemp_quarterly <- dt_unemp[, .(unrate = mean(UNRATE, na.rm = TRUE)), by = quarter]
# dt_wti_quarterly[, quarter := quarter + 0.25]
dt_unemp_quarterly[, quarter := convert_qtr_to_decimal(quarter)]



# ─────────────────────────────────────────────────────────────────────────────
# Données déjà trimestrielles
# ─────────────────────────────────────────────────────────────────────────────

# PIB réel
dt_gdp <- load_csv_dt("GDPC1.csv")
dt_gdp[, quarter := as.yearqtr(DATE)]
dt_gdp_quarterly <- dt_gdp[, .(gdp = GDPC1), by = quarter]
dt_gdp_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Investissement fixe réel
dt_inv <- load_csv_dt("GPDIC1.csv")
dt_inv[, quarter := as.yearqtr(DATE)]
dt_inv_quarterly <- dt_inv[, .(inv = GPDIC1), by = quarter]
dt_inv_quarterly[, quarter := convert_qtr_to_decimal(quarter)]

# Données GPR
# dt_gpr_quarterly <- fread("data/processed/data_gpr_quarterly.csv")
# dt_gpr_quarterly[,Date := NULL]
# setnames(dt_gpr_quarterly, old = "Date_quarter", new = "quarter")

# dt_gpr_quarterly <-  dt_gpr_quarterly[,.(quarter,GPRD,GPRD_ACT,GPRD_THREAT)]

# ─────────────────────────────────────────────────────────────────────────────
# Fusion des données trimestrielles
# ─────────────────────────────────────────────────────────────────────────────

list_tables <- list(
  dt_sp500_quarterly,
  dt_vix_quarterly,
  dt_cpi_quarterly,
  dt_pop16_quarterly,
  dt_priv_emp_quarterly,
  dt_wti_quarterly,
  dt_gs2_quarterly,
  dt_nfci_quarterly,
  dt_epu_quarterly,
  dt_gdp_quarterly,
  dt_inv_quarterly,
  dt_unemp_quarterly,
  dt_gpr_quarterly,
  dt_t10Y2Y_quarterly,
  dt_t10Y3M_quarterly,
  dt_payems_quarterly
)

data_macro <- Reduce(function(x, y) merge(x, y, by = "quarter", all = TRUE), list_tables)
setorder(data_macro, quarter)

setnames(data_macro, old = "quarter", new = "Date_quarter")
data_macro[, Date := as.Date(paste0(
  floor(Date_quarter), "-",
  sprintf("%02d", 3 * (Date_quarter %% 1 / 0.25) + 1),
  "-01"
))]

setcolorder(data_macro, c("Date_quarter", "Date"))


# Étape 1 : Identifier les colonnes à vérifier (hors "quarter")
cols_to_check <- setdiff(names(data_macro), "quarter")

# Étape 2 : Marquer les lignes sans aucun NA
data_macro[, complete_row := !Reduce(`|`, lapply(.SD, is.na)),
         .SDcols = cols_to_check]

# Étape 3 : Identifier les groupes consécutifs de lignes complètes
data_macro[, group := rleid(complete_row)]

# Étape 4 : Trouver le plus long groupe consécutif de lignes sans NA
valid_groups <- data_macro[complete_row == TRUE, .N, by = group][order(-N)]
max_valid_group <- valid_groups[1, group]

# Étape 5 : Extraire le sous-échantillon propre
data_macro <- data_macro[group == max_valid_group & complete_row == TRUE]

# Étape 6 : Supprimer les colonnes auxiliaires
data_macro[, c("complete_row", "group") := NULL]


fwrite(data_macro,"data/processed/data_macro_no_na.csv")