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

# Création d'un data set pour le comparer à celui de Caldara
data_no_transform <- copy(data_macro)
setnames(data_no_transform, old = "quarter", new = "Date_quarter")
data_no_transform[, Date := as.Date(paste0(
  floor(Date_quarter), "-",
  sprintf("%02d", 3 * (Date_quarter %% 1 / 0.25) + 1),
  "-01"
))]
setcolorder(data_no_transform, c("Date_quarter", "Date"))

data_no_transform <- data_no_transform[Date <= as.Date("2022-01-01") & Date >= as.Date("1990-01-01")]
fwrite(data_no_transform, file = "data/processed/data_raw.csv")

# ─────────────────────────────────────────────────────────────────────────────
# Transformations à la Caldara & Iacoviello (cf page 12 appendix)
# ─────────────────────────────────────────────────────────────────────────────

# Variables transformées pour le modèle VAR 

data_macro[, log_inv_pc               := 100 * log(inv / pop16)]        # log(investissement réel par tête)
data_macro[, log_gdp_pc               := 100 * log(gdp / pop16)]        # log(PIB réel par tête)
data_macro[, log_gdp                  := 100 * log(gdp)]
data_macro[, log_gdp_yoy_pct          := 100 * (log(gdp) - shift(log(gdp), 4L, type = "lag"))]    # log(croissance du PIB YOY) en %
data_macro[, log_gdp_annualized_pct   := 400 * (log(gdp) - log(shift(log(gdp), 1L, type = "lag")))]    # log de la croissance du PIB QOQ annualisée en %
data_macro[, gdp_yoy_pct              := 100 * (log(gdp) - shift(log(gdp), 4L, type = "lag"))]    # croissance du PIB YOY en %
data_macro[, gdp_annualized_pct       := 400 * (log(gdp) - shift(log(gdp), 1L, type = "lag"))]    # croissance du PIB QOQ annualisée en %
data_macro[, log_hours_pc             := 100 * log(priv_emp / pop16)]   # log(emploi privé par tête)
data_macro[, log_sp500_real           := 100 * log(sp500 / cpi)]        # log(S&P500 réel)
data_macro[, log_oil_real             := 100 * log(wti / cpi)]          # log(prix du pétrole réel)
data_macro[, infl_annualized_pct      := 400 * (log(cpi) - shift(log(cpi), 1L, type = "lag"))]  # inflation trimestrielle annualisée en %
data_macro[, infl_yoy_pct             := 100 * (log(cpi) - shift(log(cpi), 4L, type = "lag"))] # inflation glissante sur 1 an en %
data_macro[, log_payems               := 100 * log(payems)]          # log(prix du pétrole réel)


data_macro[, log_GPRD                 := 100 * log(GPRD)]
data_macro[, logdiff_GPRD             := 100 * (log(GPRD) - shift(log(GPRD)))]
data_macro[, diff_GPRD                := 100 * (GPRD - shift(GPRD))]
data_macro[, log_GPRD_ACT             := 100 * log(GPRD_ACT)]
data_macro[, logdiff_GPRD_ACT         := 100 * (log(GPRD_ACT) - shift(log(GPRD_ACT)))]
data_macro[, diff_GPRD_ACT            := 100 * (GPRD_ACT - shift(GPRD_ACT))]
data_macro[, log_GPRD_THREAT          := 100 * log(GPRD_THREAT)]
data_macro[, logdiff_GPRD_THREAT      := 100 * (log(GPRD_THREAT) - shift(log(GPRD_THREAT)))]
data_macro[, diff_GPRD_THREAT         := 100 * (GPRD_THREAT - shift(GPRD_THREAT))]                                                                                       

setnames(data_macro, old = "quarter", new = "Date_quarter")
data_macro[, Date := as.Date(paste0(
  floor(Date_quarter), "-",
  sprintf("%02d", 3 * (Date_quarter %% 1 / 0.25) + 1),
  "-01"
))]

setcolorder(data_macro, c("Date_quarter", "Date"))
# ─────────────────────────────────────────────────────────────────────────────
# Construction du jeu de données final pour le VAR
# ─────────────────────────────────────────────────────────────────────────────

data_var <- data_macro[, .(
  Date_quarter,
  Date,
  log_inv_pc,           # log(real investment per capita)
  log_gdp_pc,           # log(real GDP per capita)
  log_gdp,
  log_gdp_yoy_pct,
  log_gdp_annualized_pct,
  gdp_yoy_pct,
  gdp_annualized_pct,
  log_hours_pc,         # log(private employment per capita)
  log_sp500_real,       # log(real S&P500 index)
  log_oil_real,         # log(real oil price)
  infl_annualized_pct,  # inflation trimestrielle annualisée
  infl_yoy_pct,         # inflation glissante sur 1 an en %
  log_payems,
  vix,                  # VIX (volatilité implicite)
  gs2,                  # taux 2 ans
  t10Y2Y,               # Spread T10 ans T2ans
  t10Y3M,               # Spread T10 ans T2ans
  nfci,                 # NFCI
  epu,                  # indice d'incertitude politique
  cpi,                  # Inflation brut
  sp500,                # sp500 brut
  gdp,                  # real GDP
  pop16,                # Population des +16 aux USA
  inv,                  # Investissement fixe brut réel
  priv_emp,             # Nombre d’emplois privés
  wti,                  # Prix moyen trimestriel du baril de pétrole WTI
  unrate,
  payems,
  # GPR    
  logdiff_GPRD ,      
  diff_GPRD,          
  GPRD,
  log_GPRD,
  GPRD_ACT,           
  log_GPRD_ACT,      
  logdiff_GPRD_ACT,  
  diff_GPRD_ACT,      
  GPRD_THREAT,      
  log_GPRD_THREAT,    
  logdiff_GPRD_THREAT,
  diff_GPRD_THREAT
)]


# ─────────────────────────────────────────────────────────────────────────────
# Filtrage de la plus longue période sans NA (intersection complète)
# et ajout de la date réelle par trimestre
# ─────────────────────────────────────────────────────────────────────────────

# Étape 1 : Identifier les colonnes à vérifier (hors "quarter")
cols_to_check <- setdiff(names(data_var), "quarter")

# Étape 2 : Marquer les lignes sans aucun NA
data_var[, complete_row := !Reduce(`|`, lapply(.SD, is.na)),
         .SDcols = cols_to_check]

# Étape 3 : Identifier les groupes consécutifs de lignes complètes
data_var[, group := rleid(complete_row)]

# Étape 4 : Trouver le plus long groupe consécutif de lignes sans NA
valid_groups <- data_var[complete_row == TRUE, .N, by = group][order(-N)]
max_valid_group <- valid_groups[1, group]

# Étape 5 : Extraire le sous-échantillon propre
data_var_clean <- data_var[group == max_valid_group & complete_row == TRUE]

# Étape 6 : Supprimer les colonnes auxiliaires
data_var_clean[, c("complete_row", "group") := NULL]

# Sauvegarder le jeu complet avec toutes les variables macroéconomiques
fwrite(data_macro, file = "data/processed/data_macro_full.csv")

# Sauvegarder uniquement le jeu prêt pour estimation VAR
fwrite(data_var_clean, file = "data/processed/data_var_for_model.csv")
