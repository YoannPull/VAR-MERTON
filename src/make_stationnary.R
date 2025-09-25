# Chargement des librairies nécessaires
library(caret)
library(data.table)
library(tseries)
library(urca)

# Chargement des fonctions externes
source("src/utils/split.R")
source("src/utils/timeseries_test.R")

# Chargement des données
data_macro <- fread("data/processed/data_macro_no_na.csv")
data_macro_full <- fread("data/processed/data_macro_full.csv")

# Séparation des données en ensembles d'apprentissage, de validation et de test
result <- split_data(data_macro, train_prop = 0.8, validation_prop = 0.1, test_prop = 0.1)
df <- result$trainData
validData <- result$validData
testData <- result$testData



###### ON VERIFIE LA STATIO SEULEMENT SUR TRAIN !
# Suppression des colonnes inutiles
df[, c("Date", "GPRD_THREAT", "GPRD_ACT") := NULL]

# Vérification de la stationnarité des variables
stationarity_results <- check_stationarity(df)
list_I0 <- stationarity_results$statio_col
not_statio_col <- stationarity_results$not_statio_col

# Utilisation des variables non stationnaires pour vérifier la nécessité de différencier
results_check_I1 <- check_I1(df, not_statio_col)
list_I1 <- results_check_I1$list_I1
list_I2 <- results_check_I1$list_I2

# Affichage des résultats des vérifications de stationnarité
print(list_I1)
print(list_I2)



#### MAINTENANT QUE L'ON A NOS LISTES ON PEUX TRANSFORMER SUR TRAIN, VALIDATION
### TEST

## On se remet sur les données MACRO_FULL pour ne pas perdre deux lignes
## quant on différencie sur test. On peut utiliser les la diff sur macro pour split les données
col_names <- colnames(data_macro)
df_diff_1 <- copy(data_macro_full[, ..col_names])

# Suppression des colonnes non pertinentes pour la différenciation
df_diff_1[, c("Date", "GPRD_THREAT", "GPRD_ACT") := NULL]
df_diff_1 <- df_diff_1[, !..list_I1]
df_diff_1 <- df_diff_1[, !..list_I0]

# Appliquer la première différenciation (I(1)) sur les variables
for (col in list_I1) {
  new_name <- paste0("diff1_", col)
  df_diff_1[[new_name]] <- c(NA, diff(data_macro_full[[col]], differences = 1))
}

# Affichage des premières lignes après différenciation
head(df_diff_1)

# Créer un DataFrame pour la deuxième différenciation (I(2))
df_diff_2 <- copy(data_macro_full[, ..col_names])
df_diff_2[, c("Date", "GPRD_THREAT", "GPRD_ACT") := NULL]
df_diff_2 <- df_diff_2[, !..list_I2]

# Appliquer la deuxième différenciation (I(2)) sur les variables
for (col in list_I2) {
  new_name <- paste0("diff2_", col)
  df_diff_2[[new_name]] <- c(NA, NA, diff(data_macro_full[[col]], differences = 2))
}

# Fusionner les données différenciées en un seul DataFrame
df_statio <- merge(df_diff_1[, !..list_I2], df_diff_2[, !..list_I1], all.y = FALSE, by = "Date_quarter")
df_statio <- merge(data_macro[, c("Date_quarter")], df_statio)
results_statio <- check_stationarity(df_statio)
results_statio$not_statio_col

fwrite(df_statio, "data/processed/data_macro_statio.csv")