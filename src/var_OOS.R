library(caret)
library(data.table)
library(vars)

# Charger les données
df <- fread("data/processed/data_macro_statio.csv")

# Exemple d'utilisation de la fonction de split des données
result <- split_data(df, train_prop = 0.8, validation_prop = 0.1, test_prop = 0.1)

# Accéder aux différents jeux de données
df_train <- result$trainData
df_valid <- result$validData
df_test <- result$testData

# Résumé des tailles des jeux de données
cat("\n Data training size :", nrow(df_train), "\n",
    "Data validation size :", nrow(df_valid), "\n",
    "Data test size :", nrow(df_test), "\n")

#### Prétraitement des données d'entraînement ####
df_train <- df_train[-1, ]  # Enlever la première ligne (due à la différenciation)
df_train$Date_quarter <-  NULL

# Assurez-vous que diff1_GPRD existe dans df_train et df_valid
# Supposons que diff1_GPRD est déjà calculée dans df.

# Sélectionner toutes les colonnes sauf diff1_GPRD pour la génération de combinaisons
all_vars <- setdiff(colnames(df_train), "diff1_GPRD")

# Liste pour stocker les résultats
results <- list()

# Test de toutes les combinaisons de 3 à 6 variables
for (n_vars in 3:6) {
  
  # Générer toutes les combinaisons possibles de variables
  combs <- combn(all_vars, n_vars, simplify = FALSE)
  
  for (comb in combs) {
    # Ajouter diff1_GPRD à la combinaison de variables
    selected_vars <- c("diff1_GPRD", comb)
    
    # Préparer les données pour le modèle VAR
    df_train_selected <- df_train[, selected_vars, with = FALSE]
    df_valid_selected <- df_valid[, selected_vars, with = FALSE]
    
    # Ajuster un modèle VAR
    var_model <- tryCatch({
      VAR(df_train_selected, p = 2)  # Utiliser 2 lags par exemple, ajustez selon votre analyse
    }, error = function(e) {
      return(NULL)  # En cas d'erreur, ignorer cette combinaison
    })
    
    # Si le modèle VAR a pu être ajusté, évaluer la performance
    if (!is.null(var_model)) {
      # Prédire sur l'ensemble de validation
      forecast <- predict(var_model, n.ahead = nrow(df_valid_selected))
      
      # Afficher la structure de forecast pour voir comment extraire les prévisions
      print(str(forecast))
      
      # Extraire les prévisions pour diff1_GPRD
      # En fonction de la structure de forecast, cela pourrait changer
      if ("diff1_GPRD" %in% colnames(forecast$fcst)) {
        predicted_values <- forecast$fcst[,"diff1_GPRD", drop = FALSE]  # Ajustez en fonction de la structure exacte
      } else {
        # Si les colonnes ne correspondent pas, essayer d'extraire diff1_GPRD
        predicted_values <- forecast$fcst[,1]  # Utilisez la première colonne comme prévision pour diff1_GPRD
      }
      
      # Calculer la performance (par exemple, RMSE)
      rmse <- sqrt(mean((predicted_values - df_valid_selected$diff1_GPRD)^2))
      
      # Stocker le résultat
      results[[paste(selected_vars, collapse = "_")]] <- list(
        model = var_model,
        rmse = rmse,
        selected_vars = selected_vars
      )
    }
  }
}

# Résultats triés par RMSE
sorted_results <- results[order(sapply(results, function(x) x$rmse))]

# # Afficher les meilleurs modèles selon RMSE
# cat("\nBest models based on RMSE:\n")
# for (i in 1:length(sorted_results)) {
#   cat("Model with variables:", sorted_results[[i]]$selected_vars, 
#       "RMSE:", sorted_results[[i]]$rmse, "\n")
# }
# 
# # Le meilleur modèle est le premier
# best_model <- sorted_results[[1]]
# cat("Best model selected with variables:", best_model$selected_vars, "with RMSE:", best_model$rmse, "\n")
