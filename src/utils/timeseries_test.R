library(urca)
library(tseries)
library(caret)
library(data.table)


# Fonction pour vérifier la stationnarité
check_stationarity <- function(data) {
  
  df <- data[,-1] # do not take the key 
  results <- list()
  list_statio <- c()
  list_not_statio <- c()
  for (col in colnames(df)) {
    adf_test <- adf.test(na.omit(df[[col]]), alternative = "stationary")
    
    results[[col]] <- adf_test
    
    if (adf_test$p.value < 0.05) {
      list_statio <- append(list_statio,col)
    } else {
      list_not_statio <- append(list_not_statio,col)
    }
  }
  
  
  return(list(results = results, statio_col = list_statio,
              not_statio_col = list_not_statio))
}

check_I1 <- function(data, not_statio_col) {
  # Do a check of NA in df
  list_I1 <- c()  # Liste pour stocker les variables I(1)
  list_I2 <- c()
  df <- data[,-1] # do not take the key 
  
  for (col in not_statio_col) {
    diff_series <- diff(df[[col]])  # Première différence de la série
    adf_test_diff <- adf.test(diff_series, alternative = "stationary")  # Test ADF sur la première différence
    
    if (adf_test_diff$p.value < 0.05) { 
      list_I1 <- append(list_I1, col)  
    } else {
      list_I2 <- append(list_I2,col)
    }
  }
  
  return(list(list_I1 = list_I1,list_I2 = list_I2))  # Retourner la liste des variables I(1)
}




engle_granger_test_df <- function(data,not_statio_col) {
  # Do a check of NA in df
  df <- data[,-1] # do not take the key
  coint_pairs <- list()        # Liste des paires cointégrées
  not_coint_pairs <- list()    # Liste des paires non cointégrées
  
  # Boucles pour itérer sur toutes les paires de variables
  for (x in not_statio_col){
    for (y in not_statio_col){
      if (x != y){
        # Extraire les séries temporelles
        y_val <- df[[x]]
        x_val <- df[[y]]
        
        # Estimer la régression entre les deux séries
        model <- lm(y_val ~ x_val)
        
        # Calculer les résidus de la régression
        residuals <- model$residuals
        
        # Tester la stationnarité des résidus avec le test ADF
        adf_test <- adf.test(residuals, alternative = "stationary")
        
        # Interprétation du test ADF
        if (adf_test$p.value < 0.05) {
          # Paires cointégrées
          coint_pairs <- append(coint_pairs, list(c(x, y)))
        } else {
          # Paires non cointégrées
          not_coint_pairs <- append(not_coint_pairs, list(c(x, y)))
        }
      }
    }
  }
  
  # Retourner les résultats sous forme de liste
  return(list(coint_pairs = coint_pairs, not_coint_pairs = not_coint_pairs))
}


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