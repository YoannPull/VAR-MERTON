# ─────────────────────────────────────────────
# 0. Chargement des packages
# ─────────────────────────────────────────────
library(data.table)
library(ggplot2)
library(scales)

# ─────────────────────────────────────────────
# 1. Chargement et nettoyage des données
# ─────────────────────────────────────────────

fichier <- "data/raw/comparaison_data_caldara.csv"
dt <- fread(fichier, dec = ",")



# Nettoyage des colonnes : convertir les virgules décimales
clean_numeric_columns <- function(dt) {
  for (col in names(dt)) {
    if (is.character(dt[[col]]) || is.factor(dt[[col]])) {
      dt[[col]] <- as.numeric(gsub(",", ".", as.character(dt[[col]])))
    }
  }
  return(dt)
}
dt <- clean_numeric_columns(dt)


data_macro <- fread("data/processed/data_var_for_model.csv")

data_macro <- data_macro[,.(Date_quarter,GPRD)]
data_macro <- merge(data_macro,dt[,.(Date_quarter,C_GPRD)],
                  by.x = "Date_quarter", by.y = "Date_quarter")

data_macro_long <- melt(
  data_macro,
  id.vars = "Date_quarter",
  measure.vars = c("GPRD","C_GPRD"),
  variable.name = "Type",
  value.name = "VAR"
)

# Tracé
ggplot(data_macro_long, aes(x = Date_quarter, y = VAR, color = Type)) +
  geom_line() +
  labs(title = "Évolution des courbes VIX",
       x = "Date",
       y = "Valeur",
       color = "Type de var") +
  theme_minimal()
# ─────────────────────────────────────────────
# 2. Fonction de comparaison avancée
# ─────────────────────────────────────────────

compare_variables <- function(dt, col_c, col_ref, date_col = "Date_quarter") {
  x <- dt[[col_c]]
  y0 <- dt[[col_ref]]
  date <- dt[[date_col]]
  
  if (all(is.na(x)) || all(is.na(y0))) return(NULL)
  
  # Générer les 3 versions (non décalée, lead, lag)
  y_lead <- shift(y0, -1)
  y_lag <- shift(y0, 1)
  
  compare_align <- function(y) {
    ratio <- x / y
    factor <- median(ratio[is.finite(ratio)], na.rm = TRUE)
    y_rescaled <- y * factor
    mae <- mean(abs(x - y), na.rm = TRUE)
    cor_val <- cor(x, y, use = "complete.obs")
    mae_scaled <- mean(abs(x - y_rescaled), na.rm = TRUE)
    cor_scaled <- cor(x, y_rescaled, use = "complete.obs")
    list(mae = mae, cor = cor_val,
         mae_scaled = mae_scaled, cor_scaled = cor_scaled,
         factor = factor)
  }
  
  results <- list(
    "none" = compare_align(y0),
    "lead" = compare_align(y_lead),
    "lag"  = compare_align(y_lag)
  )
  
  # Sélection du meilleur alignement
  maes <- sapply(results, function(res) res$mae_scaled)
  best_shift <- names(which.min(maes))
  best_result <- results[[best_shift]]
  y_best <- switch(best_shift, none = y0, lead = y_lead, lag = y_lag)
  y_best_rescaled <- y_best * best_result$factor
  
  only_scaling <- best_result$mae_scaled < 0.1 * maes["none"] &&
    best_result$cor_scaled > 0.99
  
  # Graphique
  plot_data <- data.table(date, `Caldara` = x, `Local` = y_best_rescaled)
  p <- ggplot(plot_data, aes(x = date)) +
    geom_line(aes(y = Caldara, color = "Caldara")) +
    geom_line(aes(y = Local, color = "Local ajustée")) +
    labs(title = paste("Comparaison :", col_c, "vs", col_ref),
         subtitle = paste0("Décalage optimal : ", best_shift,
                           if (only_scaling) paste0(" | ⚠️ Échelle différente (~x", round(best_result$factor, 3), ")")),
         x = "Date", y = "Valeur") +
    scale_color_manual(values = c("Caldara" = "red", "Local ajustée" = "blue")) +
    theme_minimal()
  
  summary <- data.table(
    variable = gsub("^C_", "", col_c),
    correlation = round(best_result$cor, 4),
    MAE = round(best_result$mae, 4),
    RMSE = round(sqrt(mean((x - y_best_rescaled)^2, na.rm = TRUE)), 4),
    MaxAbsDiff = round(max(abs(x - y_best_rescaled), na.rm = TRUE), 4),
    ScalingFactor = round(best_result$factor, 4),
    MAE_AprèsRescaling = round(best_result$mae_scaled, 4),
    Cor_AprèsRescaling = round(best_result$cor_scaled, 4),
    ÉchelleSeulement = only_scaling,
    DécalageOptimal = best_shift
  )
  
  return(list(summary = summary, plot = p))
}

# ─────────────────────────────────────────────
# 3. Comparaison en boucle
# ─────────────────────────────────────────────

# Paires à comparer (manuellement listées ici)
pairs <- list(
  c("C_epu", "epu"), c("C_gs2", "gs2"), c("C_wti", "wti"), c("C_vix", "vix"),
  c("C_cpi", "cpi"), c("C_inv", "inv"), c("C_privemp", "priv_emp"),
  c("C_gdp", "gdp"), c("C_pop16", "pop16"), c("C_sp500", "sp500"),
  c("C_NFCI", "nfci"), c("C_GPRD", "GPRD"), c("C_GPRD_ACT", "GPRD_ACT"),
  c("C_GPRD_THREAT", "GPRD_THREAT")
)

results <- lapply(pairs, function(pair) {
  compare_variables(dt, pair[1], pair[2])
})

# ─────────────────────────────────────────────
# 4. Résumé global
# ─────────────────────────────────────────────

summary_table <- rbindlist(lapply(results, `[[`, "summary"))
summary_table <- summary_table[order(-MAE)]
print(summary_table)

# ─────────────────────────────────────────────
# 5. Affichage de tous les graphiques
# ─────────────────────────────────────────────

for (i in seq_along(results)) {
  if (!is.null(results[[i]]$plot)) {
    print(results[[i]]$plot)
    # readline(prompt = "Appuyez sur [Entrée] pour continuer...") # optionnel
  }
}
