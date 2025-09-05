# ====================================================================
# Packages
# ====================================================================
library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)

# ====================================================================
# Couleurs (Square Management)
# ====================================================================
squareblue  <- "#876b3a"
squaredark  <- "#1A1A1A"
squaregold  <- "#B68B4C"

# ====================================================================
# Thème commun
# ====================================================================
theme_square <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title        = element_text(face = "bold", color = squaredark, size = 16),
      plot.subtitle     = element_text(color = squaredark),
      plot.caption      = element_text(color = "grey40", size = 9),
      axis.title        = element_text(color = squaredark),
      axis.text         = element_text(color = "grey20"),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x= element_line(linetype = "dotted", linewidth = 0.25),
      panel.grid.major.y= element_line(linetype = "dotted", linewidth = 0.25)
    )
}

blank_titles <- theme(
  plot.title    = element_blank(),
  plot.subtitle = element_blank(),
  plot.caption  = element_blank()
)

# ====================================================================
# Données GPR (trimestrielles)
# ====================================================================
gpr_data <- fread("data/processed/data_var_for_model.csv")
gpr_data <- gpr_data[, .(Date_quarter, Date, GPRD, GPRD_ACT, GPRD_THREAT)]

# Date = Date si dispo, sinon Date_quarter
dcol <- if ("Date" %in% names(gpr_data) && !all(is.na(gpr_data$Date))) "Date" else "Date_quarter"
gpr_data[, date := as.Date(get(dcol))]
setorder(gpr_data, date)

# Quarter string "YYYYQn" pour joindre avec events_df
gpr_data[, quarter := paste0(format(date, "%Y"), "Q", ((as.integer(format(date, "%m")) - 1L) %/% 3L) + 1L)]

# ====================================================================
# Événements fournis (au format quarter = "YYYYQn")
# ====================================================================
events_df <- data.frame(
  quarter = c("2001Q4", "1990Q3", "2005Q3", "2014Q3", "2006Q3", "2023Q4", "1991Q1", "2022Q1", "2003Q1"),
  label   = c("11 September", "Kuwait Invasion", "London Bombings", "Escalation Ukraine/Russia",
              "Transatlantic Aircraft Plot", "Israel-Hamas War", "Gulf War",
              "Russia invades Ukraine", "Iraq Invasion")
)
events_dt <- as.data.table(events_df)

# Convertir les quarters évènement -> dates via la série GPR (join par quarter)
evt_gpr <- merge(events_dt, unique(gpr_data[, .(quarter, date)]),
                 by = "quarter", all.x = TRUE)
evt_gpr <- evt_gpr[!is.na(date)]  # garder ceux présents dans la série

# Position des labels au ras de la courbe GPR
setkey(gpr_data, date)
yrange_gpr   <- range(gpr_data$GPRD, na.rm = TRUE)
y_offset_gpr <- diff(yrange_gpr) * 0.03
events_pos_gpr <- gpr_data[evt_gpr, on = .(date), roll = "nearest"][
  , .(x = date, y = pmin(GPRD + y_offset_gpr, yrange_gpr[2]), label)]

# ====================================================================
# GRAPHIQUE 1 — GPR avec événements (sans titre/description)
# ====================================================================
graph_gpr_beau <-
  ggplot(gpr_data, aes(x = date, y = GPRD)) +
  geom_line(linewidth = 0.9, alpha = 0.95, color = squareblue) +
  geom_vline(data = evt_gpr, aes(xintercept = date),
             linetype = "dashed", alpha = 0.6, color = squaregold) +
  geom_label_repel(
    data = events_pos_gpr,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    direction = "y",
    size = 3.2,
    fill = "white", label.size = 0,
    box.padding = 0.2, point.padding = 0.15,
    min.segment.length = 0,
    seed = 1,
    segment.color = squaregold, segment.size = 0.2, segment.alpha = 0.5
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(big.mark = " ")) +
  labs(x = "Date", y = "GPR") +
  theme_square() + blank_titles

print(graph_gpr_beau)

# ====================================================================
# GRAPHIQUE 2 — Choc structurel e_{1,t} (sans titre/description)
# ====================================================================
e1_csv <- "output/e1_median.csv"
e1_dt <- fread(e1_csv)
e1_dt[, date := as.Date(date)]
setnames(e1_dt,"median_e1","e1")
e1_dt[, quarter := paste0(format(date, "%Y"), "Q", ((as.integer(format(date, "%m")) - 1L) %/% 3L) + 1L)]

evt_e1 <- merge(events_dt, unique(e1_dt[, .(quarter, date)]),
                by = "quarter", all.x = TRUE)
evt_e1 <- evt_e1[!is.na(date)]

setkey(e1_dt, date)
yrange_e1   <- range(e1_dt$e1, na.rm = TRUE)
y_offset_e1 <- diff(yrange_e1) * 0.05
events_pos_e1 <- e1_dt[evt_e1, on = .(date), roll = "nearest"][
  , .(x = date, y = pmax(pmin(e1 + y_offset_e1, yrange_e1[2]), yrange_e1[1]), label)]

p_e1_events <-
  ggplot(e1_dt, aes(x = date, y = e1)) +
  geom_line(linewidth = 0.9, alpha = 0.95, color = squareblue) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50",
             linetype = "dashed") +
  geom_vline(data = evt_e1, aes(xintercept = date),
             linetype = "dashed", alpha = 0.6, color = squaregold) +
  geom_label_repel(
    data = events_pos_e1,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    direction = "y",
    size = 3.2,
    fill = "white", label.size = 0,
    box.padding = 0.2, point.padding = 0.15,
    min.segment.length = 0,
    seed = 123,
    segment.color = squaregold, segment.size = 0.2, segment.alpha = 0.5
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(big.mark = " ")) +
  labs(x = "Date", y = expression(e[1])) +
  theme_square() + blank_titles

print(p_e1_events)

# ====================================================================
# Fonction générique IRF (sans titre) + option d'expansion Y
# ====================================================================
plot_bands_smart <- function(bands_df, ylab_txt = "Response", y_expand_mult = 0) {
  p <- ggplot(bands_df, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = squaregold, alpha = 0.30) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), fill = squaregold, alpha = 0.55) +
    geom_line(linewidth = 1.05, color = squareblue) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50") +
    scale_y_continuous(expand = expansion(mult = c(y_expand_mult, y_expand_mult))) +
    labs(x = "Horizon (quarters)", y = ylab_txt) +
    theme_square() + blank_titles
  
  if ("variable" %in% names(bands_df) && length(unique(bands_df$variable)) > 1) {
    p <- p + facet_wrap(~ variable, scales = "free_y")
  }
  p
}

# ====================================================================
# Autres graphiques (IRF, Z, PD) — sans titres ; Z/PD avec y un peu élargie
# ====================================================================
var_bands <- fread("output/irf/VAR/oirf_bands.csv")
z_bands   <- fread("output/irf/Z/psiZ_bands.csv")
pd_bands  <- fread("output/irf/PD/girf_pd_bands.csv")

p_var <- plot_bands_smart(var_bands, ylab_txt = "Response", y_expand_mult = 0.05)
print(p_var)

p_z   <- plot_bands_smart(z_bands,  ylab_txt = expression(psi[Z](h)), y_expand_mult = 0.12)
print(p_z)

p_pd  <- plot_bands_smart(pd_bands, ylab_txt = expression(Delta~PD~"(pp)"), y_expand_mult = 0.12)
print(p_pd)

# ====================================================================
# SCENARIO : (IRF, Z, PD) — mêmes réglages
# ====================================================================
var_bands_scenario <- fread("output_scenario/irf/VAR/oirf_bands.csv")
z_bands_scenario   <- fread("output_scenario/irf/Z/psiZ_bands.csv")
pd_bands_scenario  <- fread("output_scenario/irf/PD/girf_pd_bands.csv")

p_var_scenario <- plot_bands_smart(var_bands_scenario, ylab_txt = "Response", y_expand_mult = 0.05)
print(p_var_scenario)

p_z_scenario  <- plot_bands_smart(z_bands_scenario,  ylab_txt = expression(psi[Z](h)), y_expand_mult = 0.12)
print(p_z_scenario)

p_pd_scenario <- plot_bands_smart(pd_bands_scenario, ylab_txt = expression(Delta~PD~"(pp)"), y_expand_mult = 0.12)
print(p_pd_scenario)

# ====================================================================
# SAUVEGARDE DES FIGURES
# ====================================================================
dir.create("output/img", showWarnings = FALSE, recursive = TRUE)

ggsave("output/img/graph_gpr.png", graph_gpr_beau, width = 9, height = 5, dpi = 300)
ggsave("output/img/graph_e1.png", p_e1_events, width = 9, height = 5, dpi = 300)
ggsave("output/img/graph_var.png", p_var, width = 9, height = 5, dpi = 300)
ggsave("output/img/graph_z.png", p_z, width = 9, height = 5, dpi = 300)
ggsave("output/img/graph_pd.png", p_pd, width = 9, height = 5, dpi = 300)
ggsave("output/img/graph_var_scenario.png", p_var_scenario, width = 9, height = 5, dpi = 300)
ggsave("output/img/graph_z_scenario.png", p_z_scenario, width = 9, height = 5, dpi = 300)
ggsave("output/img/graph_pd_scenario.png", p_pd_scenario, width = 9, height = 5, dpi = 300)








# Conversion en data.table
setDT(events_df)

# Conversion des "1990Q3" etc. en dates (fin de trimestre)
events_df[, date := as.Date(as.yearqtr(quarter, format = "%YQ%q"), frac = 1)]

# Jointure avec les valeurs e1
events_with_e1 <- merge(events_df, e1_dt[, .(quarter, e1)], by = "quarter", all.x = TRUE)

events_with_e1



# Conversion en data.table
setDT(events_df)

# Conversion des "1990Q3" etc. en dates (fin de trimestre)
events_df[, date := as.Date(as.yearqtr(quarter, format = "%YQ%q"), frac = 1)]

# Jointure avec les valeurs e1
events_with_gpr <- merge(events_df, gpr_data[, .(quarter, GPRD)],
                         by = "quarter", all.x = TRUE)

events_with_gpr

max(gpr_data$GPRD)


# Comparer les max
max_gpr <- gpr_data[, .(date = quarter, gpr = GPRD)][which.max(gpr)]
max_e1  <- e1_dt[which.max(e1)]

max_gpr
max_e1

# ====================================================================
# GRAPHIQUE COMBINÉ — GPRD + e_{1,t} (axe secondaire) AVEC POINTS
# ====================================================================
# Alignement sur les dates communes (e1 commence à t = p+1)
combo <- merge(
  gpr_data[, .(date, GPRD)],
  e1_dt[, .(date, e1)],
  by = "date", all = FALSE
)

# Mapping linéaire : e1 → échelle GPR (par l'écart-type)
m_gpr <- mean(combo$GPRD, na.rm = TRUE)
sd_gpr <- sd(combo$GPRD,   na.rm = TRUE)
m_e1  <- mean(combo$e1,    na.rm = TRUE)
sd_e1 <- sd(combo$e1,      na.rm = TRUE)
scale_e1_to_gpr <- sd_gpr / sd_e1

combo[, e1_on_gpr := (e1 - m_e1) * scale_e1_to_gpr + m_gpr]

p_gpr_e1 <-
  ggplot(combo, aes(x = date)) +
  # Lignes
  geom_line(aes(y = GPRD), linewidth = 0.9, color = "red",  alpha = 0.95) +
  geom_line(aes(y = e1_on_gpr), linewidth = 0.8, color = "blue", alpha = 0.95) +
  # Points (mêmes couleurs)
  geom_point(aes(y = GPRD),      size = 1.6, color = "black",  alpha = 0.75) +
  geom_point(aes(y = e1_on_gpr), size = 1.6, color = "black", alpha = 0.75) +
  # Repères d'événements
  geom_vline(data = evt_gpr, aes(xintercept = date),
             linetype = "dashed", alpha = 0.35, color = squaregold) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(
    name = "GPR",
    labels = scales::label_number(big.mark = " "),
    sec.axis = sec_axis(~ (. - m_gpr) / scale_e1_to_gpr + m_e1,
                        name = expression(e[1]))
  ) +
  labs(x = "Date") +
  theme_square() + blank_titles

print(p_gpr_e1)
ggsave("output/img/graph_gpr_e1_combo.png", p_gpr_e1, width = 9, height = 5, dpi = 300)






# ====================================================================
# GRAPHIQUE — IRF PD : Baseline vs Scénario (même figure, clair)
# ====================================================================

# --- IRF PD : ruban Baseline + bornes Scénario uniquement ---
# --- IRF PD : ruban Baseline + bornes Scénario uniquement ---
p_pd_bands_outline <- ggplot() +
  geom_ribbon(data = pd_bands,
              aes(horizon, ymin = lower90, ymax = upper90),
              fill = squareblue, alpha = 0.18) +
  geom_ribbon(data = pd_bands,
              aes(horizon, ymin = lower68, ymax = upper68),
              fill = squareblue, alpha = 0.34) +
  # Scénario : uniquement les lignes des bornes
  geom_line(data = pd_bands_scenario, aes(horizon, upper90),
            color = squaregold, linetype = "dotted", linewidth = 0.8) +
  geom_line(data = pd_bands_scenario, aes(horizon, lower90),
            color = squaregold, linetype = "dotted", linewidth = 0.8) +
  geom_line(data = pd_bands_scenario, aes(horizon, upper68),
            color = squaregold, linetype = "longdash", linewidth = 1.0) +
  geom_line(data = pd_bands_scenario, aes(horizon, lower68),
            color = squaregold, linetype = "longdash", linewidth = 1.0) +
  # Médianes
  geom_line(data = pd_bands, aes(horizon, median, color = "Baseline"), linewidth = 1.1) +
  geom_line(data = pd_bands_scenario, aes(horizon, median, color = "Scénario"),
            linewidth = 1.1, linetype = "longdash") +
  scale_color_manual(values = c("Baseline" = squareblue, "Scénario" = squaregold), name = NULL) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50", linetype = "dashed") +
  scale_y_continuous(expand = expansion(mult = c(0.14, 0.14))) +
  labs(x = "Horizon (trimestres)", y = expression(Delta~PD~"(pp)")) +
  theme_square() + blank_titles

print(p_pd_bands_outline)
ggsave("output/img/graph_pd_bands_outline.png", p_pd_bands_outline, width = 9, height = 5, dpi = 300)
