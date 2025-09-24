# ================================================================
#   BVAR (NIW) + OIRF (Cholesky) — multi-spécifications (VAR only)
#   - Sortie de TOUS les graphiques IRF dans un même dossier
#   - Thème "theme_square" + titres/sous-titres/captions masqués
#   - Mesures de stabilité (part tirages stables, |λ|max)
#   - Couleurs doré/bleu explicites sur les rubans et la ligne
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mniw)
  library(parallel)
  library(ggplot2)
  library(lmtest)     # DIAG (optionnel)
  library(tseries)    # DIAG (optionnel)
})

options(stringsAsFactors = FALSE)
set.seed(123)  # reproductibilité globale

# ===================== COULEURS & THEME ========================= #
# Couleurs Square Management
squareblue  <- "#876b3a"  # ligne médiane (bleu)
squaredark  <- "#1A1A1A"
squaregold  <- "#B68B4C"  # rubans (doré)

theme_square <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", color = squaredark, size = 16),
      plot.subtitle      = element_text(color = squaredark),
      plot.caption       = element_text(color = "grey40", size = 9),
      axis.title         = element_text(color = squaredark),
      axis.text          = element_text(color = "grey20"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.25),
      panel.grid.major.y = element_line(linetype = "dotted", linewidth = 0.25),
      strip.text         = element_text(face = "bold", color = squaredark)
    )
}

blank_titles <- theme(
  plot.title    = element_blank(),
  plot.subtitle = element_blank(),
  plot.caption  = element_blank()
)

theme_set(theme_square())

# ----------------------- PARAMÈTRES GÉNÉRAUX ----------------------- #
data_var_path <- "data/processed/data_var_for_model.csv"

# Dossier racine des sorties
out_dir <- "output_robustness/var_specification/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Dossier UNIQUE pour TOUTES les images IRF
out_img_dir <- file.path(out_dir, "IRF_ALL")
dir.create(out_img_dir, showWarnings = FALSE, recursive = TRUE)

# Paramètres VAR + IRF
p_lags       <- 2           # ordre VAR
nrep_draws   <- 20000       # nb de tirages postérieurs B|Sigma ~ MNIW
h            <- 12          # horizons IRF (0..h)
impulse_ix   <- 1           # choc sur la 1ère variable (log_GPRD)
chol_jitter  <- 1e-10       # robustesse Cholesky
qs_irf       <- c(0.05,0.16,0.50,0.84,0.95)
min_stab_pct <- 0.50        # avertissement si < 50 %
M_target     <- 20000L      # nombre max de tirages stables conservés

# ----------------------- SPÉCIFICATIONS TESTÉES ----------------------- #
# IMPORTANT : 'log_GPRD' est toujours en 1er dans le VAR
spec_list <- list(
  S1_main      = c("vix","log_sp500_real","log_oil_real","log_hours_pc","log_gdp_pc"),
  S2_appendix  = c("vix","log_sp500_real","gs2","log_inv_pc","log_gdp_pc","log_hours_pc"),
  S3           = c("vix","log_sp500_real","log_oil_real","log_inv_pc","log_gdp_pc","log_hours_pc"),
  S4           = c("vix","log_sp500_real","t10Y2Y","t10Y3M","log_gdp_pc","log_hours_pc"),
  S5           = c("vix","nfci","epu","log_sp500_real","gs2","log_gdp_pc"),
  S6           = c("log_inv_pc","log_gdp_pc","log_hours_pc","log_oil_real","infl_yoy_pct"),
  S7           = c("vix","log_sp500_real","gs2","t10Y3M","nfci","epu"),
  S8           = c("gdp_yoy_pct","infl_annualized_pct","log_oil_real","vix","log_sp500_real"),
  S9           = c("vix","log_sp500_real","log_oil_real","log_gdp_pc"),
  S10          = c("t10Y2Y","gs2","log_gdp_pc","log_hours_pc","log_sp500_real")
)

# Traçabilité : table des specs testées
spec_tbl <- data.table(
  spec = names(spec_list),
  k_minus1 = vapply(spec_list, length, integer(1)),
  variables = vapply(spec_list, function(v) paste(v, collapse = ", "), character(1))
)
fwrite(spec_tbl, file.path(out_dir, "00_specs_tested.csv"))
cat("==> Spécifications testées écrites dans 00_specs_tested.csv\n")

# ----------------------- OUTILS VAR/BVAR ----------------------- #
sanitize_name <- function(x) gsub("[^A-Za-z0-9_\\-]+","_", x)

simulate_bvar_niw <- function(Y, p = 2, nrep = 5000, seed = 123) {
  set.seed(seed)
  Y <- as.matrix(Y)
  T <- nrow(Y); k <- ncol(Y)
  Y_t <- Y[(p+1):T, , drop = FALSE]
  X_t <- do.call(cbind, lapply(1:p, function(l) Y[(p+1-l):(T-l), , drop = FALSE]))
  X_t <- cbind(1, X_t)
  T_eff <- nrow(X_t); m <- ncol(X_t)
  nu_post <- T_eff - m
  if (nu_post <= (k - 1)) stop(sprintf("IW non définie: besoin T_eff - m > k - 1 (ici %d <= %d).", nu_post, k - 1))
  XtX <- crossprod(X_t); XtY <- crossprod(X_t, Y_t)
  B_ols <- qr.coef(qr(XtX), XtY)
  U_ols <- Y_t - X_t %*% B_ols
  S     <- crossprod(U_ols)
  XtX_inv <- chol2inv(chol(XtX))
  draws_B     <- array(NA_real_, dim = c(m, k, nrep))
  draws_Sigma <- array(NA_real_, dim = c(k, k, nrep))
  for (i in 1:nrep) {
    d <- mniw::rMNIW(1, Lambda = B_ols, Sigma = XtX_inv, Psi = S, nu = nu_post)
    draws_B[,, i]     <- d$X
    draws_Sigma[,, i] <- d$V
  }
  list(B = draws_B, S = draws_Sigma, p = p, k = k, m = m, T_eff = T_eff, nu_post = nu_post)
}

build_companion <- function(B, k, p) {
  A_comp <- matrix(0, nrow = k * p, ncol = k * p)
  A_stack <- t(B[-1, , drop = FALSE])                 # (k x (k*p)) = [A1 A2 ... Ap]
  A_comp[1:k, 1:(k*p)] <- A_stack
  if (p > 1) A_comp[(k+1):(k*p), 1:(k*(p-1))] <- diag(k*(p-1))
  A_comp
}

is_stable <- function(B, k, p) {
  rho <- max(Mod(eigen(build_companion(B, k, p), only.values = TRUE)$values))
  rho < 1 - 1e-10
}

compute_ma_coefficients <- function(B, p, H) {
  k <- ncol(B)
  A_comp <- build_companion(B, k, p)
  Psi <- array(0.0, dim = c(H + 1, k, k),
               dimnames = list(horizon = 0:H, variable = NULL, shock_col = NULL))
  for (j in 1:k) {
    delta <- rep(0, k); delta[j] <- 1
    state <- c(delta, rep(0, k*(p-1)))
    Psi[1, , j] <- delta
    for (hh in 1:H) {
      state <- A_comp %*% state
      Psi[hh + 1, , j] <- state[1:k]
    }
  }
  Psi
}

compute_oirf <- function(psi_draws, Sigma_kept, impulse_idx = 1,
                         qs = c(0.05,0.16,0.50,0.84,0.95),
                         jitter = 1e-10) {
  H <- dim(psi_draws)[1] - 1
  k <- dim(psi_draws)[2]
  M <- dim(psi_draws)[4]
  variables <- dimnames(psi_draws)$variable
  oirf_draws <- array(NA_real_, dim = c(H + 1, k, M),
                      dimnames = list(horizon = 0:H, variable = variables, draw = NULL))
  for (m in 1:M) {
    Psi_m   <- psi_draws[,,, m, drop=FALSE][,,,1]
    Sigma_m <- Sigma_kept[,, m]
    L <- tryCatch(t(chol(Sigma_m)),
                  error = function(e) t(chol(Sigma_m + diag(jitter, k))))
    delta_o <- L[, impulse_idx]  # choc d’1 écart-type structurel
    for (hh in 0:H) oirf_draws[hh+1, , m] <- Psi_m[hh+1, , ] %*% delta_o
  }
  out <- vector("list", k)
  for (r in 1:k) {
    mat <- oirf_draws[, r, , drop = FALSE]
    qh  <- t(apply(mat[,1,], 1, quantile, probs = qs, na.rm = TRUE))
    out[[r]] <- data.frame(horizon = 0:H,
                           lower90 = qh[,1], lower68 = qh[,2], median = qh[,3],
                           upper68 = qh[,4], upper90 = qh[,5],
                           variable = variables[r],
                           method = "OIRF (Cholesky)")
  }
  list(draws = oirf_draws, bands = do.call(rbind, out))
}

# --- Fonction de plot : couleurs doré/bleu + pas de titres -------------
plot_oirf_panels <- function(bands_df, var_order, impulse_name,
                             out_path_png, out_path_pdf = NULL,
                             title_prefix = "Orthogonalized IRFs (Cholesky)",
                             use_blank_titles = TRUE) {
  
  # copie + colonnes numériques (si CSV a importé en character)
  df <- as.data.frame(bands_df)
  must_num <- c("horizon","lower90","lower68","median","upper68","upper90")
  for (cc in must_num) if (!is.numeric(df[[cc]])) df[[cc]] <- as.numeric(df[[cc]])
  
  # ordre des facettes
  vars_in <- unique(df$variable)
  stopifnot(all(var_order %in% vars_in))
  df$variable <- factor(df$variable, levels = var_order)
  
  # nb de colonnes de facettes
  n_panels <- length(var_order)
  ncol_facets <- min(3L, max(1L, ceiling(n_panels / 2)))
  
  plot_bands_smart <- function(bands_df, ylab_txt = "Response", y_expand_mult = 0) {
    p <- ggplot(bands_df, aes(x = horizon, y = median)) +
      geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = squaregold, alpha = 0.30) +
      geom_ribbon(aes(ymin = lower68, ymax = upper68), fill = squaregold, alpha = 0.55) +
      geom_line(linewidth = 1.05, color = squareblue) +
      geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50") +
      scale_y_continuous(expand = expansion(mult = c(y_expand_mult, y_expand_mult))) +
      # >>> force des ticks entiers sur l’axe des x
      scale_x_continuous(
        breaks = function(x) seq(max(1, ceiling(x[1])), floor(x[2]), by = 1),
        minor_breaks = NULL
        # optionnel : pour démarrer à 1 exactement
        #, limits = c(1, NA)
      ) +
      labs(x = "Horizon (quarters)", y = ylab_txt) +
      theme_square() + blank_titles
    
    if ("variable" %in% names(bands_df) && length(unique(bands_df$variable)) > 1) {
      p <- p + facet_wrap(~ variable, scales = "free_y")
    }
    p
  }
  
  gg <- plot_bands_smart(df, ylab_txt = "Response", y_expand_mult = 0.05)
  
  ggsave(out_path_png, gg, width = 12, height = 8, dpi = 300)
  if (!is.null(out_path_pdf)) ggsave(out_path_pdf, gg, width = 12, height = 8)
  invisible(gg)
}

max_root <- function(Bmat, p, k) {
  max(Mod(eigen(build_companion(Bmat, k, p), only.values = TRUE)$values))
}

# ----------------------- CHARGEMENT DONNÉES ----------------------- #
stopifnot(file.exists(data_var_path))
DT_raw <- fread(data_var_path)

# ----------------------- TRAITEMENT PAR SPÉC ----------------------- #
process_one_spec <- function(sname, allowed_vect, DT_raw) {
  cat("\n==================== SPEC:", sname, "====================\n")
  i_var_str <- c("log_GPRD", allowed_vect)
  
  # Dossiers par spéc pour artefacts bruts (RDS/CSV), mais IMAGES dans out_img_dir
  spec_dir  <- file.path(out_dir, paste0("spec_", sanitize_name(sname)))
  var_dir   <- file.path(spec_dir, "irf", "VAR")
  dir.create(var_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Vérif colonnes
  missing_cols <- setdiff(i_var_str, names(DT_raw))
  if (length(missing_cols) > 0) {
    stop(sprintf("Colonnes manquantes pour %s : %s", sname, paste(missing_cols, collapse = ", ")))
  }
  
  # Sous-ensemble données (ordre exact)
  DT <- copy(DT_raw[, ..i_var_str])
  setcolorder(DT, i_var_str)
  if (anyNA(DT)) {
    n_before <- nrow(DT); DT <- na.omit(DT); n_after <- nrow(DT)
    cat(sprintf("NA détectés, lignes supprimées: %d -> %d\n", n_before, n_after))
  }
  
  # Série ts trimestrielle (dates facultatives)
  Y <- ts(as.matrix(DT), start = c(1990, 2), frequency = 4)
  variables <- colnames(Y)
  k <- ncol(Y)
  
  # Estimation BVAR (tirages postérieurs)
  res <- simulate_bvar_niw(Y, p = p_lags, nrep = nrep_draws,
                           seed = 123 + which(names(spec_list) == sname))
  
  # Stabilité des tirages
  keep <- vapply(1:dim(res$B)[3], function(i) is_stable(res$B[,, i], k, p_lags), logical(1))
  share <- mean(keep)
  idx_kept <- which(keep)
  cat(sprintf("Tirages stables: %d / %d (%.1f%%).\n", sum(keep), length(keep), 100*share))
  if (share < min_stab_pct) {
    warning(sprintf("[%s] Part de tirages stables < %.0f%% (%.1f%%).", sname, 100*min_stab_pct, 100*share))
  }
  
  # Échantillonnage (limite M_target)
  if (length(idx_kept) > M_target) {
    set.seed(42)
    idx_kept <- sample(idx_kept, M_target)
    cat(sprintf("Échantillonnage de %d tirages stables (sur %d).\n", M_target, sum(keep)))
  }
  
  # Conserve uniquement les tirages stables
  B_kept     <- res$B[,, idx_kept, drop = FALSE]
  Sigma_kept <- res$S[,, idx_kept, drop = FALSE]
  M_kept     <- dim(B_kept)[3]
  if (M_kept == 0) {
    warning(sprintf("[%s] Aucun tirage stable — spec ignorée.", sname))
    return(data.table(
      spec = sname, k = k, p = p_lags, nrep = nrep_draws, stable_share = 0,
      stable_n = 0, total_n = length(keep), max_root_median = NA_real_, max_root_p95 = NA_real_,
      variables = paste(i_var_str, collapse = ", "),
      irf_png = NA_character_, irf_pdf = NA_character_
    ))
  }
  
  # Diagnostics stabilité (|λ|max)
  roots <- vapply(1:M_kept, function(i) max_root(B_kept[,,i], p_lags, k), numeric(1))
  roots_med <- median(roots)
  roots_p95 <- as.numeric(quantile(roots, .95))
  
  # Ψ_h pour chaque tirage stable (parallèle)
  num_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(num_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  parallel::clusterExport(
    cl,
    varlist = c("B_kept","p_lags","h","compute_ma_coefficients","build_companion"),
    envir = environment()
  )
  parallel::clusterSetRNGStream(cl, 1234)
  
  psi_list <- tryCatch(
    parallel::parLapply(cl, 1:M_kept, function(ii) {
      compute_ma_coefficients(B_kept[,, ii], p = p_lags, H = h)
    }),
    error = function(e) {
      message("parLapply a échoué (", conditionMessage(e), "). Fallback en séquentiel…")
      lapply(seq_len(M_kept), function(ii) {
        compute_ma_coefficients(B_kept[,, ii], p = p_lags, H = h)
      })
    }
  )
  
  psi_draws <- array(
    NA_real_,
    dim = c(h + 1, k, k, M_kept),
    dimnames = list(
      horizon   = 0:h,
      variable  = variables,
      shock_col = variables,
      draw      = NULL
    )
  )
  for (jj in seq_len(M_kept)) psi_draws[ , , , jj] <- psi_list[[jj]]
  
  # OIRF (choc sur log_GPRD)
  oirf_res   <- compute_oirf(psi_draws, Sigma_kept, impulse_idx = impulse_ix,
                             qs = qs_irf, jitter = chol_jitter)
  oirf_draws <- oirf_res$draws
  oirf_bands <- oirf_res$bands
  
  # Sauvegardes brutes (RDS/CSV) PAR SPEC (pour traçabilité)
  saveRDS(list(Psi_draws = psi_draws, B_kept = B_kept, Sigma_kept = Sigma_kept,
               p = p_lags, k = k, h = h, variables = variables),
          file = file.path(spec_dir, "var_kernel.rds"))
  saveRDS(oirf_draws, file = file.path(var_dir, "oirf_draws.rds"))
  fwrite(oirf_bands,  file = file.path(var_dir, "oirf_bands.csv"))
  
  # ---------- Chemins IMAGES (UN SEUL DOSSIER POUR TOUTES LES FIGS) ----------
  impulse_name <- variables[impulse_ix]
  basefile <- sprintf("IRF_%s_oirf_panels_%s", sanitize_name(sname), sanitize_name(impulse_name))
  png_path  <- file.path(out_img_dir, paste0(basefile, ".png"))
  pdf_path  <- file.path(out_img_dir, paste0(basefile, ".pdf"))
  
  # Plot IRF (ordre = i_var_str, impulse = log_GPRD) avec thème + titres masqués
  plot_oirf_panels(oirf_bands,
                   var_order    = i_var_str,
                   impulse_name = impulse_name,
                   out_path_png = png_path,
                   out_path_pdf = pdf_path,
                   title_prefix = "Orthogonalized IRFs (Cholesky)",
                   use_blank_titles = TRUE)
  
  cat("Graphiques IRF -> ", png_path, " ; ", pdf_path, "\n", sep = "")
  
  # Résumé pour le catalog
  data.table(
    spec = sname,
    k = k,
    p = p_lags,
    nrep = nrep_draws,
    stable_share = round(share, 4),
    stable_n = sum(keep),
    total_n = length(keep),
    max_root_median = round(roots_med, 4),
    max_root_p95 = round(roots_p95, 4),
    variables = paste(i_var_str, collapse = ", "),
    irf_png = png_path,
    irf_pdf = pdf_path
  )
}

# ----------------------- BOUCLE SUR SPECS ----------------------- #
catalog <- rbindlist(lapply(names(spec_list), function(sname) {
  allowed_vect <- spec_list[[sname]]
  process_one_spec(sname, allowed_vect, DT_raw)
}), fill = TRUE)

# ----------------------- SORTIES SYNTHÈSE ----------------------- #
fwrite(catalog, file.path(out_dir, "99_catalog_summary.csv"))
print(catalog)

cat("\n==> Résumés écrits dans:\n",
    "- 00_specs_tested.csv (liste des specs)\n",
    "- 99_catalog_summary.csv (stabilité, |λ|max, chemins des graphes)\n",
    "Les IRFs (PNG/PDF) sont dans: ", normalizePath(out_img_dir), "\n", sep = "")