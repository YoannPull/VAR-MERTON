# ================================================================
#   BVAR (NIW) + OIRF (Cholesky) + Injection vers Z (Koop shock)
#   + μ_{t+h} exact + s^2(h), s^2_δ(h) + GIRF PD (Koop)
#   + DIAG (1) stabilité & (2) résidus
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mniw)
  library(parallel)
  library(ggplot2)
  library(lmtest)    # DIAG (2)
  library(tseries)   # DIAG (2)
  library(matrixStats) 
})

theme_set(theme_minimal(base_size = 14))
options(stringsAsFactors = FALSE)
set.seed(123)  # reproductibilité globale

# ----------------------- PARAMÈTRES GÉNÉRAUX ----------------------- #
data_var_path   <- "data/processed/data_var_for_model.csv"
best_model_path <- "scripts/model/best_model.rds"   # produit par ton pipeline Z existant
out_dir         <- "output/scenario/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dt_e1_median <- fread("output/e1_median.csv")

# Variables du VAR (ordre fixé, GPR en premier pour identification récursive)
allowed_vect <-  c("vix","log_sp500_real","log_oil_real",
                   "log_hours_pc","log_gdp_pc","nfci")

i_var_str <- c("log_GPRD",allowed_vect)

# Paramètres VAR + IRF
p_lags     <- 2        # ordre VAR
nrep       <- 20000    # nb de tirages postérieurs (B, Σ)
seed       <- 123
h          <- 12       # horizons (0..h)
impulse_ix <- 1        # choc sur GPR (1ère variable)
chol_jitter <- 1e-10   # robustesse Cholesky
# shock_scale <- max(dt_e1_median$median_e1)
shock_scale <- 1

# (Optionnel) standardiser Z sur l'échelle du training pour éviter saturation PD
scale_Z <- FALSE  # passe à TRUE si besoin

# ----------------------- OUTILS VAR/BVAR ----------------------- #
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
  
  # OLS stable
  B_ols <- qr.coef(qr(XtX), XtY)
  U_ols <- Y_t - X_t %*% B_ols
  S     <- crossprod(U_ols)
  
  # (XtX)^(-1) stable
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

# stabilité Hamilton TSA, cercle unitaire
is_stable <- function(B, k, p) {
  rho <- max(Mod(eigen(build_companion(B, k, p), only.values = TRUE)$values))
  rho < 1 - 1e-10
}

# Décompose B (m x k) en intercept c (k) et liste {A_1,..,A_p}
extract_c_A_list <- function(B, k, p) {
  c_vec <- as.numeric(B[1, ])
  A_stack <- t(B[-1, , drop = FALSE])                # k x (k*p)
  A_list <- vector("list", p)
  for (i in 1:p) A_list[[i]] <- A_stack[, ((i-1)*k + 1):(i*k), drop = FALSE]
  list(c = c_vec, A = A_list)
}

# Ψ_h (MA) : (h+1) x k x k
compute_ma_coefficients <- function(B, p, H) {
  k <- ncol(B)
  A_comp <- build_companion(B, k, p)
  Psi <- array(0.0, dim = c(H + 1, k, k), dimnames = list(horizon = 0:H, variable = NULL, shock_col = NULL))
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

# Forecast baseline E[Y_{t+h}|Ω] (zéro-chocs)
forecast_baseline_path <- function(B, Y_hist, p, H) {
  k <- ncol(Y_hist)
  decomp <- extract_c_A_list(B, k, p)
  c_vec <- decomp$c; A_list <- decomp$A
  lag_list <- lapply(1:p, function(ell) Y_hist[nrow(Y_hist) - ell + 1, , drop = FALSE])
  Y_fore <- matrix(NA_real_, nrow = H + 1, ncol = k)
  for (hh in 0:H) {
    y_next <- c_vec
    for (ell in 1:p) y_next <- y_next + A_list[[ell]] %*% as.numeric(lag_list[[ell]])
    Y_fore[hh + 1, ] <- y_next
    if (p > 1) for (ell in seq(p, 2)) lag_list[[ell]] <- lag_list[[ell-1]]
    lag_list[[1]] <- matrix(y_next, nrow = 1)
  }
  Y_fore
}

# ----------------------- CHARGEMENT DES DONNÉES VAR ----------------------- #
stopifnot(file.exists(data_var_path))
DT <- fread(data_var_path)
dates <- DT$Date
if ("Date" %in% names(DT))         DT[, Date := NULL]
if ("Date_quarter" %in% names(DT)) DT[, Date_quarter := NULL]

missing_cols <- setdiff(i_var_str, names(DT))
if (length(missing_cols) > 0) stop("Colonnes manquantes dans data_var_for_model: ",
                                   paste(missing_cols, collapse = ", "))
setcolorder(DT, i_var_str)
if (anyNA(DT)) DT <- na.omit(DT)

Y <- ts(as.matrix(DT[, ..i_var_str]), start = c(1990, 2), frequency = 4)
variables <- colnames(Y)

# ----------------------- ESTIMATION + Ψ (PARALLÈLE) ----------------------- #
res <- simulate_bvar_niw(Y, p = p_lags, nrep = nrep, seed = seed)
k <- res$k; p <- res$p

keep <- vapply(1:dim(res$B)[3], function(i) is_stable(res$B[,, i], k, p), logical(1))
message(sprintf("Tirages stables: %d / %d (%.1f%%).", sum(keep), length(keep), 100*mean(keep)))
idx_kept <- which(keep)
if (length(idx_kept) < 50) warning("Peu de tirages stables — bandes potentiellement bruitées.")

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

clusterExport(cl, varlist = c("res","idx_kept","p","h","build_companion","compute_ma_coefficients"),
              envir = environment())

psi_list <- parLapply(cl, idx_kept, function(i) {
  compute_ma_coefficients(res$B[,, i], p = p, H = h)
})

M_kept <- length(idx_kept)
psi_draws <- array(NA_real_, dim = c(h + 1, k, k, M_kept),
                   dimnames = list(horizon = 0:h, variable = variables, shock_col = variables, draw = NULL))
for (jj in seq_len(M_kept)) psi_draws[ , , , jj] <- psi_list[[jj]]

B_kept     <- res$B[,, idx_kept]   # (m x k x M_kept)
Sigma_kept <- res$S[,, idx_kept]   # (k x k x M_kept)
saveRDS(list(Psi_draws = psi_draws, B_kept = B_kept, Sigma_kept = Sigma_kept,
             p = p, k = k, h = h, variables = variables),
        file = file.path(out_dir, "var_kernel.rds"))

cat("Noyau VAR (Ψ, B, Σ) sauvegardé dans:\n",
    "-", file.path(out_dir, "var_kernel.rds"), "\n")

# ==================== DIAG (1) STABILITÉ ==================== #
max_root <- function(Bmat, p, k) {
  max(Mod(eigen(build_companion(Bmat, k, p), only.values = TRUE)$values))
}
roots <- vapply(1:dim(B_kept)[3], function(i) max_root(B_kept[,,i], p_lags, k), numeric(1))
stable_share <- mean(roots < 1)
cat(sprintf("DIAG(1) Stabilité — part tirages stables: %.1f%%; médiane |λ|max: %.3f; pct95: %.3f\n",
            100*stable_share, median(roots), quantile(roots, .95)))
png(file.path(out_dir,"diag1_roots_hist.png"), width=900, height=500)
hist(roots, breaks=30, main="Max root (post draws)", xlab="|λ|max"); abline(v=1, col="red", lwd=2)
dev.off()

# ==================== OIRF (Cholesky) ==================== #
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
    
    # Cholesky robuste (Sigma = L %*% t(L), L triangulaire inférieure)
    L <- tryCatch(t(chol(Sigma_m)),
                  error = function(e) t(chol(Sigma_m + diag(jitter, k))))
    delta_o <- shock_scale * L[, impulse_idx]  # ici c'est notre choc d'un écart-type
    cat("AFFICHAGE DE DELTA_O : \n")
    print(delta_o)
    
    for (hh in 0:H) {
      oirf_draws[hh+1, , m] <- Psi_m[hh+1, , ] %*% delta_o
    }
  }
  
  # Bandes
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

oirf_res   <- compute_oirf(psi_draws, Sigma_kept, impulse_idx = impulse_ix, jitter = chol_jitter)
oirf_draws <- oirf_res$draws
oirf_bands <- oirf_res$bands

saveRDS(oirf_draws, file = file.path(out_dir, "oirf_draws.rds"))
fwrite(oirf_bands,  file = file.path(out_dir, "oirf_bands.csv"))
cat("Fichiers OIRF écrits :\n",
    "-", file.path(out_dir, "oirf_draws.rds"), "\n",
    "-", file.path(out_dir, "oirf_bands.csv"), "\n")

plot_oirf <- function(df, title_txt = "Orthogonalized IRFs (Cholesky)") {
  ggplot(df, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), alpha = 0.30) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.55) +
    geom_line(linewidth = 1.05) +
    geom_hline(yintercept = 0) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(x = "Horizon (quarters)", y = "Response", title = title_txt) +
    theme(strip.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))
}
print(plot_oirf(oirf_bands))

# ----------------------- OUTILS MODELE Z (lm) & TERMES ----------------------- #
extract_terms_from_best_model <- function(best_model_path,
                                          allowed_bases = allowed_vect,
                                          gpr_name = "log_GPRD",
                                          vars_in_var = variables) {
  stopifnot(file.exists(best_model_path))
  best_model <- readRDS(best_model_path)
  betas <- coef(best_model)
  beta0 <- unname(betas["(Intercept)"])
  betas <- betas[names(betas) != "(Intercept)"]
  if (length(betas) == 0) stop("Le modèle ne contient pas de prédicteurs (seulement l'intercept).")
  parse_term <- function(name) {
    m <- regexec("^(.*)_lag([0-9]+)$", name)
    g <- regmatches(name, m)[[1]]
    if (length(g) != 3) return(list(base = NA_character_, lag = NA_integer_))
    list(base = g[2], lag = as.integer(g[3]))
  }
  terms_df <- data.frame(term = names(betas),
                         base = NA_character_, lag = NA_integer_, beta = as.numeric(betas),
                         stringsAsFactors = FALSE)
  for (i in seq_len(nrow(terms_df))) {
    p <- parse_term(terms_df$term[i])
    terms_df$base[i] <- p$base
    terms_df$lag[i]  <- p$lag
  }
  in_var   <- terms_df$base %in% vars_in_var
  not_gpr  <- terms_df$base != gpr_name
  allowed  <- terms_df$base %in% allowed_bases
  keep_row <- in_var & not_gpr & allowed
  if (any(!keep_row)) {
    dropped <- terms_df$term[!keep_row]
    if (length(dropped) > 0) {
      warning(sprintf("Termes ignorés (hors VAR / GPR / non-autorisés): %s",
                      paste(dropped, collapse = ", ")))
    }
    terms_df <- terms_df[keep_row, , drop = FALSE]
  }
  if (nrow(terms_df) == 0) stop("Aucun terme valide pour l'injection vers Z.")
  terms_df$col <- match(terms_df$base, vars_in_var)
  list(terms = terms_df, beta0 = beta0, best_model = best_model)
}

em <- extract_terms_from_best_model(best_model_path)
terms_df   <- em$terms
beta0      <- em$beta0
best_model <- em$best_model

# ----------------------- INJECTION OIRF -> Z (même choc) ----------------------- #
inject_oirf_into_Z <- function(var_kernel_rds, terms_df, impulse_idx = 1, jitter = 1e-10) {
  kern <- readRDS(var_kernel_rds)
  Psi_draws  <- kern$Psi_draws     # (H+1) x k x k x M
  Sigma_kept <- kern$Sigma_kept    # (k x k x M)
  variables  <- kern$variables
  H <- kern$h; k <- length(variables); M <- dim(Psi_draws)[4]
  
  psiZ_draws <- matrix(0.0, nrow = H + 1, ncol = M,
                       dimnames = list(horizon = 0:H, draw = NULL))
  
  for (m in 1:M) {
    Sigma <- Sigma_kept[,,m]
    L <- tryCatch(t(chol(Sigma)),
                  error = function(e) t(chol(Sigma + diag(jitter, k))))
    delta_o <- shock_scale * L[, impulse_idx]  # même choc que l’OIRF
    
    Psi_m <- Psi_draws[,,,m, drop=FALSE][,,,1]
    irfY <- matrix(0.0, nrow = H + 1, ncol = k)
    for (hh in 0:H) irfY[hh+1, ] <- Psi_m[hh+1, , ] %*% delta_o
    
    acc <- numeric(H + 1)
    for (r in 1:nrow(terms_df)) {
      lag_r  <- terms_df$lag[r]
      col_r  <- terms_df$col[r]
      beta_r <- terms_df$beta[r]
      contrib <- if (lag_r == 0) irfY[, col_r] else c(rep(0, lag_r), irfY[1:(H + 1 - lag_r), col_r])
      acc <- acc + beta_r * contrib
    }
    psiZ_draws[, m] <- acc
  }
  
  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(psiZ_draws, 1, quantile, probs = qs, na.rm = TRUE))
  psiZ_bands <- data.frame(
    horizon = 0:H,
    lower90 = qmat[,1], lower68 = qmat[,2], median = qmat[,3],
    upper68 = qmat[,4], upper90 = qmat[,5]
  )
  list(psiZ_draws = psiZ_draws, psiZ_bands = psiZ_bands, terms_used = terms_df)
}

inj_oirf <- inject_oirf_into_Z(
  var_kernel_rds = file.path(out_dir, "var_kernel.rds"),
  terms_df = terms_df,
  impulse_idx = impulse_ix,
  jitter = chol_jitter
)

saveRDS(inj_oirf$psiZ_draws, file = file.path(out_dir, "psiZ_draws_OIRF.rds"))
fwrite(inj_oirf$psiZ_bands,  file = file.path(out_dir, "psiZ_bands_OIRF.csv"))
fwrite(inj_oirf$terms_used,  file = file.path(out_dir, "psiZ_terms_used_OIRF.csv"))

cat("Injection OIRF -> Z sauvegardée dans:\n",
    "-", file.path(out_dir, "psiZ_draws_OIRF.rds"), "\n",
    "-", file.path(out_dir, "psiZ_bands_OIRF.csv"), "\n",
    "-", file.path(out_dir, "psiZ_terms_used_OIRF.csv"), "\n")

plot_psiZ_OIRF <- function(df, title_txt = "OIRF → Z (Koop shock)") {
  ggplot(df, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), alpha = 0.30) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.55) +
    geom_line(linewidth = 1.05) +
    geom_hline(yintercept = 0) +
    labs(x = "Horizon (quarters)", y = "ψ_Z(h)", title = title_txt) +
    theme(strip.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))
}
print(plot_psiZ_OIRF(inj_oirf$psiZ_bands))

# ----------------------- OUTILS POUR VARIANCES (Ψ_h -> Σ) ----------------------- #
build_selection_from_terms <- function(terms_df, variables) {
  bases <- terms_df$base; lags <- terms_df$lag
  m <- length(bases); k <- length(variables)
  Lmax <- max(lags)
  S_list <- vector("list", Lmax + 1)
  for (ell in 0:Lmax) {
    S <- matrix(0L, nrow = m, ncol = k)
    rows <- which(lags == ell)
    if (length(rows) > 0) for (r in rows) S[r, match(bases[r], variables)] <- 1L
    S_list[[ell + 1]] <- S
  }
  list(S_list = S_list, Lmax = Lmax)
}

build_G_hq <- function(h, q, S_list, Lmax, Psi_h_all) {
  m <- nrow(S_list[[1]]); k <- ncol(S_list[[1]])
  G <- matrix(0.0, nrow = m, ncol = k)
  for (a in 0:Lmax) {
    idx <- (h - a) - q
    if (idx >= 0) {
      Psi_idx <- Psi_h_all[idx + 1, , ]
      G <- G + S_list[[a + 1]] %*% Psi_idx
    }
  }
  G
}

# --- Baseline μ_{t+h} EXACTE : E[Z_{t+h}|Ω] --- #
compute_mu_baseline_draws <- function(DT_var, variables, terms_df, beta0, B_kept, p, H) {
  Tn <- nrow(DT_var); k <- length(variables); M <- dim(B_kept)[3]
  Lmax <- max(terms_df$lag)
  if (Tn <= Lmax) stop("Pas assez d'historique pour couvrir le Lmax des lags.")
  Y_all <- as.matrix(DT_var[, ..variables])
  Y_hist_for_forecast <- Y_all[(Tn - p + 1):Tn, , drop = FALSE]
  
  # FIX off-by-one : prendre Lmax+1 lignes (t-Lmax ... t)
  hist_tail <- Y_all[(Tn - Lmax):Tn, , drop = FALSE]
  
  mu_draws <- matrix(NA_real_, nrow = H + 1, ncol = M,
                     dimnames = list(horizon = 0:H, draw = NULL))
  for (m in 1:M) {
    Bm <- B_kept[,, m]
    Y_fore <- forecast_baseline_path(Bm, Y_hist_for_forecast, p, H)
    mu_h <- rep(beta0, H + 1)
    for (r in 1:nrow(terms_df)) {
      base_j <- terms_df$col[r]; lag_r <- terms_df$lag[r]; beta_r <- terms_df$beta[r]
      for (hh in 0:H) {
        val <- if (hh >= lag_r) {
          Y_fore[hh - lag_r + 1, base_j]
        } else {
          j <- lag_r - hh # j ∈ [1, Lmax]
          hist_tail[(Lmax + 1) - j, base_j]  # ligne cohérente avec t-j
        }
        mu_h[hh + 1] <- mu_h[hh + 1] + beta_r * val
      }
    }
    mu_draws[, m] <- mu_h
  }
  mu_draws
}

# --- s^2(h) et s^2_δ(h) (linéaire gaussien) --- #
compute_s2_one_draw <- function(Psi_draw_one, Sigma_one, terms_df, beta_vec, sigma_eta2) {
  H <- dim(Psi_draw_one)[1] - 1
  variables <- dimnames(Psi_draw_one)$variable
  sel <- build_selection_from_terms(terms_df, variables)
  S_list <- sel$S_list; Lmax <- sel$Lmax
  m <- nrow(S_list[[1]])
  s2      <- numeric(H + 1)
  s2_cond <- numeric(H + 1)
  for (hh in 0:H) {
    Sigma_s      <- matrix(0.0, m, m)
    Sigma_s_cond <- matrix(0.0, m, m)
    for (q in 0:hh) {
      G_hq <- build_G_hq(hh, q, S_list, Lmax, Psi_draw_one)
      Sigma_s      <- Sigma_s      + G_hq %*% Sigma_one %*% t(G_hq)
      if (q >= 1)  Sigma_s_cond <- Sigma_s_cond + G_hq %*% Sigma_one %*% t(G_hq) # Koop: exclure q=0
    }
    s2[hh + 1]      <- as.numeric(t(beta_vec) %*% Sigma_s      %*% beta_vec) + sigma_eta2
    s2_cond[hh + 1] <- as.numeric(t(beta_vec) %*% Sigma_s_cond %*% beta_vec) + sigma_eta2
  }
  list(s2 = s2, s2_cond = s2_cond)
}

# ===================== CALCUL μ_t+h, s^2, s^2_δ (par tirage) ===================== #
kernel_path <- file.path(out_dir, "var_kernel.rds")
stopifnot(file.exists(kernel_path))
kern <- readRDS(kernel_path)
Psi_draws   <- kern$Psi_draws
Sigma_kept  <- kern$Sigma_kept
B_kept      <- kern$B_kept
variables   <- kern$variables
H           <- kern$h
M           <- dim(Psi_draws)[4]

co <- coef(best_model)
beta0 <- unname(co["(Intercept)"])
betas <- co[names(co) != "(Intercept)"]
terms <- terms_df
beta_vec <- terms$beta
sigma_eta2 <- (summary(best_model)$sigma)^2   # variance de l'epsilon du modèle Z

mu_draws <- compute_mu_baseline_draws(DT_var = DT, variables = variables,
                                      terms_df = terms, beta0 = beta0,
                                      B_kept = B_kept, p = p, H = H)
saveRDS(mu_draws, file = file.path(out_dir, "z_mu_draws.rds"))

s2_draws        <- matrix(NA_real_, nrow = H + 1, ncol = M, dimnames = list(horizon = 0:H, draw = NULL))
s2_delta_draws  <- matrix(NA_real_, nrow = H + 1, ncol = M, dimnames = list(horizon = 0:H, draw = NULL))
for (m in 1:M) {
  Psi_m <- Psi_draws[,,, m, drop = FALSE][,,,1]
  Sigma_m <- Sigma_kept[,, m]
  res_m <- compute_s2_one_draw(Psi_m, Sigma_m, terms, beta_vec, sigma_eta2)
  s2_draws[ , m]       <- pmax(res_m$s2, 0)        # garde-fou
  s2_delta_draws[ , m] <- pmax(res_m$s2_cond, 0)   # garde-fou
}
saveRDS(list(mu_draws = mu_draws, s2_draws = s2_draws, s2_delta_draws = s2_delta_draws, terms_used = terms),
        file = file.path(out_dir, "z_mu_and_var_components_OIRF.rds"))

# ==================== DIAG (2) RESIDUS (moyenne postérieure) ==================== #
Y_mat <- as.matrix(as.data.frame(DT)[, i_var_str])
Y_t <- Y_mat[(p_lags+1):nrow(Y_mat), , drop=FALSE]
X_t <- do.call(cbind, lapply(1:p_lags, function(l) Y_mat[(p_lags+1-l):(nrow(Y_mat)-l), , drop=FALSE]))
X_t <- cbind(1, X_t)
B_bar <- apply(B_kept, c(1,2), mean)
U_bar <- Y_t - X_t %*% B_bar

lb_p <- apply(U_bar, 2, function(u) Box.test(as.numeric(u), lag = 12, type = "Ljung-Box")$p.value)
bp_p <- apply(U_bar, 2, function(u) bptest(lm(u ~ X_t[,-1]))$p.value)
jb_p <- apply(U_bar, 2, function(u) tseries::jarque.bera.test(as.numeric(u))$p.value)

resid_diag <- data.frame(variable = colnames(U_bar), p_LjungBox = lb_p, p_BreuschPagan = bp_p, p_JarqueBera = jb_p, row.names = NULL)
fwrite(resid_diag, file.path(out_dir,"diag2_residual_tests.csv"))
cat("DIAG(2) Résidus — p-values sauvegardées dans diag2_residual_tests.csv\n")

# ===================== (Option) STANDARDISATION DE Z ===================== #
# Utile si tu observes une saturation de la PD (ΔPD ~ 0 partout).
# NOTE: cette étape suppose que Z_train (réponse du modèle lm best_model) reflète l'échelle cible.
if (scale_Z) {
  Z_train <- model.response(model.frame(best_model))
  mZ <- mean(Z_train); sZ <- sd(Z_train)
  if (!is.finite(sZ) || sZ <= 0) stop("Echec standardisation: sd(Z_train) non positif.")
  # On ne touche PAS à mu_draws en niveau si ton mapping PD est défini sur Z standardisé.
  # Ici, on convertit tout vers l'échelle standardisée :
  mu_draws        <- (mu_draws - mZ) / sZ
  inj_oirf$psiZ_draws <- inj_oirf$psiZ_draws / sZ
  s2_draws        <- s2_draws       / (sZ^2)
  s2_delta_draws  <- s2_delta_draws / (sZ^2)
}

# ===================== GIRF de la PD (forme fermée) — SINGLE (p,ρ) ===================== #
compute_girf_pd_single <- function(psiZ_draws, mu_draws, s2_draws, s2_delta_draws,
                                   p, rho, qs = c(0.05, 0.16, 0.50, 0.84, 0.95)) {
  stopifnot(is.finite(p), p > 0, p < 1, is.finite(rho), rho > 0, rho < 1)
  H <- nrow(psiZ_draws) - 1
  M <- ncol(psiZ_draws)
  qpi <- qnorm(p)
  s2_draws       <- pmax(s2_draws, 0)
  s2_delta_draws <- pmax(s2_delta_draws, 0)
  pd_diff_draws <- matrix(NA_real_, nrow = H + 1, ncol = M,
                          dimnames = list(horizon = 0:H, draw = NULL))
  for (m in 1:M) {
    mu_base   <- mu_draws[, m]
    mu_delta  <- mu_draws[, m] + psiZ_draws[, m]
    s_h       <- sqrt( (1 - rho) + rho * s2_draws[, m] )
    s_h_delta <- sqrt( (1 - rho) + rho * s2_delta_draws[, m] )
    pd_base  <- pnorm( (qpi - sqrt(rho) * mu_base ) / s_h )
    pd_delta <- pnorm( (qpi - sqrt(rho) * mu_delta) / s_h_delta )
    pd_diff_draws[, m] <- pd_delta - pd_base
  }
  # Bandes (quantiles)
  qmat <- t(apply(pd_diff_draws, 1, quantile, probs = qs, na.rm = TRUE))
  bands <- data.frame(horizon = 0:H,
                      lower90 = qmat[,1], lower68 = qmat[,2], median = qmat[,3],
                      upper68 = qmat[,4], upper90 = qmat[,5],
                      p = p, rho = rho,
                      row.names = NULL)
  list(draws = pd_diff_draws, bands = bands)
}

# ---- TES paramètres ----
p0   <- 0.007954558   # PD inconditionnelle
rho0 <- 0.02383827    # corrélation d'actif
stopifnot(p0 > 0 && p0 < 1, rho0 > 0 && rho0 < 1)

girf_pd_single_OIRF <- compute_girf_pd_single(
  psiZ_draws     = readRDS(file.path(out_dir, "psiZ_draws_OIRF.rds")),
  mu_draws       = mu_draws,
  s2_draws       = s2_draws,
  s2_delta_draws = s2_delta_draws,
  p = p0, rho = rho0
)

# Sauvegardes
fwrite(girf_pd_single_OIRF$bands, file = file.path(out_dir, "girf_pd_bands_OIRF_SINGLE.csv"))
saveRDS(girf_pd_single_OIRF,       file = file.path(out_dir, "girf_pd_single_OIRF.rds"))
cat("GIRF-PD (Koop, choc OIRF) sauvegardé dans:\n",
    "-", file.path(out_dir, "girf_pd_bands_OIRF_SINGLE.csv"), "\n",
    "-", file.path(out_dir, "girf_pd_single_OIRF.rds"), "\n")

# ----------------------- Plot GIRF-PD (Koop) ----------------------- #
plot_girf_pd_single_OIRF <- function(bands_df,
                                     title_txt = "GIRF de la PD (Koop, choc OIRF)") {
  p_sel <- unique(bands_df$p); rho_sel <- unique(bands_df$rho)
  ggplot(bands_df, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), alpha = 0.30) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.55) +
    geom_line(linewidth = 1.05) +
    geom_hline(yintercept = 0) +
    labs(x = "Horizon (quarters)", y = "Δ PD (pts de prob.)", title = title_txt,
         subtitle = paste0("p = ", round(p_sel*100, 3), "% ; ρ = ", signif(rho_sel, 4))) +
    theme(strip.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))
}
print(plot_girf_pd_single_OIRF(girf_pd_single_OIRF$bands))

# ----------------------- Sanity check signe (optionnel) ----------------------- #
psiZ_med <- inj_oirf$psiZ_bands$median
dPD_med  <- girf_pd_single_OIRF$bands$median
if (length(psiZ_med) >= 1 && length(dPD_med) >= 1) {
  s0 <- sign(psiZ_med[1]); t0 <- sign(dPD_med[1])
  if (s0 < 0 && t0 <= 0) warning("Incohérence possible: ψ_Z(0)<0 mais ΔPD(0) ≤ 0 (avec ρ>0, PD devrait ↑).")
  if (s0 > 0 && t0 >= 0) warning("Incohérence possible: ψ_Z(0)>0 mais ΔPD(0) ≥ 0 (avec ρ>0, PD devrait ↓).")
}



message(sprintf("Tirages stables: %d / %d (%.1f%%).", sum(keep), length(keep), 100*mean(keep)))





# ===================== EXPORTS POUR LES GRAPHIQUES =====================

irf_dir <- file.path(out_dir, "irf")
var_dir <- file.path(irf_dir, "VAR")
z_dir   <- file.path(irf_dir, "Z")
pd_dir  <- file.path(irf_dir, "PD")
dir.create(var_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(z_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(pd_dir,  showWarnings = FALSE, recursive = TRUE)

# -- 1) VAR (OIRF des variables Y)
# Nécessite: oirf_draws (array) et oirf_bands (data.frame) déjà créés plus haut
data.table::fwrite(oirf_bands, file.path(var_dir, "oirf_bands.csv"))
saveRDS(oirf_bands,            file = file.path(var_dir, "oirf_bands.rds"))
saveRDS(oirf_draws,            file = file.path(var_dir, "oirf_draws.rds"))

# -- 2) Z (injection ψZ et composantes μ, s², s²δ)
# Nécessite: inj_oirf$psiZ_bands, inj_oirf$psiZ_draws, mu_draws, s2_draws, s2_delta_draws
data.table::fwrite(inj_oirf$psiZ_bands, file.path(z_dir, "psiZ_bands.csv"))
saveRDS(inj_oirf$psiZ_bands,            file = file.path(z_dir, "psiZ_bands.rds"))
saveRDS(inj_oirf$psiZ_draws,            file = file.path(z_dir, "psiZ_draws.rds"))
saveRDS(list(mu_draws = mu_draws,
             s2_draws = s2_draws,
             s2_delta_draws = s2_delta_draws,
             terms_used = terms_df),
        file = file.path(z_dir, "z_mu_s2_components.rds"))

# (petit résumé compact utile pour plots rapides)
med_mu   <- apply(mu_draws, 1, median, na.rm = TRUE)
med_s2   <- apply(s2_draws, 1, median, na.rm = TRUE)
med_s2_d <- apply(s2_delta_draws, 1, median, na.rm = TRUE)
data.table::fwrite(
  data.table::data.table(horizon = 0:H, mu_med = med_mu,
                         s2_med = med_s2, s2_delta_med = med_s2_d),
  file.path(z_dir, "z_mu_s2_components_median.csv")
)

# -- 3) PD (GIRF ΔPD)
# Nécessite: girf_pd_single_OIRF$bands (df) et girf_pd_single_OIRF$draws (matrix)
data.table::fwrite(girf_pd_single_OIRF$bands, file.path(pd_dir, "girf_pd_bands.csv"))
saveRDS(girf_pd_single_OIRF$bands,            file = file.path(pd_dir, "girf_pd_bands.rds"))
saveRDS(girf_pd_single_OIRF$draws,            file = file.path(pd_dir, "girf_pd_draws.rds"))

# (résumé médian)
pd_med <- apply(girf_pd_single_OIRF$draws, 1, median, na.rm = TRUE)
data.table::fwrite(
  data.table::data.table(horizon = 0:H, dPD_median = pd_med, p = p0, rho = rho0),
  file.path(pd_dir, "girf_pd_median.csv")
)

cat("\n==> FICHIERS ÉCRITS pour les graphes :\n",
    "- VAR : ", var_dir, "\n",
    "- Z   : ", z_dir,  "\n",
    "- PD  : ", pd_dir,  "\n")