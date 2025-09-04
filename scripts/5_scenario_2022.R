# ================================================================
#   BVAR (NIW) + OIRF (Cholesky) + Injection vers Z (Koop shock)
#   + μ_{t+h} exact + s^2(h), s^2_δ(h) + GIRF PD (Koop)
#   + DIAG (1) stabilité & (2) résidus
#   + REPLAY 2022: (A) choc unique 2022Q1 ; (B) séquence 2022Q1-Q4
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
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
out_dir         <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Variables du VAR (ordre fixé, GPR en premier pour identification récursive)
allowed_vect <- c("vix","log_sp500_real","log_oil_real","log_hours_pc",
                  "log_gdp_pc","nfci")
i_var_str <- c("log_GPRD", allowed_vect)

# Paramètres VAR + IRF
p_lags     <- 2        # ordre VAR
nrep       <- 20000    # nb de tirages postérieurs (B, Σ)
seed       <- 123
h          <- 12       # horizons (0..h)
impulse_ix <- 1        # choc sur GPR (1ère variable)
chol_jitter <- 1e-10   # robustesse Cholesky

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
dates <- as.Date(DT$Date)
if ("Date" %in% names(DT))         DT[, Date := NULL]
if ("Date_quarter" %in% names(DT)) DT[, Date_quarter := NULL]

missing_cols <- setdiff(i_var_str, names(DT))
if (length(missing_cols) > 0) stop("Colonnes manquantes dans data_var_for_model: ",
                                   paste(missing_cols, collapse = ", "))
setcolorder(DT, i_var_str)
if (anyNA(DT)) DT <- na.omit(DT)

Y <- ts(as.matrix(DT[, ..i_var_str]), frequency = 4)
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

# ==================== OIRF (Cholesky, choc unitaire) ==================== #
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
    delta_o <- L[, impulse_idx]
    for (hh in 0:H) {
      oirf_draws[hh+1, , m] <- Psi_m[hh+1, , ] %*% delta_o
    }
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

oirf_res   <- compute_oirf(psi_draws, Sigma_kept, impulse_idx = impulse_ix, jitter = chol_jitter)
oirf_draws <- oirf_res$draws
oirf_bands <- oirf_res$bands
saveRDS(oirf_draws, file = file.path(out_dir, "oirf_draws.rds"))
fwrite(oirf_bands,  file = file.path(out_dir, "oirf_bands.csv"))

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

# ----------------------- INJECTION OIRF -> Z (choc unitaire) ----------------------- #
inject_oirf_into_Z <- function(var_kernel_rds, terms_df, impulse_idx = 1, jitter = 1e-10) {
  kern <- readRDS(var_kernel_rds)
  Psi_draws  <- kern$Psi_draws
  Sigma_kept <- kern$Sigma_kept
  variables  <- kern$variables
  H <- kern$h; k <- length(variables); M <- dim(Psi_draws)[4]
  psiZ_draws <- matrix(0.0, nrow = H + 1, ncol = M,
                       dimnames = list(horizon = 0:H, draw = NULL))
  for (m in 1:M) {
    Sigma <- Sigma_kept[,,m]
    L <- tryCatch(t(chol(Sigma)),
                  error = function(e) t(chol(Sigma + diag(jitter, k))))
    delta_basis <- L[, impulse_idx]  # P e1
    Psi_m <- Psi_draws[,,,m, drop=FALSE][,,,1]
    irfY <- matrix(0.0, nrow = H + 1, ncol = k)
    for (hh in 0:H) irfY[hh+1, ] <- Psi_m[hh+1, , ] %*% delta_basis
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
forecast_mu_baseline_draws <- function(DT_var, variables, terms_df, beta0, B_kept, p, H) {
  Tn <- nrow(DT_var); k <- length(variables); M <- dim(B_kept)[3]
  Lmax <- max(terms_df$lag)
  if (Tn <= Lmax) stop("Pas assez d'historique pour couvrir le Lmax des lags.")
  Y_all <- as.matrix(DT_var[, ..variables])
  Y_hist_for_forecast <- Y_all[(Tn - p + 1):Tn, , drop = FALSE]
  hist_tail <- Y_all[(Tn - Lmax):Tn, , drop = FALSE]  # t-Lmax ... t
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
          j <- lag_r - hh
          hist_tail[(nrow(hist_tail) - j), base_j]
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
sigma_eta2 <- (summary(best_model)$sigma)^2

mu_draws <- forecast_mu_baseline_draws(DT_var = DT, variables = variables,
                                       terms_df = terms, beta0 = beta0,
                                       B_kept = B_kept, p = p_lags, H = H)
saveRDS(mu_draws, file = file.path(out_dir, "z_mu_draws.rds"))

s2_draws        <- matrix(NA_real_, nrow = H + 1, ncol = M, dimnames = list(horizon = 0:H, draw = NULL))
s2_delta_draws  <- matrix(NA_real_, nrow = H + 1, ncol = M, dimnames = list(horizon = 0:H, draw = NULL))
for (m in 1:M) {
  Psi_m <- Psi_draws[,,, m, drop = FALSE][,,,1]
  Sigma_m <- Sigma_kept[,, m]
  res_m <- compute_s2_one_draw(Psi_m, Sigma_m, terms, beta_vec, sigma_eta2)
  s2_draws[ , m]       <- pmax(res_m$s2, 0)
  s2_delta_draws[ , m] <- pmax(res_m$s2_cond, 0)
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

# ===================== (Option) STANDARDISATION DE Z ===================== #
if (scale_Z) {
  Z_train <- model.response(model.frame(best_model))
  mZ <- mean(Z_train); sZ <- sd(Z_train)
  if (!is.finite(sZ) || sZ <= 0) stop("Echec standardisation: sd(Z_train) non positif.")
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

# ---- GIRF-PD pour le choc OIRF unitaire ----
girf_pd_single_OIRF <- compute_girf_pd_single(
  psiZ_draws     = readRDS(file.path(out_dir, "psiZ_draws_OIRF.rds")),
  mu_draws       = mu_draws,
  s2_draws       = s2_draws,
  s2_delta_draws = s2_delta_draws,
  p = p0, rho = rho0
)
fwrite(girf_pd_single_OIRF$bands, file = file.path(out_dir, "girf_pd_bands_OIRF_SINGLE.csv"))
saveRDS(girf_pd_single_OIRF,       file = file.path(out_dir, "girf_pd_single_OIRF.rds"))

# ----------------------- PLOTS RAPIDES ----------------------- #
plot_bands_generic <- function(df, ylab, title_txt) {
  ggplot(df, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), alpha = 0.30) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.55) +
    geom_line(linewidth = 1.05) +
    geom_hline(yintercept = 0) +
    labs(x = "Horizon (quarters)", y = ylab, title = title_txt) +
    theme(strip.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))
}
# (Exemple d'affichage si interactif)
# print(plot_bands_generic(oirf_bands, "Response", "Orthogonalized IRFs (Cholesky)"))
# print(plot_bands_generic(inj_oirf$psiZ_bands, "ψ_Z(h)", "OIRF → Z (Koop shock)"))
# print(plot_bands_generic(girf_pd_single_OIRF$bands, "Δ PD (pp)", "GIRF de la PD (Koop, choc OIRF)"))

# ===================== REPLAY 2022 : CHOC CALIBRÉ & SÉQUENCE ===================== #
# 1) Récupère les chocs structurels historiques (posterior mean)
compute_structural_shocks <- function(Y_mat, B_bar, Sigma_bar, p_lags) {
  # résidus réduits (échantillon à partir de t=p+1)
  Tn <- nrow(Y_mat)
  Y_t <- Y_mat[(p_lags+1):Tn, , drop=FALSE]
  X_t <- do.call(cbind, lapply(1:p_lags, function(l) Y_mat[(p_lags+1-l):(Tn-l), , drop=FALSE]))
  X_t <- cbind(1, X_t)
  U_bar <- Y_t - X_t %*% B_bar    # réduits
  # chocs structurels ε_t = P^{-1} u_t ; P = t(chol(Sigma))'
  L_bar <- t(chol(Sigma_bar + diag(1e-12, ncol(Sigma_bar))))
  P_bar <- L_bar                   # ici P = L (par construction ci-dessus)
  eps <- t(solve(P_bar, t(U_bar))) # (T_eff x k)
  list(eps = eps)                  # même alignement que Y_t
}

# 2) Sélectionne 2022Q1..Q4 et calibre δ (GPR = 1ère colonne)
get_2022_eps_gpr <- function(eps, dates, p_lags) {
  # dates alignées sur Y_t: indices (p+1):T
  dt_eff <- dates[(p_lags+1):length(dates)]
  yr <- year(dt_eff); qr <- quarter(dt_eff)
  idx <- which(yr == 2022 & qr %in% 1:4)
  if (length(idx) == 0) stop("Aucune observation 2022Q1..Q4 trouvée dans l'échantillon.")
  list(idx = idx, eps_gpr = eps[idx, 1])
}

# 3) Injection d'une SEQUENCE δ_h (générique) vers Z
inject_sequence_into_Z <- function(var_kernel_rds, terms_df, delta_seq, impulse_idx = 1, jitter = 1e-10) {
  kern <- readRDS(var_kernel_rds)
  Psi_draws  <- kern$Psi_draws
  Sigma_kept <- kern$Sigma_kept
  variables  <- kern$variables
  H <- kern$h; k <- length(variables); M <- dim(Psi_draws)[4]
  stopifnot(length(delta_seq) >= H + 1)
  psiZ_draws <- matrix(0.0, nrow = H + 1, ncol = M,
                       dimnames = list(horizon = 0:H, draw = NULL))
  for (m in 1:M) {
    Sigma <- Sigma_kept[,,m]
    L <- tryCatch(t(chol(Sigma)),
                  error = function(e) t(chol(Sigma + diag(jitter, k))))
    delta_basis <- L[, impulse_idx]  # P e1
    Psi_m <- Psi_draws[,,,m, drop=FALSE][,,,1]
    irfY <- matrix(0.0, nrow = H + 1, ncol = k)
    for (hh in 0:H) {
      acc <- rep(0.0, k)
      for (q in 0:hh) {
        acc <- acc + (Psi_m[hh+1 - q, , ] %*% delta_basis) * delta_seq[q+1]
      }
      irfY[hh+1, ] <- acc
    }
    accZ <- numeric(H + 1)
    for (r in 1:nrow(terms_df)) {
      lag_r  <- terms_df$lag[r]
      col_r  <- terms_df$col[r]
      beta_r <- terms_df$beta[r]
      contrib <- if (lag_r == 0) irfY[, col_r] else c(rep(0, lag_r), irfY[1:(H + 1 - lag_r), col_r])
      accZ <- accZ + beta_r * contrib
    }
    psiZ_draws[, m] <- accZ
  }
  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(psiZ_draws, 1, quantile, probs = qs, na.rm = TRUE))
  psiZ_bands <- data.frame(
    horizon = 0:H,
    lower90 = qmat[,1], lower68 = qmat[,2], median = qmat[,3],
    upper68 = qmat[,4], upper90 = qmat[,5]
  )
  list(psiZ_draws = psiZ_draws, psiZ_bands = psiZ_bands)
}

# 4) Variances conditionnelles pour une SEQUENCE (exclut tous les q où δ_q est imposé)
compute_s2_sequence_one_draw <- function(Psi_draw_one, Sigma_one, terms_df, beta_vec, sigma_eta2, delta_seq) {
  H <- dim(Psi_draw_one)[1] - 1
  variables <- dimnames(Psi_draw_one)$variable
  sel <- build_selection_from_terms(terms_df, variables)
  S_list <- sel$S_list; Lmax <- sel$Lmax
  m <- nrow(S_list[[1]])
  s2      <- numeric(H + 1)
  s2_cond <- numeric(H + 1)
  shock_q <- which(abs(delta_seq) > 0) - 1L  # positions q imposées
  for (hh in 0:H) {
    Sigma_s      <- matrix(0.0, m, m)
    Sigma_s_cond <- matrix(0.0, m, m)
    for (q in 0:hh) {
      G_hq <- build_G_hq(hh, q, S_list, Lmax, Psi_draw_one)
      Sigma_s <- Sigma_s + G_hq %*% Sigma_one %*% t(G_hq)
      if (!(q %in% shock_q)) {
        Sigma_s_cond <- Sigma_s_cond + G_hq %*% Sigma_one %*% t(G_hq)
      }
    }
    s2[hh + 1]      <- as.numeric(t(beta_vec) %*% Sigma_s      %*% beta_vec) + sigma_eta2
    s2_cond[hh + 1] <- as.numeric(t(beta_vec) %*% Sigma_s_cond %*% beta_vec) + sigma_eta2
  }
  list(s2 = s2, s2_cond = s2_cond)
}

# ===================== CONSTRUIT LES δ POUR 2022 ===================== #
# Posterior means pour B et Σ (pour extraire eps historiques)
Sigma_bar <- apply(Sigma_kept, c(1,2), mean)
B_bar     <- apply(B_kept,   c(1,2), mean)

shocks <- compute_structural_shocks(as.matrix(DT[, ..i_var_str]), B_bar, Sigma_bar, p_lags)
sel2022 <- get_2022_eps_gpr(shocks$eps, dates, p_lags)
eps_2022_gpr <- as.numeric(sel2022$eps_gpr)

# (A) Choc unique 2022Q1 (calibré) -> δ_0 = s * ε^{GPR}_{2022Q1}
s_scale_single <- 1.0
delta_single <- rep(0, H + 1)
if (length(eps_2022_gpr) >= 1) delta_single[1] <- s_scale_single * eps_2022_gpr[1]

# (B) Séquence 2022Q1-Q4 -> δ_q = s * ε^{GPR}_{2022Q1+q}, q=0..3
s_scale_seq <- 1.0
delta_seq <- rep(0, H + 1)
Lseq <- min(4, length(eps_2022_gpr))
if (Lseq > 0) delta_seq[1:Lseq] <- s_scale_seq * eps_2022_gpr[1:Lseq]

# ===================== INJECTION SEQUENCE -> Z & VARIANCES ===================== #
# (A) Choc unique calibré 2022Q1
inj_single <- inject_sequence_into_Z(file.path(out_dir, "var_kernel.rds"),
                                     terms_df, delta_single, impulse_idx = impulse_ix, jitter = chol_jitter)
saveRDS(inj_single$psiZ_draws, file = file.path(out_dir, "psiZ_draws_REPLAY2022_single.rds"))

# (B) Séquence 2022Q1-Q4
inj_seq <- inject_sequence_into_Z(file.path(out_dir, "var_kernel.rds"),
                                  terms_df, delta_seq, impulse_idx = impulse_ix, jitter = chol_jitter)
saveRDS(inj_seq$psiZ_draws, file = file.path(out_dir, "psiZ_draws_REPLAY2022_sequence.rds"))

# Variances conditionnelles pour (A) et (B)
s2_delta_single <- matrix(NA_real_, nrow = H + 1, ncol = M, dimnames = list(horizon = 0:H, draw = NULL))
s2_delta_seq    <- matrix(NA_real_, nrow = H + 1, ncol = M, dimnames = list(horizon = 0:H, draw = NULL))
for (m in 1:M) {
  Psi_m <- Psi_draws[,,, m, drop = FALSE][,,,1]
  Sigma_m <- Sigma_kept[,, m]
  r1 <- compute_s2_sequence_one_draw(Psi_m, Sigma_m, terms, beta_vec, sigma_eta2, delta_single)
  r2 <- compute_s2_sequence_one_draw(Psi_m, Sigma_m, terms, beta_vec, sigma_eta2, delta_seq)
  s2_delta_single[, m] <- pmax(r1$s2_cond, 0)
  s2_delta_seq[, m]    <- pmax(r2$s2_cond, 0)
}

# ===================== GIRF-PD pour REPLAY 2022 ===================== #
girf_pd_single_REPLAY2022 <- compute_girf_pd_single(
  psiZ_draws     = readRDS(file.path(out_dir, "psiZ_draws_REPLAY2022_single.rds")),
  mu_draws       = mu_draws,
  s2_draws       = s2_draws,
  s2_delta_draws = s2_delta_single,
  p = p0, rho = rho0
)
fwrite(girf_pd_single_REPLAY2022$bands, file = file.path(out_dir, "girf_pd_bands_REPLAY2022_SINGLE.csv"))
saveRDS(girf_pd_single_REPLAY2022,       file = file.path(out_dir, "girf_pd_single_REPLAY2022.rds"))

girf_pd_seq_REPLAY2022 <- compute_girf_pd_single(
  psiZ_draws     = readRDS(file.path(out_dir, "psiZ_draws_REPLAY2022_sequence.rds")),
  mu_draws       = mu_draws,
  s2_draws       = s2_draws,
  s2_delta_draws = s2_delta_seq,
  p = p0, rho = rho0
)
fwrite(girf_pd_seq_REPLAY2022$bands, file = file.path(out_dir, "girf_pd_bands_REPLAY2022_SEQUENCE.csv"))
saveRDS(girf_pd_seq_REPLAY2022,       file = file.path(out_dir, "girf_pd_sequence_REPLAY2022.rds"))

# ----------------------- PLOTS RAPIDES (replay) ----------------------- #
# print(plot_bands_generic(inj_single$psiZ_bands, "ψ_Z(h)", "REPLAY 2022 — Single shock (2022Q1 calibrated) → Z"))
# print(plot_bands_generic(girf_pd_single_REPLAY2022$bands, "Δ PD (pp)", "REPLAY 2022 — Single shock → PD"))
# print(plot_bands_generic(inj_seq$psiZ_bands, "ψ_Z(h)", "REPLAY 2022 — Sequence (Q1..Q4) → Z"))
# print(plot_bands_generic(girf_pd_seq_REPLAY2022$bands, "Δ PD (pp)", "REPLAY 2022 — Sequence → PD"))

# ----------------------- LOG FINAL ----------------------- #
message(sprintf("Tirages stables: %d / %d (%.1f%%).", sum(keep), length(keep), 100*mean(keep)))
cat("Outputs écrits dans le dossier:", normalizePath(out_dir), "\n")
