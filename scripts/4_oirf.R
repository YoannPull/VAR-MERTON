# ================================================================
#   BVAR (NIW) + OIRF (Cholesky) + ψZ (Koop) + μ_{t+h} exacte
#   + s^2(h), s^2_δ(h) + GIRF-PD + e_{1,t}
#   -- Script unique & optimisé --
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mniw)
  library(parallel)
  library(ggplot2)
  library(lmtest)     # DIAG (résidus)
  library(tseries)    # DIAG (résidus)
  library(matrixStats)
})

theme_set(theme_minimal(base_size = 14))
options(stringsAsFactors = FALSE)
set.seed(123)  # reproductibilité globale

# ----------------------- PARAMÈTRES GÉNÉRAUX ----------------------- #
data_var_path   <- "data/processed/data_var_for_model.csv"
best_model_path <- "scripts/model/best_model.rds"   # modèle lm de Z
out_dir         <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Variables du VAR (ordre fixé, GPR en premier pour identification récursive)
allowed_vect <- c("vix","log_sp500_real","log_oil_real","log_hours_pc",
                  "log_gdp_pc","nfci")
i_var_str <- c("log_GPRD", allowed_vect)

# Paramètres VAR + IRF
p_lags      <- 2          # ordre VAR
nrep        <- 20000      # nb de tirages postérieurs
seed        <- 123
h           <- 12         # horizons (0..h)
impulse_ix  <- 1          # choc sur GPR (1ère variable)
chol_jitter <- 1e-10      # robustesse Cholesky

# (Optionnel) standardiser Z sur l'échelle du training
scale_Z <- FALSE

# ---- Paramètres GIRF-PD (forme fermée, single (p, rho)) ---- #
p0   <- 0.007954558   # PD inconditionnelle
rho0 <- 0.02383827    # corrélation d'actif
stopifnot(p0 > 0 && p0 < 1, rho0 > 0 && rho0 < 1)

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

extract_c_A_list <- function(B, k, p) {
  c_vec <- as.numeric(B[1, ])
  A_stack <- t(B[-1, , drop = FALSE])                # k x (k*p)
  A_list <- vector("list", p)
  for (i in 1:p) A_list[[i]] <- A_stack[, ((i-1)*k + 1):(i*k), drop = FALSE]
  list(c = c_vec, A = A_list)
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
DT_raw <- fread(data_var_path)
# Conserve la colonne de date si elle existe
if ("Date" %in% names(DT_raw)) {
  dates <- as.Date(DT_raw$Date)
} else {
  dates <- NULL
}

# Ne garder QUE les colonnes du VAR dans le bon ordre (évite l'erreur de dimnames)
missing_cols <- setdiff(i_var_str, names(DT_raw))
if (length(missing_cols) > 0) stop("Colonnes manquantes: ", paste(missing_cols, collapse = ", "))
DT <- copy(DT_raw[, ..i_var_str])
setcolorder(DT, i_var_str)

if (anyNA(DT)) DT <- na.omit(DT)

Y <- ts(as.matrix(DT), start = c(1990, 2), frequency = 4)
variables <- colnames(Y)  # référence unique des noms
k <- ncol(Y)

# ----------------------- ESTIMATION + Ψ (PARALLÈLE) ----------------------- #
res <- simulate_bvar_niw(Y, p = p_lags, nrep = nrep, seed = seed)

keep <- vapply(1:dim(res$B)[3], function(i) is_stable(res$B[,, i], k, p_lags), logical(1))
message(sprintf("Tirages stables: %d / %d (%.1f%%).", sum(keep), length(keep), 100*mean(keep)))
idx_kept <- which(keep)

# --- Sélectionne exactement 10 000 tirages stables (ou tous si < 10k) ---
M_target <- 10000L
if (length(idx_kept) >= M_target) {
  set.seed(seed)
  idx_kept <- sample(idx_kept, M_target)
  message(sprintf("Échantillonnage de %d tirages stables (sur %d).", M_target, sum(keep)))
} else {
  warning(sprintf("Seulement %d tirages stables disponibles — on gardera tout.", length(idx_kept)))
}

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

clusterExport(cl, varlist = c("res","idx_kept","p_lags","h","build_companion","compute_ma_coefficients"),
              envir = environment())

psi_list <- parLapply(cl, idx_kept, function(i) {
  compute_ma_coefficients(res$B[,, i], p = p_lags, H = h)
})

M_kept <- length(idx_kept)
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

B_kept     <- res$B[,, idx_kept]   # (m x k x M_kept)
Sigma_kept <- res$S[,, idx_kept]   # (k x k x M_kept)
saveRDS(list(Psi_draws = psi_draws, B_kept = B_kept, Sigma_kept = Sigma_kept,
             p = p_lags, k = k, h = h, variables = variables),
        file = file.path(out_dir, "var_kernel.rds"))

cat("Noyau VAR (Ψ, B, Σ) -> ", file.path(out_dir, "var_kernel.rds"), "
")

# ==================== DIAG STABILITÉ (rapide) ==================== #
max_root <- function(Bmat, p, k) {
  max(Mod(eigen(build_companion(Bmat, k, p), only.values = TRUE)$values))
}
roots <- vapply(1:dim(B_kept)[3], function(i) max_root(B_kept[,,i], p_lags, k), numeric(1))
stable_share <- mean(roots < 1)
cat(sprintf("Stabilité: %.1f%% stables; médiane |λ|max: %.3f; pct95: %.3f
",
            100*stable_share, median(roots), quantile(roots, .95)))

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
    L <- tryCatch(t(chol(Sigma_m)),
                  error = function(e) t(chol(Sigma_m + diag(jitter, k))))
    delta_o <- L[, impulse_idx]  # choc d’1 écart-type structurel
    for (hh in 0:H) oirf_draws[hh+1, , m] <- Psi_m[hh+1, , ] %*% delta_o
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
cat("OIRF -> ",
    file.path(out_dir, "oirf_draws.rds"), " ; ",
    file.path(out_dir, "oirf_bands.csv"), "
", sep = "")

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
  if (length(betas) == 0) stop("Modèle Z sans prédicteurs.")
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
    if (length(dropped) > 0) warning("Termes ignorés: ", paste(dropped, collapse = ", "))
    terms_df <- terms_df[keep_row, , drop = FALSE]
  }
  if (nrow(terms_df) == 0) stop("Aucun terme Z valide pour l'injection.")
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
    delta_o <- L[, impulse_idx]  # choc unitaire (structurel)
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
cat("Injection OIRF -> Z -> ",
    file.path(out_dir, "psiZ_draws_OIRF.rds"), " ; ",
    file.path(out_dir, "psiZ_bands_OIRF.csv"), "
", sep = "")

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

# --- μ_{t+h} baseline (exacte) --- #
compute_mu_baseline_draws <- function(DT_var, variables, terms_df, beta0, B_kept, p, H) {
  Tn <- nrow(DT_var); k <- length(variables); M <- dim(B_kept)[3]
  Lmax <- max(terms_df$lag)
  if (Tn <= Lmax) stop("Pas assez d'historique pour couvrir Lmax.")
  Y_all <- as.matrix(DT_var[, ..variables])
  Y_hist_for_forecast <- Y_all[(Tn - p + 1):Tn, , drop = FALSE]
  hist_tail <- Y_all[(Tn - Lmax):Tn, , drop = FALSE]  # (Lmax+1) lignes
  
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
          hist_tail[(Lmax + 1) - j, base_j]
        }
        mu_h[hh + 1] <- mu_h[hh + 1] + beta_r * val
      }
    }
    mu_draws[, m] <- mu_h
  }
  mu_draws
}

# --- s^2(h) et s^2_δ(h) --- #
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

# ===================== μ, s^2, s^2_δ (par tirage) ===================== #
kernel_path <- file.path(out_dir, "var_kernel.rds")
stopifnot(file.exists(kernel_path))
kern <- readRDS(kernel_path)
Psi_draws   <- kern$Psi_draws
Sigma_kept  <- kern$Sigma_kept
B_kept      <- kern$B_kept
variables_k <- kern$variables  # identique à 'variables'
H           <- kern$h
M           <- dim(Psi_draws)[4]

co <- coef(best_model)
beta0 <- unname(co["(Intercept)"])
betas <- co[names(co) != "(Intercept)"]
terms <- terms_df
beta_vec <- terms$beta
sigma_eta2 <- (summary(best_model)$sigma)^2

mu_draws <- compute_mu_baseline_draws(DT_var = as.data.table(DT), variables = variables,
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

# ==================== DIAG RÉSIDUS (B médian) ==================== #
Y_mat <- as.matrix(as.data.frame(DT)[, i_var_str])
Y_t <- Y_mat[(p_lags+1):nrow(Y_mat), , drop=FALSE]
X_t <- do.call(cbind, lapply(1:p_lags, function(l) Y_mat[(p_lags+1-l):(nrow(Y_mat)-l), , drop=FALSE]))
X_t <- cbind(1, X_t)
B_bar <- apply(B_kept, c(1,2), median)
U_bar <- Y_t - X_t %*% B_bar

lb_p <- apply(U_bar, 2, function(u) Box.test(as.numeric(u), lag = 12, type = "Ljung-Box")$p.value)
bp_p <- apply(U_bar, 2, function(u) bptest(lm(u ~ X_t[,-1]))$p.value)
jb_p <- apply(U_bar, 2, function(u) tseries::jarque.bera.test(as.numeric(u))$p.value)
resid_diag <- data.frame(variable = colnames(U_bar), p_LjungBox = lb_p, p_BreuschPagan = bp_p, p_JarqueBera = jb_p, row.names = NULL)
fwrite(resid_diag, file.path(out_dir,"diag_residual_tests.csv"))
cat("Résidus -> ", file.path(out_dir,"diag_residual_tests.csv"), "
", sep = "")

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

# ===================== GIRF de la PD (forme fermée) ===================== #
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

girf_pd_single_OIRF <- compute_girf_pd_single(
  psiZ_draws     = readRDS(file.path(out_dir, "psiZ_draws_OIRF.rds")),
  mu_draws       = mu_draws,
  s2_draws       = s2_draws,
  s2_delta_draws = s2_delta_draws,
  p = p0, rho = rho0
)
fwrite(girf_pd_single_OIRF$bands, file = file.path(out_dir, "girf_pd_bands_OIRF_SINGLE.csv"))
saveRDS(girf_pd_single_OIRF,       file = file.path(out_dir, "girf_pd_single_OIRF.rds"))
cat("GIRF-PD (Koop, choc OIRF) -> ",
    file.path(out_dir, "girf_pd_bands_OIRF_SINGLE.csv"), "
", sep = "")

# ==================== Innovations structurelles e_{j,t} ==================== #
compute_structural_innovations <- function(B_kept, Sigma_kept, Y_mat, p_lags, dates_vec = NULL,
                                           agg = c("median","mean")) {
  agg <- match.arg(agg)
  B_typ <- apply(B_kept, c(1,2), if (agg=="median") median else mean)
  Sigma_typ <- apply(Sigma_kept, c(1,2), if (agg=="median") median else mean)
  Sigma_typ <- (Sigma_typ + t(Sigma_typ)) / 2
  Y_t <- Y_mat[(p_lags+1):nrow(Y_mat), , drop = FALSE]
  X_t <- do.call(cbind, lapply(1:p_lags, function(l) Y_mat[(p_lags+1-l):(nrow(Y_mat)-l), , drop = FALSE]))
  X_t <- cbind(1, X_t)
  U_t <- Y_t - X_t %*% B_typ
  L_typ <- tryCatch(t(chol(Sigma_typ)),
                    error = function(e) t(chol(Sigma_typ + diag(1e-10, ncol(Sigma_typ)))))
  E_t <- t(forwardsolve(L_typ, t(U_t)))  # (T_eff x k), Var ≈ I
  if (!is.null(dates_vec)) {
    dates_eff <- as.Date(dates_vec[(p_lags+1):length(dates_vec)])
  } else {
    dates_eff <- seq_len(nrow(E_t))
  }
  list(E = E_t,
       e1_dt = data.table(date = dates_eff, e1 = as.numeric(E_t[,1])))
}

E_out <- compute_structural_innovations(B_kept, Sigma_kept, Y_mat, p_lags, dates_vec = dates, agg = "median")
E_all <- E_out$E
e1_dt <- E_out$e1_dt
fwrite(e1_dt, file.path(out_dir, "e1_series.csv"))
saveRDS(list(E = E_all, e1_dt = e1_dt), file = file.path(out_dir, "struct_innovations.rds"))
cat("e1_t -> ", file.path(out_dir, "e1_series.csv"), "
", sep = "")

# ==================== (Option) Bandes postérieures de e_{1,t} ==================== #
compute_e1_posterior_bands <- function(B_kept, Sigma_kept, Y_mat, p_lags, dates_vec,
                                       qs = c(0.16, 0.50, 0.84), jitter = 1e-10) {
  Y_t <- Y_mat[(p_lags+1):nrow(Y_mat), , drop = FALSE]
  X_t <- do.call(cbind, lapply(1:p_lags, function(l) Y_mat[(p_lags+1-l):(nrow(Y_mat)-l), , drop = FALSE]))
  X_t <- cbind(1, X_t)
  T_eff <- nrow(Y_t); k <- ncol(Y_t); M <- dim(B_kept)[3]
  e1_draws <- matrix(NA_real_, nrow = T_eff, ncol = M)
  for (m in 1:M) {
    U_m <- Y_t - X_t %*% B_kept[,, m]
    L_m <- tryCatch(t(chol(Sigma_kept[,, m])),
                    error = function(e) t(chol(Sigma_kept[,, m] + diag(jitter, k))))
    E_m <- t(forwardsolve(L_m, t(U_m)))
    e1_draws[, m] <- E_m[, 1]
  }
  Q <- t(apply(e1_draws, 1, quantile, probs = qs, na.rm = TRUE))
  data.table(
    date = if (!is.null(dates_vec)) as.Date(dates_vec[(p_lags+1):length(dates_vec)]) else seq_len(T_eff),
    q16 = Q[, 1], q50 = Q[, 2], q84 = Q[, 3]
  )
}

# ==================== (NOUVEAU) e1_t pour TOUS les tirages (matrice T_eff x M) ==================== #
compute_e1_all_draws <- function(B_kept, Sigma_kept, Y_mat, p_lags, dates_vec = NULL, jitter = 1e-10) {
  Y_t <- Y_mat[(p_lags+1):nrow(Y_mat), , drop = FALSE]
  X_t <- do.call(cbind, lapply(1:p_lags, function(l) Y_mat[(p_lags+1-l):(nrow(Y_mat)-l), , drop = FALSE]))
  X_t <- cbind(1, X_t)
  T_eff <- nrow(Y_t); k <- ncol(Y_t); M <- dim(B_kept)[3]
  e1_draws <- matrix(NA_real_, nrow = T_eff, ncol = M)
  for (m in 1:M) {
    U_m <- Y_t - X_t %*% B_kept[,, m]
    L_m <- tryCatch(t(chol(Sigma_kept[,, m])),
                    error = function(e) t(chol(Sigma_kept[,, m] + diag(jitter, k))))
    E_m <- t(forwardsolve(L_m, t(U_m)))
    e1_draws[, m] <- E_m[, 1]
  }
  dates_eff <- if (!is.null(dates_vec)) as.Date(dates_vec[(p_lags+1):length(dates_vec)]) else seq_len(T_eff)
  colnames(e1_draws) <- sprintf("draw_%05d", seq_len(ncol(e1_draws)))
  rownames(e1_draws) <- as.character(dates_eff)
  list(e1_draws = e1_draws, dates_eff = dates_eff)
}

# Calcul + sauvegardes (grand fichier, ~ quelques dizaines de Mo)
all_e1 <- compute_e1_all_draws(B_kept, Sigma_kept, Y_mat, p_lags, dates)
e1_draws_mat <- all_e1$e1_draws
saveRDS(list(e1_draws = e1_draws_mat, dates = all_e1$dates_eff),
        file = file.path(out_dir, "e1_draws_all.rds"), compress = "xz")

# Export "wide" (1 ligne = date, 10k colonnes = tirages)
DT_e1_wide <- as.data.table(e1_draws_mat)
DT_e1_wide[, date := all_e1$dates_eff]
setcolorder(DT_e1_wide, c("date", setdiff(colnames(DT_e1_wide), "date")))
fwrite(DT_e1_wide, file.path(out_dir, "e1_draws_all_wide.csv.gz"))

# Export "long" (optionnel) : 1.2M+ lignes si 120 dates x 10k tirages
# DT_e1_long <- melt(DT_e1_wide, id.vars = "date", variable.name = "draw", value.name = "e1")
# fwrite(DT_e1_long, file.path(out_dir, "e1_draws_all_long.csv.gz"))



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


# Colonnes de tirages
draw_cols <- grep("^draw_", names(DT_e1_wide), value = TRUE)

# Médiane par date (ligne) sur les tirages
e1_median_by_date <- DT_e1_wide[, .(
  median_e1 = rowMedians(as.matrix(.SD), na.rm = TRUE)
), .SDcols = draw_cols]
e1_median_by_date[, date := DT_e1_wide$date]
setcolorder(e1_median_by_date, c("date", "median_e1"))

fwrite(e1_median_by_date, "output/e1_median.csv")
