# --------------------------------------------------------------------
# Build the Section 5 result tables
# Run after fit_models.R.
# --------------------------------------------------------------------

library(rSPDE)
library(Matrix)
library(fmesher)
library(numDeriv)

# delta-method CIs for (beta_pw, beta_star, g_X, beta_eff) on an a=2 hybrid fit
delta_ci <- function(fit, mesh, Xn, fem, x_name = "elev") {
  th <- fit$mle_par_orig
  fx <- fit$coeff$fixed_effects
  beta_pw <- if (x_name %in% names(fx)) as.numeric(fx[[x_name]]) else NA_real_
  se_pw   <- if (x_name %in% names(fx)) as.numeric(fit$std_errors$std_fixed[[x_name]]) else NA_real_
  gfun <- function(theta) {
    kappa <- exp(theta[["kappa"]])
    bfc <- theta[["beta_x1"]]
    st <- standardize_forcing(mesh, Xn, kappa=kappa, beta_fc=bfc, alpha=2, fem=fem)
    c(st$beta_star, st$g_X, st$beta_eff)
  }
  est <- gfun(th)

  H <- numDeriv::hessian(fit$lik_fun, th)
  J <- numDeriv::jacobian(gfun, th)
  se <- sqrt(diag(J %*% solve(H) %*% t(J)))
  data.frame(quantity  = c("beta_pw","beta_star","g_X","beta_eff"),
             estimate  = c(beta_pw, est),
             se        = c(se_pw, se),
             lo        = c(beta_pw, est) - 1.959964*c(se_pw, se),
             hi        = c(beta_pw, est) + 1.959964*c(se_pw, se))
}

# Extract parameters
sc1 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1) NA_real_ else x
}

hyper_row <- function(name, fit) {
  re <- fit$coeff$random_effects
  a  <- sc1(alpha_from_fit(fit))
  k  <- sc1(kappa_from_fit(fit))
  tau <- sc1(tryCatch(re[["tau"]], error = function(e) NA_real_))
  sig <- sc1(fit$coeff$measurement_error)
  data.frame(model=name, alpha=a, kappa=k, kappa_inv_km=1/k,
             practical_range_km = if (is.na(a) || is.na(k)) NA_real_ else sqrt(8*(a-1))/k,
             tau=tau, sigma_e=sig,
             loglik=sc1(fit$loglik))
}

source(file.path("..", "common", "standardize.R"))

regs <- c(COLORADO = "Colorado MAM", NORWAYANNUAL = "Norway annual",
          NORWAYDJF = "Norway DJF")

# Accumulators for the three tables (one entry per region).
cv_rows     <- list()
decomp_rows <- list()
hyper_rows  <- list()

# model display order for the CV table
cv_order <- c("OLS",
              "Additive_est","Additive_a2","Forced_est","Forced_a2",
              "Hybrid_est","Hybrid_a2","TwoScale_a2","Coast_est","Coast_a2",
              "BW_est","BW_a2","SpatPlus_a2")

for(tag in names(regs)){
  load(sprintf("results/fits_%s.RData", tag))
  mesh <- p$mesh; Xn <- p$X_nodes
  fem  <- fmesher::fm_fem(mesh)

  scores <- as.data.frame(pr$scores)
  write.csv(scores, sprintf("results/cv_scores_%s.csv", tag), row.names = FALSE)

  # Decomposition + delta-method CIs (a=2 hybrid-family fits)
  decomp_fits <- intersect(c("Forced_a2","Hybrid_a2","Coast_a2"), names(fits))
  decomp <- do.call(rbind, lapply(decomp_fits, function(m) {
    dc <- delta_ci(fits[[m]], mesh, Xn, fem)
    cbind(model = m, dc)
  }))
  write.csv(decomp, sprintf("%s/decomp_%s.csv", "results/", tag), row.names = FALSE)

  # Hyperparameter table
  hyper <- do.call(rbind, Map(hyper_row, names(fits), fits))
  rownames(hyper) <- NULL
  write.csv(hyper, sprintf("%s/hyper_%s.csv", "results/", tag), row.names = FALSE)

  # accumulate rows for the three printed tables

  ## Table 1
  cv <- scores
  cv <- cv[order(match(cv$Model, cv_order)), ]
  cv_rows[[tag]] <- data.frame(region = regs[[tag]], cv, row.names = NULL)

  ## Table 2
  b_ols <- as.numeric(fits$OLS$coeff$fixed_effects["elev"])
  b_add <- as.numeric(fits$Additive_a2$coeff$fixed_effects["elev"])
  hyb   <- fits$Hybrid_a2
  st    <- standardized_summary(hyb, mesh, Xn, fem = fem)
  th    <- hyb$mle_par_orig
  g <- function(theta) { k <- exp(theta[["kappa"]]); b <- theta[["beta_x1"]]
    s <- standardize_forcing(mesh, Xn, kappa = k, beta_fc = b, alpha = 2, fem = fem)
    c(s$beta_star, s$beta_eff) }
  V  <- solve(numDeriv::hessian(hyb$lik_fun, th)); J <- numDeriv::jacobian(g, th)
  se <- sqrt(diag(J %*% V %*% t(J)))
  se_pw <- as.numeric(hyb$std_errors$std_fixed["elev"])
  lr    <- 2 * (hyb$loglik - fits$Additive_a2$loglik)
  decomp_rows[[tag]] <- data.frame(
    region       = regs[[tag]],
    beta_ols     = b_ols,
    beta_add     = b_add,
    beta_pw      = st$beta_pw,   beta_pw_se   = se_pw,
    beta_star    = st$beta_star, beta_star_se = se[1],
    g_X          = st$g_X,
    beta_eff     = st$beta_eff,  beta_eff_se  = se[2],
    LR           = lr,
    row.names = NULL)

  # hyperparameters + CV RMSE at alpha=2 vs estimated alpha.
  # Two-scale has only an alpha=2 fit, so its alpha_est / rmse_est are NA.
  rownames(hyper) <- hyper$model
  rownames(cv)    <- cv$Model
  mods    <- c("Additive", "Forced", "Hybrid", "TwoScale", "Coast")
  mod_lab <- c(Additive = "Additive", Forced = "Forced", Hybrid = "Hybrid",
               TwoScale = "Two-scale", Coast = "Coast")
  present <- mods[paste0(mods, "_a2") %in% rownames(hyper)]
  hyper_rows[[tag]] <- do.call(rbind, lapply(present, function(m) {
    a2 <- paste0(m, "_a2"); es <- paste0(m, "_est")
    data.frame(
      region       = regs[[tag]],
      model        = mod_lab[[m]],
      alpha_est    = if (es %in% rownames(hyper)) hyper[es, "alpha"] else NA_real_,
      kappa_inv_km = hyper[a2, "kappa_inv_km"],
      range_km     = hyper[a2, "practical_range_km"],
      sigma_e      = hyper[a2, "sigma_e"],
      rmse_a2      = cv[a2, "rmse"],
      rmse_est     = if (es %in% rownames(cv)) cv[es, "rmse"] else NA_real_,
      row.names = NULL)
  }))
}


# Assemble and print the three tables as data frames
cv_table     <- do.call(rbind, cv_rows)
rownames(cv_table)     <- NULL
decomp_table <- do.call(rbind, decomp_rows)
rownames(decomp_table) <- NULL
hyper_table  <- do.call(rbind, hyper_rows)
rownames(hyper_table)  <- NULL

# round numeric columns for display
round_num <- function(df, d = 3) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = d)
  df
}

hyper_disp <- transform(hyper_table,
                        alpha_est    = round(alpha_est, 2),
                        kappa_inv_km = round(kappa_inv_km, 0),
                        range_km     = round(range_km, 0),
                        sigma_e      = round(sigma_e, 2),
                        rmse_a2      = round(rmse_a2, 3),
                        rmse_est     = round(rmse_est, 3))
cat("\n=== Table 1: fitted hyperparameters and effect of estimating alpha ===\n")
cat("(rmse_a2, rmse_est = 10-fold CV RMSE at alpha=2 and at the estimated alpha)\n")
print(hyper_disp, row.names = FALSE)

cat("\n=== Table 2: scale-free decomposition (delta-method SEs) ===\n")
print(round_num(decomp_table, 3), row.names = FALSE)

# Full per-model cross-validation scores (10-fold).
cat("\n=== Cross-validation scores (10-fold, per model) ===\n")
print(round_num(cv_table, 3), row.names = FALSE)
