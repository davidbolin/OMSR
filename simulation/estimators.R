## ---------------------------------------------------------------------
## The estimators compared in the simulation study
##
## All operator-matched / additive fits (and the Forced-fixed oracle) are
## rspde_lme fits.  Only the RSR (Reich-Hodges-Zadnik
## restricted spatial regression) estimator needs extra machinery: its
## point estimate of the fixed effect equals OLS, but the reported
## standard error is the RSR model-based variance sigma_e_hat^2 (X'X)^{-1}
## with the nugget sigma_e_hat^2 estimated by restricted maximum likelihood.  
## ---------------------------------------------------------------------

library(rSPDE)
library(mgcv)


## Safe wrapper: time a fit, swallow errors/warnings, return NULL on fail.
.try_time <- function(expr) {
  t <- system.time(
    val <- tryCatch(suppressWarnings(force(expr)), error = function(e) NULL)
  )["elapsed"]
  list(val = val, time = as.numeric(t))
}

.row <- function(estimator, param, estimate, se, truth, time, ok = TRUE) {
  data.frame(estimator = estimator, param = param,
             estimate = estimate, se = se, truth = truth,
             time = time, ok = ok, stringsAsFactors = FALSE)
}

.na_row <- function(estimator, param, truth, time) {
  .row(estimator, param, NA_real_, NA_real_, truth, time, ok = FALSE)
}

## Helpers to pull a fixed / forcing coefficient out of rspde_lme 
.fixed_coef <- function(fit, name) {
  est <- tryCatch(unname(fit$coeff$fixed_effects[name]), error = function(e) NA)
  se  <- tryCatch(unname(fit$std_errors$std_fixed[name]),  error = function(e) NA)
  c(est = est, se = se)
}

.forcing_coef <- function(fit) {            
  est <- tryCatch(unname(fit$coeff$random_effects["beta_x1"]),
                  error = function(e) NA)
  se  <- tryCatch(unname(fit$std_errors$std_random["beta_x1"]),
                  error = function(e) NA)
  c(est = est, se = se)
}

## ============================ estimators =============================

## 1. Additive-GLS: y ~ X with a Whittle-Matern random effect (alpha = 2).
fit_additive <- function(ds, world, models) {
  r <- .try_time(rspde_lme(y ~ X, data = ds$data, model = models$model1,
                           loc = c("x1", "x2"),
                           model_options = list(fix_alpha = 2),
                           optim_method = "BFGS"))
  if (is.null(r$val)) return(.na_row("Additive-GLS", "beta", ds$truth$beta0, r$time))
  cc <- .fixed_coef(r$val, "X")
  .row("Additive-GLS", "beta", cc["est"], cc["se"], ds$truth$beta0, r$time)
}

## 2. BW-empirical: pre-smooth X with a Gaussian kernel of fixed bandwidth equal 
##    to the field's practical range, then additive GLS.  
fit_bw <- function(ds, world, models) {
  d <- ds$data
  d$Sx <- gaussian_kernel_smooth(ds$loc, d$X, ds$loc, world$range_u)
  r <- .try_time(rspde_lme(y ~ Sx, data = d, model = models$model1,
                           loc = c("x1", "x2"),
                           model_options = list(fix_alpha = 2),
                           optim_method = "BFGS"))
  if (is.null(r$val)) return(.na_row("BW-empirical", "beta", ds$truth$beta0, r$time))
  cc <- .fixed_coef(r$val, "Sx")
  .row("BW-empirical", "beta", cc["est"], cc["se"], ds$truth$beta0, r$time)
}

## 3. Spatial+: regress X on a thin-plate spline of location, use the
##    residual covariate, then additive GLS (Dupont-Wood-Augustin 2022).
fit_splus <- function(ds, world, models, k = 30) {
  d <- ds$data
  gx <- tryCatch(mgcv::gam(X ~ s(x1, x2, k = k), data = d, method = "REML"),
                 error = function(e) NULL)
  if (is.null(gx)) return(.na_row("Spatial+", "beta", ds$truth$beta0, 0))
  d$Rx <- d$X - as.numeric(predict(gx, d))
  r <- .try_time(rspde_lme(y ~ Rx, data = d, model = models$model1,
                           loc = c("x1", "x2"),
                           model_options = list(fix_alpha = 2),
                           optim_method = "BFGS"))
  if (is.null(r$val)) return(.na_row("Spatial+", "beta", ds$truth$beta0, r$time))
  cc <- .fixed_coef(r$val, "Rx")
  .row("Spatial+", "beta", cc["est"], cc["se"], ds$truth$beta0, r$time)
}

## 4. RSR: restricted spatial regression (Reich-Hodges-Zadnik).  The
##    spatial random effect is restricted to the orthogonal complement of
##    the fixed-effect design, P_perp = I - X(X'X)^{-1}X'. 
fit_rsr <- function(ds, world, models) {
  tt <- system.time({
    val <- tryCatch({
      Xd <- cbind(1, ds$data$X); y <- ds$data$y; n <- length(y)
      XtXi <- solve(crossprod(Xd))
      bhat <- as.numeric(XtXi %*% crossprod(Xd, y))      # = OLS
      Pp   <- diag(n) - Xd %*% XtXi %*% t(Xd)            # restriction P_perp
      r    <- as.numeric(y - Xd %*% bhat)                # = P_perp y
      C0 <- world$C0; G <- world$G; A <- ds$A
      ## profile over kappa; inner profile over (1/tau^2, sigma_e^2) is
      ## cheap once B = P_perp (A K^{-1} C0 K^{-1} A') P_perp is diagonalised.
      prof <- function(lk) {
        K  <- exp(lk)^2 * C0 + G
        Wt <- Matrix::solve(K, Matrix::t(A))             # K^{-1} A'  (m x n)
        Su <- as.matrix(Matrix::crossprod(Wt, C0 %*% Wt))# A K^{-1} C0 K^{-1} A'
        B  <- Pp %*% Su %*% Pp
        eg <- eigen(B, symmetric = TRUE)
        d  <- pmax(eg$values, 0); z <- as.numeric(crossprod(eg$vectors, r))
        inner <- function(p) {                           # p = (log 1/tau^2, log sigma_e^2)
          ev <- d * exp(p[1]) + exp(p[2])
          0.5 * (sum(log(ev)) + sum(z^2 / ev))
        }
        o <- optim(c(log(1 / world$tau^2), log(0.05)), inner, method = "BFGS")
        list(nll = o$value, s2 = exp(o$par[2]))
      }
      ok <- optimize(function(lk) prof(lk)$nll, log(c(0.05, 20)))
      s2e <- prof(ok$minimum)$s2
      list(est = bhat[2], se = sqrt(s2e * XtXi[2, 2]))
    }, error = function(e) NULL)
  })["elapsed"]
  if (is.null(val)) return(.na_row("RSR", "beta", ds$truth$beta0, as.numeric(tt)))
  .row("RSR", "beta", val$est, val$se, ds$truth$beta0, as.numeric(tt))
}

## 5. Forced-fixed: forced model with the TRUE operator (oracle).  
fit_forced_fixed <- function(ds, world, models) {
  r <- .try_time(rspde_lme(y ~ 1, data = ds$data, model = models$model2,
                           loc = c("x1", "x2"),
                           model_options = list(fix_alpha = 2,
                                                fix_kappa = world$kappa,
                                                fix_tau   = world$tau)))
  if (is.null(r$val)) return(.na_row("Forced-fixed", "beta_fc", ds$truth$beta0, r$time))
  cc <- .forcing_coef(r$val)
  .row("Forced-fixed", "beta_fc", cc["est"], cc["se"], ds$truth$beta0, r$time)
}

## 6. Forced: forced model with the operator (kappa, tau) and beta_fc estimated
fit_forced <- function(ds, world, models) {
  r <- .try_time(rspde_lme(y ~ 1, data = ds$data, model = models$model2,
                           loc = c("x1", "x2"),
                           model_options = list(fix_alpha = 2)))
  if (is.null(r$val)) return(.na_row("Forced", "beta_fc", ds$truth$beta0, r$time))
  cc <- .forcing_coef(r$val)
  .row("Forced", "beta_fc", cc["est"], cc["se"], ds$truth$beta0, r$time)
}

## 7. Hybrid: covariate enters both pointwise and as forcing 
fit_hybrid <- function(ds, world, models) {
  r <- .try_time(rspde_lme(y ~ X, data = ds$data, model = models$model2,
                           loc = c("x1", "x2"),
                           model_options = list(fix_alpha = 2),
                           optim_method = "BFGS"))
  if (is.null(r$val)) {
    return(rbind(
      .na_row("Hybrid (pw)", "beta_pw", ds$truth$beta_pw, r$time),
      .na_row("Hybrid (fc)", "beta_fc", ds$truth$beta_fc, r$time)))
  }
  pw <- .fixed_coef(r$val, "X")
  fc <- .forcing_coef(r$val)
  rbind(
    .row("Hybrid (pw)", "beta_pw", pw["est"], pw["se"], ds$truth$beta_pw, r$time),
    .row("Hybrid (fc)", "beta_fc", fc["est"], fc["se"], ds$truth$beta_fc, r$time))
}

## ---------------------------------------------------------------------
## fit_all(): build the data-dependent models once and run every
## estimator on a single dataset.  Returns a tidy data.frame.
## ---------------------------------------------------------------------
fit_all <- function(ds, world) {
  models <- list(
    model1 = spde.matern.operators(mesh = world$mesh, alpha = 2),
    model2 = hybrid.spde(mesh = world$mesh, alpha = 2, X = ds$Xnodes)
  )
  out <- rbind(
    fit_additive(ds, world, models),
    fit_bw(ds, world, models),
    fit_splus(ds, world, models),
    fit_rsr(ds, world, models),
    fit_forced_fixed(ds, world, models),
    fit_forced(ds, world, models),
    fit_hybrid(ds, world, models)
  )
  out$scenario <- ds$scenario
  out$n        <- ds$n
  rownames(out) <- NULL
  out
}
