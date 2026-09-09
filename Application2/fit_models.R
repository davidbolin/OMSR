# ---------------------------------------------------------------------
# Fit the NO2 models
# ---------------------------------------------------------------------
library(Matrix)
library(fmesher)
library(INLA)
library(rSPDE)

set.seed(1)
NSTART    <- 2       # multi-start restarts for the operator models
PREC_BETA <- 1e-5    # weakly-informative prior precision on fixed effects

d <- readRDS("data/no2_data.rds")
invisible(list2env(d, environment())) # copy the data to the environment


# LUR supervised forward selection
# Add the candidate that most improves the adjusted R^2 of an OLS regression
# of y on the covariates chosen so far, but only if its coefficient
# has the physically expected sign and the gain exceeds 0.01.
fwd_select <- function(cols){
  sel<-character(0)                          # covariates selected so far
  aR2<- -Inf                                 # their adjusted R^2
  rem<-cols                                  # remaining candidates
  repeat {
    bv<-NA                                   # best covariate to add this round
    ba<- if(length(sel)) aR2 + 0.01 else -Inf
    for (v in rem){
      f<-lm(y~., data=as.data.frame(cbind(y = y, cand[, c(sel,v), drop=FALSE])))
      cf<-coef(f)[v]
      a<-summary(f)$adj.r.squared
      if(!is.na(cf) && sign(cf) == signs[[v]] && a>ba){   # right sign and best improvement
        ba<-a
        bv<-v
      }
    }
    if(is.na(bv)) break
    sel<-c(sel, bv)
    rem<-setdiff(rem, bv)
    aR2<-ba
  }
  sel
}
bg_cols  <- grep("urban_|natural_|dist_coast", colnames(cand), value=TRUE)
sel_full <- fwd_select(colnames(cand))       # covariates chosen for the LUR models
lur_rhs  <- paste(sel_full, collapse = " + ")
N_CAND   <- ncol(cand)                        # candidates offered to the LUR selection
N_LUR    <- length(sel_full)
cat(sprintf("LUR forward selection kept %d/%d covariates: %s\n",
            N_LUR, N_CAND, paste(sel_full, collapse = ", ")))

# INLA building blocks
A_obs <- fm_basis(mesh, loc = loc_km)
matern_model <- rspde.matern(mesh = mesh, nu = NU, parameterization = "matern")
A_field <- rspde.make.A(mesh = mesh, loc = loc_km, nu = NU)
field_idx <- rspde.make.index(name="field", mesh = mesh, nu = NU)
cdf <- as.data.frame(cand[, unique(c(sel_full, bg_cols, "own_traffic")), drop = FALSE])
cb_all <- c(list(Intercept = 1, own_ind = own_ind, own_pop = own_pop), as.list(cdf))

# prior calibration from a field-only fit
cat("Calibrating field priors (field-only INLA fit)...\n")
stk0 <- inla.stack(data=list(y = y), A = list(A_field, 1),
                   effects = list(field_idx, cb_all), tag = "s")
f0 <- inla(y ~ -1 + Intercept + f(field, model = matern_model),
           family="gaussian", data = inla.stack.data(stk0),
           control.predictor=list(A = inla.stack.A(stk0)),
           control.fixed=list(prec = PREC_BETA, prec.intercept = 0),
           control.inla=list(int.strategy = "eb"))

s0 <- summary(rspde.result(f0, "field", matern_model))

# posterior mean of the first summary row whose name matches pattern `p`
g1<-function(s,p){
  i<-grep(p,rownames(s), ignore.case = TRUE)
  if(length(i)) s[i[1], "mean"] else NA_real_
}
range0 <- g1(s0,"range")           # calibrated field range and marginal sd, used
sd0 <- g1(s0,"std|sigma|dev")      # to set the priors/starts for the operator fields
kappa0 <- sqrt(8*NU)/range0        # SPDE kappa
tau0 <- 1/(sd0*sqrt(4*pi)*kappa0)  # SPDE tau
th_mat <- f0$mode$theta            # warm start for the plain-field (matern) models
Xall <- cbind(Xtraf, Xpt, Xpop)    # three forcing channels

cat("Building operator-field model objects...\n")
hybM <- rspde.hybrid.matern(mesh = mesh, X = Xall,
                            prior.tau = list(mean = log(tau0), prec = 5),
                            prior.kappa = list(mean = log(kappa0),prec = 5),
                            prior.beta_x = list(mean = rep(0,3),prec = rep(0.001, 3)),
                            start.ltau = log(tau0), start.lkappa = log(kappa0),
                            start.beta_x = rep(0,3))
km0 <- sqrt(8*NU)/15        # starting kappa_mu (forcing range ~15 km)

# the two-scale model
hyb2 <- rspde.hybrid.matern(mesh = mesh, X = Xall, separate_kappa_mu = TRUE,
                            prior.tau = list(mean = log(tau0), prec = 5),
                            prior.kappa = list(mean = log(kappa0),prec=5),
                            prior.kappa_mu = list(mean = log(km0),prec = 2),
                            prior.beta_x = list(mean = rep(0, 3), prec = rep(0.001, 3)),
                            start.ltau = log(tau0), start.lkappa = log(kappa0),
                            start.lkappa_mu = log(km0), start.beta_x = rep(0,3))

# build a single-source operator field for covariate Xn with prior range r0
mkh <- function(Xn, r0){
  k0<-sqrt(8*NU)/r0
  t0<-1/(sd0*sqrt(4*pi)*k0)
  rspde.hybrid.matern(mesh = mesh, X = as.matrix(Xn),
                      prior.tau = list(mean = log(t0), prec = 5),
                      prior.kappa = list(mean = log(k0), prec = 3),
                      prior.beta_x = list(mean = 0, prec = 0.001),
                      start.ltau = log(t0), start.lkappa = log(k0),
                      start.beta_x = 0)
}
hyb_t <- mkh(Xtraf, 5)
hyb_i <- mkh(Xpt, 15)
hyb_p <- mkh(Xpop, 8)

# like mkh but for a multi-column X (several sources sharing one field)
mkh_multi <- function(Xn, r0){
  Xn <- as.matrix(Xn)
  p <- ncol(Xn)
  k0 <- sqrt(8*NU)/r0
  t0 <- 1/(sd0*sqrt(4*pi)*k0)
  rspde.hybrid.matern(mesh = mesh, X = Xn, prior.tau = list(mean = log(t0), prec = 5),
                      prior.kappa = list(mean = log(k0),prec = 3),
                      prior.beta_x = list(mean = rep(0,p),prec = rep(0.001, p)),
                      start.ltau = log(t0), start.lkappa = log(k0),
                      start.beta_x = rep(0,p))
}
hyb_loc <- mkh_multi(cbind(Xtraf, Xpop), 2)
hyb_ind <- mkh_multi(Xpt, 20)

# fitting helpers
PW <- "own_traffic + own_ind + own_pop"    # the three pointwise source covariates

# one INLA fit with the shared control options; `th` is an optional warm start
fit_inla <- function(form, stk, th = NULL)
  inla(form, family = "gaussian", data = inla.stack.data(stk),
       control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
       control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
       control.mode      = if (is.null(th)) NULL else list(theta = th, restart = TRUE),
       control.fixed     = list(prec = PREC_BETA, prec.intercept = 0),
       control.inla      = list(int.strategy = "eb"))

mlik <- function(f) as.numeric(f$mlik)[1]

# fit once, then NSTART perturbed restarts, keeping the best marginal likelihood.
fit_multistart <- function(form, stk) {
  best <- fit_inla(form, stk)
  ml <- mlik(best)
  th0 <- best$mode$theta
  if (NSTART > 0 && length(th0)) for (s in 1:NSTART) {
    set.seed(7 * s)
    fs <- tryCatch(fit_inla(form, stk, th = th0 + rnorm(length(th0), sd = 1.5)),
                   error = function(e) NULL)
    if (!is.null(fs)) {
      ml <- c(ml, mlik(fs))
      if (mlik(fs) > mlik(best) + 1e-3) best <- fs
    }
  }
  list(fit = best, spread = max(ml) - min(ml))
}

# build the fixed-effects formula: Intercept [+ covariates] [+ field terms]
mk_form <- function(cov, field = "") {
  rhs <- if (nzchar(cov) && cov != "1") paste("Intercept +", cov) else "Intercept"
  as.formula(paste("y ~ -1 +", rhs, field))
}

# leave-one-out group-CV RMSE and log-score
gcv <- function(f) {
  g <- inla.group.cv(f, num.level.sets = M_LGOCV)
  if (is.null(g$cv)) return(c(rmse = NA, logs = NA))
  c(rmse = sqrt(mean((y - g$mean)^2, na.rm = TRUE)), logs = -mean(log(g$cv), na.rm = TRUE))
}

matern_range <- function(f) {
  rr <- rspde.result(f, "field", matern_model)
  if (!is.null(rr$summary.range)) rr$summary.range$mean else NA
}

cgen_range <- function(f, pattern, which = 1) {
  i <- grep(pattern, rownames(f$summary.hyperpar))
  if (length(i) >= which) sqrt(8 * NU) / exp(f$summary.hyperpar[i[which], "mean"]) else NA
}

# record a fitted model: store it and add its comparison-table row.
fits <- list()
results <- list()
record <- function(name, fit, field_km, spread = 0, preset = 0, n_cov = 0) {
  fits[[name]] <<- fit
  cv <- gcv(fit)
  results[[name]] <<- data.frame(model = name,
    preset_candidates = preset,
    n_param = length(fit$names.fixed) + nrow(fit$summary.hyperpar),
    n_cov = n_cov,
    mlik = mlik(fit), field_km = field_km,
    DIC = fit$dic$dic, WAIC = fit$waic$waic,
    CV_logscore = as.numeric(cv["logs"]), CV_rmse = as.numeric(cv["rmse"]), row.names = NULL)
  cat(sprintf("  %-13s | npar %2d (cov %d) | mlik %.1f (spread %.2f) | field %.1f km | DIC %.0f | WAIC %.0f | logs %.3f | CV %.3f\n",
              name, results[[name]]$n_param, n_cov, mlik(fit), spread, field_km,
              fit$dic$dic, fit$waic$waic, results[[name]]$CV_logscore, results[[name]]$CV_rmse))
}

# fit each model
cat("Fitting (INLA); multi-start on the operator models:\n")

# 1. field-only
stk <- inla.stack(data = list(y = y), A = list(A_field, 1), effects = list(field_idx, cb_all), tag = "s")
fit <- fit_inla(mk_form("", "+ f(field, model = matern_model)"), stk, th = th_mat)
record("field-only", fit, field_km = matern_range(fit))

# 2. LUR
stk <- inla.stack(data = list(y = y), A = list(1), effects = list(cb_all), tag = "s")
fit <- fit_inla(mk_form(lur_rhs), stk)
record("LUR (OLS)", fit, field_km = NA, preset = N_CAND, n_cov = N_LUR)

# 3. LUR+field
stk <- inla.stack(data = list(y = y), A = list(A_field, 1), effects = list(field_idx, cb_all), tag = "s")
fit <- fit_inla(mk_form(lur_rhs, "+ f(field, model = matern_model)"), stk, th = th_mat)
record("LUR+field", fit, field_km = matern_range(fit), preset = N_CAND, n_cov = N_LUR)

# 4. forced
stk <- inla.stack(data = list(y = y), A = list(A_obs, 1),
                  effects = list(list(field = seq_len(mesh$n)), cb_all), tag = "s")
m <- fit_multistart(mk_form("", "+ f(field, model = hybM)"), stk)
record("forced(bg+op)", m$fit, field_km = cgen_range(m$fit, "^Theta[0-9]+ for field$", 2), spread = m$spread, n_cov = 3)

# 5. hybrid
stk <- inla.stack(data = list(y = y), A = list(A_obs, 1),
                  effects = list(list(field = seq_len(mesh$n)), cb_all), tag = "s")
m <- fit_multistart(mk_form(PW, "+ f(field, model = hybM)"), stk)
record("hybrid(bg+op)", m$fit, field_km = cgen_range(m$fit, "^Theta[0-9]+ for field$", 2), spread = m$spread, n_cov = 3)

# 6. two-scale
stk <- inla.stack(data = list(y = y), A = list(A_obs, 1),
                  effects = list(list(field = seq_len(mesh$n)), cb_all), tag = "s")
m <- fit_multistart(mk_form(PW, "+ f(field, model = hyb2)"), stk)
record("two-scale", m$fit, field_km = cgen_range(m$fit, "^Theta[0-9]+ for field$", 2), spread = m$spread, n_cov = 3)

# 7. three-field
stk <- inla.stack(data = list(y = y), A = list(A_obs, A_obs, A_obs, 1),
                  effects = list(list(field_t = seq_len(mesh$n)), list(field_i = seq_len(mesh$n)),
                                 list(field_p = seq_len(mesh$n)), cb_all), tag = "s")
m <- fit_multistart(mk_form(PW,
       "+ f(field_t, model = hyb_t) + f(field_i, model = hyb_i) + f(field_p, model = hyb_p)"), stk)
record("three-field", m$fit, field_km = cgen_range(m$fit, "^Theta2 for field_t$"), spread = m$spread, n_cov = 3)

# 8. split
stk <- inla.stack(data = list(y = y), A = list(A_obs, A_obs, 1),
                  effects = list(list(field_loc = seq_len(mesh$n)), list(field_ind = seq_len(mesh$n)),
                                 cb_all), tag = "s")
m <- fit_multistart(mk_form(PW,
       "+ f(field_loc, model = hyb_loc) + f(field_ind, model = hyb_ind)"), stk)
record("split(ind|loc)", m$fit, field_km = cgen_range(m$fit, "^Theta2 for field_ind$"), spread = m$spread, n_cov = 3)

tab <- do.call(rbind, results)

dir.create("results", showWarnings = FALSE)
save(fits, tab, y, n, loc_km, mesh, Xtraf, Xpt, Xpop, NU, M_LGOCV,
     matern_model, hybM, hyb2, hyb_t, hyb_i, hyb_p, hyb_loc, hyb_ind,
     file = "results/be_no2_fits.RData")
