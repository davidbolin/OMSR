## =====================================================================
## Build the Section 5 covariance comparison 
##
## Run after fit_models.R and make_tables.R.  
## =====================================================================

library(fields)
library(terra)
library(geodata)
library(rSPDE)
library(Matrix)
library(fmesher)


GRID_KM  <- 8       # quadrature grid (km)
deg_km   <- 111.32
regs     <- c(COLORADO = "Colorado MAM", 
              NORWAYANNUAL = "Norway annual", 
              NORWAYDJF = "Norway DJF")

# regional data: stations (y, elev, xy in km) + elevation field on a  regular grid 
prep <- function(tag) {
  if (tag == "COLORADO") {
    data(COmonthlyMet)
    ok <- !is.na(CO.tmean.MAM.climate) & !is.na(CO.elev) & !is.na(CO.loc[,1]) & !is.na(CO.loc[,2])
    y <- as.numeric(CO.tmean.MAM.climate[ok])
    loc <- CO.loc[ok,,drop=FALSE]
    elev <- as.numeric(CO.elev[ok])/1000
    g <- CO.elevGrid
    glon <- rep(g$x, times = length(g$y))
    glat <- rep(g$y, each = length(g$x))
    gz0  <- as.numeric(g$z)/1000
    fin <- is.finite(gz0)
    glon <- glon[fin]
    glat <- glat[fin]
    gz0 <- gz0[fin]
  } else {
    no <- readRDS(file.path("data", "norway_tavg.rds"))
    no <- no[no$lat < 65, ]
    yv <- if (tag == "NORWAYDJF") no$tavg_djf else no$tavg_ann
    ok <- complete.cases(yv, no$lon, no$lat, no$elev)
    y <- yv[ok]
    loc <- cbind(no$lon[ok], no$lat[ok])
    elev <- no$elev[ok]/1000
    rs  <- lapply(c("NOR","SWE","DNK"), function(cc) elevation_30s(country=cc, path="data/dem", mask=FALSE))
    dem <- do.call(terra::mosaic, c(rs, list(fun = "mean")))
    bb  <- terra::ext(min(loc[,1])-1, max(loc[,1])+1, min(loc[,2])-1, max(loc[,2])+1)
    demc <- terra::crop(dem, bb)
    xy0 <- terra::xyFromCell(demc, 1:terra::ncell(demc))
    z0  <- terra::values(demc)[,1]
    fin <- !is.na(z0)
    glon <- xy0[fin,1]
    glat <- xy0[fin,2]
    gz0 <- z0[fin]/1000
    gz0[gz0 < 0] <- 0
  }
  lat0 <- mean(loc[,2])
  lon0 <- mean(loc[,1])
  cs <- cos(lat0*pi/180)
  xy <- cbind((loc[,1]-lon0)*deg_km*cs, (loc[,2]-lat0)*deg_km)      # stations, km 
  fx <- (glon-lon0)*deg_km*cs; fy <- (glat-lat0)*deg_km            # DEM cell centres, km
  ## bin the DEM to a regular GRID_KM grid (area-mean elevation per cell)
  bi <- floor(fx/GRID_KM)
  bj <- floor(fy/GRID_KM)
  key <- paste(bi, bj)
  ux <- tapply((bi+0.5)*GRID_KM, key, `[`, 1)
  uy <- tapply((bj+0.5)*GRID_KM, key, `[`, 1)
  gz <- tapply(gz0, key, mean)
  list(y=y, elev=elev, xy=xy, ux=as.numeric(ux), uy=as.numeric(uy),
       gz=as.numeric(gz), wq=rep(GRID_KM^2, length(gz)), n=length(y))
}

# kernels
green  <- function(r, kappa) besselK(kappa*r, 0)/(2*pi)                 
matcor <- function(D, kappa) { 
  z <- kappa*D 
  R <- z*besselK(z,1)
  diag(R) <- 1
  R 
} 

# ML objective and GLS 
make_obj <- function(P) {
  Dss <- as.matrix(dist(P$xy))
  Dsg <- sqrt(outer(P$xy[,1], P$ux, "-")^2 + outer(P$xy[,2], P$uy, "-")^2)
  Dsg[Dsg < 0.5*GRID_KM] <- 0.5*GRID_KM         # floor K0 singularity
  Selev <- function(kappa) as.numeric(green(Dsg, kappa) %*% (P$wq*P$gz))
  n <- P$n
  obj <- function(par, kap_mu = NULL) {
    kappa <- exp(par[1])
    phi <- exp(par[2])
    if (!is.finite(kappa+phi)) return(1e10)
    z <- Selev(if (is.null(kap_mu)) kappa else kap_mu)
    M <- cbind(1, P$elev, z)
    V <- matcor(Dss, kappa)
    diag(V) <- diag(V) + phi
    L <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(L)) return(1e10)
    Mi <- backsolve(L, forwardsolve(t(L), M))
    yi <- backsolve(L, forwardsolve(t(L), P$y))
    A  <- crossprod(M, Mi)
    b  <- tryCatch(solve(A, crossprod(M, yi)), error = function(e) NULL)
    if (is.null(b)) return(1e10)
    r  <- P$y - M %*% b
    ri <- backsolve(L, forwardsolve(t(L), r))
    rss <- sum(r*ri)
    logdetV <- 2*sum(log(diag(L)))
    val <- n*log(rss/n) + logdetV
    if (is.finite(val)) val else 1e10
  }
  gls <- function(kappa, phi, kap_mu = NULL) {
    z <- Selev(if (is.null(kap_mu)) kappa else kap_mu)
    M <- cbind(1, P$elev, z)
    p <- ncol(M)
    V <- matcor(Dss, kappa); diag(V) <- diag(V) + phi
    L <- chol(V)
    Mi <- backsolve(L, forwardsolve(t(L), M))
    yi <- backsolve(L, forwardsolve(t(L), P$y))
    b  <- solve(crossprod(M, Mi), crossprod(M, yi))
    r <- P$y - M %*% b
    ri <- backsolve(L, forwardsolve(t(L), r))
    list(beta = as.numeric(b), sigU2 = sum(r*ri)/(n-p))
  }
  list(obj = obj, gls = gls, Selev = Selev, Dss = Dss)
}

# scale-free diagnostics g_X and vartheta on the coarse grid 
diagnostics <- function(P, kappa, diag_km = 16) {
  ## re-bin the field to a coarser grid so the dense G x G solves stay feasible
  bi <- floor(P$ux/diag_km); bj <- floor(P$uy/diag_km); key <- paste(bi, bj)
  ux <- as.numeric(tapply((bi+0.5)*diag_km, key, `[`, 1))
  uy <- as.numeric(tapply((bj+0.5)*diag_km, key, `[`, 1))
  gz <- as.numeric(tapply(P$gz, key, mean));  wq <- rep(diag_km^2, length(gz))
  P <- list(ux = ux, uy = uy, gz = gz, wq = wq)
  X <- P$gz - sum(P$wq*P$gz)/sum(P$wq)          # centre the field
  Dgg <- sqrt(outer(P$ux, P$ux, "-")^2 + outer(P$uy, P$uy, "-")^2)
  Dgg[Dgg < 0.5*diag_km] <- 0.5*diag_km
  ## unit-DC-gain smoother S~ = kappa^2 L^{-1} normalised so S~1 = 1. 
  W   <- green(Dgg, kappa)
  SXt <- as.numeric(W %*% (P$wq*X)) / as.numeric(W %*% P$wq)         
  l2  <- function(a, b) sum(P$wq*a*b)
  gX  <- sqrt(l2(SXt,SXt)/l2(X,X))
  Sig <- matcor(Dgg, kappa)                                 
  Lc  <- chol(Sig + diag(1e-8, nrow(Sig)))
  si  <- function(v) backsolve(Lc, forwardsolve(t(Lc), v))
  q12 <- sum(X*si(SXt))
  q11 <- sum(X*si(X))
  q22 <- sum(SXt*si(SXt))
  c(g_X = gX, vartheta = 1 - q12^2/(q11*q22))
}

# 10-fold CV RMSE at the fitted (kappa, phi, sigU2) 
cv_rmse <- function(P, E, kappa, phi, sigU2) {
  set.seed(26)
  n <- P$n
  fold <- sample(rep(1:10, length.out = n))
  z <- E$Selev(kappa)
  M <- cbind(1, P$elev, z)
  S <- matcor(E$Dss, kappa)*sigU2; diag(S) <- diag(S) + phi*sigU2
  pred <- numeric(n)
  for (k in 1:10) { 
    te <- which(fold==k)
    tr <- setdiff(1:n, te)
    Lt <- chol(S[tr,tr])
    bi <- solve(crossprod(M[tr,], backsolve(Lt, forwardsolve(t(Lt), M[tr,]))),
                crossprod(M[tr,], backsolve(Lt, forwardsolve(t(Lt), P$y[tr]))))
    rt <- P$y[tr] - M[tr,] %*% bi
    pred[te] <- M[te,] %*% bi + S[te,tr] %*% backsolve(Lt, forwardsolve(t(Lt), rt))
  }
  sqrt(mean((P$y - pred)^2))
}

# fit one regime by ML, + two-scale ML LRT 
fit_region <- function(tag) {
  P <- prep(tag)
  start <- c(log(1/150), log(0.05))
  ## ML fit with the fixed effects profiled out (GLS); optimise only (kappa, phi)
  E <- make_obj(P)
  o <- optim(start, E$obj, method = "Nelder-Mead", control = list(reltol = 1e-10, maxit = 1200))
  kap <- exp(o$par[1])
  phi <- exp(o$par[2])
  f <- E$gls(kap, phi)
  dg <- diagnostics(P, kap)
  out <- list(ML = data.frame(regime = regs[tag], method = "cov ML (mesh-free)", 
                              n = P$n,
                              kappa_inv_km = 1/kap, range_km = sqrt(8)/kap,
                              beta_pw = f$beta[2], beta_star = f$beta[3]/kap^2, 
                              g_X = dg["g_X"], vartheta = dg["vartheta"], 
                              beta_eff = (f$beta[3]/kap^2)*dg["g_X"],
                              sigma_e = sqrt(phi*f$sigU2), 
                              cv_rmse = cv_rmse(P, E, kap, phi, f$sigU2)))
  
  ## two-scale (free kappa_mu) + matching LRT, also ML (matched fit = o above)
  ots <- optim(c(o$par, log(1/300)), function(pp) E$obj(pp[1:2], kap_mu = exp(pp[3])),
               method = "Nelder-Mead", control = list(reltol = 1e-10, maxit = 2000))
  lr  <- max(0, o$value - ots$value)
  attr(out, "matching") <- data.frame(regime = regs[tag],
                                      kappa_mu_inv_km = 1/exp(ots$par[3]), 
                                      lr_2dLL = lr, p = 1 - pchisq(lr, 1))
  out
}

# run all regimes and assemble the comparison 
cov_rows <- list()
match_rows <- list()
for (tag in names(regs)) {
  cat("fitting", tag, "...\n")
  flush.console()
  fr <- fit_region(tag)
  cov_rows[[tag]] <- fr$ML
  match_rows[[tag]] <- attr(fr, "matching")
}
cov_tab <- do.call(rbind, cov_rows)
rownames(cov_tab) <- NULL
match_tab <- do.call(rbind, match_rows)

## FEM/rSPDE benchmark (Hybrid_a2) from the stored results
fem_rows <- lapply(names(regs), function(tag) {
  h <- read.csv(sprintf("results/hyper_%s.csv", tag))
  rownames(h) <- h$model
  d <- read.csv(sprintf("results/decomp_%s.csv", tag))
  cv <- read.csv(sprintf("results/cv_scores_%s.csv", tag))
  rownames(cv) <- cv$Model
  g <- function(q) d$estimate[d$model=="Hybrid_a2" & d$quantity==q]
  data.frame(regime = regs[tag], method = "FEM (ML)", n = NA,
             kappa_inv_km = h["Hybrid_a2","kappa_inv_km"], 
             range_km = h["Hybrid_a2","practical_range_km"],
             beta_pw = g("beta_pw"), beta_star = g("beta_star"), g_X = g("g_X"),
             vartheta = NA, beta_eff = g("beta_eff"), 
             sigma_e = h["Hybrid_a2","sigma_e"],
             cv_rmse = cv["Hybrid_a2","rmse"])
})

fem_tab <- do.call(rbind, fem_rows)

full <- rbind(cov_tab, fem_tab)
full <- full[order(match(full$regime, regs), full$method), ]
write.csv(full, "results/covariance_compare.csv", row.names = FALSE)
write.csv(match_tab, "results/covariance_matching.csv", row.names = FALSE)

## Table 3 itself is printed at the end, once the per-evaluation timings
## (the eval-ms column) have been benchmarked below.
cat("\n===== operator-matching (two-scale) ML LRT =====\n")
print(transform(match_tab, kappa_mu_inv_km = round(kappa_mu_inv_km,0),
                lr_2dLL = round(lr_2dLL,2), p = round(p,3)), row.names = FALSE)
cat("\nWrote results/covariance_compare.csv, results/covariance_matching.csv\n")


bench_ms<-function(f, reps){ 
  f()
  1000*median(replicate(reps, system.time(f())[["elapsed"]])) 
}

res<-list()
for(tag in names(regs)){
  ## FEM: a single rSPDE negative-log-likelihood evaluation at the MLE
  e<-new.env()
  load(sprintf("results/fits_%s.RData",tag),envir=e)
  fit<-e$fits$Hybrid_a2
  par<-fit$mle_par_orig
  reps<-if(tag=="COLORADO")200 else 400
  t_fem<-bench_ms(function() fit$lik_fun(par), reps)
  ## covariance: a single profiled-ML objective evaluation (n x grid convolution + n x n chol + GLS)
  P<-prep(tag)
  E<-make_obj(P)
  t_cov<-bench_ms(function() E$obj(c(log(1/150),log(0.05))), reps)
  res[[tag]]<-data.frame(regime=tag, n=P$n, mesh_nodes=nrow(e$p$mesh$loc), 
                         grid_cells=length(P$gz),
                         fem_eval_ms=round(t_fem), cov_eval_ms=round(t_cov))
  cat(sprintf("%-14s n=%d  FEM %.0f ms  cov %.0f ms  (FEM/cov = %.2f)\n",regs[tag],
              P$n,t_fem,t_cov,t_fem/t_cov))
  flush.console()
}
tab<-do.call(rbind,res)
rownames(tab)<-NULL
write.csv(tab, "results/covariance_timing.csv", row.names=FALSE)

## ---- Table 3: hybrid (alpha = 2) fitted through FEM vs covariance --------
## Manuscript Table 3.  Merge the per-evaluation timings (eval_ms) onto the
## comparison table: FEM rows take fem_eval_ms, covariance rows cov_eval_ms.
## beta_eff = beta_star_fc * g_X is the "beta*_fc g_X" column; cv_rmse is the
## 10-fold CV RMSE; eval_ms is the median wall-clock of one log-likelihood eval.
full$tag     <- names(regs)[match(full$regime, regs)]
full$eval_ms <- ifelse(grepl("^FEM", full$method),
                       tab$fem_eval_ms[match(full$tag, tab$regime)],
                       tab$cov_eval_ms[match(full$tag, tab$regime)])
full$method  <- ifelse(grepl("^FEM", full$method), "FEM", "covariance")
table3 <- transform(full[, c("regime","method","kappa_inv_km","beta_pw",
                             "beta_eff","sigma_e","cv_rmse","eval_ms")],
                    kappa_inv_km = round(kappa_inv_km, 0),
                    beta_pw      = round(beta_pw, 2),
                    beta_eff     = round(beta_eff, 2),
                    sigma_e      = round(sigma_e, 2),
                    cv_rmse      = round(cv_rmse, 3),
                    eval_ms      = round(eval_ms, 0))
rownames(table3) <- NULL
cat("\n=== Table 3: hybrid (alpha=2) via FEM vs covariance ===\n")
cat("(beta_eff = beta_star_fc * g_X; cv_rmse = 10-fold CV; eval_ms = median wall-clock of one log-likelihood evaluation)\n")
print(table3, row.names = FALSE)
