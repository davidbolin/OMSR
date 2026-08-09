## =====================================================================
## Fit the models for Section 6 
## =====================================================================

library(Matrix)
library(fmesher)
library(sf)
library(INLA)
library(rSPDE)
library(terra)

set.seed(1)
ALPHA <- 2          # fixed SPDE smoothness (alpha=2 throughout)
NU <- ALPHA - 1     # Matern smoothness nu
NSTART  <- 2      # multi-start restarts for the operator models
M_LGOCV <- -1L    # inla.group.cv level sets: -1 = leave-one-out

cn <- read.csv("data/curieuzeneuzen_2018.csv")
cn <- cn[is.finite(cn$lon)&is.finite(cn$lat)&is.finite(cn$no2), ]
n <- nrow(cn)
y <- cn$no2
lon0 <- mean(cn$lon)
lat0 <- mean(cn$lat)
deg_km <- 111.32
## project lon/lat onto a local equirectangular grid in km, centred on the data
to_km <- function(lo,la) cbind((lo-lon0)*deg_km*cos(lat0*pi/180), (la-lat0)*deg_km)
loc_km <- to_km(cn$lon, cn$lat)
pts <- read.csv("data/be_point_sources_2018.csv")
pts_km <- to_km(pts$lon, pts$lat)

dom <- rbind(loc_km, pts_km)
mesh <- fm_mesh_2d(loc.domain=dom, boundary=fm_nonconvex_hull(dom,convex=10,concave=15),
                   max.edge=c(2,15), cutoff=1.0, offset=c(5,30))
cat(sprintf("Read %d NO2 tubes and %d facilities; mesh has %d nodes\n", n, nrow(pts), mesh$n))
nodes <- mesh$loc[,1:2]
ns2 <- rowSums(nodes^2)          # squared node norms, reused in the distance formula
chunk <- 512L

## index of the nearest mesh node for each row of P, computed in chunks so the
## full points-by-nodes distance matrix is never formed at once
assign_node <- function(P){
  out<-integer(nrow(P))
  for(i in seq.int(1L,nrow(P),by=chunk)){
    j<-min(i+chunk-1L,nrow(P))
    qi<-P[i:j,,drop=FALSE]
    D2<-outer(rowSums(qi^2),ns2,`+`)-2*tcrossprod(qi,nodes)
    out[i:j]<-max.col(-D2,ties.method="first")
  }
  out
}

node_of_obs <- assign_node(loc_km)          # nearest node for each observation station
fem <- fm_fem(mesh)
area <- as.numeric(Matrix::diag(fem$c0))    # area associated with each mesh node
## accumulate the weights `w` of points P onto their nearest node and divide by
## the nodal area to get a density (weight per km^2); `sp` is an overall scale
nodal_density <- function(P,w,sp=1){
  no<-assign_node(P)
  s<-numeric(mesh$n)
  ag<-tapply(w,no,sum)
  s[as.integer(names(ag))]<-ag
  s*sp/area
}
## the three forcing fields as nodal densities: road-traffic intensity (OSM
## roads), industrial NOx (E-PRTR point sources), and residential population
cat("Building forcing fields (traffic, industry, population)...\n")
roads <- st_segmentize(st_transform(st_zm(st_read("data/be_traffic_roads.gpkg",quiet=TRUE)),31370),
                       dfMaxLength=200)
rp <- st_cast(roads,"POINT",warn=FALSE)
rxy <- st_coordinates(st_transform(rp,4326))
tk <- to_km(rxy[,1],rxy[,2])
kr <- tk[,1]>=min(dom[,1])-10 & tk[,1]<=max(dom[,1])+10 & tk[,2]>=min(dom[,2])-10 & tk[,2]<=max(dom[,2])+10
traf_dens <- nodal_density(tk[kr,,drop=FALSE], rp$intensity[kr], sp=0.2)
pt_dens   <- nodal_density(pts_km, pts$nox_kg, sp=1)
pv <- rep(0, mesh$n)

popr <- terra::rast("data/be_population.tif")
nlon <- lon0 + nodes[,1]/(deg_km*cos(lat0*pi/180))
nlat <- lat0 + nodes[,2]/deg_km
pv <- terra::extract(popr, cbind(nlon,nlat))[,1]
pv[!is.finite(pv)] <- 0

## standardised nodal forcing covariates -- the three X channels the operator models use
Xtraf <- as.numeric(scale(traf_dens))
Xpt <- as.numeric(scale(pt_dens))
Xpop <- as.numeric(scale(pv))

## own_ind / own_pop are the pointwise source covariates (own_traffic enters
## through the `cand` column of the same name, so no separate variable is needed)
own_ind     <- as.numeric(scale(pt_dens[node_of_obs]))
own_pop     <- as.numeric(scale(pv[node_of_obs]))


## ---- candidate LUR predictors + supervised forward selection -------------
## The land-use regression (LUR) reference models use classic buffer covariates:
## for each station, the total traffic / industrial emission within several radii,
## plus distance-to- and count-of nearby industry and some land-use fractions.
## A supervised forward selection then keeps the subset that best explains NO2.
cat("Building LUR buffer predictors...\n")
Dp <- sqrt(outer(loc_km[,1],pts_km[,1],"-")^2 + outer(loc_km[,2],pts_km[,2],"-")^2)  # station-to-facility distances (km)
rt <- c(0.1,0.3,0.5,1,2)     # traffic buffer radii (km)
ri <- c(1,5,10,25)           # industry buffer radii (km)
wt <- traf_dens*area         # nodal traffic mass (density x nodal area)
wi <- pt_dens*area           # nodal industrial mass
## buffer sums: for each station add up the nodal mass within each radius,
## looping over blocks of stations to keep the station-by-node matrix small
B <- matrix(0, n, length(rt)+length(ri))
for (s in seq.int(1L,n,by=1000L)) {
  e<-min(s+999L,n)
  D2 <- outer(loc_km[s:e,1],nodes[,1],"-")^2 + outer(loc_km[s:e,2],nodes[,2],"-")^2
  for (j in seq_along(rt)) B[s:e,j] <- (D2<=rt[j]^2)%*%wt
  for (j in seq_along(ri)) B[s:e,length(rt)+j] <- (D2<=ri[j]^2)%*%wi
}

## candidate matrix: traffic buffers + pointwise traffic, industry buffers,
## distance to nearest industry, and count of industries within 10 km
cand <- cbind(traf_0.1=B[,1],traf_0.3=B[,2],traf_0.5=B[,3],traf_1=B[,4],traf_2=B[,5],
              own_traffic=traf_dens[node_of_obs],
              ind_1=B[,6],ind_5=B[,7],ind_10=B[,8],ind_25=B[,9],
              dist_ind=apply(Dp,1,min), cnt_ind10=rowSums(Dp<=10))

## physically expected sign of each covariate's effect on NO2 (more emission ->
## more NO2; distance to industry is the exception and should be negative)
signs <- c(traf_0.1=1,traf_0.3=1,traf_0.5=1,traf_1=1,traf_2=1,own_traffic=1,
           ind_1=1,ind_5=1,ind_10=1,ind_25=1,dist_ind=-1,cnt_ind10=1)

## append the extra land-use / coast covariates (urban & natural fractions at
## several radii, distance to coast), matched to the stations by `code`
ex <- read.csv("data/be_extra_covariates.csv")
ex<-ex[match(cn$code, ex$code),]
for (v in setdiff(names(ex),"code")) {
  cand<-cbind(cand, ex[[v]])
  colnames(cand)[ncol(cand)]<-v
  signs[v]<- if(grepl("natural",v)) -1 else 1   # natural land use is expected to lower NO2
}

cand <- scale(cand)          # standardise all candidates
cand[is.na(cand)] <- 0

## Supervised forward selection: greedily add the candidate that most improves
## the adjusted R^2 of an OLS regression of y on the covariates chosen so far --
## but only if its coefficient has the physically expected sign and the gain
## exceeds 0.01.  Stop when nothing qualifies; return the chosen covariate names.
fwd_select <- function(cols){
  sel<-character(0)                          # covariates selected so far
  aR2<- -Inf                                 # their adjusted R^2
  rem<-cols                                  # remaining candidates
  repeat {
    bv<-NA                                   # best covariate to add this round
    ba<- if(length(sel)) aR2+0.01 else -Inf  # adj-R^2 it must beat to be added
    for (v in rem){
      f<-lm(y~., data=as.data.frame(cbind(y=y, cand[,c(sel,v),drop=FALSE])))
      cf<-coef(f)[v]
      a<-summary(f)$adj.r.squared
      if(!is.na(cf)&&sign(cf)==signs[[v]]&&a>ba){   # right sign and best improvement
        ba<-a
        bv<-v
      }
    }
    if(is.na(bv)) break                      # no admissible improvement -> done
    sel<-c(sel,bv)
    rem<-setdiff(rem,bv)
    aR2<-ba
  }
  sel
}

bg_cols  <- grep("urban_|natural_|dist_coast", colnames(cand), value=TRUE)  # land-use "background" covariates

cat("Running LUR forward selection...\n")
sel_full <- fwd_select(colnames(cand))       # covariates chosen for the LUR models
lur_rhs  <- paste(sel_full, collapse=" + ")
N_CAND   <- ncol(cand)                        # pre-set candidates offered to the LUR selection
N_LUR    <- length(sel_full)                  # of which this many were selected (LUR / LUR+field)
cat(sprintf("LUR predictors: %s | operator background: none\n",
            paste(sel_full,collapse=", ")))

## INLA building blocks 
PREC_BETA <- 1e-5
A_obs <- fm_basis(mesh, loc=loc_km)
matern_model <- rspde.matern(mesh=mesh, nu=NU, parameterization="matern")
A_field <- rspde.make.A(mesh=mesh, loc=loc_km, nu=NU)
field_idx <- rspde.make.index(name="field", mesh=mesh, nu=NU)
## covariate columns available to every formula (pointwise + selected LUR buffers)
cdf <- as.data.frame(cand[, unique(c(sel_full, bg_cols, "own_traffic")), drop=FALSE])
cb_all <- c(list(Intercept=1, own_ind=own_ind, own_pop=own_pop), as.list(cdf))

## prior calibration from a field-only fit
cat("Calibrating field priors (field-only INLA fit)...\n")
stk0 <- inla.stack(data=list(y=y), A=list(A_field,1), effects=list(field_idx, cb_all), tag="s")
f0 <- inla(y ~ -1 + Intercept + f(field, model=matern_model), 
           family="gaussian", data=inla.stack.data(stk0),
           control.predictor=list(A=inla.stack.A(stk0)), 
           control.fixed=list(prec=PREC_BETA,prec.intercept=0),
           control.inla=list(int.strategy="eb"))

s0 <- summary(rspde.result(f0,"field",matern_model))

## posterior mean of the first summary row whose name matches pattern `p`
g1<-function(s,p){
  i<-grep(p,rownames(s),ignore.case=TRUE)
  if(length(i))s[i[1],"mean"] else NA_real_
}
range0 <- g1(s0,"range")           # calibrated field range and marginal sd, used
sd0 <- g1(s0,"std|sigma|dev")      # to set the priors/starts for the operator fields
if(is.na(range0)) range0<-6
if(is.na(sd0)) sd0<-sd(y)
kappa0 <- sqrt(8*NU)/range0        # SPDE kappa (inverse range)
tau0 <- 1/(sd0*sqrt(4*pi)*kappa0)  # SPDE tau (precision scale)
th_mat <- f0$mode$theta             # warm start for the plain-field (matern) models
Xall <- cbind(Xtraf, Xpt, Xpop)     # three forcing channels

cat("Building operator-field model objects...\n")
## hybM: one operator-forced field driven by all three sources (kappa_mu = kappa,
## so the forcing shares the field's range) -- for the forced & hybrid models
hybM <- rspde.hybrid.matern(mesh=mesh, X=Xall,
                            prior.tau=list(mean=log(tau0),prec=5), 
                            prior.kappa=list(mean=log(kappa0),prec=5),
                            prior.beta_x=list(mean=rep(0,3),prec=rep(0.001,3)),
                            start.ltau=log(tau0), start.lkappa=log(kappa0), 
                            start.beta_x=rep(0,3))
km0 <- sqrt(8*NU)/15                # starting kappa_mu (forcing range ~15 km)

## hyb2: as hybM but with a free forcing range kappa_mu != kappa -- the two-scale
## model, used to test the operator-matching assumption
hyb2 <- rspde.hybrid.matern(mesh=mesh, X=Xall, separate_kappa_mu=TRUE,
    prior.tau=list(mean=log(tau0),prec=5), prior.kappa=list(mean=log(kappa0),prec=5),
    prior.kappa_mu=list(mean=log(km0),prec=2), prior.beta_x=list(mean=rep(0,3),prec=rep(0.001,3)),
    start.ltau=log(tau0), start.lkappa=log(kappa0), start.lkappa_mu=log(km0), start.beta_x=rep(0,3))

# build a single-source operator field for covariate Xn with prior range r0
mkh <- function(Xn, r0){
  k0<-sqrt(8*NU)/r0
  t0<-1/(sd0*sqrt(4*pi)*k0)
  rspde.hybrid.matern(mesh=mesh, X=as.matrix(Xn),
                      prior.tau=list(mean=log(t0),prec=5),
                      prior.kappa=list(mean=log(k0),prec=3), 
                      prior.beta_x=list(mean=0,prec=0.001),
                      start.ltau=log(t0), start.lkappa=log(k0), 
                      start.beta_x=0) 
}
hyb_t <- mkh(Xtraf,5)               # per-source fields for the three-field model
hyb_i <- mkh(Xpt,15)
hyb_p <- mkh(Xpop,8)

## mkh_multi: like mkh but for a multi-column X (several sources sharing one field)
mkh_multi <- function(Xn, r0){
  Xn<-as.matrix(Xn)
  p<-ncol(Xn)
  k0<-sqrt(8*NU)/r0
  t0<-1/(sd0*sqrt(4*pi)*k0)
  rspde.hybrid.matern(mesh=mesh, X=Xn, prior.tau=list(mean=log(t0),prec=5), 
                      prior.kappa=list(mean=log(k0),prec=3),
                      prior.beta_x=list(mean=rep(0,p),prec=rep(0.001,p)), 
                      start.ltau=log(t0), start.lkappa=log(k0), 
                      start.beta_x=rep(0,p)) 
}
hyb_loc <- mkh_multi(cbind(Xtraf,Xpop), 2)
hyb_ind <- mkh_multi(Xpt, 20)

## ---- fitting helpers -----------------------------------------------------
PW <- "own_traffic + own_ind + own_pop"    # the three pointwise source covariates

## one INLA fit with the shared control options; `th` is an optional warm start
fit_inla <- function(form, stk, th = NULL)
  inla(form, family = "gaussian", data = inla.stack.data(stk),
       control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
       control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
       control.mode      = if (is.null(th)) NULL else list(theta = th, restart = TRUE),
       control.fixed     = list(prec = PREC_BETA, prec.intercept = 0),
       control.inla      = list(int.strategy = "eb"))

mlik <- function(f) as.numeric(f$mlik)[1]   # marginal log-likelihood of an INLA fit

## operator models are multimodal: fit once, then NSTART perturbed restarts,
## keeping the best marginal likelihood.  Returns the fit and the mlik spread.
fit_multistart <- function(form, stk) {
  best <- fit_inla(form, stk)
  ml <- mlik(best)
  th0 <- best$mode$theta
  if (NSTART > 0 && length(th0)) for (s in 1:NSTART) {
    set.seed(7 * s)
    fs <- tryCatch(fit_inla(form, stk, th = th0 + rnorm(length(th0), sd = 1.5)),
                   error = function(e) NULL)
    if (!is.null(fs)) { ml <- c(ml, mlik(fs)); if (mlik(fs) > mlik(best) + 1e-3) best <- fs }
  }
  list(fit = best, spread = max(ml) - min(ml))
}

## build the fixed-effects formula: Intercept [+ covariates] [+ field terms]
mk_form <- function(cov, field = "") {
  rhs <- if (nzchar(cov) && cov != "1") paste("Intercept +", cov) else "Intercept"
  as.formula(paste("y ~ -1 +", rhs, field))
}

## leave-one-out group-CV RMSE and log-score ($cv = predictive density;
## LS = -mean(log(cv)), lower better, as in Adin et al. 2024)
gcv <- function(f) {
  g <- inla.group.cv(f, num.level.sets = M_LGOCV)
  if (is.null(g$cv)) return(c(rmse = NA, logs = NA))
  c(rmse = sqrt(mean((y - g$mean)^2, na.rm = TRUE)), logs = -mean(log(g$cv), na.rm = TRUE))
}
## practical range kappa^-1 (km) of a Matern field, or of a cgeneric field from
## its log-kappa hyperparameter row
matern_range <- function(f) { 
  rr <- rspde.result(f, "field", matern_model)
  if (!is.null(rr$summary.range)) rr$summary.range$mean else NA 
}

cgen_range <- function(f, pattern, which = 1) {
  i <- grep(pattern, rownames(f$summary.hyperpar))
  if (length(i) >= which) sqrt(8 * NU) / exp(f$summary.hyperpar[i[which], "mean"]) else NA 
}

## record a fitted model: store it and add its comparison-table row.
## `preset` = pre-set candidates offered to selection (LUR models only);
## `n_cov`  = covariates estimated from data (the "(cov.)" count in Table 4).
fits <- list()
results <- list()
record <- function(name, fit, field_km, spread = 0, preset = 0, n_cov = 0) {
  fits[[name]] <<- fit
  cv <- gcv(fit)
  results[[name]] <<- data.frame(model = name,
    preset_candidates = preset,
    n_param = length(fit$names.fixed) + nrow(fit$summary.hyperpar),  # fixed coefs + hyperpar (incl. nugget)
    n_cov = n_cov,
    mlik = mlik(fit), field_km = field_km,
    DIC = fit$dic$dic, WAIC = fit$waic$waic,
    CV_logscore = as.numeric(cv["logs"]), CV_rmse = as.numeric(cv["rmse"]), row.names = NULL)
  cat(sprintf("  %-13s | npar %2d (cov %d) | mlik %.1f (spread %.2f) | field %.1f km | DIC %.0f | WAIC %.0f | logs %.3f | CV %.3f\n",
              name, results[[name]]$n_param, n_cov, mlik(fit), spread, field_km,
              fit$dic$dic, fit$waic$waic, results[[name]]$CV_logscore, results[[name]]$CV_rmse))
}

## fit each model
cat("Fitting (INLA); multi-start on the operator models:\n")

## 1. field-only: a centred Whittle-Matern field, no covariates
stk <- inla.stack(data = list(y = y), A = list(A_field, 1), effects = list(field_idx, cb_all), tag = "s")
fit <- fit_inla(mk_form("", "+ f(field, model = matern_model)"), stk, th = th_mat)
record("field-only", fit, field_km = matern_range(fit))

## 2. LUR (OLS): forward-selected buffer/land-use covariates, no field
stk <- inla.stack(data = list(y = y), A = list(1), effects = list(cb_all), tag = "s")
fit <- fit_inla(mk_form(lur_rhs), stk)
record("LUR (OLS)", fit, field_km = NA, preset = N_CAND, n_cov = N_LUR)

## 3. LUR+field: the LUR covariates plus a Matern field
stk <- inla.stack(data = list(y = y), A = list(A_field, 1), effects = list(field_idx, cb_all), tag = "s")
fit <- fit_inla(mk_form(lur_rhs, "+ f(field, model = matern_model)"), stk, th = th_mat)
record("LUR+field", fit, field_km = matern_range(fit), preset = N_CAND, n_cov = N_LUR)

## 4. forced(bg+op): the forcing enters only through the operator field (no pointwise terms)
stk <- inla.stack(data = list(y = y), A = list(A_obs, 1),
                  effects = list(list(field = seq_len(mesh$n)), cb_all), tag = "s")
m <- fit_multistart(mk_form("", "+ f(field, model = hybM)"), stk)
record("forced(bg+op)", m$fit, field_km = cgen_range(m$fit, "^Theta[0-9]+ for field$", 2), spread = m$spread, n_cov = 3)

## 5. hybrid(bg+op): pointwise sources + the operator-forced field
stk <- inla.stack(data = list(y = y), A = list(A_obs, 1),
                  effects = list(list(field = seq_len(mesh$n)), cb_all), tag = "s")
m <- fit_multistart(mk_form(PW, "+ f(field, model = hybM)"), stk)
record("hybrid(bg+op)", m$fit, field_km = cgen_range(m$fit, "^Theta[0-9]+ for field$", 2), spread = m$spread, n_cov = 3)

## 6. two-scale: hybrid with a free forcing range (kappa_mu != kappa), to test operator matching
stk <- inla.stack(data = list(y = y), A = list(A_obs, 1),
                  effects = list(list(field = seq_len(mesh$n)), cb_all), tag = "s")
m <- fit_multistart(mk_form(PW, "+ f(field, model = hyb2)"), stk)
record("two-scale", m$fit, field_km = cgen_range(m$fit, "^Theta[0-9]+ for field$", 2), spread = m$spread, n_cov = 3)

## 7. three-field: a separate operator field per source (traffic / industry / population)
stk <- inla.stack(data = list(y = y), A = list(A_obs, A_obs, A_obs, 1),
                  effects = list(list(field_t = seq_len(mesh$n)), list(field_i = seq_len(mesh$n)),
                                 list(field_p = seq_len(mesh$n)), cb_all), tag = "s")
m <- fit_multistart(mk_form(PW,
       "+ f(field_t, model = hyb_t) + f(field_i, model = hyb_i) + f(field_p, model = hyb_p)"), stk)
record("three-field", m$fit, field_km = cgen_range(m$fit, "^Theta2 for field_t$"), spread = m$spread, n_cov = 3)

## 8. split(ind|loc): traffic+population share one local field; industry gets its own
stk <- inla.stack(data = list(y = y), A = list(A_obs, A_obs, 1),
                  effects = list(list(field_loc = seq_len(mesh$n)), list(field_ind = seq_len(mesh$n)),
                                 cb_all), tag = "s")
m <- fit_multistart(mk_form(PW,
       "+ f(field_loc, model = hyb_loc) + f(field_ind, model = hyb_ind)"), stk)
record("split(ind|loc)", m$fit, field_km = cgen_range(m$fit, "^Theta2 for field_ind$"), spread = m$spread, n_cov = 3)

tab <- do.call(rbind, results)

## save the fitted models plus the objects the downstream scripts need
save(fits, tab, y, n, loc_km, mesh, Xtraf, Xpt, Xpop, NU, M_LGOCV,
     file = "results/be_no2_fits.RData")