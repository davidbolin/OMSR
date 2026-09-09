# --------------------------------------------------------------------
# Fit the models and run the 10-fold CV.
# --------------------------------------------------------------------

library(Matrix)
library(fmesher)
library(rSPDE)
library(fields)
library(geodata)
library(terra)
library(sf)
library(rnaturalearth)

SEED_CV  <- 26

source("../common/utils.R")

# warm-start options: reuse the a=2 fit's parameters as starting values
# but estimate alpha (convert the inherited fix_alpha into start_alpha).
warm_opts <- function(a2fit) {
  mo <- rSPDE:::extract_starting_values(a2fit)
  mo$start_alpha <- as.numeric(mo$fix_alpha)
  mo$fix_alpha <- NULL
  mo
}

# distance to coast (km) for lon/lat points, from Natural Earth coastline
dist_to_coast_km <- function(lon, lat) {
  sf::sf_use_s2(TRUE)
  coast <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf")
  pts   <- sf::st_as_sf(data.frame(lon = lon, lat = lat),
                        coords = c("lon", "lat"), crs = 4326)
  d <- sf::st_distance(pts, sf::st_union(sf::st_geometry(coast)))
  as.numeric(d) / 1000
}

# Data preparation per regime
prep_data <- function(REGION, SEASON = NA) {
  deg_to_km <- 111.32
  if (REGION == "COLORADO") {
    data(COmonthlyMet)
    # NOTE: in the fields package CO.tmean.MAM.climate is mislabelled.  Use the
    # manually computed mean (min + max)/2 as the response.
    y_all <- (CO.tmin.MAM.climate + CO.tmax.MAM.climate) / 2
    loc_all <- CO.loc
    elev_all <- CO.elev
    ok <- !is.na(y_all) & !is.na(elev_all) & !is.na(loc_all[,1]) & !is.na(loc_all[,2])
    y <- as.numeric(y_all[ok])
    loc <- loc_all[ok,,drop=FALSE]
    elev <- as.numeric(elev_all[ok])
    X_obs <- elev / 1000
    lat0 <- mean(loc[,2])
    loc_km <- cbind(x = (loc[,1]-mean(loc[,1]))*deg_to_km*cos(lat0*pi/180),
                    y = (loc[,2]-lat0)*deg_to_km)
    bnd  <- fm_nonconvex_hull(loc_km, convex = 60, concave = 80)
    mesh <- fm_mesh_2d(loc = loc_km, boundary = bnd, max.edge = c(55,160),
                       cutoff = 25, offset = c(40,180))
    mesh_lonlat <- cbind(lon = mesh$loc[,1]/(deg_to_km*cos(lat0*pi/180)) + mean(loc_all[,1]),
                         lat = mesh$loc[,2]/deg_to_km + lat0)
    elev_nodes <- interp.surface(CO.elevGrid, mesh_lonlat)
    if (any(is.na(elev_nodes))) for (i in which(is.na(elev_nodes))) {
      ix <- which.min(abs(CO.elevGrid$x - mesh_lonlat[i,1]))
      iy <- which.min(abs(CO.elevGrid$y - mesh_lonlat[i,2]))
      elev_nodes[i] <- CO.elevGrid$z[ix, iy]
    }
    X_nodes <- elev_nodes / 1000
    data <- data.frame(y=y, x1=loc_km[,1], x2=loc_km[,2], elev=elev/1000)
    has_coast <- FALSE
  } else { # NORWAY
    no <- readRDS(file.path("data", "norway_tavg.rds"))
    no <- no[no$lat < 65, ]
    y <- if (SEASON == "DJF") no$tavg_djf else no$tavg_ann
    loc <- cbind(no$lon, no$lat)
    elev <- no$elev / 1000
    ok <- complete.cases(y, loc, elev)
    y <- y[ok]
    loc <- loc[ok,]
    elev <- elev[ok]
    lat0 <- mean(loc[,2])
    lon0 <- mean(loc[,1])
    loc_km <- cbind(x=(loc[,1]-lon0)*deg_to_km*cos(lat0*pi/180), y=(loc[,2]-lat0)*deg_to_km)
    X_obs <- elev
    bnd  <- fm_nonconvex_hull(loc_km, convex = 80, concave = 120)
    mesh <- fm_mesh_2d(loc = loc_km, boundary = bnd, max.edge = c(40,120),
                       cutoff = 10, offset = c(50,180))
    # DEM at mesh nodes L2-projected onto the FEM basis: each node value is
    # the area-average of the 30 arc-sec DEM over its basis-function support
    node_lon <- mesh$loc[,1]/(deg_to_km*cos(lat0*pi/180)) + lon0
    node_lat <- mesh$loc[,2]/deg_to_km + lat0
    rs  <- lapply(c("NOR","SWE","DNK"), function(cc)
             elevation_30s(country=cc, path="data", mask=FALSE))
    dem <- do.call(terra::mosaic, c(rs, list(fun="mean")))
    bb  <- terra::ext(min(node_lon)-0.3, max(node_lon)+0.3,
                      min(node_lat)-0.3, max(node_lat)+0.3)
    demc <- terra::crop(dem, bb)
    xy <- terra::xyFromCell(demc, 1:terra::ncell(demc))
    zz <- terra::values(demc)[,1]
    okz <- !is.na(zz)
    xy <- xy[okz,]
    zz <- zz[okz]
    zz[zz < 0] <- 0
    fx <- (xy[,1]-lon0)*deg_to_km*cos(lat0*pi/180)
    fy <- (xy[,2]-lat0)*deg_to_km
    wq <- cos(xy[,2]*pi/180)                # cell-area weight
    Afine <- fm_basis(mesh, loc = cbind(fx, fy))
    ins <- Matrix::rowSums(Afine) > 0.999
    Afine <- Afine[ins,]
    zf <- zz[ins]
    wf <- wq[ins]
    Xnode_m <- as.numeric(Matrix::crossprod(Afine, wf*zf)) / as.numeric(Matrix::crossprod(Afine, wf))
    Xnode_m[is.na(Xnode_m)] <- 0
    X_nodes <- Xnode_m / 1000
    coast   <- dist_to_coast_km(loc[,1], loc[,2])
    data <- data.frame(y=y, x1=loc_km[,1], x2=loc_km[,2], elev=elev,
                       coast = as.numeric(scale(coast)))
    has_coast <- TRUE
  }
  list(REGION=REGION, SEASON=SEASON, data=data, mesh=mesh, X_nodes=X_nodes,
       X_obs=X_obs, loc_km=loc_km, has_coast=has_coast, n=length(y))
}


# Main code to run one case
run_region <- function(REGION, SEASON = NA) {
  cat("Running ", REGION, SEASON, "\n")
  tag <- if (REGION=="NORWAY") paste0("NORWAY", SEASON) else REGION

  #Prepare data for the region
  p <- prep_data(REGION, SEASON)

  #Fit all models
  data <- p$data
  mesh <- p$mesh
  model1 <- spde.matern.operators(mesh = mesh, alpha = 2)
  model2 <- hybrid.spde(mesh = mesh, alpha = 2, X = p$X_nodes)

  data$Selev <- gaussian_kernel_smooth(p$loc_km, p$X_obs, p$loc_km, 20)
  tps_X <- mgcv::gam(elev ~ s(x1, x2, k = 30), data = data, method = "REML")
  data$Relev <- data$elev - as.numeric(predict(tps_X, data))

  cat("Fit alpha=2 models\n")
  res0    <- rspde_lme(y ~ elev, data = data)
  add_a2  <- rspde_lme(y ~ elev, data=data, model=model1, loc=c("x1","x2"), model_options=list(fix_alpha = 2), optim_method="BFGS")
  for_a2  <- rspde_lme(y ~ 1,    data=data, model=model2, loc=c("x1","x2"), model_options=list(fix_alpha = 2))
  hyb_a2  <- rspde_lme(y ~ elev, data=data, model=model2, loc=c("x1","x2"), model_options=list(fix_alpha = 2), optim_method="BFGS")
  bw_a2   <- rspde_lme(y ~ Selev,data=data, model=model1, loc=c("x1","x2"), model_options=list(fix_alpha = 2))
  sp_a2   <- rspde_lme(y ~ Relev,data=data, model=model1, loc=c("x1","x2"), model_options=list(fix_alpha = 2))

  cat("Fit general alpha models\n")
  add_w <- rspde_lme(y ~ elev, data=data, model=model1, loc=c("x1","x2"), model_options=warm_opts(add_a2), optim_method="BFGS")
  for_w <- rspde_lme(y ~ 1,    data=data, model=model2, loc=c("x1","x2"), model_options=warm_opts(for_a2))
  hyb_w <- rspde_lme(y ~ elev, data=data, model=model2, loc=c("x1","x2"), model_options=warm_opts(hyb_a2), optim_method="BFGS")
  bw_w  <- rspde_lme(y ~ Selev,data=data, model=model1, loc=c("x1","x2"), model_options=warm_opts(bw_a2))

  fits <- list(OLS=res0, Additive_est=add_w, Additive_a2=add_a2,
               Forced_est=for_w, Forced_a2=for_a2,
               Hybrid_est=hyb_w, Hybrid_a2=hyb_a2,
               BW_est=bw_w, BW_a2=bw_a2, SpatPlus_a2=sp_a2)

  cat("Fit two-scale models\n")
  kap0 <- as.numeric(hyb_a2$coeff$random_effects[["kappa"]])
  model_ts <- hybrid.spde(mesh = mesh, alpha = 2, X = p$X_nodes, kappa_mu = kap0)
  fits$TwoScale_a2 <- rspde_lme(y ~ elev, data=data, model=model_ts, loc=c("x1","x2"),
                                model_options = rSPDE:::extract_starting_values(hyb_a2),
                                optim_method="BFGS")

  addcoast <- NULL
  if (p$has_coast) {
    cat("Fit coastal models\n")
    fits$Coast_a2 <- rspde_lme(y ~ elev + coast, data=data, model=model2, loc=c("x1","x2"),
                               model_options=list(fix_alpha = 2), optim_method="BFGS")
    fits$Coast_est  <- rspde_lme(y ~ elev + coast, data=data, model=model2, loc=c("x1","x2"),
                          model_options=warm_opts(fits$Coast_a2), optim_method="BFGS")
    # additive field + coast: a figure-only comparator scored by the
    # leave-group-out CV in make_lgocv.R. Fitted here and saved separately so
    # that script need not refit it; kept out of `fits` so it does not enter
    # the k-fold CV or the tables.
    addcoast <- rspde_lme(y ~ elev + coast, data=data, model=model1, loc=c("x1","x2"),
                          model_options=list(fix_alpha = 2), optim_method="BFGS")
  }
  attr(fits, "data") <- data

  # Cross-validation
  cat("Run crossvalidation\n")
  pr <- rSPDE::posterior_crossvalidation(object = fits, true_CV = TRUE,
                                         seed = SEED_CV, parallel_folds = TRUE,
                                         print = FALSE)

  save(p, fits, pr, data, addcoast, file = sprintf("results/fits_%s.RData", tag))

}

run_region("COLORADO", NA)
run_region("NORWAY", "ANNUAL")
run_region("NORWAY", "DJF")
