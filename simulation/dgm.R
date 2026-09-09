# ---------------------------------------------------------------------
# Data-generating mechanism for the simulation study
#
# Four scenarios
# (S1) Forced-SPDE truth    Y = beta_fc S X + U           (pw,fc)=(0,1)
# (S2) Rough X, additive    Y = beta_pw X   + U, X rough  (pw,fc)=(1,0)
# (S3) Smooth X, additive   Y = beta_pw X   + U, X smooth (pw,fc)=(1,0)
# (S4) Mismatched operator  Y = beta_fc S0 X + U,         (pw,fc)=(0,1)
#      with S0 a Gaussian-kernel smoother.
# ---------------------------------------------------------------------

library(Matrix)
library(fmesher)
library(rSPDE)

source("../common/utils.R")


# build_world has everything that does not change across replicates
build_world <- function(L         = 10,        # domain side length
                        range_u   = 2.5,       # practical range of U
                        var_u     = 1,         # marginal variance of U
                        nu_u      = 1,         # alpha = nu + d/2 = 2
                        max.edge  = c(0.6, 2.5),
                        cutoff    = 0.4,
                        nu_x      = c(rough = 0.5, medium = 1, smooth = 2)) {

  # Fixed mesh on [0, L]^2 with a boundary extension to limit edge
  # effects (the field is observed only on the inner square)
  dom  <- cbind(c(0, L, L, 0), c(0, 0, L, L))
  bnd  <- fm_extensions(dom, convex = c(0.2 * L, 0.4 * L))
  mesh <- fm_mesh_2d(boundary = bnd, max.edge = max.edge, cutoff = cutoff)

  fem    <- fm_fem(mesh)
  cdiag  <- as.numeric(Matrix::rowSums(fem$c0))   # mass-lumped diagonal
  C0     <- Matrix::Diagonal(x = cdiag)
  Ci     <- Matrix::Diagonal(x = 1 / cdiag)
  G      <- fem$g1

  # kappa from the practical-range convention range = sqrt(8 nu)/kappa
  kappa <- sqrt(8 * nu_u) / range_u
  K     <- kappa^2 * C0 + G

  # Calibrate tau so the FEM field U = tau^{-1} K^{-1} w, w ~ N(0, C0),
  # has interior marginal variance var_u.  Var(U_i) = tau^{-2} (P1^{-1})_ii
  # with P1 = K C0^{-1} K. Average over a cluster of central nodes.
  P1   <- Matrix::t(K) %*% Ci %*% K
  d2c  <- (mesh$loc[, 1] - L / 2)^2 + (mesh$loc[, 2] - L / 2)^2
  idx0 <- order(d2c)[seq_len(min(25, mesh$n))]
  Lf   <- Matrix::Cholesky(P1)
  v0   <- vapply(idx0, function(i) {
    e <- numeric(mesh$n); e[i] <- 1
    as.numeric(Matrix::solve(Lf, e)[i])
  }, numeric(1))
  tau  <- sqrt(mean(v0) / var_u)

  # Operators for generating the covariate field X at the requested
  # smoothness levels (unit marginal variance, same practical range)
  xops <- lapply(nu_x, function(nx) {
    matern.operators(range = range_u, sigma = 1, nu = nx,
                     mesh = mesh, parameterization = "matern")
  })
  names(xops) <- names(nu_x)

  list(mesh = mesh, cdiag = cdiag, C0 = C0, G = G, K = K,
       kappa = kappa, tau = tau,
       range_u = range_u, var_u = var_u, L = L, xops = xops)
}

# Sample U at the mesh nodes from the calibrated mass-lumped precision Q,
sample_U <- function(world) {
  m <- world$mesh$n
  w <- sqrt(world$cdiag) * rnorm(m)
  as.numeric((1 / world$tau) * Matrix::solve(world$K, w))
}

# Scenario table.
SCENARIOS <- list(
  S1 = list(label = "S1: forced-SPDE truth",  gen = "forced",
            x_smooth = "medium", beta_pw = 0, beta_fc = 1, beta0 = 1),
  S2 = list(label = "S2: rough X, additive",  gen = "additive",
            x_smooth = "rough",  beta_pw = 1, beta_fc = 0, beta0 = 1),
  S3 = list(label = "S3: smooth X, additive", gen = "additive",
            x_smooth = "smooth", beta_pw = 1, beta_fc = 0, beta0 = 1),
  S4 = list(label = "S4: mismatched operator", gen = "mismatch",
            x_smooth = "medium", beta_pw = 0, beta_fc = 1, beta0 = 1)
)


# Draw one replicate for a given scenario and n.
# Returns the observed data frame plus everything an estimator needs.
simulate_dataset <- function(world, scenario, n, sigma_e = 0.1,
                             h_kernel = 1.5) {
  sc  <- SCENARIOS[[scenario]]
  L   <- world$L

  # Uniform observation locations on the inner square and the projector
  loc <- cbind(runif(n, 0, L), runif(n, 0, L))
  A   <- fm_basis(world$mesh, loc = loc)

  # Covariate field at the mesh nodes (and at the observation locations)
  Xn   <- as.numeric(simulate(world$xops[[sc$x_smooth]]))
  Xobs <- as.numeric(A %*% Xn)

  # Random effect
  u <- sample_U(world)

  # Latent signal at the mesh nodes (the forcing channel, if any)
  if (sc$gen == "forced") {
    # beta_fc * S X with S = L^{-1} (alpha = 2) at the true kappa
    fmod   <- hybrid.spde(mesh = world$mesh, alpha = 2, X = Xn,
                          beta_X = sc$beta_fc, kappa = world$kappa,
                          tau = world$tau)
    signal <- fmod$mu
  } else if (sc$gen == "mismatch") {
    # beta_fc * S0 X with S0 a Gaussian-kernel smoother (mismatched)
    S0x    <- gaussian_kernel_smooth(world$mesh$loc[, 1:2], Xn,
                                     world$mesh$loc[, 1:2], h_kernel)
    signal <- sc$beta_fc * S0x
  } else {                       # additive
    signal <- rep(0, world$mesh$n)
  }

  # Observations: pointwise channel uses X at the observation locations
  y <- sc$beta_pw * Xobs + as.numeric(A %*% (signal + u)) + rnorm(n, sd = sigma_e)

  data <- data.frame(y = y, x1 = loc[, 1], x2 = loc[, 2], X = Xobs)

  list(data = data, Xnodes = Xn, loc = loc, A = A,
       truth = list(beta_pw = sc$beta_pw, beta_fc = sc$beta_fc,
                    beta0 = sc$beta0),
       scenario = scenario, n = n, sigma_e = sigma_e)
}
