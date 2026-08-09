## ---------------------------------------------------------------------
## standardize.R -- scale-free reporting of the forcing channel
##
## Implements the standardized quadruple of Section "A scale-free
## parametrisation of the forcing coefficient" (Proposition
## prop:scalefree) of the paper:
##
##   beta_star = beta_fc / kappa^alpha   (DC-standardised coefficient;
##               NOTE: in the rspde_lme / hybrid.spde implementation the
##               forcing enters the mean as L^{-alpha/2}(beta_fc X), with
##               no tau factor, so the implementation-level DC gain is
##               kappa^{-alpha}, not (tau kappa^alpha)^{-1} as in the
##               paper's parametrisation of the SPDE. beta_star is the
##               same quantity in both conventions.)
##   g_X       = ||S-tilde X||_{L2} / ||X||_{L2} in (0, 1]
##               (covariate-adapted effective gain; S-tilde = unit-DC-gain
##               smoother kappa^alpha L^{-alpha/2})
##   vartheta  = channel-separation index in [0, 1)
##               (Theorem thm:geometry; variance-inflation factor for the
##               channel split is 1/vartheta)
##   beta_star * g_X = variance-matched effective forcing coefficient
##               (equals the SD(SX)/SD(X) rescaling up to the empirical
##               vs L2 measure of spread)
##
## All quantities are computed with mass-lumped FEM matrices, matching
## the convention used in hybrid.spde / rspde_lme. Only alpha = 2 is
## supported: for fractional alpha the operator L^{-alpha/2} requires
## the rational approximation, and the standardized quantities should
## then be computed from the rational operator directly.
##
## The covariate is L2(D)-centred by default: the constant component of
## X is absorbed by the intercept (the constant is the Neumann
## eigenfunction of L), so the gain and separation indices are only
## meaningful for the non-constant part.
## ---------------------------------------------------------------------

## shared helpers (channel_separation_index used below); sourced relative
## to the project directory, matching how this file is itself sourced.
source("../common/utils.R")

standardize_forcing <- function(mesh, X_nodes, kappa, beta_fc,
                                beta_pw = NA_real_, alpha = 2,
                                center = TRUE, fem = NULL) {
  if (!isTRUE(all.equal(alpha, 2))) {
    warning("standardize_forcing: only alpha = 2 is supported; returning NA")
    return(list(beta_star = NA_real_, g_X = NA_real_,
                vartheta = NA_real_, beta_eff = NA_real_,
                beta_large = NA_real_))
  }
  if (is.null(fem)) fem <- fmesher::fm_fem(mesh)
  Cd <- as.numeric(Matrix::rowSums(fem$c0))   # lumped mass matrix diag
  G  <- fem$g1                                # stiffness matrix
  K  <- kappa^2 * Matrix::Diagonal(x = Cd) + G

  x <- as.numeric(X_nodes)
  if (center) x <- x - sum(Cd * x) / sum(Cd)  # remove constant component
  l2 <- function(v) sum(Cd * v * v)           # ||v||^2_{L2}, lumped

  ## Unit-DC-gain smoothed covariate: S-tilde x = kappa^2 K^{-1} C x
  Sx  <- kappa^alpha * as.numeric(Matrix::solve(K, Cd * x))
  g_X <- sqrt(l2(Sx) / l2(x))

  ## Channel-separation index vartheta(X) (x already centred above).
  vartheta <- channel_separation_index(K, Cd, x, center = FALSE)

  beta_star <- beta_fc / kappa^alpha
  list(beta_star  = as.numeric(beta_star),
       g_X        = as.numeric(g_X),
       vartheta   = as.numeric(vartheta),
       beta_eff   = as.numeric(beta_star * g_X),
       beta_large = as.numeric(beta_pw + beta_star))
}

## ---------------------------------------------------------------------
## joint_separation_forcing -- multivariate channel-separation diagnostic
## (Remark "Several forcing covariates" / rem:multivariate).
##
## With p forcing covariates the hybrid model regresses on the 2p functions
## (X_1,...,X_p, S X_1,...,S X_p). Their geometry is the 2p x 2p
## Cameron--Martin Gram G_p, assembled from the same sparse matrices as the
## fit (K = kappa^2 C + G, lumped mass C, stiffness G):
##
##   forcing--forcing      (S X_j, S X_k)_C = x_j' C x_k
##   pointwise--forcing    (X_j, S X_k)_C   = tau  x_j' K x_k
##   pointwise--pointwise  (X_j, X_k)_C     = tau^2 (K x_j)' C^{-1} (K x_k)
##
## Normalising G_p to a correlation matrix R_p, the m-th coefficient has
## variance-inflation factor VIF_m = (R_p^{-1})_mm >= 1 and generalised
## channel-separation index vartheta_m = 1/(R_p^{-1})_mm = 1 - R_m^2 in
## (0,1], the tolerance of that regressor against the other 2p-1. For p = 1
## this reduces exactly to the marginal vartheta(X) of standardize_forcing.
## (tau cancels under the correlation normalisation, so its value is
## immaterial; it is carried only to match the construction of the remark.)
##
## Returns a data.frame with one row per coefficient (p pointwise rows then
## p forcing rows), plus the correlation matrix R_p as an attribute.
## Only alpha = 2 is supported (as for standardize_forcing).
## ---------------------------------------------------------------------
joint_separation_forcing <- function(mesh, X_list, kappa, alpha = 2,
                                     tau = 1, center = TRUE, fem = NULL,
                                     names = NULL) {
  if (!isTRUE(all.equal(alpha, 2)))
    stop("joint_separation_forcing: only alpha = 2 is supported")
  if (is.matrix(X_list) || is.data.frame(X_list))
    X_list <- lapply(seq_len(ncol(X_list)), function(j) as.numeric(X_list[, j]))
  p <- length(X_list)
  if (is.null(names)) names <- if (!is.null(base::names(X_list)))
    base::names(X_list) else paste0("X", seq_len(p))

  if (is.null(fem)) fem <- fmesher::fm_fem(mesh)
  Cd <- as.numeric(Matrix::rowSums(fem$c0))          # lumped mass diag
  G  <- fem$g1                                        # stiffness
  K  <- kappa^2 * Matrix::Diagonal(x = Cd) + G

  ## centred covariate node vectors and their K-images
  X  <- lapply(X_list, function(v) { v <- as.numeric(v)
    if (center) v <- v - sum(Cd * v) / sum(Cd); v })
  KX <- lapply(X, function(v) as.numeric(K %*% v))

  ## channel order: pointwise 1..p, then forcing 1..p
  Gm <- matrix(0, 2 * p, 2 * p)
  for (j in seq_len(p)) for (k in seq_len(p)) {
    pp <- tau^2 * sum(KX[[j]] * KX[[k]] / Cd)         # (X_j, X_k)_C
    ff <- sum(Cd * X[[j]] * X[[k]])                   # (S X_j, S X_k)_C
    pf <- tau * sum(X[[j]] * KX[[k]])                 # (X_j, S X_k)_C
    Gm[j,     k]     <- pp
    Gm[p + j, p + k] <- ff
    Gm[j,     p + k] <- pf
    Gm[p + k, j]     <- pf
  }
  d  <- sqrt(diag(Gm))
  Rp <- Gm / outer(d, d)                              # correlation matrix
  Ri <- solve(Rp)
  vif <- diag(Ri)                                     # (R_p^{-1})_mm

  data.frame(
    coef      = c(paste0("pw_",  names), paste0("fc_", names)),
    covariate = rep(names, 2),
    type      = rep(c("pointwise", "forcing"), each = p),
    vif       = as.numeric(vif),
    vartheta  = as.numeric(1 / vif),
    row.names = NULL,
    check.names = FALSE
  ) -> out
  attr(out, "R_p") <- Rp
  out
}

## Robust extractors for rspde_lme fits (random_effects ordering:
## [1] alpha, [2] tau/sigma, [3] kappa, [4] beta_fc for hybrid.spde fits;
## prefer name-based lookup when names are present).
re_par <- function(res, name, pos) {
  if (is.null(res)) return(NA_real_)
  re <- res$coeff$random_effects
  if (!is.null(names(re)) && name %in% names(re)) return(as.numeric(re[[name]]))
  ln <- paste0("log_", name)
  if (!is.null(names(re)) && ln %in% names(re)) return(exp(as.numeric(re[[ln]])))
  if (length(re) >= pos) return(as.numeric(re[[pos]])) else return(NA_real_)
}
alpha_from_fit <- function(res) re_par(res, "alpha", 1)
kappa_from_fit <- function(res) re_par(res, "kappa", 3)
betafc_from_fit <- function(res) re_par(res, "beta", 4)

## One-line summary for a hybrid rspde_lme fit
standardized_summary <- function(res, mesh, X_nodes, x_name = "elev",
                                 fem = NULL) {
  fx <- res$coeff$fixed_effects
  beta_pw <- if (x_name %in% names(fx)) as.numeric(fx[[x_name]]) else NA_real_
  kap <- kappa_from_fit(res)
  alf <- alpha_from_fit(res)
  bfc <- betafc_from_fit(res)
  st  <- standardize_forcing(mesh, X_nodes, kappa = kap, beta_fc = bfc,
                             beta_pw = beta_pw, alpha = round(alf, 8),
                             fem = fem)
  data.frame(alpha = alf, kappa = kap, range = sqrt(8 * (alf - 1)) / kap,
             beta_pw = beta_pw, beta_fc = bfc,
             beta_star = st$beta_star, g_X = st$g_X,
             vartheta = st$vartheta, beta_eff = st$beta_eff,
             loglik = res$loglik)
}
