# ---------------------------------------------------------------------
# Utilities shared across the two applications and the simulation study
# ---------------------------------------------------------------------

# Gaussian-kernel smoother.  Smooth values `y` observed at the rows of
# `x` onto the query points `q` with bandwidth `h` (row-normalised
# Gaussian weights, so a constant field is reproduced exactly).  `x` and
# `q` are coordinate matrices with matching columns.
gaussian_kernel_smooth <- function(x, y, q, h) {
  D2 <- outer(rowSums(q^2), rowSums(x^2), `+`) - 2 * tcrossprod(q, x)
  W  <- exp(-D2 / (2 * h^2))
  rs <- rowSums(W)
  rs[rs == 0] <- 1
  as.numeric((W %*% y) / rs)
}

# Channel-separation index vartheta(X) in [0,1).
channel_separation_index <- function(K, Cd, x, center = TRUE) {
  x <- as.numeric(x)
  if (center) x <- x - sum(Cd * x) / sum(Cd)   # remove constant component
  Kx  <- as.numeric(K %*% x)
  q11 <- sum(Kx * Kx / Cd)     # (K x)' C^{-1} (K x), lumped C
  q12 <- sum(x * Kx)           # x' K x
  q22 <- sum(Cd * x * x)       # x' C x
  1 - q12^2 / (q11 * q22)
}
