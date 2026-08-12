## ---------------------------------------------------------------------
## The simulation study on the empirical size of the two likelihood-ratio tests
## used in the applications.
## ---------------------------------------------------------------------

source("dgm.R")                        
library(rSPDE)
library(parallel)

R         <- 200                                  # replicates per (test, truth, n)
n_grid    <- c(64, 100, 250)                      # sample sizes (cover Section 5)
n_cores   <- max(1, parallel::detectCores() - 1)  # parallel workers
sigma_e   <- 0.1                                  # observation noise sd
base_seed <- 2026


## the three (test, truth) cells and their printed labels
cells <- list(
  list(test = "matching",   scen = "S1", label = "Matching (kappa_mu = kappa)"),
  list(test = "hyb_vs_add", scen = "S2", label = "Hybrid vs. additive"),
  list(test = "hyb_vs_add", scen = "S3", label = "Hybrid vs. additive"))

cat(sprintf("LRT size check: R=%d, n={%s}, cores=%d\n",
            R, paste(n_grid, collapse = ","), n_cores))

world  <- build_world()
model1 <- spde.matern.operators(mesh = world$mesh, alpha = 2)  
cat(sprintf("  mesh nodes = %d, kappa = %.4f, tau = %.4f\n",
            world$mesh$n, world$kappa, world$tau))

## Safe fit: swallow errors/warnings, return NULL on failure.
fit_safe <- function(expr) tryCatch(suppressWarnings(force(expr)), error = function(e) NULL)

## One replicate -> the LR statistic 2*Delta-loglik for the given test
## (NA if any fit fails).  Seed depends only on (test, scen, n, rep_id).
run_one <- function(test, scen, n, rep_id) {
  set.seed(base_seed + 1e5 * match(test, c("matching", "hyb_vs_add")) +
           1000 * match(scen, names(SCENARIOS)) +
           97 * which(n_grid == n) + rep_id)
  ds <- simulate_dataset(world, scen, n, sigma_e = sigma_e)

  model2 <- hybrid.spde(mesh = world$mesh, alpha = 2, X = ds$Xnodes)
  hy <- fit_safe(rspde_lme(y ~ X, data = ds$data, model = model2,
                           loc = c("x1", "x2"),
                           model_options = list(fix_alpha = 2),
                           optim_method = "BFGS"))
  if (is.null(hy)) return(NA_real_)

  if (test == "matching") {
    kap0     <- as.numeric(hy$coeff$random_effects[["kappa"]])
    model_ts <- hybrid.spde(mesh = world$mesh, alpha = 2, X = ds$Xnodes, kappa_mu = kap0)
    ts <- fit_safe(rspde_lme(y ~ X, data = ds$data, model = model_ts,
                             loc = c("x1", "x2"),
                             model_options = rSPDE:::extract_starting_values(hy),
                             optim_method = "BFGS"))
    if (is.null(ts)) return(NA_real_)
    return(max(0, 2 * (ts$loglik - hy$loglik)))
  }

  ad <- fit_safe(rspde_lme(y ~ X, data = ds$data, model = model1,
                           loc = c("x1", "x2"),
                           model_options = list(fix_alpha = 2),
                           optim_method = "BFGS"))
  if (is.null(ad)) return(NA_real_)
  max(0, 2 * (hy$loglik - ad$loglik))
}

## run each (test, truth, n) cell
crit  <- qchisq(c(0.95, 0.99), df = 1)    # chi^2_1 critical values
ncell <- length(cells) * length(n_grid)
t0    <- Sys.time()
raw   <- list(); k <- 0L
for (cl in cells) for (n in n_grid) {
  k  <- k + 1L; b0 <- Sys.time()
  lr <- unlist(mclapply(seq_len(R), function(r)
                 tryCatch(run_one(cl$test, cl$scen, n, r),
                          error = function(e) NA_real_),
                 mc.cores = n_cores, mc.preschedule = FALSE))
  raw[[k]] <- data.frame(test = cl$test, scenario = cl$scen, n = n,
                         rep = seq_len(R), lr = lr, stringsAsFactors = FALSE)
  ok <- lr[is.finite(lr)]
  cat(sprintf("[%d/%d] %-26s %s n=%4d : size5=%.3f size1=%.3f  (%d/%d ok, %4.1fs, elapsed %.1f min)\n",
              k, ncell, cl$label, cl$scen, n,
              mean(ok > crit[1]), mean(ok > crit[2]), length(ok), R,
              as.numeric(difftime(Sys.time(), b0, units = "secs")),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, raw)
attr(res, "config") <- list(R = R, n_grid = n_grid, sigma_e = sigma_e,
                            base_seed = base_seed)
saveRDS(res, "lrt_size.rds")

## Print table as in Appendix B
rate <- function(test, scen, level) {
  vapply(n_grid, function(nn) {
    d <- res[res$test == test & res$scenario == scen & res$n == nn, "lr"]
    mean(d[is.finite(d)] > crit[level])
  }, numeric(1))
}
row_wide <- function(test, scen, label) {
  p5 <- rate(test, scen, 1); p1 <- rate(test, scen, 2)
  data.frame(Test = label, Truth = scen,
             `5%:64` = p5[1], `5%:100` = p5[2], `5%:250` = p5[3],
             `1%:64` = p1[1], `1%:100` = p1[2], `1%:250` = p1[3],
             check.names = FALSE, stringsAsFactors = FALSE)
}
tab <- rbind(
  row_wide("matching",   "S1", "Matching (kappa_mu=kappa)"),
  row_wide("hyb_vs_add", "S2", "Hybrid vs. additive"),
  row_wide("hyb_vs_add", "S3", "Hybrid vs. additive"))
tab[3:8] <- lapply(tab[3:8], function(z) sprintf("%.3f", z))

cat(sprintf("\nEmpirical size of the LRTs vs chi^2_1; R=%d per cell, MC SE ~ %.3f at 5%%:\n",
            R, sqrt(0.05 * 0.95 / R)))
print(tab, row.names = FALSE)
