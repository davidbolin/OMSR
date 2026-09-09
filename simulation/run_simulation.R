# ---------------------------------------------------------------------
# Main file for the simulation study
#
# Builds the fixed world (mesh + true operator), then loops over the
# grid (scenario x n x replicate), fitting all estimators of
# to each replicate.  Replicates are run in parallel.
# ---------------------------------------------------------------------

source("dgm.R")
source("estimators.R")
library(parallel)

R         <- 100                                  # Monte Carlo replicates per (scenario, n)
n_grid    <- c(50, 100, 250, 500, 1000)           # sample sizes
n_cores   <- max(1, parallel::detectCores() - 1)  # parallel workers
sigma_e   <- 0.1                                  # observation noise sd
base_seed <- 2026
out_path  <- "sim_results.rds"
scenarios <- names(SCENARIOS)

cat(sprintf("Simulation study: R=%d, n={%s}, scenarios={%s}, sigma_e=%g, cores=%d\n",
            R, paste(n_grid, collapse = ","), paste(scenarios, collapse = ","),
            sigma_e, n_cores))

# build the fixed world once
cat("Building mesh and calibrating the true operator ...\n")
world <- build_world()
cat(sprintf("  mesh nodes = %d, kappa = %.4f, tau = %.4f (range %.2f, var %.2f)\n",
            world$mesh$n, world$kappa, world$tau, world$range_u, world$var_u))

# one replicate: fit every estimator on a (scenario, n, rep)
# The seed depends only on (scenario, n, rep), so the run is reproducible
# regardless of how the grid is scheduled.
run_one <- function(scenario, n, rep) {
  set.seed(base_seed + 1000 * match(scenario, scenarios) +
           97 * which(n_grid == n) + rep)
  ds  <- simulate_dataset(world, scenario, n, sigma_e = sigma_e)
  res <- fit_all(ds, world)
  res$rep      <- rep
  res$vartheta <- channel_separation_index(world$K, world$cdiag, ds$Xnodes)
  res
}

# run the grid block by block, printing progress
n_cells <- length(n_grid) * length(scenarios)
cat(sprintf("Total fits: %d datasets x 7 estimators, in %d (scenario x n) blocks\n",
            n_cells * R, n_cells))

t0 <- Sys.time()
blocks    <- vector("list", n_cells)
cell      <- 0L
succeeded <- 0L
for (n in n_grid) for (sc in scenarios) {
  cell <- cell + 1L
  b0   <- Sys.time()
  part <- mclapply(seq_len(R), function(rep) {
    tryCatch(run_one(sc, n, rep), error = function(e) {
      message(sprintf("  [n=%d %s rep=%d] failed: %s", n, sc, rep,
                      conditionMessage(e))); NULL
    })
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  ok             <- sum(!vapply(part, is.null, logical(1)))
  succeeded      <- succeeded + ok
  blocks[[cell]] <- do.call(rbind, part)
  cat(sprintf("[%2d/%2d] n=%4d %-3s : %d/%d replicates in %5.1fs  (total elapsed %.1f min)\n",
              cell, n_cells, n, sc, ok, R,
              as.numeric(difftime(Sys.time(), b0, units = "secs")),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
elapsed <- difftime(Sys.time(), t0, units = "mins")

results <- do.call(rbind, blocks)
cat(sprintf("Done in %.1f min. Collected %d rows (%d/%d replicates succeeded).\n",
            as.numeric(elapsed), nrow(results), succeeded, n_cells * R))

attr(results, "config") <- list(R = R, n_grid = n_grid, sigma_e = sigma_e,
                                 base_seed = base_seed, kappa = world$kappa,
                                 tau = world$tau, mesh_n = world$mesh$n)
saveRDS(results, out_path)
cat(sprintf("Saved raw results to %s\n", out_path))
