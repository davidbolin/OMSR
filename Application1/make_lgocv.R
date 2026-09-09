# --------------------------------------------------------------------
# Build the Section 5 CV figure
# Run after fit_models.R and make_tables.R.
# --------------------------------------------------------------------

library(rSPDE)
library(Matrix)
library(fmesher)
library(ggplot2)
library(patchwork)

REGIONS   <- c("COLORADO", "NORWAYANNUAL", "NORWAYDJF")
D_GRID    <- c(0, 25, 50, 75, 100, 150, 200)   # buffer radii (km)
MIN_TRAIN <- 8                                # skip a station if too few training points remain
DM_NBLOCK <- 4      # DM block bootstrap: NBLOCK x NBLOCK spatial grid of stations
DM_NBOOT  <- 2000   # number of block-bootstrap resamples
DM_SEED   <- 1      # seed for block-bootstrap reproducibility

# per-observation Gaussian predictive scores
obs_scores <- function(mu, v, y) list(logs  = 0.5*log(2*pi*v) + 0.5*(y-mu)^2/v,
                                      sqerr = (y-mu)^2)
# paired mean loss difference (model - reference) with a spatial block-
# bootstrap standard error: the station domain is split into a
# DM_NBLOCK x DM_NBLOCK grid and whole blocks are resampled with replacement,
# so the SE reflects spatial dependence between neighbouring stations' scores.
# `blocks` is a list of station-index vectors, one per non-empty block.
paired_dm <- function(loss_m, loss_ref, blocks, B) {
  d  <- loss_m - loss_ref
  md <- mean(d, na.rm = TRUE)
  K  <- length(blocks)
  bm <- replicate(B, mean(d[unlist(blocks[sample.int(K, K, replace = TRUE)],
                                   use.names = FALSE)], na.rm = TRUE))
  se <- sd(bm)
  c(diff = md, se = se, t = md/se, p = 2*pnorm(-abs(md/se)))
}

compute <- function(tag) {
  e <- new.env()
  load(sprintf("results/fits_%s.RData", tag), envir = e)
  data <- e$p$data
  fits <- list(Additive = e$fits$Additive_a2, Hybrid = e$fits$Hybrid_a2)
  if ("Coast_a2" %in% names(e$fits)) {           # Norway: add coast variants
    fits <- list(Additive = fits$Additive, AddCoast = e$addcoast,
                 Hybrid = fits$Hybrid, HybCoast = e$fits$Coast_a2)
  }
  fits <- c(fits, list(OLS = e$fits$OLS, Forced = e$fits$Forced_a2,
                       BW = e$fits$BW_a2, SpatPlus = e$fits$SpatPlus_a2))
  loc <- e$p$loc_km
  n <- nrow(loc)
  Dm <- as.matrix(dist(loc))
  y <- data$y
  MU <- VAR <- array(NA, dim = c(n, length(fits), length(D_GRID)),
                     dimnames = list(NULL, names(fits), as.character(D_GRID)))
  for (di in seq_along(D_GRID)) {
    d <- D_GRID[di]
    for (i in 1:n) {
      grp <- which(Dm[i, ] <= d)
      tr <- setdiff(1:n, grp)   # buffer group masked; condition on rest
      if (length(tr) < MIN_TRAIN) next
      pr <- rSPDE::posterior_crossvalidation(object = fits,
              train_test_indices = list(list(train = tr, test = grp)),
              true_CV = FALSE, print = FALSE)
      for (mi in seq_along(fits)) {
        MU[i, mi, di] <- pr$mu[[mi]][i]
        VAR[i, mi, di] <- pr$var[[mi]][i]
      }
    }
    cat(sprintf("  %s d=%3d km\n", tag, d))
  }
  list(MU = MU, VAR = VAR, y = y, D = D_GRID, models = names(fits), loc = loc)
}

rows <- list()
dm_rows <- list()
for (tag in REGIONS) {
  o <- compute(tag)
  saveRDS(o, sprintf("results/lgocv_perstation_%s.rds", tag))
  R <- apply(o$MU, c(2, 3), function(mu) sqrt(mean((o$y - mu)^2, na.rm = TRUE)))
  L <- sapply(seq_along(o$D), function(di) sapply(seq_along(o$models), function(mi) {
        mu <- o$MU[, mi, di]; v <- o$VAR[, mi, di]
        mean(0.5*log(2*pi*v) + 0.5*(o$y - mu)^2/v, na.rm = TRUE) }))
  rows[[tag]] <- data.frame(region = tag,
                            d_km     = rep(o$D, each = length(o$models)),
                            model    = rep(o$models, length(o$D)),
                            rmse     = as.vector(R),
                            logscore = as.vector(L))

  # paired Diebold-Mariano at d = 0 (leave-one-out), reference = Hybrid, with a
  # spatial block-bootstrap SE over a DM_NBLOCK x DM_NBLOCK grid.
  i0     <- which(o$D == 0)
  bx     <- cut(o$loc[, 1], DM_NBLOCK, labels = FALSE)
  by     <- cut(o$loc[, 2], DM_NBLOCK, labels = FALSE)
  blocks <- split(seq_along(o$y), paste(bx, by))   # non-empty blocks only
  ref    <- obs_scores(o$MU[, "Hybrid", i0], o$VAR[, "Hybrid", i0], o$y)
  set.seed(DM_SEED)
  paired <- do.call(rbind, lapply(setdiff(o$models, "Hybrid"), function(m) {  # skip self-comparison
    s  <- obs_scores(o$MU[, m, i0], o$VAR[, m, i0], o$y)
    dl <- paired_dm(s$logs,  ref$logs,  blocks, DM_NBOOT)
    dq <- paired_dm(s$sqerr, ref$sqerr, blocks, DM_NBOOT)
    data.frame(region = tag, model = m,
               logs_diff = dl["diff"], logs_se = dl["se"], logs_p = dl["p"],
               mse_diff = dq["diff"], mse_se = dq["se"], mse_p = dq["p"],
               rmse = sqrt(mean(s$sqerr, na.rm = TRUE)), row.names = NULL)
  }))
  write.csv(paired, sprintf("results/paired_%s.csv", tag), row.names = FALSE)
  dm_rows[[tag]] <- paired
}
write.csv(do.call(rbind, rows), "results/lgocv_curves.csv", row.names = FALSE)

# LOO paired DM vs Hybrid
cat("\n== Leave-one-out (d=0) paired Diebold-Mariano vs Hybrid ==\n")
dm_all <- do.call(rbind, dm_rows); rownames(dm_all) <- NULL
print(within(dm_all, {
  logs_p <- round(logs_p, 3)
  mse_p <- round(mse_p, 3)
  logs_diff <- round(logs_diff, 4)
  mse_diff <- round(mse_diff, 4)
})[, c("region","model","logs_diff","logs_p","mse_diff","mse_p")], row.names = FALSE)

regions <- c(COLORADO = "Colorado MAM", NORWAYANNUAL = "Norway annual", NORWAYDJF = "Norway DJF")
mod_lev <- c("Additive", "AddCoast", "Hybrid", "HybCoast", "OLS", "Forced", "BW", "SpatPlus")
mod_lab <- c(Additive = "Additive", AddCoast = "Additive + coast", Hybrid = "Hybrid",
             HybCoast = "Hybrid + coast", OLS = "OLS", Forced = "Forced", BW = "BW",
             SpatPlus = "Spatial+")
COL <- c(Additive = "grey45", AddCoast = "#377eb8", Hybrid = "#e41a1c", HybCoast = "#4daf4a",
         OLS = "black", Forced = "#984ea3", BW = "#ff7f00", SpatPlus = "#a65628")
LTY <- c(Additive = "solid", AddCoast = "solid", Hybrid = "solid", HybCoast = "solid",
         OLS = "dashed", Forced = "dashed", BW = "dashed", SpatPlus = "dashed")

rows <- list()
for (tg in names(regions)) {
  O <- readRDS(sprintf("results/lgocv_perstation_%s.rds", tg))
  for (mi in seq_along(O$models)) for (di in seq_along(O$D)) {
    mu <- O$MU[, mi, di]; v <- O$VAR[, mi, di]
    rows[[length(rows) + 1]] <- data.frame(
      region = regions[tg], model = O$models[mi], d = O$D[di],
      rmse = sqrt(mean((O$y - mu)^2, na.rm = TRUE)),
      logs = mean(0.5*log(2*pi*v) + 0.5*(O$y - mu)^2/v, na.rm = TRUE))
  }
}
d <- do.call(rbind, rows)
d$region <- factor(d$region, levels = unname(regions))
d$model  <- factor(d$model, levels = mod_lev)

row_plot <- function(yvar, ylab, show_strip, show_x) {
  ggplot(d, aes(.data[["d"]], .data[[yvar]], colour = model, linetype = model)) +
    geom_line(linewidth = 0.6) + geom_point(size = 0.8) +
    facet_wrap(~ region, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = COL, breaks = mod_lev, labels = mod_lab) +
    scale_linetype_manual(values = LTY, breaks = mod_lev, labels = mod_lab) +
    guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
    labs(x = if (show_x) expression("prediction buffer radius  " * italic(d) * "  (km)") else NULL,
         y = ylab) +
    theme_bw(base_size = 11) +
    theme(legend.title = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = if (show_strip) element_text() else element_blank())
}

fig <- row_plot("rmse", expression("LGO-CV RMSE (" * degree * "C)"), TRUE,  FALSE) /
  row_plot("logs", "LGO-CV log score",                          FALSE, TRUE) +
  plot_layout(guides = "collect") & theme(legend.position = "top")

print(fig)

