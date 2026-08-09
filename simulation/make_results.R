## ---------------------------------------------------------------------
## Aggregate the raw replicates into summary metrics, and produce figure
## ---------------------------------------------------------------------
source("dgm.R")                        # SCENARIOS labels / ordering
suppressMessages({ library(ggplot2); library(ggh4x) })

results <- readRDS("sim_results.rds")

Z975       <- qnorm(0.975)
DIV_THRESH <- 4                        # |beta_hat - truth| > DIV_THRESH counts as diverged
EST_ORDER  <- c("Additive-GLS", "BW-empirical", "Spatial+", "RSR",
                "Forced-fixed", "Forced", "Hybrid (pw)", "Hybrid (fc)")
SCEN_ORDER <- names(SCENARIOS)

## ---- aggregate the replicate-level results into summary metrics ------
## Mean-based bias/RMSE/coverage/width (as reported in the paper) plus the
## robust median bias; `fail` = share of non-converged or diverged fits.
aggregate_metrics <- function(df) {
  keys <- with(df, paste(scenario, n, estimator, param, sep = "\r"))
  do.call(rbind, lapply(split(df, keys), function(d_all) {
    truth <- d_all$truth[1]; ntot <- nrow(d_all)
    d  <- d_all[d_all$ok & is.finite(d_all$estimate), ]
    lo <- d$estimate - Z975 * d$se; hi <- d$estimate + Z975 * d$se
    have_se  <- is.finite(d$se)
    diverged <- abs(d$estimate - truth) > DIV_THRESH
    data.frame(
      scenario = d_all$scenario[1], n = d_all$n[1],
      estimator = d_all$estimator[1], param = d_all$param[1],
      bias  = mean(d$estimate) - truth,
      rmse  = sqrt(mean((d$estimate - truth)^2)),
      cov95 = if (any(have_se)) mean((lo <= truth & truth <= hi)[have_se]) else NA,
      width = if (any(have_se)) mean((hi - lo)[have_se]) else NA,
      med_bias = median(d$estimate) - truth,
      fail  = (ntot - nrow(d) + sum(diverged)) / ntot,   # non-conv + diverged
      time  = mean(d$time),
      stringsAsFactors = FALSE)
  }))
}
metrics <- aggregate_metrics(results)
metrics$estimator <- factor(metrics$estimator, levels = EST_ORDER)
metrics$scenario  <- factor(metrics$scenario,  levels = SCEN_ORDER)
metrics <- metrics[order(metrics$n, metrics$scenario, metrics$estimator), ]

## ---- channel-separation index vartheta(X) per scenario (quoted in the
## appendix as 0.46, 0.51, 0.36, 0.46 for S1-S4) -----------------------
if ("vartheta" %in% names(results)) {
  vth <- aggregate(vartheta ~ scenario, data = results, FUN = mean)
  vth <- vth[match(SCEN_ORDER, vth$scenario), ]
  cat("Channel-separation index vartheta(X) (mean over replicates):\n")
  print(data.frame(Scenario = vth$scenario, vartheta = sprintf("%.3f", vth$vartheta)),
        row.names = FALSE)
  cat("\n")
}

## ---- summary metrics per sample size (the numbers behind the figure) -
for (nn in sort(unique(metrics$n))) {
  cat(sprintf("================ n = %d ================\n", nn))
  sub <- metrics[metrics$n == nn, ]
  print(data.frame(
    Scenario  = as.character(sub$scenario),
    Estimator = as.character(sub$estimator),
    Bias    = sprintf("% .2f", sub$bias),
    RMSE    = sprintf("%.2f", sub$rmse),
    Cov95   = sprintf("%.2f", sub$cov95),
    Width   = sprintf("%.2f", sub$width),
    MedBias = sprintf("% .2f", sub$med_bias),
    Fail    = sprintf("%.0f", 100 * sub$fail),
    Time    = sprintf("%.2f", sub$time)), row.names = FALSE)
  cat("\n")
}

## ---- the appendix curve figure (Figure fig:sim-curves) --------------
## bias / RMSE / coverage vs n, metric x scenario facets, one line per
## estimator; competitors dashed, the operator-matched estimators solid.
SCEN_LAB <- c(S1 = "S1: forced truth", S2 = "S2: rough X, additive",
              S3 = "S3: smooth X, additive", S4 = "S4: mismatched op.")
## Okabe-Ito colour-blind-safe palette (8 estimators).
PAL <- c("Additive-GLS" = "#000000", "BW-empirical" = "#E69F00",
         "Spatial+" = "#56B4E9", "RSR" = "#009E73",
         "Forced-fixed" = "#F0E442", "Forced" = "#0072B2",
         "Hybrid (pw)" = "#D55E00", "Hybrid (fc)" = "#CC79A7")
COMPET     <- c("Additive-GLS", "BW-empirical", "Spatial+", "RSR")   # drawn dashed
MET_LEVELS <- c("Bias", "RMSE", "Coverage")

base <- metrics[, c("scenario", "n", "estimator")]
long <- do.call(rbind, list(
  cbind(base, metric = "Bias",     value = metrics$bias),
  cbind(base, metric = "RMSE",     value = metrics$rmse),
  cbind(base, metric = "Coverage", value = metrics$cov95)))
long <- long[long$estimator %in% EST_ORDER, ]
long$estimator <- factor(long$estimator, levels = EST_ORDER)
long$metric    <- factor(long$metric, levels = MET_LEVELS)
long$scenario  <- factor(long$scenario, levels = names(SCEN_LAB))
long$lt <- ifelse(as.character(long$estimator) %in% COMPET, "competitor", "proposed")

## per-row reference lines: bias -> 0, coverage -> 0.95
refs <- data.frame(metric = factor(c("Bias", "Coverage"), levels = levels(long$metric)),
                   y = c(0, 0.95))

fig <- ggplot(long, aes(n, value, colour = estimator, linetype = lt)) +
  geom_hline(data = refs, aes(yintercept = y), colour = "grey60",
             linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.1) +
  ## metric rows share nothing (independent y per panel), scenario columns on top
  facet_grid2(metric ~ scenario, scales = "free_y", independent = "y", switch = "y") +
  scale_x_log10(breaks = sort(unique(long$n))) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_linetype_manual(values = c(competitor = "22", proposed = "solid"), guide = "none") +
  labs(x = "sample size n", y = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.background.x = element_rect(fill = "grey95", colour = NA),
        strip.background.y = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 90, face = "bold"),
        legend.key.width = unit(1.2, "lines")) +
  guides(colour = guide_legend(nrow = 1))

print(fig)
