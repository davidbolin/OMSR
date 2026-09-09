# --------------------------------------------------------------------
# Make the leave-distance-out cross-validation figure.
# Run after fit_models.R.
# --------------------------------------------------------------------

library(INLA)
library(rSPDE)
library(ggplot2)
library(patchwork)
library(RANN)

load("results/be_no2_fits.RData")

Kmax <- 600L
nn <- RANN::nn2(loc_km, k = min(Kmax, n))   # self is the 1st neighbour (dist 0)
D_GRID <- c(0,0.25,0.5,1,2,3,5)
scoreset <- function(f, grp){
  g <- inla.group.cv(f, groups = grp)
  if (is.null(g$cv)) return(c(rmse=NA, logs=NA))
  c(rmse = sqrt(mean((y-g$mean)^2, na.rm=TRUE)), logs = -mean(log(g$cv), na.rm=TRUE))
}

cur <- list()
cat("--- leave-distance-out CV curves ---\n")
for (db in D_GRID) {
  grp <- lapply(seq_len(n), function(i){ dd<-nn$nn.dists[i,]; nn$nn.idx[i, dd<=db] })
  for (spec in names(fits)) {
    sc <- scoreset(fits[[spec]], grp)
    cur[[length(cur)+1]] <- data.frame(model=spec, buffer=db,
                                       rmse=sc["rmse"], logs=sc["logs"], row.names=NULL)
  }
  cat("  dist d=",db,"km\n")
}
d <- do.call(rbind, cur)

lev <- c("field-only","LUR (OLS)","LUR+field",
         "forced(bg+op)","hybrid(bg+op)","two-scale","three-field","split(ind|loc)")
lev <- lev[lev %in% unique(d$model)]
relab <- c("field-only"="field only","LUR (OLS)"="LUR (OLS)",
           "LUR+field"="LUR + field","forced(bg+op)"="operator forced",
           "hybrid(bg+op)"="operator hybrid","two-scale"="two-scale",
           "three-field"="three-field",
           "split(ind|loc)"="two-field (industry split)")

d <- d[d$model %in% lev, ]
d$model <- factor(relab[as.character(d$model)], levels = unname(relab[lev]))

xlab_d   <- "removal radius around test point (km)"
ylab_for <- list(logs = "CV log-score (lower better)",
                 rmse = expression("CV RMSE ("*mu*"g m"^-3*", lower better)"))

COL <- c("field only"="grey55", "LUR (OLS)"="black", "LUR + field"="#1b9e77",
         "operator forced"="#d95f02", "operator hybrid"="#7570b3", "two-scale"="#e7298a",
         "three-field"="#a6761d", "two-field (industry split)"="#e31a1c")
LTY <- c("field only"="dotted", "LUR (OLS)"="solid", "LUR + field"="twodash",
         "operator forced"="dashed", "operator hybrid"="dotdash", "two-scale"="longdash",
         "three-field"="dashed", "two-field (industry split)"="solid")

panel <- function(score) {
  dd <- d
  dd$score <- dd[[score]]
  ggplot(dd, aes(buffer, score, colour = model, linetype = model)) +
    geom_line(linewidth = 0.6) + geom_point(size = 0.9) +
    scale_colour_manual(values = COL) + scale_linetype_manual(values = LTY) +
    theme_bw(base_size = 11) +
    theme(legend.title = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = xlab_d, y = ylab_for[[score]])
}

print((panel("logs") | panel("rmse")) + plot_layout(guides = "collect"))
