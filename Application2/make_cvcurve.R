## =====================================================================
## Reproduces Figure 5 
## Run after fit_models.R.
## =====================================================================

library(INLA)
library(rSPDE)
library(ggplot2)
library(patchwork) 
library(RANN)
load("results/be_no2_fits.RData")   

## compute the buffered CV curves 
Kmax <- 600L
nn <- RANN::nn2(loc_km, k = min(Kmax, n))   # self is the 1st neighbour (dist 0)
N_GRID <- c(0,1,2,5,10,20,50,100)
N_GRID <- N_GRID[N_GRID < n]
D_GRID <- c(0,0.25,0.5,1,2,3,5)
scoreset <- function(f, grp){ 
  g <- inla.group.cv(f, groups = grp)
  if (is.null(g$cv)) return(c(rmse=NA, logs=NA))
  c(rmse = sqrt(mean((y-g$mean)^2, na.rm=TRUE)), logs = -mean(log(g$cv), na.rm=TRUE)) 
}  

cur <- list()
add <- function(typ, b, grp){ 
  for (spec in names(fits)) {
    sc <- scoreset(fits[[spec]], grp)
    cur[[length(cur)+1]] <<- data.frame(model=spec, type=typ, buffer=b, 
                                        rmse=sc["rmse"], logs=sc["logs"], 
                                        row.names=NULL) 
  } 
}

cat("--- buffered CV curves (count + distance) ---\n"); flush.console()
for (nb in N_GRID) { 
  grp <- lapply(seq_len(n), function(i) nn$nn.idx[i, 1:(nb+1)])
  add("count", nb, grp)
  cat("  count n=",nb,"\n")
}

for (db in D_GRID) { 
  grp <- lapply(seq_len(n), function(i){ 
    dd<-nn$nn.dists[i,]; nn$nn.idx[i, dd<=db] })
  add("distance", db, grp)
  cat("  dist d=",db,"km\n")
}
write.csv(do.call(rbind, cur), "results/be_no2_cvcurve.csv", row.names = FALSE)

d <- read.csv("results/be_no2_cvcurve.csv")

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

xlab_for <- list(count = "neighbours removed around test point (n)",
                 distance = "removal radius around test point (km)")
ylab_for <- list(logs = "CV log-score (lower better)",
                 rmse = expression("CV RMSE ("*mu*"g m"^-3*", lower better)"))
ttl_for  <- list(count = "Leave-out by count", distance = "Leave-out by distance")

COL <- c("field only"="grey55", "LUR (OLS)"="black", "LUR + field"="#1b9e77",
         "operator forced"="#d95f02", "operator hybrid"="#7570b3", "two-scale"="#e7298a",
         "three-field"="#a6761d", "two-field (industry split)"="#e31a1c")
LTY <- c("field only"="dotted", "LUR (OLS)"="solid", "LUR + field"="twodash",
         "operator forced"="dashed", "operator hybrid"="dotdash", "two-scale"="longdash",
         "three-field"="dashed", "two-field (industry split)"="solid")

## one panel: a given score (logs|rmse) against a given buffer type 
panel <- function(score, typ, show_title, show_xlab) {
  dd <- d[d$type == typ, ]; dd$score <- dd[[score]]
  ggplot(dd, aes(buffer, score, colour = model, linetype = model)) +
    geom_line(linewidth = 0.6) + geom_point(size = 0.9) +
    scale_colour_manual(values = COL) + scale_linetype_manual(values = LTY) +
    theme_bw(base_size = 11) +
    theme(legend.title = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 11, hjust = 0.5)) +
    labs(x = if (show_xlab) xlab_for[[typ]] else NULL,
         y = ylab_for[[score]],
         title = if (show_title) ttl_for[[typ]] else NULL)
}

## rows = metric (log-score, then RMSE); columns = count, distance
fig <- (panel("logs", "count", TRUE,  FALSE) | panel("logs", "distance", TRUE,  FALSE)) /
       (panel("rmse", "count", FALSE, TRUE ) | panel("rmse", "distance", FALSE, TRUE )) +
       plot_layout(guides = "collect")

print(fig)
