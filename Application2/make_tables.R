## =====================================================================
## Build Tables 4 and 5
## Run after fit_models.R.
## =====================================================================

library(Matrix)
library(fmesher)
library(INLA)
library(rSPDE)
source("../common/standardize.R")

load("results/be_no2_fits.RData")   # fits, tab, y, n, loc_km, mesh, Xtraf, Xpt, Xpop, NU, M_LGOCV

## ---- CV-score uncertainty (Table 4's last column): paired (Diebold--Mariano)
## difference in log-score vs the best model, with a spatial block-bootstrap SE.
perobs <- lapply(fits, function(f){ g <- inla.group.cv(f, num.level.sets=M_LGOCV)
  if (is.null(g$cv)) return(NULL); -log(g$cv) })                 # per-observation log-score
perobs <- perobs[!vapply(perobs, is.null, logical(1))]
refname <- names(which.min(vapply(perobs, function(p) mean(p, na.rm=TRUE), numeric(1))))  # best (lowest) log-score
reflogs <- perobs[[refname]]
## coarse spatial blocks (12x12 quantile grid) for the block bootstrap
nb <- 12L
bx <- cut(loc_km[,1], quantile(loc_km[,1], seq(0,1,length.out=nb+1)), include.lowest=TRUE, labels=FALSE)
by <- cut(loc_km[,2], quantile(loc_km[,2], seq(0,1,length.out=nb+1)), include.lowest=TRUE, labels=FALSE)
blk <- paste(bx, by); byblk <- split(seq_len(n), blk); ublk <- names(byblk); B <- 2000L
bootse <- function(d){ 
  ms <- numeric(B)
  for (b in seq_len(B)) { 
    set.seed(b)
    s <- sample(ublk, length(ublk), replace=TRUE)
    ms[b] <- mean(d[unlist(byblk[s], use.names=FALSE)], na.rm=TRUE) 
  }
  sd(ms) 
}
cvse <- do.call(rbind, lapply(names(perobs), function(nm){ dlog <- perobs[[nm]] - reflogs
  data.frame(model=nm, ref=refname, dlogs = mean(dlog, na.rm=TRUE), dlogs_se = bootse(dlog),
             row.names=NULL) }))
cvse$dlogs_z <- cvse$dlogs / cvse$dlogs_se  
write.csv(cvse, "results/be_no2_inla_cvse.csv", row.names=FALSE)

## Table 4: model comparison
mi <- match(tab$model, cvse$model)
t4 <- data.frame(model = tab$model, preset_candidates = tab$preset_candidates,
                 n_param = tab$n_param, n_cov = tab$n_cov,
                 CV_logscore = tab$CV_logscore, CV_rmse = tab$CV_rmse,
                 dlogs = cvse$dlogs[mi], dlogs_se = cvse$dlogs_se[mi])
cat(sprintf("\n=== Table 4: model comparison (dlogs = Delta log-score vs best model, %s) ===\n", refname))
print(t4, digits = 4, row.names = FALSE)
write.csv(tab, "results/be_no2_inla_compare.csv", row.names=FALSE)

## three-field: the three per-source ranges 
f3f <- fits[["three-field"]]
hp <- f3f$summary.hyperpar
f3_rng <- function(nm){ i <- grep(paste0("^Theta2 for ", nm, "$"), rownames(hp))
  if (length(i)) c(mean = sqrt(8*NU)/exp(hp[i,"mean"]),
                   lo   = sqrt(8*NU)/exp(hp[i,"0.975quant"]),
                   hi   = sqrt(8*NU)/exp(hp[i,"0.025quant"])) else c(mean=NA,lo=NA,hi=NA) }
r3f <- rbind(traffic=f3_rng("field_t"), industry=f3_rng("field_i"), population=f3_rng("field_p"))
cat("\n=== three-field: per-source ranges kappa^-1 (km) [mean, 95% CrI] ===\n"); print(round(r3f,2))
write.csv(data.frame(source=rownames(r3f), r3f, row.names=NULL), "results/be_no2_threefield_ranges.csv", row.names=FALSE)

## industry-split two-field: the two ranges 
fsp <- fits[["split(ind|loc)"]]
hp <- fsp$summary.hyperpar
sp_rng <- function(nm){ i <- grep(paste0("^Theta2 for ", nm, "$"), rownames(hp))
  if (length(i)) c(mean=sqrt(8*NU)/exp(hp[i,"mean"]), lo=sqrt(8*NU)/exp(hp[i,"0.975quant"]), hi=sqrt(8*NU)/exp(hp[i,"0.025quant"])) else c(mean=NA,lo=NA,hi=NA) }
spr <- rbind(`local (traffic+pop)`=sp_rng("field_loc"), industry=sp_rng("field_ind"))
cat("\n=== split two-field: ranges kappa^-1 (km) [mean, 95% CrI] ===\n"); print(round(spr,2))
write.csv(data.frame(field=rownames(spr), spr, row.names=NULL), "results/be_no2_split_ranges.csv", row.names=FALSE)

## credible intervals for the hybrid parameters 
fh <- fits[["hybrid(bg+op)"]]
hp <- fh$summary.hyperpar; tr <- grep("^Theta[0-9]+ for field$", rownames(hp))
kap <- exp(hp[tr[2],"mean"]); kci <- exp(as.numeric(hp[tr[2],c("0.975quant","0.025quant")]))  # note: invert for range
fx <- fh$summary.fixed[, c("mean","0.025quant","0.975quant")]
pr <- rbind(
  setNames(data.frame(param=rownames(fx), fx), c("param","mean","lo","hi")),
  data.frame(param=c("beta_fc traffic","beta_fc industry","beta_fc population"),
             mean=hp[tr[3:5],"mean"], lo=hp[tr[3:5],"0.025quant"], hi=hp[tr[3:5],"0.975quant"]),
  data.frame(param="range kappa^-1 (km)", mean=sqrt(8*NU)/kap, lo=sqrt(8*NU)/kci[1], hi=sqrt(8*NU)/kci[2]))
cat("\n=== hybrid: posterior mean + 95% CrI ===\n"); print(pr, digits=4, row.names=FALSE)
write.csv(pr, "results/be_no2_inla_params.csv", row.names=FALSE)

## scale-free decomposition (3 channels) with CrI 
bpw <- function(nm) if(nm %in% rownames(fh$summary.fixed)) as.numeric(fh$summary.fixed[nm,c("mean","0.025quant","0.975quant")]) else rep(NA_real_,3)
chans <- list(Traffic=list(X=Xtraf, j=tr[3], pw="own_traffic"),
              `Industry (point)`=list(X=Xpt, j=tr[4], pw="own_ind"),
              Population=list(X=Xpop, j=tr[5], pw="own_pop"))
drows <- lapply(names(chans), function(nm){ 
  ch<-chans[[nm]]
  bfc<-hp[ch$j,"mean"]
  bfc_ci<-as.numeric(hp[ch$j,c("0.025quant","0.975quant")])
  pw<-bpw(ch$pw)
  st<-standardize_forcing(mesh, ch$X, kappa=kap, beta_fc=bfc, beta_pw=pw[1])  
  data.frame(channel=nm, beta_pw=pw[1], beta_pw_lo=pw[2], beta_pw_hi=pw[3],
    beta_star=st$beta_star, beta_star_lo=bfc_ci[1]/kap^2, 
    beta_star_hi=bfc_ci[2]/kap^2, g_X=st$g_X, beta_eff=st$beta_eff, 
    beta_eff_lo=bfc_ci[1]/kap^2*st$g_X, 
    beta_eff_hi=bfc_ci[2]/kap^2*st$g_X,
    vartheta=st$vartheta) 
  })
dtab <- do.call(rbind, drows)

## joint channel-separation diagnostic 
tau_h <- exp(hp[tr[1],"mean"])
js <- joint_separation_forcing(mesh, list(Traffic=Xtraf, Industry=Xpt, Population=Xpop),
                               kappa=kap, tau=tau_h)  # alpha=2 (default)
jfc <- js[js$type=="forcing", ]
jkey <- c(Traffic="Traffic", `Industry (point)`="Industry", Population="Population")
dtab$vartheta_m <- jfc$vartheta[match(jkey[dtab$channel], jfc$covariate)]

cat("\n=== scale-free decomposition (hybrid; per-SD; 95% CrI) ===\n"); print(dtab, digits=3, row.names=FALSE)
write.csv(dtab, "results/be_no2_inla_decomp.csv", row.names=FALSE)
## Cameron--Martin correlations between the forcing channels (fc-fc block of
## R_p), reported in the text: the traffic/population pair is the collinear one.
Rp <- attr(js, "R_p"); dimnames(Rp) <- list(js$coef, js$coef)
cat(sprintf("\nforcing-channel CM correlations: traffic-pop %.3f | traffic-ind %.3f | pop-ind %.3f\n",
            Rp["fc_Traffic","fc_Population"], Rp["fc_Traffic","fc_Industry"], Rp["fc_Population","fc_Industry"]))
write.csv(round(Rp,4), "results/be_no2_joint_Rp.csv")