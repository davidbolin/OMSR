## =====================================================================
## Data figure for the Flanders NO2 application 
## =====================================================================
library(fmesher)
library(ggplot2)
library(viridis)
library(maps)
library(patchwork) 
library(sf)
library(terra)

cn <- read.csv("data/curieuzeneuzen_2018.csv")
cn <- cn[is.finite(cn$lon)&is.finite(cn$lat)&is.finite(cn$no2), ]
lon0 <- mean(cn$lon)
lat0 <- mean(cn$lat)
deg_km <- 111.32
to_km <- function(lo,la) cbind(x=(lo-lon0)*deg_km*cos(lat0*pi/180), y=(la-lat0)*deg_km)
pts <- read.csv("data/be_point_sources_2018.csv")
loc_km <- to_km(cn$lon,cn$lat)
pts_km <- to_km(pts$lon,pts$lat)
dom <- rbind(loc_km, pts_km)
mesh <- fm_as_mesh_2d(fm_mesh_2d(loc.domain=dom, 
                                 boundary=fm_nonconvex_hull(dom,convex=10,concave=15),
                                 max.edge=c(2,15), cutoff=1.0, offset=c(5,30)))
xlim <- range(mesh$loc[,1]); ylim <- range(mesh$loc[,2])

## mesh edges as plain segments in km coords 
tv <- mesh$graph$tv
mseg <- rbind(tv[,c(1,2)], tv[,c(2,3)], tv[,c(3,1)])
medf <- data.frame(x=mesh$loc[mseg[,1],1], y=mesh$loc[mseg[,1],2],
                   xend=mesh$loc[mseg[,2],1], yend=mesh$loc[mseg[,2],2])

## regional outline, projected to km
mp <- map("world", region=c("Belgium","Netherlands","France","Germany","Luxembourg"), plot=FALSE)
xy <- to_km(mp$x, mp$y)
ok <- xy[,1]>=xlim[1]-100 & xy[,1]<=xlim[2]+100 & xy[,2]>=ylim[1]-100 & xy[,2]<=ylim[2]+100
xy[!ok,] <- NA
outline <- data.frame(x=xy[,1], y=xy[,2])

## road network points (subsampled) for the traffic panel
roads <- st_zm(st_read("data/be_traffic_roads.gpkg", quiet=TRUE))
rp <- st_cast(roads,"POINT",warn=FALSE)
rxy <- st_coordinates(st_transform(rp,4326))
rkm <- to_km(rxy[,1],rxy[,2])
rd <- data.frame(x=rkm[,1], y=rkm[,2], intensity=rp$intensity)
rd <- rd[order(rd$intensity), ]
if (nrow(rd) > 12000) rd <- rd[round(seq(1,nrow(rd),length.out=12000)), ]

## WorldPop population -> km grid for the population panel
popr <- rast("data/be_population.tif")
popr <- crop(popr, ext(min(cn$lon)-0.15, max(cn$lon)+0.15, min(cn$lat)-0.15, max(cn$lat)+0.15))
popr <- aggregate(popr, fact=10, fun="sum", na.rm=TRUE)   # ~1 km cells
pxy <- xyFromCell(popr, 1:ncell(popr))
pv <- terra::values(popr)[,1]
pkm <- to_km(pxy[,1],pxy[,2])
popdf <- data.frame(x=pkm[,1], y=pkm[,2], pop=pv)   

stn <- data.frame(x=loc_km[,1], y=loc_km[,2], no2=cn$no2)
stn <- stn[order(stn$no2), ]
src <- data.frame(x=pts_km[,1], y=pts_km[,2], nox_kt=pts$nox_kg/1e6)
src <- src[order(src$nox_kt), ]  # small first, large plotted on top

base <- function(p, ttl) p +
  geom_path(data=outline, aes(x,y), colour="grey25", linewidth=0.3, 
            inherit.aes=FALSE, na.rm=FALSE) +
  coord_fixed(ratio=1, xlim=xlim, ylim=ylim, expand=TRUE) +
  labs(x="x (km)", y="y (km)", title=ttl) + theme_bw(base_size=10) +
  theme(plot.title=element_text(size=10, hjust=0.5), 
        panel.grid.minor=element_blank(),
        legend.position="bottom", legend.title=element_text(size=8), 
        legend.text=element_text(size=7),
        legend.key.height=unit(0.25,"cm"), legend.key.width=unit(1.0,"cm"),
        legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"))

p_no2 <- base(ggplot() +
  geom_segment(data=medf, aes(x=x, y=y, xend=xend, yend=yend), colour="grey85", 
               linewidth=0.08, inherit.aes=FALSE) +
  geom_point(data=stn, aes(x,y,colour=pmin(no2,50)), size=0.3) +
  scale_colour_viridis(name=expression(NO[2]~(mu*g~m^{-3}))), 
  "NO2 (tubes over analysis mesh)")

p_pop <- base(ggplot() +
  geom_raster(data=popdf, aes(x,y,fill=log10(pop+1))) +
  scale_fill_viridis(name=expression(log[10]*"(persons / km"^2*")"), 
                     option="mako", na.value="white"),
  "Residential population")

p_traf <- base(ggplot() +
  geom_point(data=rd, aes(x,y,colour=log10(intensity)), size=0.25) +
  scale_colour_viridis(name=expression(log[10]~veh/day), 
                       option="rocket", direction=-1), "Road traffic")

p_ind <- base(ggplot() +
  geom_point(data=src, aes(x,y,fill=nox_kt), shape=21, colour="grey20", 
             stroke=0.2, size=1.8, alpha=0.9) +
  scale_fill_viridis(name=expression(NO[x]~(kt/yr)), option="plasma", 
                     direction=-1, trans="log10",
                     breaks=c(0.2,0.5,1,2,5), labels=c("0.2","0.5","1","2","5")), 
  "Industrial emissions (E-PRTR)")

fig <- (p_no2 | p_pop) / (p_traf | p_ind) & theme(plot.margin = margin(2, 3, 2, 3))
print(fig)
