## =====================================================================
## Data figure for the Flanders NO2 application 
## Reads the pre-computed layers in data/no2_figure_data.rds 
## =====================================================================
library(ggplot2)
library(viridis)
library(patchwork)

fd <- readRDS("data/no2_figure_data.rds")
invisible(list2env(fd, environment())) #copy the fields in fd to the environment

cbar <- function(ord) guide_colourbar(title.position="top", title.hjust=0.5,
                                      barwidth=unit(3.0, "cm"), barheight=unit(0.26, "cm"),
                                      ticks.colour="grey40", frame.colour=NA,
                                      order=ord)

base <- function(p, ttl) p +
  geom_path(data=outline, aes(x,y), colour="grey25", linewidth=0.3,
            inherit.aes=FALSE, na.rm=FALSE) +
  coord_fixed(ratio=1, xlim=xlim, ylim=ylim, expand=TRUE) +
  labs(x="x (km)", y="y (km)", title=ttl) + theme_bw(base_size=10) +
  theme(plot.title=element_text(size=10, hjust=0.5),
        panel.grid.minor=element_blank(),
        legend.title=element_text(size=10), legend.text=element_text(size=9),
        legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"))

p_no2 <- base(ggplot() +
  geom_segment(data=medf, aes(x=x, y=y, xend=xend, yend=yend), colour="grey85",
               linewidth=0.08, inherit.aes=FALSE) +
  geom_point(data=stn, aes(x,y,colour=pmin(no2,50)), size=0.3) +
  scale_colour_viridis(name=expression(NO[2]~(mu*g~m^{-3})), guide=cbar(1)),
  "NO2 (tubes over analysis mesh)")

p_pop <- base(ggplot() +
  geom_raster(data=popdf, aes(x,y,fill=log10(pop+1))) +
  scale_fill_viridis(name=expression(log[10]*"(persons / km"^2*")"),
                     option="mako", na.value="white", guide=cbar(2)),
  "Residential population")

## merged emission sources: traffic as a red road network (lines),
## industry as cividis point sources -- distinct by geometry and hue.
p_src <- base(ggplot() +
  geom_path(data=rd, aes(x,y,group=grp,colour=log10(intensity)), linewidth=0.28) +
  scale_colour_gradient(name=expression(log[10]~veh/day),
                        low="#fcbba1", high="#67000d", guide=cbar(3)) +
  geom_point(data=src, aes(x,y,fill=nox_kt), shape=21, colour="grey20",
             stroke=0.2, size=1.8, alpha=0.9) +
  scale_fill_viridis(name=expression(NO[x]~(kt/yr)), option="cividis",
                     direction=-1, trans="log10", breaks=c(0.01,0.1,1,10),
                     labels=c("0.01","0.1","1","10"), guide=cbar(4)),
  "Emission sources (traffic + industry)")

## shared y axis: keep it only on the leftmost panel
noy <- theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
             axis.ticks.y=element_blank())
p_pop <- p_pop + noy
p_src <- p_src + noy

fig <- (p_no2 | p_pop | p_src) + plot_layout(guides="collect") &
  theme(legend.position="bottom", legend.box="horizontal",
        legend.justification="center", plot.margin=margin(2,3,2,3))

print(fig)