## =====================================================================
## Data + triangulation figure for the three temperature regimes
## Run after fit_models.R.
## =====================================================================

library(fmesher)
library(fields)   
library(ggplot2)
library(viridis)
library(maps)     
library(patchwork)


DEG_KM <- 111.32
project <- function(lon, lat, lon0, lat0)
  cbind(x = (lon - lon0) * DEG_KM * cos(lat0 * pi/180), y = (lat - lat0) * DEG_KM)

## Reference (lon0, lat0) for each regime, re-derived from the raw station
## locations exactly as prep_data() in fit_models.R, so the projected outline
## lines up with the saved km coordinates.
ref_lonlat <- function(tag) {
  if (tag == "COLORADO") {
    data(COmonthlyMet)
    ## use (min + max)/2; fields' CO.tmean is mislabelled (= CO.tmin). The
    ## plotted temperatures come from the saved fit (e$data$y); this mask
    ## only sets the projection centre.
    tmean <- (CO.tmin.MAM.climate + CO.tmax.MAM.climate) / 2
    ok  <- !is.na(tmean) & !is.na(CO.elev) & !is.na(CO.loc[,1]) & !is.na(CO.loc[,2])
    loc <- CO.loc[ok, , drop = FALSE]
  } else {
    no  <- readRDS(file.path("data", "norway_tavg.rds")); no <- no[no$lat < 65, ]
    yv  <- if (tag == "NORWAYDJF") no$tavg_djf else no$tavg_ann
    loc <- cbind(no$lon, no$lat)
    elev <- no$elev / 1000
    ok  <- complete.cases(yv, loc, elev)
    loc <- loc[ok, ]
  }
  c(lon0 = mean(loc[,1]), lat0 = mean(loc[,2]))
}

## outline (NA-separated polygons) projected to the panel's km grid and clipped
outline_df <- function(map_data, ref, xlim, ylim, pad = 200) {
  xy <- project(map_data$x, map_data$y, ref["lon0"], ref["lat0"])
  ok <- xy[,1] >= xlim[1]-pad & xy[,1] <= xlim[2]+pad & xy[,2] >= ylim[1]-pad & xy[,2] <= ylim[2]+pad
  xy[!ok, ] <- NA
  data.frame(x = xy[,1], y = xy[,2])
}

## one panel: mesh + outline + stations coloured by temperature
make_panel <- function(tag, title, map_data, clim) {
  e <- new.env(); load(sprintf("results/fits_%s.RData", tag), envir = e)
  mesh <- fm_as_mesh_2d(e$p$mesh)
  dat  <- data.frame(x1 = e$p$loc_km[,1], x2 = e$p$loc_km[,2], y = e$data$y)
  xlim <- range(mesh$loc[,1])
  ylim <- range(mesh$loc[,2])
  odf  <- outline_df(map_data, ref_lonlat(tag), xlim, ylim)
  ggplot() +
    geom_fm(data = mesh) +
    geom_path(data = odf, aes(x = x, y = y), colour = "grey25", linewidth = 0.4, na.rm = FALSE) +
    geom_point(data = dat, aes(x = x1, y = x2, colour = y), size = 2) +
    scale_colour_viridis(name = expression("temp " * degree * "C"), limits = clim) +
    coord_sf(xlim = xlim, ylim = ylim, expand = TRUE) +
    labs(title = title, x = "x (km)", y = "y (km)") +
    guides(colour = guide_colourbar(barheight = unit(4.2, "cm"), barwidth = unit(0.55, "cm"),
                                    ticks.colour = "black", frame.colour = "black")) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 11, hjust = 0.5), legend.position = "right",
          legend.title = element_text(size = 10), legend.text = element_text(size = 9),
          panel.grid.minor = element_blank())
}

m_co <- map("state", region = c("colorado","wyoming","utah","new mexico",
                                "arizona","kansas","nebraska","oklahoma"), plot = FALSE)
m_no <- map("world", region = c("Norway","Sweden","Finland","Denmark",
                                "UK","Ireland","Russia"), plot = FALSE)

## shared colour scale across the three panels
yr <- range(vapply(c("COLORADO","NORWAYANNUAL","NORWAYDJF"), function(tg) {
  e <- new.env(); load(sprintf("results/fits_%s.RData", tg), envir = e); range(e$data$y)
}, numeric(2)))

p_co  <- make_panel("COLORADO",     "Colorado MAM",     m_co, yr) + theme(axis.title.x = element_blank())
p_ann <- make_panel("NORWAYANNUAL", "S. Norway annual", m_no, yr) + theme(axis.title.y = element_blank())
p_djf <- make_panel("NORWAYDJF",    "S. Norway DJF",    m_no, yr) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

fig <- p_co + p_ann + p_djf + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "right")

print(fig)
