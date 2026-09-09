# Flanders NO2 application

To reproduce

1. Run `fit_models.R` to fit all eight models (writes `results/be_no2_fits.RData`).
2. Run `make_tables.R` to reproduce the two tables.
3. Run `make_cvcurve.R` to reproduce the cross-validation figures.
4. Run `data_figure.R` to make the data figure.

Requires `rSPDE` (>= 2.6.0) and `R-INLA`.

## Data

- `data/no2_data.rds` — the analysis mesh, the NO2 response, the node-level forcing
  covariates (traffic / industry / population) and the LUR station covariates.
- `data/no2_figure_data.rds` — the lightweight map layers for the data figure.

These are produced from the following public sources: 

- **NO2 tubes** — CurieuzeNeuzen Flanders 2018 open dataset
  (`https://github.com/plovercode/course-python-geospatial`, `CN_Flanders_open_dataset.csv`).
- **Population** — WorldPop 2020 ~100 m rasters, native for BE and NL with
  FR/DE fill (`https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/<ISO3>/`).
- **Road traffic** — OpenStreetMap major roads over the mesh region.
- **Industrial NOx point sources** — EEA E-PRTR v15 air releases
  (`F1_4_Air_Releases_Facilities.csv`,
  `https://sdi.eea.europa.eu/catalogue/srv/api/records/3461f4ab-a3ee-4af2-bc11-95e651a8d0ba`);
- **Land use / distance-to-coast** — CORINE Land Cover 2018 and coastline via
  `rnaturalearth`, used for the LUR baseline covariates.
