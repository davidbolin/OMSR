# Flanders NO2 application 

Reproduces the NO2 application.

## To reproduce

1. Run fit_models.R to fit all models.
2. Run make_tables.R to reproduce the two tables. 
3. Run make_cvcurve.R to reproduce the CV plot. 
4. Run data_figure.R to make the data figure. 

## Data 

The data folder contains all required data. The CurieuzeNeuzen Flanders open dataset was 
downloaded from 

 https://github.com/plovercode/course-python-geospatial/notebooks/data/CN_Flanders_open_dataset.csv

The population covariates were computed from from WorldPop BEL 2020, 

 https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/BEL/bel_ppp_2020.tif
 
 Land use covariates are taken from 
 
  https://raw.githubusercontent.com/plovercode/course-python-geospatial/main/notebooks/data/CLC2018_V2020_20u1_flanders.tif
 
Distance to coast was computed using the `rnaturalearth` package and traffic intensity was computed using the 
`osmdata` package. Finally, industrial emissions were taken from 

https://sdi.eea.europa.eu/catalogue/srv/api/records/3461f4ab-a3ee-4af2-bc11-95e651a8d0ba


