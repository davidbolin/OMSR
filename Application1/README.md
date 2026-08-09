# Temperature–altitude application 

Reproduces the temperature–altitude application across three climate regimes.

## To reproduce

1. Run fit_models.R to fit all models
2. Run make_tables.R to produce the result tables
3. Run make_lgocv.R to produce the CV figure
4. Run compare_covariance.R to produce the comparison with the covariance-based approach
5. Run make_data_figure.R to produce the data + triangulation figure for the three regimes

## Data

The data folder contains the monthly mean-temperature climate normals (1991–2020) 
at Norwegian weather stations, paired with station elevation. The data is taken 
from NOAA/NCEI **Global Historical Climatology Network – Monthly (GHCN-M) version 4**, 
element **TAVG** (monthly mean temperature), **quality-controlled unadjusted (QCU)**.

**URL:** https://www.ncei.noaa.gov/pub/data/ghcn/v4/ (file `ghcnm.tavg.latest.qcu.tar.gz`). 
Note that GHCN-M v4 QCU is updated continuously, and `ghcnm.tavg.latest.qcu.*`
always points at the newest snapshot, so the exact station set can shift
slightly over time. The committed files were filtered from the
`ghcnm.v4.0.1.20260519` snapshot.

