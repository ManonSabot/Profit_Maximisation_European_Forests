# TractLSM (a Tractable simplified Land Surface Model)

The **TractLSM** Compares a plant profit maximisation (ProfitMax) approach
with the standard photosynthesis - transpiration coupling (Control) approach
used in most LSMs.

The model is organised as follows:

```bash
TractLSM
├── run_homogeneous_surface.py
├── run_utils.py
├── CH2OCoupler
│   ├── ProfitMax.py
│   ├── USO.py
├── SPAC
│   ├── canatm.py
│   ├── hydraulics.py
│   ├── leaf.py
│   ├── soil.py
├── TraitCoordination
│   ├── optimise_kmax.py
└── Utils
    ├── build_final_forcings.py
    ├── built_in_plots.py
    ├── calculate_solar_geometry.py
    ├── constants_and_conversions.py
    ├── cru_climate_lat_lon.py
    ├── default_params.py
    ├── general_utils.py
    ├── met_flux_LAI_site_level.py
    ├── modis_lai_lat_lon.py
    └── weather_generator.py
```

&nbsp;

`run_homogeneous_surface.py` is where the forcings are read, the main routines
called, and the output writen. `run_utils.py` contains support functions for
these actions.

&nbsp;

The `CH2OCoupler/` is where you can find the `ProfitMax.py` approach, which is
derived/adapted from the work of
[Sperry et al. (2017)](https://doi.org/10.1111/pce.12852).
The Control flux coupling method in `USO.py` and it uses the
[Medlyn et al. (2011)](https://doi.org/10.1111/j.1365-2486.2010.02375.x) model.

&nbsp;

The model's biogeophysical routines can be found in the `SPAC/` repository,
ranging from micrometeorology (`canatm.py`) to plant hydraulics
(`hydraulics.py`).

&nbsp;

`Trait Coordination/` contains a routine which calculates the maximum hydraulic
conductance (k<sub>max</sub>) based on the long-term background climate at
site-level, and assuming coordination between the photosynthetic and hydraulic
traits. This routine is not called on an instantaneous basis, but can rather be
used to calculate the input parameter k<sub>max</sub>, needed to run the
**Profit<sub>max</sub>**, when it is not provided.

&nbsp;

All support routines (automating the format of the input files, etc.) can be
found in `Utils/`.

&nbsp;

Manon Sabot: [m.e.b.sabot@gmail.com](mailto:m.e.b.sabot@gmail.com?subject=[ProfitMax_Europe_Code]%20Source%20Han%20Sans)

