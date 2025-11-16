JADE: Julia-code for Astrophysical Dust Evolution

JADE is a Julia module for modeling dust evolution in protoplanetary disks, including coagulation, fragmentation, sublimation, and advection. This code was developed for the study presented in [arXiv:2509.11036](https://arxiv.org/abs/2509.11036).

## Features

- **Dust Coagulation & Fragmentation**: Solves the Smoluchowski coagulation equation with physically-motivated fragmentation models based on turbulent velocities and material properties
- **Turbulent Velocity Calculations**: Computes relative velocities between dust grains accounting for Brownian motion, turbulence, and Epstein/Stokes drag regimes
- **Sublimation**: Includes temperature-dependent dust sublimation with user-provaporization data
- **Radiative Transfer** (only simplest for this moment): Solves ray tracing for stellar irradiation to compute dust temperatures
- **Advection**: Explicit advection scheme for gas and dust in 2D spherical coordinates (R, θ)
- **Parallelization**: Multi-threaded implementation for performance

## Requirements

- Julia ≥1.10
- Required packages: `OrdinaryDiffEq`, `SpecialFunctions`, `DelimitedFiles`, `Interpolations`

## Quick Start

```julia
using jade

# Initialize disk model
ppd = jade.ppd_t()
f_init!(ppd,
    n_mbin=20,      # Number of dust size bins
    r_min_au=0.1,   # Inner radius [AU]
    r_max_au=2.0,   # Outer radius [AU]
    n_r=32,         # Radial grid cells
    th_min=1.35,    # Min polar angle [rad]
    th_max=pi/2,    # Max polar angle [rad]
    n_th=30,        # Theta grid cells
    order=1         # Ghost zone order
)

# Run simulation
jade.f_major_loop!(ppd,
    tlim=1e4 * jade.yr,  # Simulation time [s]
    dt0 =1e2 * jade.yr   # Initial timestep
)
```

## Key Physics

- **Dust size distribution**: Logarithmic bins from 0.1 μm to 10 cm (default)
- **Fragmentation model**: Based on energy threshold with power-law fragment distribution
- **Turbulence**: α-prescription with Reynolds number dependence
- **Sublimation**: Reads comma-separated table (T, μ, η) and interpolates vapor pressure

## Output Functions

```julia
# Extract gas density, temperature, dust surface density
r_edges, th_edges, rho_g = f_get_heatmap_stat!(ppd, f_rho_g)
r_edges, th_edges, T     = f_get_heatmap_stat!(ppd, f_T)
r_edges, th_edges, sigma_dust = f_get_heatmap_stat!(ppd, f_rho_m)

# Get size-integrated dust size distribution
r_centers, th_centers, rho_m_3d = f_extract_rho_m!(ppd)
```

## Citation

If you use this code, please cite the paper associated with arXiv:2509.11036
