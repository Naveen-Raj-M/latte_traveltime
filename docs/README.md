# The spherical branch of LATTE

This branch computes first-arrival traveltimes (`x_eikonal3`), first-arrival
traveltime tomography (`x_fatt3`) and source location with optional joint
tomography (`x_tloc3`) on a spherical Earth. Inputs and outputs are in
longitude, WGS84 geodetic latitude and depth below mean sea level. Acoustic
and elastic (P and S) media are supported; computation is on CPUs with MPI
and OpenMP.

## Documents

Each step of the conversion has its own document, with the files changed,
the equations and the verification:

| Step | Document | Subject |
|---|---|---|
| 0 | [spherical_00_scope_reduction.md](spherical_00_scope_reduction.md) | What was removed: 2D, TRTT, CUDA, ray paths and other out-of-scope features |
| 1 | [spherical_01_coordinates_and_parameters.md](spherical_01_coordinates_and_parameters.md) | The spherical geometry, parameter names, input checks, offsets, adaptive range |
| 2 | [spherical_02_forward_eikonal.md](spherical_02_forward_eikonal.md) | The spherical eikonal solver |
| 3 | [spherical_03_adjoint_and_gradient.md](spherical_03_adjoint_and_gradient.md) | The adjoint solver, receiver injection, volume-weighted gradients |
| 4 | [spherical_04_processing_and_regularization.md](spherical_04_processing_and_regularization.md) | Gradient filters and model regularization with exact distances |
| 5 | [spherical_05_tloc_event_steps.md](spherical_05_tloc_event_steps.md) | TLOC event steps in meters, required starting positions |
| 6 | [spherical_06_final_verification.md](spherical_06_final_verification.md) | Clean build and the full test suite |

The tests are in [`test/spherical`](../test/spherical).

## Geometry in one paragraph

The Earth is a sphere of radius `earth_radius` (default 6371000 m). A point
(lon, lat, depth) lies at longitude lon, at the geocentric latitude that
corresponds to the WGS84 geodetic latitude lat on the ellipsoid surface, and
at radius `earth_radius - depth`. The model grid is regular in longitude,
geodetic latitude and depth. Every distance, area and volume of the grid is
the exact one of this sphere; the sphere itself is the only approximation,
apart from the window filters, which use a reference cell size.

## Units

| Quantity | Unit |
|---|---|
| Longitude, latitude (grid, geometry, `slon`, `slat`) | degrees |
| Depth (grid, geometry, `sdepth`); negative above sea level | model length unit, meters by default |
| Velocity | length unit per second |
| Every length parameter of processing and regularization, adaptive margins and tapers, step limits of `slon` and `slat` | length unit |
| `reg_tikhonov_lambda` | length unit squared |

## Minimal parameter file

```
nlon = 101
nlat = 81
ndepth = 41
dlon = 0.05
dlat = 0.05
ddepth = 1000
olon = 100
olat = 20
odepth = 0

ns = 10
file_geometry = geometry/geometry.txt
model_update = vp
file_vp = model/vp_init.bin
dir_record = data
dir_working = fatt
niter_max = 20
step_max_vp = 100
```

`nlon`, `nlat`, `ndepth`, `dlon`, `dlat`, `ddepth`, `olon` and `olat` are
required. Each axis needs at least 2 nodes, and the grid must lie between
-90 and 90 degrees latitude.

## Files

- **Models** keep LATTE's layout: a 3D array of 32-bit floats of shape
  (`ndepth`, `nlat`, `nlon`) in Fortran order, so depth varies fastest.
- **Geometry files** keep their layout. The position columns are longitude,
  geodetic latitude and depth below mean sea level.
- **TLOC starting values**: `file_slon`, `file_slat`, `file_sdepth` and
  `file_st0` are required for the parameters being updated. They hold one
  value per event, in geometry order.
- **Outputs**: models, gradients, traveltimes and TLOC event positions are
  on the same grid and in the same units as the inputs.

## Parameter names

Every Cartesian parameter that referred to x, y or z now refers to lon,
lat or depth. For example, `nx` became `nlon`, `grad_smoothx` became
`grad_smooth_lon`, and `sx` became `slon`. A parameter file with a
Cartesian name stops the program at startup and names the replacement. The
full table is in the step 1 document.

## Not supported on this branch

The program stops at startup if a parameter file asks for any of these:

- 2D, TRTT and reflections, ray paths, traveltime-gradient output.
- ML and clustering-fit source regularization.
- The dip filter and the `structure` (ANDF) model regularization.
- Wavenumber filtering along latitude.
- `adjoint_receiver = nearest`.

The build stops if asked for CUDA or the legacy solver.
