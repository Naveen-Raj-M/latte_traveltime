# Step 1: Spherical coordinates, parameter names and input checks

This step introduces the spherical-Earth geometry, renames every
parameter that referred to x, y or z, and makes the programs stop on
parameters of the Cartesian code or of removed features. It also makes
the geometry input, offsets and adaptive ranges spherical.

**Intermediate state.** The traveltime solver, the adjoint solver, the
gradients and the gradient filters are still the Cartesian ones after this
step. They treat a longitude or latitude interval in degrees as if it were
a length, so traveltimes and inversion results are not yet meaningful.
Steps 2 to 5 convert them. Use this step only to check inputs.

## Files changed

| File | What changed |
|---|---|
| `src/module_spherical.f90` | New. The spherical-Earth geometry: WGS84 latitude conversion, Earth-centred coordinates, and the exact cell volumes, arc lengths and reference spacings of the grid. |
| `src/module_parameters.f90` | Reads `earth_radius`. Checks the parameter file for old and unsupported parameters. Reads the grid, selection and adaptive-range parameters under their new names. Sets up the spherical grid. Computes offsets from Earth-centred positions. Converts adaptive-range margins from meters to degrees. Relabels printouts to lon, lat and depth. |
| `src/module_utility.f90` | The exchange-geometry file of the source-receiver exchange is labelled lon, lat and depth and written at full single precision. |
| `src/3d/module_gradient.f90` | Gradient-processing parameters and option values renamed. |
| `src/3d/module_regularization.f90` | Regularization and MPI-decomposition parameters renamed. |
| `src/3d/module_vars.f90`, `src/3d/submodule_gradient_tloc.f90`, `src/module_inversion_*.f90` | TLOC source parameters renamed from `sx`, `sy`, `sz` to `slon`, `slat`, `sdepth`. |
| `src/3d/Makefile` | Compiles `module_spherical.o` first. |
| `src/README.md` | Build notes for the spherical branch. |

## The coordinate system

The model grid is regular in three axes:

| Axis | Grid index | Unit | Increases |
|---|---|---|---|
| Longitude | `i`, 1 to `nlon` | degrees | eastward |
| Geodetic latitude, WGS84 | `j`, 1 to `nlat` | degrees | northward |
| Depth below mean sea level | `k`, 1 to `ndepth` | model length unit, by default meters | downward |

The model files keep their layout: a 3D array of shape (`ndepth`, `nlat`,
`nlon`) in Fortran order, so depth varies fastest. Geometry files keep their
columns; the first three columns are now longitude, geodetic latitude and
depth below mean sea level. Stations above sea level have negative depth.

### From a point to the sphere

The Earth is a sphere of radius `R` (`earth_radius`, default 6371000 m). A
point with longitude `lon`, geodetic latitude `lat` and depth `d` lies at

```
lambda = lon
theta  = atan2( (1 - e^2) sin(lat), cos(lat) )        geocentric latitude
r      = R - d
```

where `e^2 = f (2 - f)` and `f = 1/298.257223563` are the WGS84 values. The
latitude relation is the exact geodetic-to-geocentric relation on the WGS84
ellipsoid surface. It is applied at every depth, so a grid column of
constant (lon, lat) is a radial line of the sphere. This, and the sphere
itself, are the only approximations of the geometry.

The inverse relation, used to report positions, is

```
lat = atan2( sin(theta), (1 - e^2) cos(theta) )
```

The Earth-centred Cartesian position of the point is

```
X = r ( cos(theta) cos(lambda),  cos(theta) sin(lambda),  sin(theta) )
```

### The cells of the grid

Every grid node owns a full cell centred on it. The faces of the cells lie
half a grid interval from the nodes, including half an interval outside the
grid at the boundary nodes, as in the Cartesian code. In geodetic latitude
and depth, face `j` and face `k` lie at

```
lat_face(j) = olat + (j - 1.5) dlat,     j = 1, ..., nlat + 1
d_face(k)   = odepth + (k - 1.5) ddepth, k = 1, ..., ndepth + 1
```

so that node `j` lies between faces `j` and `j + 1`, and node `k` between
faces `k` (shallower) and `k + 1` (deeper). The face latitudes are converted
to geocentric latitudes `theta_face(j)` with the relation above, and the face
depths to radii `r_face(k) = R - d_face(k)`.

The exact volume of the cell of node (i, j, k) is the integral of
`r^2 cos(theta) dr dtheta dlambda` over the cell:

```
V(j, k) = (r_face(k)^3 - r_face(k+1)^3)/3 * (sin(theta_face(j+1)) - sin(theta_face(j))) * dlambda
```

with `dlambda = dlon * pi/180`. It does not depend on `i`.

The arc lengths between neighbouring nodes, which later steps use as grid
spacings, are

```
along a parallel:  h_lon(j, k)       = r_k cos(theta_j) dlambda
along a meridian:  h_lat(j1, j2, k)  = r_k | theta_j2 - theta_j1 |
along a radius:    h_depth           = ddepth
```

Because the grid is regular in geodetic latitude, the geocentric spacing
`theta_{j+1} - theta_j` varies slightly from row to row; `h_lat` uses the
actual geocentric latitudes of the two rows.

All these quantities are computed once, in double precision, by
`init_spherical_grid` when the grid is set up, and stored in 1D arrays of
longitude, latitude and depth. A node's Earth-centred position is the
product of these arrays, so no 3D coordinate arrays are allocated.

### Reference spacing and length of a degree

Two derived quantities are used where a physical length must become a
number of cells or degrees.

The reference spacing is the size of one grid interval at the centre of the
grid in longitude and latitude, at the top level:

```
h_ref_lon   = r_1 cos(theta_c) dlambda
h_ref_lat   = r_1 (dtheta/dlat)(lat_c) dlat_rad
h_ref_depth = ddepth
```

where `lat_c = olat + (nlat - 1) dlat/2`, `theta_c` is its geocentric
latitude, and

```
dtheta/dlat = (1 - e^2) / ( cos^2(lat) + (1 - e^2)^2 sin^2(lat) )
```

is the exact derivative of the geocentric latitude with respect to the
geodetic latitude. Step 4 uses it for the window filters.

The smallest length of one degree anywhere in the grid is

```
m_lon = r_min min( cos(theta_south), cos(theta_north) ) pi/180
m_lat = r_min min over the latitude range of (dtheta/dlat) pi/180
```

where `r_min` is the radius of the deepest level. `dtheta/dlat` is smallest
at the equator and grows toward the poles, so its minimum is at the equator
if the grid spans it, and otherwise at the row nearest to it.

### Checks on the grid

The program stops with a message unless:

- `nlon`, `nlat` and `ndepth` are each at least 2.
- `dlon`, `dlat` and `ddepth` are positive.
- `nlon * dlon` is at most 360 degrees, so cells do not overlap.
- The latitude cells, including the half cells beyond the first and last
  rows, lie strictly between -90 and 90 degrees.
- The deepest cell lies above the centre of the Earth.

Longitudes are used as given. The grid and the geometry must use the same
longitude convention, for example both 0 to 360 or both -180 to 180.
Points outside the grid are excluded, as in the Cartesian code.

## Offsets

The offset of a receiver, used by `offset_min` and `offset_max`, is the
straight-line distance between Earth-centred positions, from the receiver
to the centroid of the source points of its shot:

```
X_c    = (1/n_s) sum over source points s of X(s)
offset = | X(receiver) - X_c |
```

This keeps the meaning of the Cartesian offset, a 3D distance that includes
the depth difference. The centroid is averaged in Earth-centred coordinates,
not in degrees.

## Adaptive range

`adp_extra_lon` and `adp_extra_lat` are now lengths in meters. They are
converted to degrees with the smallest length of a degree anywhere in the
grid, so the margin is never shorter than asked for:

```
extra_lon_deg = adp_extra_lon / m_lon
extra_lat_deg = adp_extra_lat / m_lat
```

`adp_extra_depth` stays a length. The range is then rounded outward to grid
nodes, as before.

## Parameter names

Every parameter name that referred to x, y or z now refers to lon, lat or
depth. The Cartesian name of any of these in a parameter file stops the
program with a message naming the replacement.

| Cartesian | Spherical | Unit |
|---|---|---|
| `nx`, `ny`, `nz` | `nlon`, `nlat`, `ndepth` | |
| `dx`, `dy`, `dz` | `dlon`, `dlat`, `ddepth` | degrees, degrees, length |
| `ox`, `oy`, `oz` | `olon`, `olat`, `odepth` | degrees, degrees, length |
| `nnx`, `nny`, `nnz` | `nnlon`, `nnlat`, `nndepth` | |
| `ddx`, `ddy`, `ddz` | `ddlon`, `ddlat`, `dddepth` | degrees, degrees, length |
| `oox`, `ooy`, `ooz` | `oolon`, `oolat`, `oodepth` | degrees, degrees, length |
| `xmin`, `xmax`, `ymin`, ... | `lonmin`, `lonmax`, `latmin`, `latmax`, `depthmin`, `depthmax` | degrees, length |
| `sxmin`, `sxmax`, ... | `slon_min`, `slon_max`, `slat_min`, `slat_max`, `sdepth_min`, `sdepth_max` | degrees, length |
| `rxmin`, `rxmax`, ... | `rlon_min`, `rlon_max`, `rlat_min`, `rlat_max`, `rdepth_min`, `rdepth_max` | degrees, length |
| `yn_adpx`, `yn_adpy`, `yn_adpz` | `yn_adp_lon`, `yn_adp_lat`, `yn_adp_depth` | |
| `adp_extrax`, ... | `adp_extra_lon`, `adp_extra_lat`, `adp_extra_depth` | length |
| `adp_taperx`, ... | `adp_taper_lon`, `adp_taper_lat`, `adp_taper_depth` | length |
| `vpvsratio_smoothx`, ... | `vpvsratio_smooth_lon`, `_lat`, `_depth` | length |
| `rankx`, `ranky`, `rankz` | `rank_lon`, `rank_lat`, `rank_depth` | |
| `reg_smoothx`, `reg_smoothx_<model>` | `reg_smooth_lon`, `reg_smooth_lon_<model>`, and lat, depth | length |
| `reg_similarity_smoothx`, ... | `reg_similarity_smooth_lon`, `_lat`, `_depth` | length |
| `<grad>_smoothx` | `<grad>_smooth_lon`, `_lat`, `_depth` | length |
| `<grad>_taperx` | `<grad>_taper_lon`, `_lat`, `_depth` | length |
| `<grad>_movingbalx` | `<grad>_movingbal_lon`, `_lat`, `_depth` | length |
| `<grad>_medianfiltx` | `<grad>_medianfilt_lon`, `_lat`, `_depth` | length |
| `<grad>_rmsbalx` | `<grad>_rmsbal_lon`, `_lat`, `_depth` | length |
| `<grad>_andf_smoothx` | `<grad>_andf_smooth_lon`, `_lat`, `_depth` | length |
| `<grad>_andf_rankx` | `<grad>_andf_rank_lon`, `_lat`, `_depth` | |
| `shot_<grad>_wavenumx`, `..._amps` | `shot_<grad>_wavenum_lon`, `shot_<grad>_wavenum_lon_amps`, and lat, depth | |
| `shot_<grad>_conemutez` | `shot_<grad>_conemute_depth` | length |
| model names `sx`, `sy`, `sz` | `slon`, `slat`, `sdepth` | degrees, degrees, length |
| `file_sx`, `min_sx`, `step_max_sx`, ... | `file_slon`, `min_slon`, `step_max_slon`, ... | |
| options `rmsbalx`, `rmsbaly`, `rmsbalxy`, `rmsbalz` | `rmsbal_lon`, `rmsbal_lat`, `rmsbal_lonlat`, `rmsbal_depth` | |

`<grad>` is `grad`, `grad_<model>`, `shot_grad` or `shot_grad_<model>`, as
in the Cartesian code.

New parameter:

| Name | Meaning | Default |
|---|---|---|
| `earth_radius` | Radius of the spherical Earth, in the model length unit | 6371000 |

Now required: `nlat`, `dlat`, `olon` and `olat`. They used to default to 1,
1, 0 and 0. A default origin at the equator and prime meridian would
silently place a model in the wrong place.

The units in the table are the units the parameters have once steps 2 to 5
are done. In this step, the gradient filters still divide these lengths by
`dlon` and `dlat` in degrees.

## Unsupported parameters

The program stops, naming each offending line, if the parameter file has:

- `source_regularization_method`, any `clustering_fit_*`, any `reg_ml_*`.
- Any `raypath_*`, `yn_save_traveltime_gradients`, `dir_gradient`.
- Any reflector or TRTT parameter: `nrefl`, `incident_wave`,
  `offset_min_refl`, `offset_max_refl`, `offset_max_refl_depth`,
  `step_max_refl`, `yn_refl_image`, any `refl_*`, `reflector_*`, or any name
  with a part `refl`, such as `file_refl`.
- Any `reg_andf_*`, the parameters of the removed `structure` regularization.
- Any name containing `dipfilt`.
- The model name `refl` in `model_name`, `model_update` or `model_aux`.
- The option `dipfilt` in any `process_*` list.
- The option `structure` in `model_regularization_method`.

The parameter file is read with the rules of the FLIT readpar routines: one
`name = value` per line, names are case-insensitive, a line with an empty
value is not a parameter, and reading stops at a line `exit`. Every problem
is reported before the program stops, so one run lists them all.

## Mathematical changes

The new formulas are the geometry above: the latitude conversion, the
Earth-centred position, cell volumes, arc lengths, the reference spacing,
the length of a degree, the offset and the adaptive margin. The solvers,
gradients and inversion are mathematically unchanged in this step.

## Execution changes

- Startup prints the Earth radius and the grid in lon, lat and depth.
- A parameter file with Cartesian names or unsupported features stops the
  program at startup, listing every problem.
- `init_spherical_grid` allocates the 1D geometry arrays once.
- Offsets and adaptive ranges are spherical.
- `exchange_geometry.txt`, written by `x_tloc3`, has columns id, lon, lat,
  depth in exponential format with 9 decimals.

## Verification

The scripts of these tests are in `test/spherical/step1`; `test/spherical/README.md` explains how to run them.

All tests ran in a scratch directory.

A test program checked `module_spherical`:

| Check | Result |
|---|---|
| Geocentric latitude of geodetic 45 degrees | 44.80757678 degrees, the WGS84 value |
| Round trip geodetic, geocentric, geodetic at 45 and 89.9 degrees | error 0 |
| Earth-centred position of (90 E, 0, 1000 m) | (0, 6370000, 0) m |
| Sum of all cell volumes of a 101 x 81 x 41 grid vs the exact volume of the region they tile | relative difference 1.3e-15 |
| Arc length vs straight chord between neighbouring nodes, 10 km apart | arc longer by about 1 mm, as expected |

`x_eikonal3` was run with a parameter file holding 26 Cartesian or
unsupported entries. It reported all 26, each with its replacement or
reason, and stopped. With the new names, it loaded the grid and geometry
and ran. With `offset_max = 20000`, it kept only the receiver 10 km above
the source and dropped the two receivers more than 100 km away. With 50 km
adaptive margins, the shot range was 1 degree in longitude and latitude
around the stations, as the conversion predicts. With `nlat = 1`, and with a
grid reaching past 90 degrees latitude, it stopped with the grid messages
above.
