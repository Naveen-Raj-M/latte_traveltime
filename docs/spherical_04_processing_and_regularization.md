# Step 4: Gradient processing and model regularization

This step converts the gradient filters, the adaptive-range taper and the
model regularizers to the spherical Earth. Every length parameter is now in
the model length unit (meters).

- **Exact:** Gaussian smoothing, tapers, the cone mute, wavenumber
  filtering along longitude and depth, and Tikhonov regularization use
  exact distances on the spherical grid.
- **Approximate, as agreed:** the window filters, which work on a fixed
  number of cells, convert lengths with one reference cell size.
- **Index space:** the Laplace filter and TGpV are unchanged.

**Intermediate state.** TLOC event positions still step in degrees with
degree limits; step 5 converts them. Everything else in `x_fatt3` and
`x_tloc3` is now spherical.

## Files changed

| File | What changed |
|---|---|
| `src/module_spherical_filter.f90` | New. Exact Gaussian smoothing along grid lines, exact taper, and wavenumber filtering along longitude lines with each line's spacing and along depth. |
| `src/3d/module_gradient.f90` | Per-shot processing: exact smoothing, taper and cone mute, wavenumber filtering along longitude and depth, window filters with the reference cell size, adaptive-range taper with exact distances. Whole-model processing: exact smoothing and taper, window filters and ANDF with the reference cell size. Defaults in meters. |
| `src/3d/module_regularization.f90` | Tikhonov with the exact spherical Laplacian and `reg_tikhonov_lambda` in m^2. Smoothing regularization and the Vp/Vs-similarity Gaussian are exact. Defaults in meters. |
| `src/3d/module_vars.f90` | The smoothing of the Vp/Vs ratio (`vpvsratio_smooth_*`) is exact. |
| `src/module_parameters.f90` | `shot_<grad>_wavenum_lat` and `shot_<grad>_wavenum_lat_amps` stop the program. |
| `src/3d/Makefile` | Compiles `module_spherical_filter.o`. |

## Distances along grid lines

The exact filters measure distance along the grid line itself, as arc
length. For a block of the grid whose first row is the global row `jbeg`:

```
along a parallel (row j, level k):  s_i = (i - 1) r_k cos(theta_j) dlambda
along a meridian (level k):         s_j = r_k (theta_j - theta_jbeg)
along a radius:                     s_k = (k - 1) ddepth
```

`theta` is the geocentric latitude of a row (step 1). Along a meridian the
spacing changes slightly from row to row, because the grid is regular in
geodetic latitude; these positions account for that exactly.

## Gaussian smoothing

All Gaussian smoothing goes through one routine:

- `smooth` in per-shot and whole-model gradient processing;
- `smooth` in model regularization, and its reciprocal form for velocities;
- the Gaussian of `vpvs_similarity`;
- `vpvsratio_smooth_*`.

The lengths are standard deviations, as before. Each axis is filtered in
turn, longitude then latitude then depth, by direct convolution along every
grid line:

```
f_smoothed(s_m) = sum_n G(s_n - s_m) f(s_n) / sum_n G(s_n - s_m)
G(d) = exp( -d^2 / (2 sigma^2) ),   |d| <= 6 sigma
```

Beyond the ends, the line is extended by repeated mirror reflection about
its end samples, as FLIT's Gaussian filter does. The 6-sigma truncation is
also FLIT's. Because the weights are renormalized, a constant field stays
constant.

The east-west cell width shrinks toward the poles and with depth, so a
fixed length is fewer cells at low latitude and more at high latitude. The
smoothing keeps the same physical width everywhere.

## Taper

The tapers are:

- `taper` in per-shot and whole-model processing (`<grad>_taper_lon`,
  `_lat`, `_depth`, each with one value or two values for the two sides);
- the taper of the sides of an adaptive shot range (`adp_taper_lon`,
  `_lat`, `_depth`).

A node at distance `d` from a side, measured along its grid line from the
edge node, is multiplied by

```
B(x) = 0.42 - 0.5 cos(pi x) + 0.08 cos(2 pi x),   x = min(d / L, 1)
```

the Blackman ramp of FLIT's taper, where `L` is the taper length. The
weight is 0 at the edge and reaches 1 at exactly `L`. FLIT's cell-based
taper reached 1 one cell earlier. The weights of different sides multiply.
The adaptive-range taper is no longer capped at half the shot range; where
two ramps overlap, their product applies.

## Cone mute

The cone mute keeps full weight inside a radius that grows with depth below
the origin, then ramps linearly to zero:

```
R(d)   = conemuter * ( min(d - d0, conemute_depth) / conemute_depth )^(1/conemutepower)
weight = 1                                   if dist <= R(d)
       = 1 - (dist - R(d)) / conemutetaper   if R(d) < dist < R(d) + conemutetaper
       = 0                                   otherwise, and for d < d0
```

The origin is a direction and a depth:

| `conemuteorigin` | Direction | Depth `d0` |
|---|---|---|
| `source` | Earth-centred centroid of the source points | their mean depth |
| `receiver` | Earth-centred centroid of the receivers | the shallowest receiver depth |

The horizontal distance of a node is the great-circle angle between its
column and the origin direction, times the node's own radius:

```
dist = r_k * atan2( |u x u_c|, u . u_c )
```

where `u` and `u_c` are unit vectors. The Cartesian code measured distance
on a sampling of a quarter cell and offset the radius by one sample; the
spherical code evaluates the formula above directly.

## Wavenumber filtering

`wavenumfilt` filters along longitude lines with FLIT's wavenumber filter,
each line with its own exact spacing `r_k cos(theta_j) dlambda`. Along depth
it uses `ddepth`. The wavenumbers are in cycles per model length unit, as
before. Along latitude, the spacing of a line is not even, so an FFT cannot
be exact there. The keys `shot_<grad>_wavenum_lat` and
`shot_<grad>_wavenum_lat_amps` stop the program.

## Window filters

These filters work on a fixed number of cells:

- `movingbal`, per shot and whole model;
- `medianfilt`, per shot and whole model;
- `rmsbal_lon`, `rmsbal_lat`, `rmsbal_lonlat` and `rmsbal_depth`;
- `andf`, through its smoothing lengths and its `sigma`.

A cell's size varies over the grid, so a length can't be an exact cell count
everywhere. They convert lengths with the reference cell size of step 1, at
the centre of the grid in longitude and latitude, at the top level:

```
cells = nint( L / h_ref )
h_ref = ( r_1 cos(theta_c) dlambda,   r_1 (dtheta/dlat)(lat_c) dlat_rad,   ddepth )
```

For `andf`, the smoothing lengths are divided by the matching `h_ref`, and
`sigma` by the largest of the three. This is the agreed approximation.

`laplacefilt` and TGpV work in grid-index space, unchanged. The 3 x 3 x 3
median of `vpvs_similarity` is unchanged.

## Tikhonov regularization

Tikhonov regularization solves

```
min over u of   lambda * integral |grad u|^2 dV  +  integral (u - m)^2 dV
```

whose optimality condition, `u - lambda div(grad u) = m`, is solved by
Gauss-Seidel iterations, 200 as before. With the exact finite volumes of
the spherical cells, one update is

```
u_c = ( V m_c + lambda sum_n c_n u_n ) / ( V + lambda sum_n c_n ),   c_n = A_n / h_n
```

where `V` is the node's cell volume, and for each neighbour `n` across a face
of area `A_n`, `h_n` is the arc length to it. There is no neighbour across
the edge of the grid. The face areas are the exact ones of step 3:

```
longitude faces:  A = (r_face(k)^2 - r_face(k+1)^2)/2 * (theta_face(j+1) - theta_face(j))
latitude face jf: A = (r_face(k)^2 - r_face(k+1)^2)/2 * cos(theta_face(jf)) * dlambda
radial face kf:   A = r_face(kf)^2 * (sin(theta_face(j+1)) - sin(theta_face(j))) * dlambda
```

`lambda` (`reg_tikhonov_lambda`) is now a squared length in m^2. In a
Cartesian grid of unit cells this update is the Cartesian one, where
`lambda` was effectively in squared cells.

## Defaults

Every length default keeps its Cartesian number of cells, times the
reference cell size `h_ref` along the matching axis:

| Parameter | Default |
|---|---|
| `shot_<grad>_smooth_*`, `<grad>_smooth_*` | 3 h_ref |
| `shot_<grad>_movingbal_*` | 6 h_ref |
| `<grad>_movingbal_*` | 3 h_ref |
| `<grad>_medianfilt_*`, `shot_<grad>_medianfilt_*` | 1 h_ref |
| `<grad>_rmsbal_*` | 1 h_ref |
| `<grad>_andf_smooth_lon`, `_lat`, `_depth` | 2, 2, 8 h_ref |
| `<grad>_andf_sigma` | 6 max(h_ref) |
| `shot_<grad>_conemutetaper` | 10 mean(h_ref) |
| `reg_smooth_*` | 1 h_ref |
| `reg_similarity_smooth_*` | 5 h_ref |
| `reg_tikhonov_lambda` | 10 (mean(h_ref))^2, printed at use |

## Execution changes

- All length parameters of processing and regularization are in meters.
- `reg_tikhonov_lambda` is in m^2.
- `shot_<grad>_wavenum_lat` stops the program.
- Exact smoothing is a direct convolution, so it costs more than FLIT's
  recursive filter for long smoothing lengths: about 12 sigma over the cell
  size operations per node and axis.

## Verification

The scripts of these tests are in `test/spherical/step4`;
`test/spherical/README.md` explains how to run them.

### Exact filters

`step4/test_filters.f90`, on a grid from 10° to 70° latitude with
0.2° x 0.2° x 5 km cells, where an east-west cell at 60° is about half the
width of one at 20°.

| Check | Result |
|---|---|
| Spike smoothed with sigma = 30 km east-west, at 20° and 60° | 30000.0 m at both (1.45 and 2.71 cells) |
| Same spike, north-south (30 km) and vertical (10 km) | 30000.0 m and 10000.2 m |
| Constant field smoothed | unchanged, error 0 |
| Taper weights vs the Blackman ramp at exact distances | max difference 3.0e-8 |
| Tikhonov, lambda = (20 km)^2, spike spread east-west at 20° and 60° | 15239 m and 15271 m |
| Tikhonov of a constant field | unchanged, error 0 |
| Longitude low pass (pass to 1/100 km, stop from 1/60 km), kept RMS fraction at 20° and 60° | 40 km wavelength: 0.0001 and 0.0003; 300 km wavelength: 1.0000 at both |

### Cone mute in x_fatt3

One shot of the step-3 FATT inputs, one iteration, with and without
`process_shot_grad = conemute` (origin at the receivers, radius 60 km
reached at 20 km depth, linear growth, 20 km ramp). The ratio of the two
gradients equals the exact weight computed independently from great-circle
distances, within 1.4e-6, at 50,518 nodes: 5,424 with full weight, 5,125 on
the ramp and 39,969 muted.

### FATT with processing and regularization

The step-3 FATT problem with per-shot smoothing (8 km horizontal, 2 km
vertical), whole-model smoothing (10 km, 3 km), a 6 km depth taper at the
bottom, and Tikhonov (lambda = 1e8 m^2) plus Gaussian regularization from
iteration 2 on.

| Iteration | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| Misfit, relative | 1.000 | 0.573 | 0.361 | 0.175 | 0.115 | 0.080 | 0.062 |

The recovered perturbation correlates with the true one at 0.968 between 4
and 20 km depth, against 0.944 without processing in step 3. A run with
`shot_grad_wavenum_lat` stopped with the error message.
