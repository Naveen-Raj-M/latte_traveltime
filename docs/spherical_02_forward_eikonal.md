# Step 2: Forward eikonal solver in spherical coordinates

This step converts the factorized fast-sweeping eikonal solver,
`forward_iso`, to the spherical Earth of step 1. After this step,
`x_eikonal3` computes correct first-arrival traveltimes in spherical
coordinates.

**Intermediate state.** The adjoint solver, gradients and gradient
filters are still Cartesian, so `x_fatt3` and `x_tloc3` are not yet
correct. Steps 3 to 5 convert them.

## Files changed

| File | What changed |
|---|---|
| `src/3d/module_traveltime_iso.f90` | The sweep uses the arc length to each upwind neighbour instead of a fixed `dx`, `dy`, `dz`. The background time is the straight-line distance between Earth-centred positions, with its gradient along the local east, north and down directions. The same-cell receiver time uses that distance too. A new routine, `set_grid_offsets`, locates each shot's sub-grid within the global grid. |

No other file changed. The local solvers `local_solver_2d` and
`local_solver_3d` are unchanged: they take spacings as arguments, and they
now receive the spherical ones.

## The equation

In spherical coordinates (radius `r`, geocentric latitude `theta`,
longitude `lambda`) the isotropic eikonal equation is

```
(dT/dr)^2 + (1/r^2)(dT/dtheta)^2 + (1/(r cos(theta)))^2 (dT/dlambda)^2 = 1/v^2
```

Written in the local orthonormal frame of each point, east (`e`), north
(`n`) and down (`d`), it has the Cartesian form

```
(dT/de)^2 + (dT/dn)^2 + (dT/dd)^2 = 1/v^2

dT/de = (1/(r cos(theta))) dT/dlambda
dT/dn = (1/r) dT/dtheta
dT/dd = dT/d(depth) = -dT/dr
```

This is why the Cartesian local solver can be reused: only the spacings of
its one-sided differences change.

## Discretization

The grid axes are longitude (`i`, eastward), geodetic latitude (`j`,
northward) and depth (`k`, downward). The one-sided difference toward an
upwind neighbour divides by the arc length to that neighbour:

```
toward i - 1 or i + 1:  h_lon   = r_k cos(theta_j) dlambda
toward j - 1:           h_lat-  = r_k (theta_j - theta_{j-1})
toward j + 1:           h_lat+  = r_k (theta_{j+1} - theta_j)
toward k - 1 or k + 1:  h_depth = ddepth
```

`theta_j` is the geocentric latitude of row `j`. Because the grid is
regular in geodetic latitude, `h_lat-` and `h_lat+` differ slightly, and
each octant of the solver uses the one on its side. For example, the
octant whose upwind neighbours are (i - 1, j + 1, k - 1) uses `h_lon`,
`h_lat+` and `h_depth`.

These spacings are computed for each node from the 1D geometry arrays of
step 1, using global node indices. A shot's adaptive sub-grid starts at
global node (`ioff + 1`, `joff + 1`, `koff + 1`), where

```
ioff = nint((lon_first - olon)/dlon)
joff = nint((lat_first - olat)/dlat)
koff = nint((depth_first - odepth)/ddepth)
```

## Factorization

The traveltime is factored as `T = T0 * tau`, as before. The background
time `T0` of a source point `s` with velocity `v_s` and start time `t_s`
is now the exact traveltime of a homogeneous sphere:

```
T0(x) = |X(x) - X(s)| / v_s + t_s
```

where `X` is the Earth-centred Cartesian position of step 1. For a source
of many points, `T0` is the minimum over the points, as before. Its gradient
is

```
grad T0 = (X(x) - X(s)) / (v_s |X(x) - X(s)|)
```

and the solver takes its components along the node's unit vectors:

```
e = (-sin(lambda),               cos(lambda),               0)
n = (-sin(theta) cos(lambda),   -sin(theta) sin(lambda),    cos(theta))
d = -(cos(theta) cos(lambda),    cos(theta) sin(lambda),    sin(theta))

px0 = grad T0 . e,   py0 = grad T0 . n,   pz0 = grad T0 . d
```

With the spacings above, the local solver's quadratic for `tau` at a node
is the Cartesian one:

```
sum over the 3 axes of ( a_m tau + b_m )^2 = 1/v^2
a_m = p0_m + T0_c sign_m / h_m
b_m = -tau_m T0_c sign_m / h_m
```

where `m` runs over east, north and down, `T0_c` is `T0` at the node,
`tau_m` is `tau` at the upwind neighbour, `h_m` is the spacing to it, and
`sign_m` is +1 for a neighbour at the lower index and -1 for the higher.
The causality checks and the 2D and 1D fallbacks are unchanged, with the
same spacings.

For many source points, the bins that speed up the minimum are now built
on Earth-centred positions. The distance from a node to a bin's box is
still a lower bound of its distance to every point in the bin, so the
pruning stays exact.

## Receivers

Receiver traveltimes are still interpolated trilinearly in grid-index space.
A receiver in the same grid cell as a source point still takes the analytic
time from it, now

```
t = |X(receiver) - X(source)| / v_source + t_source
```

## Precision

Earth-centred positions are always computed in double precision. Grid
positions come from the 1D geometry arrays allocated once in step 1, so no
3D coordinate arrays are allocated. Source positions are converted once per
shot into small double-precision arrays, as the Cartesian code already did
for its source coordinates.

## Execution changes

- `x_eikonal3` now gives spherical traveltimes.
- `fast_sweep_forward` no longer takes `dx`, `dy`, `dz`. It computes the
  spacings per node.
- The CPU cost per sweep rises slightly for these per-node spacings. The
  number of sweeps is unchanged; the tests below converged in 2 iterations.

## Verification

All tests ran in a scratch directory with `x_eikonal3`.

### Homogeneous Earth

A 10° by 8° by 80 km grid with velocity 6000 m/s and an off-node source
at 15.3 km depth.

| Run | Cells | Max error | Mean error |
|---|---|---|---|
| Spherical, this step | 0.1° x 0.1° x 10 km | 0.051 s | 0.013 s |
| Cartesian, main branch | 10 x 10 x 10 km | 0.069 s | 0.020 s |
| Spherical, this step | 0.1° x 0.1° x 2 km | 0.139 s | |
| Cartesian, main branch | 10 x 10 x 2 km | 0.132 s | |

Neither code is exact for a homogeneous model, even though `T0` is exact
there. The largest errors are directly below the source. At such a node,
the neighbours on both sides in longitude and latitude arrive later. The
3D update fails the causality check in every octant, and the solver falls
back to its 1D and 2D updates, which are not exact. Flat cells make this
worse. This behaviour comes from the main branch's local solver and is
unchanged here; the spherical errors match the Cartesian ones.

### Linear velocity gradient, grid refinement

Velocity increases linearly along the downward direction at the centre of
a 4° by 4° by 60 km grid: `v = 5000 + 0.02 s` m/s, where `s` is the
distance below the surface tangent plane there. This model has the exact
solution

```
T = (1/G) arccosh( 1 + G^2 |X - X_s|^2 / (2 v(X) v(X_s)) )
```

Its rays are circular arcs. For distant pairs of points, the exact ray dips
below the 60 km floor of the model, which the solver cannot follow. The
table counts only nodes whose exact ray stays above 50 km, about two thirds
of the nodes.

| Grid | Nodes | Mean error | Max error | Order of max error |
|---|---|---|---|---|
| 0.1° x 10 km | 11,767 | 0.0278 s | 0.238 s | |
| 0.05° x 5 km | 85,293 | 0.0186 s | 0.135 s | 0.82 |
| 0.025° x 2.5 km | 648,025 | 0.0079 s | 0.062 s | 1.12 |
| 0.0125° x 1.25 km | 5,049,009 | 0.0021 s | 0.016 s | 1.97 |

The error converges at first order or better. At the nodes whose ray dips
below the model, the error stays at about 1.15 s on every grid. That is the
expected effect of the model's depth limit, not a solver error.

### Velocity depending on depth, against exact spherical traveltimes

Velocity `v = 5800 + 0.025 d` m/s for depth `d` up to 200 km, on a
0.05° x 0.05° x 5 km grid. The source is 10 km deep, with surface receivers
at 0.5° steps along its latitude. Exact spherical first arrivals come from
a ray-parameter integral, which is exact for a velocity that depends only on
depth. The last two columns are the times of the same velocity profile in a
flat Earth, as the Cartesian code would compute it.

| Distance | Solver | Exact spherical | Solver error | Flat Earth | Flat-Earth error |
|---|---|---|---|---|---|
| 104 km | 17.514 s | 17.513 s | 0.001 s | 17.537 s | 0.024 s |
| 209 km | 34.097 s | 34.084 s | 0.013 s | 34.185 s | 0.100 s |
| 313 km | 49.363 s | 49.324 s | 0.038 s | 49.581 s | 0.257 s |
| 417 km | 63.102 s | 63.027 s | 0.075 s | 63.520 s | 0.493 s |
| 521 km | 75.350 s | 75.233 s | 0.118 s | 76.027 s | 0.794 s |
| 625 km | 86.257 s | 86.097 s | 0.159 s | 87.240 s | 1.143 s |
| 678 km | 91.277 s | 91.083 s | 0.194 s | 92.416 s | 1.333 s |

The solver follows the exact spherical times. Its error is about one
seventh of the flat-Earth error, and shrinks with a finer grid. From 730 km
on, the exact first-arrival ray bottoms below the 200 km floor of this test
model, at 210 km and deeper, so those distances test the model's depth
rather than the solver and are left out.
