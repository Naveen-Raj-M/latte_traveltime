# Step 3: Adjoint solver, receiver injection and volume-weighted gradients

This step converts the adjoint-state solver to the spherical Earth, makes
linear receiver injection the only option, and weights the gradients, the
DWS and the regularization term by cell volume. After this step, `x_fatt3`
and `x_tloc3` compute correct spherical gradients.

**Intermediate state.** The gradient filters and the model regularizers
(smoothing, tapers, cone mute, Tikhonov and the others) still convert
lengths with `dlon` and `dlat` in degrees; step 4 converts them. TLOC event
positions still step in degrees with degree limits; step 5 converts them.
Runs without gradient filters or regularization are already correct.

## Files changed

| File | What changed |
|---|---|
| `src/3d/module_traveltime_iso.f90` | `fast_sweep_adjoint` weights each face flux by the exact face area over the exact cell volume, and differences traveltimes over arc lengths. `adjoint_iso` has no receiver option: it always injects residuals as point sources, divided by each node's exact cell volume. The Dirichlet ("nearest") injection and its receiver flags were removed. |
| `src/module_spherical.f90` | New: the mean cell volume of the grid, relative cell volumes, and the exact face weights of the three axes. |
| `src/3d/submodule_gradient_fatt.f90` | The gradient, and the preconditioner when it is on, are weighted by relative cell volume. The receiver option is no longer passed. |
| `src/3d/submodule_gradient_tloc.f90` | Same weighting for the Vp and Vs gradients of TLOC, including the double-difference version, and for the DWS. TLOC now uses linear injection; before, it always used "nearest". |
| `src/module_inversion_regularization.f90` | The regularization term of the gradient is weighted by relative cell volume. |
| `src/module_parameters.f90` | `adjoint_receiver` defaults to `linear`; any other value stops the program. |

## The adjoint-state equation

For a shot with traveltime field `T`, the adjoint field `lambda` solves

```
-div( lambda grad T ) = sum over receivers r of  dt_r delta(x - x_r)
```

where `dt_r` is the weighted traveltime residual at receiver `r`.

### Discretization

The Cartesian LATTE solves this with a conservative upwind scheme (Leung and
Qian, 2006; equation A-9 of Taillandier et al., 2009). That scheme is a
balance of the fluxes through the six faces of a node's cell, which this
step writes out with the faces of the spherical cell.

For each axis, the traveltime gradients on the two faces of the cell are

```
a+ = (T_+ - T_c) / h_+        a- = (T_c - T_-) / h_-
```

with `T_c` at the node, `T_+` and `T_-` at the neighbours on the + and -
side, and `h` the arc length to the neighbour (step 2):

| Axis | + side | h_+ | h_- |
|---|---|---|---|
| Longitude | east, i + 1 | r_k cos(theta_j) dlambda | r_k cos(theta_j) dlambda |
| Latitude | north, j + 1 | r_k (theta_{j+1} - theta_j) | r_k (theta_j - theta_{j-1}) |
| Depth | down, k + 1 | ddepth | ddepth |

With `pos(a) = (a + |a|)/2`, `neg(a) = (a - |a|)/2`, and `w+`, `w-` the
weights of the two faces, the update of the node is

```
lhs = sum over axes of ( w+ neg(a+) - w- pos(a-) )
rhs = sum over axes of ( w- neg(a-) lambda_-  -  w+ pos(a+) lambda_+ )
lambda_c = (rhs - s_c) / lhs          (lambda_c = 0 if lhs = 0)
```

where `s_c` is the point-source density at the node. In the Cartesian code
every weight is `1/h`, and this is exactly the original update. As before,
the sweeps keep the smaller of the old and new value at each node, and a
face on the edge of the grid has `a = 0`, so no flux crosses it.

### Face weights

Each weight is the exact area of a face of the node's cell over the exact
cell volume `V(j, k)` of step 1. With the cell faces of step 1 (latitudes
`theta_face(j)` and `theta_face(j+1)`, radii `r_face(k)` above and
`r_face(k+1)` below):

```
longitude faces (both):  w = (r_face(k)^2 - r_face(k+1)^2)/2 * (theta_face(j+1) - theta_face(j)) / V
southern face:           w- = (r_face(k)^2 - r_face(k+1)^2)/2 * cos(theta_face(j))   * dlambda / V
northern face:           w+ = (r_face(k)^2 - r_face(k+1)^2)/2 * cos(theta_face(j+1)) * dlambda / V
upper face:              w- = r_face(k)^2   * (sin(theta_face(j+1)) - sin(theta_face(j))) * dlambda / V
lower face:              w+ = r_face(k+1)^2 * (sin(theta_face(j+1)) - sin(theta_face(j))) * dlambda / V
```

The differences between the two latitude weights and between the two depth
weights carry the spherical spreading: a ray tube widens or narrows as it
crosses latitudes and radii, and these weights make the flux account for it.

### Receiver injection

The residual of a receiver goes to the eight nodes of the receiver's cell
with the trilinear weights `w_n` by which `forward_iso` interpolates the
receiver's traveltime, divided by each node's own cell volume:

```
s_n += dt_r w_n / V(n)
```

This is the transpose of the interpolation, as a density. It was the
`linear` option of the Cartesian code, where every cell had the volume
`dx dy dz`. It is now the only option, for FATT and TLOC.
`adjoint_receiver = nearest` stops the program.

## Gradients

The gradient of the misfit with respect to the velocity at a node is the
adjoint field times the volume the node represents. With the relative cell
volume

```
V_rel(j, k) = V(j, k) / V_mean
```

where `V_mean` is the mean cell volume of the grid, computed once at
startup, the per-shot gradient is

```
without preconditioner:  g = -lambda V_rel / v^3
with preconditioner:     g = -lambda V_rel / ( E V_rel + eps max|E V_rel| )
```

where `E` is the preconditioner's adjoint field, which has residual -1 at
every receiver, and `eps` is `precond_eps`. With the preconditioner, `V_rel`
cancels except in the damping term.

Relative rather than absolute volume keeps the gradient magnitudes close to
those of Cartesian LATTE. The two differ by the constant factor `V_mean`,
to which step sizes, conjugate gradients, L-BFGS and the automatic
regularization scaling are insensitive.

### DWS

In TLOC, the DWS output accumulates the preconditioner field weighted by
relative cell volume, `|E V_rel|`, so it measures coverage per cell. FATT
does not accumulate a DWS; that is unchanged from the main branch.

### Regularization term

With model regularization, the gradient of a velocity model gets the term

```
g += c * V_rel * (m - m_reg)
```

the gradient of the volume integral of the squared difference. `c` is
`reg_lambda` with `const_reg`. Otherwise `c` is computed as before, as the
mean of `g` over the mean of `V_rel (m - m_reg)`, times `reg_scale`.

## Execution changes

- `adjoint_receiver` defaults to `linear`; `nearest` stops the program.
- TLOC now injects residuals as point sources, like FATT. On the main
  branch it always used the Dirichlet ("nearest") condition.
- Gradient files, `grad_vp.bin` and `grad_vs.bin`, and the DWS are weighted
  by relative cell volume.
- `adjoint_iso` no longer takes a receiver argument.

## Verification

The scripts of these tests are in `test/spherical/step3`;
`test/spherical/README.md` explains how to run them.

### Finite-difference check of the gradient

Grid: longitude 0–20°, latitude 10–60°, depth 0–600 km; cells of 0.25° and
10 km. Velocity `6000 + 0.01 d` m/s. One source at (10°, 35°, 300 km) and
3,861 surface receivers. The functional is the sum of all receiver
traveltimes, so every residual is 1. For each Gaussian velocity bump of
150 km radius and 100 m/s amplitude, the change of the functional predicted
by the adjoint gradient is compared with a central finite difference of two
forward runs.

| Bump centre | Adjoint / FD with V | Adjoint / FD without V |
|---|---|---|
| 10°, 15°, 80 km | 1.0022 | 0.7784 |
| 10°, 25°, 80 km | 1.0037 | 0.8299 |
| 10°, 45°, 80 km | 1.0037 | 1.0566 |
| 10°, 55°, 80 km | 1.0022 | 1.2966 |
| 10°, 29°, 220 km | 1.0029 | 0.8952 |
| 12°, 30°, 30 km, at the receivers | 1.0224 | 0.8798 |
| 10°, 35°, 150 km, reaching the source | 1.1989 | 1.1573 |

With the cell volume, the adjoint gradient matches the finite difference
within 0.4% wherever the bump is away from the source and receivers. Without
it, the ratio drifts by the variation of `r^2 cos(theta)` across the grid,
from 0.78 to 1.30.

The adjoint field is singular at the source and at the receivers, so the
gradient is less accurate where a bump overlaps them. The Cartesian code of
the main branch shows the same, run on a comparable flat grid with linear
injection:

| Main branch, Cartesian, bump | Adjoint / FD |
|---|---|
| Reaching the source | 1.1304 |
| At the receivers | 1.0682 |
| Away from both | 1.0464 |

### End-to-end FATT

Grid of 3° by 3° by 40 km with 0.05° by 2 km cells. True model: velocity
increasing with depth, plus a +8% and a -8% Gaussian anomaly at 12 km
depth. Synthetic data from 9 sources at 15 km depth and 49 surface stations,
computed with `x_eikonal3`. `x_fatt3` started from the model without
anomalies, with no gradient filters.

| Iteration | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| Misfit, relative | 1.000 | 0.551 | 0.310 | 0.097 | 0.031 | 0.010 | 0.008 |

The recovered perturbation correlates with the true one at 0.944 between 4
and 20 km depth. Both anomalies come back within 0.05° of their true
centres, with amplitudes of +5.6% and -5.8%. A run with
`adjoint_receiver = nearest` stopped with the error message.

### End-to-end TLOC

The same model and data, with the 9 events at their true positions but with
starting origin times of 0.3 s instead of 0, given in `file_st0`. `x_tloc3`
updated Vp and origin times from the starting velocity model.

| Iteration | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| Misfit, relative | 1.000 | 0.385 | 0.227 | 0.118 | 0.077 | 0.054 | 0.046 |

After 6 iterations the origin times were between 0.068 and 0.084 s, against
the true 0 and the starting 0.3 s, and still converging. Their remaining
offsets are nearly equal; with Vp updated at the same time, part of a common
origin-time shift is absorbed by the velocity model.

**Correction.** An earlier version of this test, and of this section, gave
the 0.3 s starting times only in the geometry, without `file_st0`. When
`st0` is updated, TLOC takes the starting origin times from the `st0` model,
not from the geometry, and without the file that model was 0, the true
value. That run therefore started at the true origin times; its reported
end values of 0.001 to 0.029 s were drift from 0, not recovery from 0.3 s.
Since step 5, `file_st0` is required when `st0` is updated, and the test
gives it.
