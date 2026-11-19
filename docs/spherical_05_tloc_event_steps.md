# Step 5: TLOC event steps in meters

This step makes the step limits of TLOC event positions physical.
`step_max_slon` and `step_max_slat` are now in meters, enforced exactly for
every event at its own latitude and depth. With it, all of `x_eikonal3`,
`x_fatt3` and `x_tloc3` are spherical.

## Files changed

| File | What changed |
|---|---|
| `src/module_inversion_step_size.f90` | Event longitude and latitude steps are limited in meters: a first step from each event's local length of a degree, then an exact arc-length check of every event's move. Defaults of `step_max_slon` and `step_max_slat` in meters. The largest event move in meters is printed. |

## Event positions and their gradients

Events, the virtual receivers of TLOC, are stored and updated as before:

| Parameter | Unit |
|---|---|
| `slon` | degrees of longitude |
| `slat` | degrees of geodetic latitude |
| `sdepth` | meters below mean sea level |
| `st0` | seconds |

Their gradients need no change. TLOC differentiates the traveltime field of
each station across the grid nodes around the event:

```
dT/dslon   = ( T(i+1) - T(i-1) ) / (2 dlon)      seconds per degree
dT/dslat   = ( T(j+1) - T(j-1) ) / (2 dlat)      seconds per degree
dT/dsdepth = ( T(k+1) - T(k-1) ) / (2 ddepth)    seconds per meter
```

The grid is regular in longitude, geodetic latitude and depth, so these are
the derivatives with respect to exactly the stored variables. The chain rule
from physical distance to degrees is built in, with no metric needed.

## Step limits

LATTE scales each parameter's step separately: the search direction `p` of
a parameter is multiplied by a step `a`, and `a` is chosen so that the
largest change stays below `step_max`. For `sdepth` and `st0` this is
unchanged.

For `slon` and `slat`, `step_max` is a distance in meters. An event `e` at
geodetic latitude `lat_e`, geocentric latitude `theta_e = theta(lat_e)` and
radius `r_e = R - depth_e` moves by the exact arc

```
longitude change dlon:  m_e = r_e cos(theta_e) |dlon| pi/180          along its parallel
latitude change dlat:   m_e = r_e | theta(lat_e + dlat) - theta_e |    along its meridian
```

where `theta()` is the WGS84 geodetic-to-geocentric relation of step 1.

**First step.** With the local length of a degree at each event,

```
longitude:  L_e = r_e cos(theta_e) pi/180
latitude:   L_e = r_e (dtheta/dlat)(lat_e) pi/180

a = step_max / max over events of ( |p_e| L_e )
```

so the event with the largest move moves about `step_max` meters.

**Check.** Each time LATTE checks a trial step `a s` (with its step scaling
`s`), it computes every event's exact arc `m_e` for `dlon` or
`dlat = a s p_e` and accepts only if

```
max over events of m_e < step_max
```

LATTE halves `s` until the check passes, so the limit holds exactly, even
though the first step used a linearized length for latitude.

The events' positions are those being updated, from the `slon`, `slat` and
`sdepth` models, or else from the geometry for a coordinate not being
updated.

## Defaults

| Parameter | Default |
|---|---|
| `step_max_slon` | 0.1 (nlon - 1) h_ref_lon meters, a tenth of the east-west extent at the reference cell size |
| `step_max_slat` | 0.1 (nlat - 1) h_ref_lat meters |
| `step_max_sdepth` | 0.1 (ndepth - 1) ddepth meters, unchanged |
| `step_max_st0` | 0.1 s, unchanged |

`h_ref` is the reference cell size of step 1, at the centre of the grid, top
level. This follows the rule of step 4: the Cartesian default in cells
times the reference cell size.

`min_slon`, `max_slon`, `min_slat` and `max_slat` stay in degrees, by
default the grid's range.

## Starting positions

As in the Cartesian examples, the starting event positions come from the
model files `file_slon`, `file_slat` and `file_sdepth`, one value per event
in geometry order. Without such a file, LATTE sets that coordinate of every
event to 0, as on the main branch. In spherical coordinates that means
longitude 0 or latitude 0, usually outside the grid. This step does not
change that behaviour.

## Execution changes

- `step_max_slon` and `step_max_slat` are in meters, and so are their
  defaults.
- The step-size output prints the largest event move by `slon` and by `slat`
  in meters.

## Verification

The scripts of this test are in `test/spherical/step5`;
`test/spherical/README.md` explains how to run it.

Relocation at 60° latitude, where a degree of longitude is about 55 km, half
its length at the equator, so a limit in degrees would differ clearly from a
limit in meters. Grid: longitude 10–14°, latitude 58–62°, depth 0–40 km;
cells 0.04° x 0.02° x 1 km, about 2.2 x 2.2 x 1 km. Velocity
`5800 + 0.03 d` m/s. Synthetic data from 16 events at 8–25 km depth and 49
surface stations, computed with `x_eikonal3`. `x_tloc3` started from the
events moved 3–6 km east and north and 3 km in depth, with the true
velocity, updating `slon`, `slat`, `sdepth` and `st0`, with
`step_max_slon = step_max_slat = step_max_sdepth = 2000` m.

| Iteration | Misfit, relative | Largest move, lon / lat / depth | Horizontal error, median / max | Depth error, median / max |
|---|---|---|---|---|
| 0 | 1.0000 | | 6780 / 7924 m | 3000 / 3000 m |
| 1 | 0.5320 | 1497 / 1497 / 1497 m | 5113 / 6493 m | 2634 / 3992 m |
| 3 | 0.0860 | 1738 / 1738 / 1738 m | 2186 / 3075 m | 1550 / 5579 m |
| 6 | 0.0047 | 1038 / 1037 / 1038 m | 318 / 1398 m | 537 / 2000 m |
| 9 | 0.0007 | 384 / 384 / 384 m | 47 / 582 m | 147 / 1282 m |
| 12 | 0.0001 | 214 / 214 / 214 m | 39 / 124 m | 68 / 554 m |

The moves were measured independently from the saved positions, as exact
arcs. No event moved more than 2000 m in longitude, latitude or depth in any
iteration; the largest move was 1738 m. The largest moves printed by
`x_tloc3` agree with the independent ones, for example 1496.7 m at
iteration 1. The three moves are equal within an iteration because the
three limits are equal and one step scaling applies to all parameters. The
run was stopped after 12 iterations; the events were still converging.
