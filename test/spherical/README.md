# Tests of the spherical branch

These scripts check each step of the spherical conversion. The expected
results, and what they show, are in the step documents in `docs/`. The tests
write their outputs into `work*` folders next to the scripts, which git
ignores.

## Requirements

- The spherical build: `x_eikonal3` and the object directory holding
  `module_spherical.o` and `spherical.mod`.
- The Intel compilers and MPI runtime on `PATH`, e.g.
  `source /path/to/oneapi/setvars.sh`.
- Python 3 with numpy. Set `PYTHON` to choose the interpreter.
- Optional: `x_eikonal3` built from the main branch, for the Cartesian
  comparison in `step2/test_homogeneous.py`.

## Run everything

```
bash run_all.sh <bindir> <objdir> [main_bindir]
```

With the default Makefile paths, `bindir` is `bin` at the top of the repo
and `objdir` is `$HOME/bin/obj/latte3`.

## Tests

| Script | Step | What it checks | Expected |
|---|---|---|---|
| `step1/run_unit_test.sh [objdir] [flitdir]` | 1 | Latitude conversion, Earth-centred positions, cell volumes, arc lengths, reference spacing, length of a degree | Geodetic 45° gives geocentric 44.80757678°. Round trips are exact. Cell volumes sum to the exact region volume within about 1e-15. Arcs exceed chords by about 1 mm. |
| `step1/run_param_checks.sh <bindir>` | 1 | Rejection of Cartesian and unsupported parameters; offsets and adaptive range | 26 errors for `param_bad.rb`. `param_good.rb` runs with `ns, nr = 1, 1` and a shot range of 101.5–102.5° by 22–23°. Grid errors for `nlat = 1` and for a grid past 90°. |
| `step2/test_homogeneous.py --bin <bindir> [--main-bin <dir>]` | 2 | Homogeneous Earth against exact straight-line times, optionally against the main branch | Spherical: max error 0.139 s with 2 km cells and 0.051 s with 10 km cells. Main branch, flat grids: 0.132 s and 0.069 s. |
| `step2/test_gradient.py --bin <bindir> [--levels N]` | 2 | Linear velocity gradient with an exact solution, over 4 grid refinements | Max error, over nodes whose ray stays in the grid, falls from 0.238 s to 0.016 s. |
| `step2/test_radial.py --bin <bindir>` | 2 | Velocity depending only on depth, against exact spherical times | Solver error about 0.19 s or less out to 678 km, while flat-Earth times are off by up to 1.33 s. Rows from 730 km are marked: their exact ray bottoms below the grid. |
| `step3/run_gradient_fd.sh [objdir] [flitdir]` | 3 | Adjoint gradient with cell volume against central finite differences, for Gaussian velocity bumps at 15–55° latitude and 30–220 km depth | Ratio of adjoint to finite difference 1.002–1.004 away from the source and receivers. Without the cell volume it drifts from 0.78 to 1.30. Near the source (1.20) and the surface receivers (1.02) the adjoint method is less accurate, as in the Cartesian code. |
| `step3/run_fatt.sh <bindir>` | 3 | Synthetic data with `x_eikonal3`, then 6 iterations of `x_fatt3`; `adjoint_receiver = nearest` | Misfit falls to 0.8% of its starting value. `nearest` stops with an error. |
| `step3/run_tloc.sh <bindir>` | 3 | `x_tloc3` updating Vp and origin times that start 0.3 s off (`file_st0`); run after `run_fatt.sh`; a run without `file_st0` | Misfit falls to 4.6% of its starting value; origin times move from 0.3 s to 0.068–0.084 s. The run without `file_st0` stops with an error. |
| `step4/run_filters.sh [objdir] [flitdir]` | 4 | Gaussian smoothing, taper, Tikhonov and longitude wavenumber filtering, at 20° and 60° latitude | Smoothing widths exactly as requested in meters at both latitudes, constant fields unchanged, taper equal to the Blackman ramp within 3e-8, Tikhonov spread 15.2 km at both latitudes, 40 km wavelength removed and 300 km kept. |
| `step4/run_conemute.sh <bindir>` | 4 | Cone mute applied by `x_fatt3` against the exact weight from great-circle distances | Ratio of muted to unmuted gradient equals the exact weight within 1.4e-6 at 50,518 nodes. |
| `step4/run_fatt_processing.sh <bindir>` | 4 | `x_fatt3` with smoothing, taper, Tikhonov and Gaussian regularization; run after `run_conemute.sh`; `shot_grad_wavenum_lat` | Misfit falls to 6% of its starting value, correlation with the true model 0.968; `wavenum_lat` stops with an error. |
| `step5/run_relocation.sh <bindir>` | 5 | `x_tloc3` relocating 16 events moved 3–6 km horizontally and 3 km in depth, at 60° latitude, with 2000 m step limits; about 15 minutes with 16 threads | After 12 iterations: horizontal error median 39 m (max 124 m), depth error median 68 m (max 554 m). No event moves more than 2000 m in longitude, latitude or depth in one iteration (largest 1738 m). A preliminary run without `file_slat` stops with an error. |
