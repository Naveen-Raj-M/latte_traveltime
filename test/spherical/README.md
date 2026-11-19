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
