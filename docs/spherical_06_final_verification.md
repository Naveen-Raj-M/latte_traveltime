# Step 6: Final verification

The spherical branch was built from scratch and the tests of steps 1 to 4
were run on that build. They reproduced every result recorded in the step
documents. The relocation test of step 5 was not rerun in this step; its
results, from the step 5 code, are in the step 5 document.

## Files changed

| File | What changed |
|---|---|
| `test/spherical/step3/run_gradient_fd.sh` | Links `module_spherical_filter.o`, which `module_vars` needs since step 4. Without it the test no longer linked. |
| `src/3d/module_vars.f90` | `file_st0` is required when `st0` is updated, like the position files (see step 5). |
| `test/spherical/step3/make_tloc_inputs.py`, `param_tloc.rb`, `run_tloc.sh` | The TLOC test gives its 0.3 s starting origin times in `file_st0`, and checks that a run without it stops. Before, the starting times were only in the geometry, which TLOC does not use for an updated `st0`; see the correction in the step 3 document. |
| `test/spherical/step5/make_inputs.py`, `param_tloc.rb` | The relocation test gives its starting origin times, the true 0, in `file_st0`, now required. The test starts from the same values as before, so its step 5 results stand. A short start of `x_tloc3` confirmed it passes the file check; the relocation itself was not rerun. |
| `docs/README.md` | New. Index of the step documents and a short guide: geometry, units, a minimal parameter file, file layouts, names, unsupported features. |

The only source change is the `file_st0` requirement.

## Build

`make clean`, then `make` in `src/3d`, with the Intel oneAPI compilers and
FLIT. The build made `x_eikonal3`, `x_fatt3` and `x_tloc3` with no errors
and no warnings.

## Results

| Test | Result | As documented |
|---|---|---|
| Step 1: geometry unit test | Geodetic 45° to geocentric 44.80757678°; cell volumes tile the region within 1.3e-15 | yes |
| Step 1: parameter checks | 26 errors for Cartesian and unsupported entries; correct offsets and adaptive range; grid errors for `nlat = 1` and a grid past 90° | yes |
| Step 2: homogeneous Earth | Max error 0.139 s (2 km cells) and 0.051 s (10 km cells); main branch 0.132 s and 0.069 s | yes |
| Step 2: gradient refinement | Max error 0.238, 0.135, 0.062, 0.016 s over 4 grids | yes |
| Step 2: velocity depending on depth | Solver error 0.195 s at 678 km; flat Earth 1.333 s | yes |
| Step 3: finite-difference gradient | Ratio 1.002 to 1.004 away from source and receivers; 0.78 to 1.30 without the cell volume | yes |
| Step 3: FATT | Misfit 0.008 after 6 iterations; `nearest` rejected | yes |
| Step 3: TLOC origin times, rerun with `file_st0` | Misfit 0.046 after 6 iterations; origin times from 0.3 s to 0.068–0.084 s; a run without `file_st0` stops | corrected in step 3 |
| Step 4: exact filters | Smoothing widths exact at 20° and 60°; taper within 3e-8; Tikhonov spread 15.2 and 15.3 km; wavenumber low pass correct | yes |
| Step 4: cone mute | Within 1.4e-6 of the exact weight at 50,518 nodes | yes |
| Step 4: FATT with processing | Misfit 0.062 after 6 iterations, correlation 0.968; `wavenum_lat` rejected | yes |

## Running the tests

`test/spherical/run_all.sh` runs every test, including the relocation of
step 5, which takes about 15 minutes with 16 OpenMP threads. To leave it
out, run the step scripts listed in `test/spherical/README.md` one by one.
