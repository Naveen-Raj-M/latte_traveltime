# Step 0: Reduce the branch to the spherical scope

This step removes every part of LATTE that the spherical branch will not
support, before any Cartesian-to-spherical conversion. Nothing in this step
changes the mathematics of the code that remains. The remaining code is
still the original Cartesian code; later steps convert it.

The scope of the branch is the 3D eikonal solver (`x_eikonal3`), 3D
first-arrival traveltime tomography (`x_fatt3`) and 3D source location with
optional joint tomography (`x_tloc3`), in acoustic or elastic media, on CPUs.

The main branch keeps everything that was removed here.

## Removed directories

| Directory | What it held |
|---|---|
| `example/` | All 18 Cartesian examples, including 2D, TRTT and topography cases |
| `test/` | Cartesian tests: 2D/3D adjoint-receiver kernels, CUDA vs CPU comparison, 3D sweep test |
| `ml/` | Python ML fault-detection models used only by ML source regularization |
| `src/2d/` | The whole 2D code |

`misc/` and `doc/` are kept unchanged. `doc/parameters.md` still describes
the Cartesian parameter names.

## Removed source files

| File | Reason |
|---|---|
| `src/main_trtt.f90` | TRTT, joint transmission-reflection tomography, is out of scope |
| `src/3d/submodule_gradient_trtt.f90` | TRTT gradients |
| `src/3d/submodule_reflector_imaging.f90` | Reflector imaging for TRTT |
| `src/3d/module_traveltime_iso_reflection.f90` | Reflection traveltimes, used by TRTT and by `x_eikonal3` with a `refl` model |
| `src/3d/module_traveltime_iso_cuda.f90` | Fortran interface to the CUDA solvers |
| `src/3d/traveltime_iso_cuda.cu` | CUDA solvers |
| `src/module_raypath.f90` | Ray tracing from traveltime gradients, used only by `x_eikonal3` |

## Changed source files

| File | What changed |
|---|---|
| `src/3d/Makefile` | Builds only `x_eikonal3`, `x_fatt3` and `x_tloc3`. All CUDA rules are gone. `use_cuda=on` or `-Dlegacy_solver` in `fflags` now stops the build with an error message. |
| `src/install.rb` | Builds only the 3D directory, in the foreground. Previously it built 2D and 3D in the background. |
| `src/main_eikonal.f90` | Rewritten for 3D first arrivals only. Removed: 2D branches, reflection traveltimes, ray-path tracing and output, traveltime-gradient output and its metadata file. The temporary debug print of the velocity profile is kept. |
| `src/3d/module_traveltime_iso.f90` | Removed the legacy fast-sweeping solver and both CUDA call sites. `forward_iso` no longer has the optional outputs `pdx`, `pdy`, `pdz` (traveltime gradients), whose only users were removed. |
| `src/module_parameters.f90` | Removed TRTT, reflector, ray-path and gradient-output variables and their reading: `dir_gradient`, `yn_save_traveltime_gradients`, `incident_wave`, `nrefl`, `yn_update_reflector`, `offset_min_refl`, `offset_max_refl`, `offset_max_refl_depth`, `yn_refl_depth_offsets`, the `refl` model limits, the function `source_receiver_coincide` (unused) and the subroutine `reflection_offset_range`. Removed source regularization switches (`source_regularization_method`, `yn_regularize_source`). |
| `src/module_utility.f90` | Removed the 2D versions of `point_in_domain`, `get_point_value_inside` and `within_the_same_grid`, the TRTT branch of `compute_shot_misfit`, and `write_gradient_metadata`. Observed data are read as one column, first arrivals only. |
| `src/3d/module_gradient.f90` | Removed reflector building, refinement and joint-update routines, the TRTT and reflector-imaging interfaces, and the dip filter (`dipfilt`). All other gradient processing options are kept. |
| `src/3d/module_regularization.f90` | Removed the `structure` (ANDF) model regularization and the whole source regularization (`clustering-fit` and `ml`). Kept: Tikhonov, smooth, TGpV and Vp/Vs similarity. |
| `src/module_inversion_gradient.f90` | Removed reflector imaging and refinement, and `misfit_of_current_model`, which only they used. |
| `src/module_inversion_regularization.f90` | Removed the source-regularization scheduling. Only Vp and Vs are regularized. |
| `src/module_inversion_step_size.f90` | Removed reflector step handling and reflection-data weights. |
| `src/module_inversion_search_direction.f90` | Removed the 2D branches. |
| `src/3d/module_vars.f90` | Removed nearest-neighbour interpolation of reflector models and the TRTT program case. |

## Kept on purpose

These are kept and will be converted to spherical coordinates in later steps:

- Gradient processing: `smooth`, `maxbal`, `rmsbal`, `movingbal`, `medianfilt`, `removenan`, `laplacefilt`, `wavenumfilt`, `mask`, `taper`, `conemute`, `scale`, `signed_power`, `andf`, `rmsbalx`, `rmsbaly`, `rmsbalxy`, `rmsbalz`.
- Model regularization: `Tikhonov`, `smooth`, `TGpV`, `vpvs_similarity`.
- Elastic media (`elastic-iso`), with P and S data.
- Both TLOC misfits (absolute and double difference) and DWS output.

## Mathematical changes

None. The forward solver, adjoint solver, gradients and inversion are
mathematically identical to the main branch for the cases that remain.

## Execution changes

- Only three executables are built: `x_eikonal3`, `x_fatt3`, `x_tloc3`.
- `x_eikonal3` computes first-arrival traveltimes only. It no longer models
  reflections, traces rays or writes traveltime-gradient files.
- Parameters of removed features are now ignored. Step 1 turns them into
  errors.
- The build refuses CUDA and the legacy solver.

## Verification

The trimmed branch compiles with the Intel compiler and FLIT with no errors
or warnings.
