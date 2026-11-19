We currently only provide traditional Makefile-based approach to compile the spherical branch of `LATTE`:

* To compile, set the path to `FLIT` with `flitdir = ...` in [Makefile](3d/Makefile), or run `ruby install.rb flitdir=...`; the default path to `FLIT` is `$HOME/flit`. `ruby install.rb` passes each argument `name=value` to `make`, which overrides the variable `name` in the Makefile.
* The build makes `x_eikonal3`, `x_fatt3` and `x_tloc3` in spherical coordinates. The spherical branch has no 2D executables, no TRTT, no CUDA executables and no legacy solver; `use_cuda=on` or `-Dlegacy_solver` stops the build with an error.
* The design and the changes of the spherical branch are documented step by step in [docs](../docs).
