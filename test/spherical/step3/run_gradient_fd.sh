#!/bin/bash
#
# Finite-difference check of the spherical adjoint-state gradient; see the
# header of test_gradient_fd.f90.
#
# Usage: bash run_gradient_fd.sh [objdir] [flitdir]
#
#   objdir   build directory holding the objects and modules of the spherical
#            branch (default: $HOME/bin/obj/latte3, the Makefile default)
#   flitdir  FLIT directory (default: $HOME/flit)
#
# The Intel compilers must be on PATH (source setvars.sh).

set -e
here=$(cd "$(dirname "$0")" && pwd)
objdir=${1:-$HOME/bin/obj/latte3}
flitdir=${2:-$HOME/flit}

# The link flags of FLIT, as its Makefile.in sets them
lflags=$(printf 'include $(flitdir)/src/Makefile.in\nprint:\n\t@echo $(base_lflags)\n' | make -s -f - flitdir="$flitdir" print)

cd "$here"
mpiifx -O2 -fpp -qmkl -qopenmp -heap-arrays -I"$objdir" -I"$flitdir/lib" test_gradient_fd.f90 \
    "$objdir/module_spherical.o" "$objdir/module_parameters.o" "$objdir/module_vars.o" \
    "$objdir/module_utility.o" "$objdir/module_traveltime_iso.o" \
    "$flitdir/lib/libflit.a" $lflags -o test_gradient_fd
./test_gradient_fd 2>&1 | grep -v "Fast sweeping\|Wall-clock"
