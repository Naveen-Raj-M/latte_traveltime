#!/bin/bash
#
# Unit test of src/module_spherical.f90: latitude conversion, Earth-centred
# positions, exact cell volumes, arc lengths, reference spacing and the
# smallest length of a degree.
#
# Usage: bash run_unit_test.sh [objdir] [flitdir]
#
#   objdir   directory holding module_spherical.o and spherical.mod from the
#            build (default: $HOME/bin/obj/latte3, the Makefile default)
#   flitdir  FLIT directory (default: $HOME/flit)
#
# The Intel compilers must be on PATH, e.g., after
#   source /path/to/oneapi/setvars.sh

set -e
here=$(cd "$(dirname "$0")" && pwd)
objdir=${1:-$HOME/bin/obj/latte3}
flitdir=${2:-$HOME/flit}

cd "$here"
mpiifx -O2 -fpp -qmkl -qopenmp -I"$objdir" -I"$flitdir/lib" test_spherical.f90 \
    "$objdir/module_spherical.o" "$flitdir/lib/libflit.a" -lstdc++ -o test_spherical
./test_spherical
