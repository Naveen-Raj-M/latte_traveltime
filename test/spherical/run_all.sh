#!/bin/bash
#
# Run all spherical-branch tests.
#
# Usage: bash run_all.sh <bindir> <objdir> [main_bindir]
#
#   bindir       directory holding x_eikonal3 of the spherical branch
#   objdir       build directory holding module_spherical.o and spherical.mod
#   main_bindir  optional: directory holding x_eikonal3 of the main branch,
#                for the Cartesian comparison of the homogeneous test
#
# Needs the Intel compilers and MPI runtime on PATH (source setvars.sh) and
# Python with numpy (set PYTHON to choose the interpreter).

set -e
here=$(cd "$(dirname "$0")" && pwd)
bindir=$(cd "$1" && pwd)
objdir=$(cd "$2" && pwd)
python=${PYTHON:-python}

echo "################ Step 1: geometry unit test ################"
bash "$here/step1/run_unit_test.sh" "$objdir"

echo "################ Step 1: parameter checks ################"
PYTHON=$python bash "$here/step1/run_param_checks.sh" "$bindir"

echo "################ Step 2: homogeneous Earth ################"
if [ -n "$3" ]; then
    $python "$here/step2/test_homogeneous.py" --bin "$bindir" --main-bin "$3"
else
    $python "$here/step2/test_homogeneous.py" --bin "$bindir"
fi

echo "################ Step 2: linear gradient, grid refinement ################"
$python "$here/step2/test_gradient.py" --bin "$bindir"

echo "################ Step 2: velocity depending on depth ################"
$python "$here/step2/test_radial.py" --bin "$bindir"

echo "################ Step 3: finite-difference check of the adjoint gradient ################"
bash "$here/step3/run_gradient_fd.sh" "$objdir"

echo "################ Step 3: end-to-end FATT and TLOC ################"
PYTHON=$python bash "$here/step3/run_fatt.sh" "$bindir"
PYTHON=$python bash "$here/step3/run_tloc.sh" "$bindir"

echo "################ Step 4: exact-distance filters ################"
bash "$here/step4/run_filters.sh" "$objdir"

echo "################ Step 4: cone mute and FATT with processing ################"
PYTHON=$python bash "$here/step4/run_conemute.sh" "$bindir"
PYTHON=$python bash "$here/step4/run_fatt_processing.sh" "$bindir"

echo "################ Step 5: TLOC relocation at 60 degrees latitude ################"
PYTHON=$python bash "$here/step5/run_relocation.sh" "$bindir"
