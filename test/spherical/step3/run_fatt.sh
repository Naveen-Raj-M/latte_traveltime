#!/bin/bash
#
# End-to-end FATT run of step 3: synthetic data from the true model with
# x_eikonal3, then 6 iterations of x_fatt3 from the starting model. Prints
# the data misfit of each iteration, which must decrease, and the result of
# a run with adjoint_receiver = nearest, which must stop with an error.
#
# Usage: bash run_fatt.sh <bindir> [workdir]
#
# Needs the Intel MPI runtime (source setvars.sh) and Python with numpy
# (set PYTHON to choose the interpreter).

set -e
here=$(cd "$(dirname "$0")" && pwd)
bindir=$(cd "$1" && pwd)
work=${2:-$here/work_fatt}
python=${PYTHON:-python}

mkdir -p "$work"
cd "$work"
cp "$here"/param_eikonal.rb "$here"/param_fatt.rb .
$python "$here/make_fatt_inputs.py"
"$bindir/x_eikonal3" param_eikonal.rb > log_eikonal.txt 2>&1
"$bindir/x_fatt3" param_fatt.rb > log_fatt.txt 2>&1
echo "Data misfit by iteration (fatt/data_misfit.txt):"
cat fatt/data_misfit.txt
echo "Run with adjoint_receiver = nearest (expect an error):"
sed 's/^dir_working = fatt/dir_working = fatt_nearest/' param_fatt.rb > param_nearest.rb
echo "adjoint_receiver = nearest" >> param_nearest.rb
"$bindir/x_fatt3" param_nearest.rb 2>&1 | grep -i "error" || true
