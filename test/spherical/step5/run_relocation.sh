#!/bin/bash
#
# TLOC relocation test of step 5 at 60 degrees latitude: synthetic data from
# the true events with x_eikonal3, then 12 iterations of x_tloc3 from moved
# events, with step limits of 2000 m in longitude, latitude and depth.
# check_relocation.py prints the moves and the location errors per
# iteration.
#
# Usage: bash run_relocation.sh <bindir> [workdir]
#
# Takes about 15 minutes with 16 OpenMP threads (OMP_NUM_THREADS=16).
#
# Needs the Intel MPI runtime (source setvars.sh) and Python with numpy
# (set PYTHON to choose the interpreter).

set -e
here=$(cd "$(dirname "$0")" && pwd)
bindir=$(cd "$1" && pwd)
work=${2:-$here/work_step5}
python=${PYTHON:-python}

mkdir -p "$work"
cd "$work"
cp "$here"/param_eikonal.rb "$here"/param_tloc.rb .
$python "$here/make_inputs.py"
"$bindir/x_eikonal3" param_eikonal.rb > log_eikonal.txt 2>&1
"$bindir/x_tloc3" param_tloc.rb > log_tloc.txt 2>&1
$python "$here/check_relocation.py"
