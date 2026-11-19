#!/bin/bash
#
# TLOC run of step 3, after run_fatt.sh in the same work directory: joint
# update of Vp and the origin times of 9 events whose origin times start at
# 0.3 s instead of 0. Prints the data misfit by iteration, which must
# decrease, and the final origin times.
#
# Usage: bash run_tloc.sh <bindir> [workdir]

set -e
here=$(cd "$(dirname "$0")" && pwd)
bindir=$(cd "$1" && pwd)
work=${2:-$here/work_fatt}
python=${PYTHON:-python}

cd "$work"
cp "$here"/param_tloc.rb .
$python "$here/make_tloc_inputs.py"
"$bindir/x_tloc3" param_tloc.rb > log_tloc.txt 2>&1
echo "Data misfit by iteration (tloc/data_misfit.txt):"
cat tloc/data_misfit.txt
last=$(dirname $(dirname $(ls tloc/iteration_*/model/updated_st0.bin | sort -t_ -k2 -n | tail -1)))
echo "Origin times after the last iteration ($last/model/updated_st0.bin):"
$python -c "import numpy as np; print(np.round(np.fromfile('$last/model/updated_st0.bin', np.float32), 3))"
