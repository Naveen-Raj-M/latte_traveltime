#!/bin/bash
#
# End-to-end FATT of step 4, after run_conemute.sh in the same work
# directory (which makes the step-3 inputs and data): 6 iterations of
# x_fatt3 with per-shot and whole-model smoothing, a depth taper, and
# Tikhonov plus Gaussian model regularization. Prints the data misfit by
# iteration, which must decrease, the correlation of the recovered and true
# perturbations, and the result of a run with shot_grad_wavenum_lat, which
# must stop with an error.
#
# Usage: bash run_fatt_processing.sh <bindir> [workdir]

set -e
here=$(cd "$(dirname "$0")" && pwd)
bindir=$(cd "$1" && pwd)
work=${2:-$here/work_step4}
python=${PYTHON:-python}

cd "$work"
cp "$here"/param_fatt_processing.rb .
"$bindir/x_fatt3" param_fatt_processing.rb > log_fatt_processing.txt 2>&1
echo "Data misfit by iteration (fatt_processing/data_misfit.txt):"
cat fatt_processing/data_misfit.txt
grep "Tikhonov lambda" log_fatt_processing.txt | head -1 || true
$python - <<'PY'
import numpy as np
shape = (61, 61, 21)
vt = np.fromfile('model/vp_true.bin', np.float32).reshape(shape)
v0 = np.fromfile('model/vp_init.bin', np.float32).reshape(shape)
vu = np.fromfile('fatt_processing/iteration_6/model/updated_vp.bin', np.float32).reshape(shape)
dt = (vt - v0)/v0
du = (vu - v0)/v0
print('correlation of recovered and true perturbation, 4-20 km depth: %.3f'
      % np.corrcoef(dt[:, :, 2:11].ravel(), du[:, :, 2:11].ravel())[0, 1])
PY
echo "Run with shot_grad_wavenum_lat (expect an error):"
sed 's/^dir_working = fatt_processing/dir_working = fatt_wavenum_lat/' param_fatt_processing.rb > param_wavenum_lat.rb
printf 'process_shot_grad = wavenumfilt\nshot_grad_wavenum_lat = 0, 0.0001\nshot_grad_wavenum_lat_amps = 1, 0\n' >> param_wavenum_lat.rb
"$bindir/x_fatt3" param_wavenum_lat.rb 2>&1 | grep "Error" || true
