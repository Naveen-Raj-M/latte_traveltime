#!/bin/bash
#
# Cone mute of step 4: one shot of the step-3 FATT inputs, one iteration of
# x_fatt3 with and without process_shot_grad = conemute; check_conemute.py
# compares the ratio of the two gradients with the exact mute weight.
#
# Usage: bash run_conemute.sh <bindir> [workdir]
#
# Needs the Intel MPI runtime (source setvars.sh) and Python with numpy
# (set PYTHON to choose the interpreter).

set -e
here=$(cd "$(dirname "$0")" && pwd)
bindir=$(cd "$1" && pwd)
work=${2:-$here/work_step4}
python=${PYTHON:-python}

mkdir -p "$work"
cd "$work"
cp "$here"/../step3/param_eikonal.rb .
$python "$here/../step3/make_fatt_inputs.py"
"$bindir/x_eikonal3" param_eikonal.rb > log_eikonal.txt 2>&1
common="nlon = 61
nlat = 61
ndepth = 21
dlon = 0.05
dlat = 0.05
ddepth = 2000
olon = 100
olat = 20
ns = 1
file_geometry = geometry/geometry.txt
model_update = vp
file_vp = model/vp_init.bin
dir_record = data
niter_max = 1"
printf '%s\ndir_working = fatt_nomute\n' "$common" > param_nomute.rb
printf '%s\ndir_working = fatt_mute\nprocess_shot_grad = conemute\nshot_grad_conemuteorigin = receiver\nshot_grad_conemuter = 60000\nshot_grad_conemute_depth = 20000\nshot_grad_conemutepower = 1\nshot_grad_conemutetaper = 20000\n' "$common" > param_mute.rb
"$bindir/x_fatt3" param_nomute.rb > log_nomute.txt 2>&1
"$bindir/x_fatt3" param_mute.rb > log_mute.txt 2>&1
$python "$here/check_conemute.py"
