#!/bin/bash
#
# Parameter checks of step 1, with x_eikonal3:
#   1. param_bad.rb has 26 Cartesian or unsupported entries; the run must
#      report all 26 and stop.
#   2. param_good.rb uses the spherical names; the run must load the grid
#      and geometry, keep only the receiver 10 km above the source
#      (offset_max = 20 km), and set a shot range of about 1 degree around
#      the stations (50 km adaptive margins).
#   3. nlat = 1 must stop with the grid-size message.
#   4. A grid reaching past 90 degrees latitude must stop with the latitude
#      message.
#
# Usage: bash run_param_checks.sh <bindir> [workdir]
#
#   bindir   directory holding x_eikonal3 of the spherical branch
#   workdir  where the runs write (default: ./work_step1 next to this script)
#
# The runtime libraries of the Intel compilers and MPI must be available,
# e.g., after source /path/to/oneapi/setvars.sh. Needs Python with numpy.

set -e
here=$(cd "$(dirname "$0")" && pwd)
bindir=$(cd "$1" && pwd)
work=${2:-$here/work_step1}
python=${PYTHON:-python}

mkdir -p "$work"
cd "$work"
cp "$here"/param_bad.rb "$here"/param_good.rb .
$python "$here/make_inputs.py"

echo "=== 1. Cartesian and unsupported parameters (expect 26 errors) ==="
"$bindir/x_eikonal3" param_bad.rb 2>&1 | grep "Error" || true

echo "=== 2. Spherical parameters ==="
"$bindir/x_eikonal3" param_good.rb 2>&1 | grep "Earth radius\|lonmin\|latmin\|depthmin\|ns, nr\|Completed" || true

echo "=== 3. nlat = 1 (expect a grid-size error) ==="
sed 's/^nlat = 31/nlat = 1/' param_good.rb > param_nlat1.rb
"$bindir/x_eikonal3" param_nlat1.rb 2>&1 | grep -i "error" || true

echo "=== 4. Grid past 90 degrees latitude (expect a latitude error) ==="
sed 's/^olat = 20/olat = 87.2/' param_good.rb > param_pole.rb
"$bindir/x_eikonal3" param_pole.rb 2>&1 | grep -i "error" || true
