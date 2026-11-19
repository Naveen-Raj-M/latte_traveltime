"""Model and geometry for the parameter-check runs of step 1.

Writes, in the current directory:
  model/vp.bin             41 x 31 x 21 homogeneous model, 6000 m/s
  geometry/geometry.txt    one shot
  geometry/shot_1.txt      a source at (102, 22.5, 10 km) and three surface
                           receivers: one 10 km above the source and two more
                           than 100 km away
"""
import os
import numpy as np

nlon, nlat, ndep = 41, 31, 21
os.makedirs('model', exist_ok=True)
os.makedirs('geometry', exist_ok=True)
# C order (nlon, nlat, ndepth) is the Fortran order (ndepth, nlat, nlon) of LATTE
np.full((nlon, nlat, ndep), 6000.0, dtype=np.float32).tofile('model/vp.bin')
with open('geometry/geometry.txt', 'w') as f:
    f.write('shot_1.txt\n')
with open('geometry/shot_1.txt', 'w') as f:
    f.write('1\n1\n')
    f.write('102.0 22.5 10000.0 0.0\n')
    f.write('3\n')
    f.write('101.0 22.0 0.0 1.0\n103.0 23.0 0.0 1.0\n102.0 22.5 0.0 1.0\n')
