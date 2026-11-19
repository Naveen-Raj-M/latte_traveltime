"""Inputs of the small end-to-end FATT run of step 3.

Grid: longitude 100 to 103, geodetic latitude 20 to 23 degrees, depth 0 to
40 km; 0.05 degree x 2 km cells. True model: 5000 + 0.06 d m/s at depth d,
with a 8% fast and a 8% slow Gaussian anomaly of 30 km radius at 12 km depth.
Starting model: the same without the anomalies. 9 sources at 15 km depth,
each recorded by 49 surface stations.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'step2'))
from common import grid, ecef

g = (100.0, 20.0, 0.0, 0.05, 0.05, 2000.0, 61, 61, 21)
L, A, D = grid(*g)
v0 = 5000.0 + 0.06*D
X = ecef(L, A, D)
anom = np.zeros_like(v0)
for lon, lat, sign in [(101.0, 21.0, 1.0), (102.0, 22.0, -1.0)]:
    d2 = np.sum((X - ecef(lon, lat, 12000.0))**2, -1)
    anom += sign*0.08*np.exp(-d2/30000.0**2)
os.makedirs('model', exist_ok=True)
(v0*(1 + anom)).astype(np.float32).tofile('model/vp_true.bin')
v0.astype(np.float32).tofile('model/vp_init.bin')

os.makedirs('geometry', exist_ok=True)
srcs = [(100.5 + 1.0*i, 20.5 + 1.0*j, 15000.0) for i in range(3) for j in range(3)]
sta = [(100.25 + 0.4167*i, 20.25 + 0.4167*j, 0.0) for i in range(7) for j in range(7)]
with open('geometry/geometry.txt', 'w') as f:
    for s in range(len(srcs)):
        f.write('shot_%d.txt\n' % (s + 1))
for s, src in enumerate(srcs):
    with open('geometry/shot_%d.txt' % (s + 1), 'w') as f:
        f.write('%d\n1\n%.4f %.4f %.1f 0.0\n%d\n' % ((s + 1,) + src + (len(sta),)))
        for r in sta:
            f.write('%.4f %.4f %.1f 1.0\n' % r)
