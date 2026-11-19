"""Inputs of the TLOC relocation test of step 5.

Grid: longitude 10 to 14, geodetic latitude 58 to 62 degrees, depth 0 to
40 km; 0.04 degree x 0.02 degree x 1 km cells (about 2.2 km x 2.2 km x 1 km
at 60 degrees). Velocity 5800 + 0.03 d m/s at depth d.

16 true events at 8 to 25 km depth, origin time 0; 49 surface stations.
geometry_true: the events at their true positions, for the synthetic data.
geometry_init: the same events moved by 3 to 6 km east and north (in
opposite senses for alternate events) and 3 km in depth, origin time 0.
model/slon_init.bin, slat_init.bin and sdepth_init.bin hold these moved
positions, the starting positions of TLOC, and model/st0_init.bin the
starting origin times, the true 0. events_true.txt lists the true
positions.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'step2'))
from common import grid, geoc, R

g = (10.0, 58.0, 0.0, 0.04, 0.02, 1000.0, 101, 201, 41)
L, A, D = grid(*g)
os.makedirs('model', exist_ok=True)
(5800.0 + 0.03*D).astype(np.float32).tofile('model/vp.bin')

rng = np.random.default_rng(5)
ev = np.column_stack([rng.uniform(11.0, 13.0, 16), rng.uniform(59.0, 61.0, 16), rng.uniform(8000, 25000, 16)])
sta = [(10.4 + 0.533*i, 58.4 + 0.533*j, 0.0) for i in range(7) for j in range(7)]
np.savetxt('events_true.txt', ev, fmt='%.6f %.6f %.1f')

# Start: move east by de and north by dn meters, and down by dz
de = np.where(np.arange(16) % 2 == 0, 1, -1)*rng.uniform(3000, 6000, 16)
dn = np.where(np.arange(16) % 3 == 0, -1, 1)*rng.uniform(3000, 6000, 16)
dz = np.where(np.arange(16) % 2 == 0, -1, 1)*3000.0
r = R - ev[:, 2]
init = ev.copy()
init[:, 0] += np.degrees(de/(r*np.cos(geoc(ev[:, 1]))))
# north: move the geocentric latitude by the arc, then back to geodetic
F = 1/298.257223563
E2 = F*(2 - F)
t = geoc(ev[:, 1]) + dn/r
init[:, 1] = np.degrees(np.arctan2(np.sin(t), (1 - E2)*np.cos(t)))
init[:, 2] += dz


def write(dirname, events):
    os.makedirs(dirname, exist_ok=True)
    with open(dirname + '/geometry.txt', 'w') as f:
        for e in range(len(events)):
            f.write('shot_%d.txt\n' % (e + 1))
    for e, p in enumerate(events):
        with open(dirname + '/shot_%d.txt' % (e + 1), 'w') as f:
            f.write('%d\n1\n%.6f %.6f %.1f 0.0\n%d\n' % ((e + 1,) + tuple(p) + (len(sta),)))
            for s in sta:
                f.write('%.6f %.6f %.1f 1.0\n' % s)


write('geometry_true', ev)
write('geometry_init', init)

# Starting positions for TLOC, as in the Cartesian examples (file_sx, ...):
# one value per event, in the order of the geometry
init[:, 0].astype(np.float32).tofile('model/slon_init.bin')
init[:, 1].astype(np.float32).tofile('model/slat_init.bin')
init[:, 2].astype(np.float32).tofile('model/sdepth_init.bin')
# Starting origin times: the true 0
np.zeros(len(init), dtype=np.float32).tofile('model/st0_init.bin')
