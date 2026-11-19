"""Compare the cone-mute weight that x_fatt3 applied with the exact one.

Reads iteration_1/model/grad_vp.bin of the runs fatt_nomute and fatt_mute
(one shot, one iteration, the second with process_shot_grad = conemute),
takes their ratio where the unmuted gradient is not negligible, and
compares it with the exact weight: 1 within the radius
R(d) = conemuter * (min(d - d0, conemute_depth)/conemute_depth)^(1/power)
of the origin column, a linear ramp to 0 over conemutetaper, 0 beyond, and
0 above the origin depth d0. The horizontal distance is the great-circle
angle between the node's column and the Earth-centred centroid of the
receivers, times the node's radius.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'step2'))
from common import grid, ecef, geoc, R

g = (100.0, 20.0, 0.0, 0.05, 0.05, 2000.0, 61, 61, 21)
L, A, D = grid(*g)
shape = L.shape
g0 = np.fromfile('fatt_nomute/iteration_1/model/grad_vp.bin', np.float32).reshape(shape)
g1 = np.fromfile('fatt_mute/iteration_1/model/grad_vp.bin', np.float32).reshape(shape)

# Receivers of shot 1 (origin = receiver: centroid of all receivers, the shallowest depth)
lines = open('geometry/shot_1.txt').read().split('\n')
nr = int(lines[3])
rec = np.array([[float(x) for x in lines[4 + i].split()[:3]] for i in range(nr)])
xc = ecef(rec[:, 0], rec[:, 1], rec[:, 2]).sum(0)
uc = xc/np.linalg.norm(xc)
d0 = rec[:, 2].min()
conemuter, conemutez, power, taper = 60000.0, 20000.0, 1.0, 20000.0

t = geoc(A)
l = np.radians(L)
un = np.stack([np.cos(t)*np.cos(l), np.cos(t)*np.sin(l), np.sin(t)], -1)
ang = np.arctan2(np.linalg.norm(np.cross(un, uc), axis=-1), un @ uc)
dist = (R - D)*ang
radius = conemuter*(np.minimum(D - d0, conemutez)/conemutez)**(1/power)
w = np.where(dist <= radius, 1.0, np.where(dist >= radius + taper, 0.0, 1 - (dist - radius)/taper))
w = np.where(D < d0, 0.0, w)

mask = np.abs(g0) > 1e-3*np.abs(g0).max()
ratio = g1[mask]/g0[mask]
print('nodes compared: %d' % mask.sum())
print('max |applied weight - exact weight| = %.2e' % np.abs(ratio - w[mask]).max())
print('nodes with weight 1: %d, with weight 0: %d, on the ramp: %d'
      % ((w[mask] == 1).sum(), (w[mask] == 0).sum(), ((w[mask] > 0) & (w[mask] < 1)).sum()))
