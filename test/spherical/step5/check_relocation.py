"""Check the TLOC relocation of step 5.

For each iteration with updated event positions, prints the misfit, the
largest exact move of any event in longitude (along its parallel) and in
latitude (along its meridian) in meters, and the median and largest distance
of the events from their true positions (horizontal great-circle arc at the
event depth, and depth). Moves must not exceed the 2000 m step limits.
"""
import glob
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'step2'))
from common import geoc, ecef, R

true = np.loadtxt('events_true.txt')
init = np.array([[float(x) for x in open('geometry_init/shot_%d.txt' % (e + 1)).read().split('\n')[2].split()[:3]]
                 for e in range(len(true))])


def horizontal_error(p):
    u = ecef(p[:, 0], p[:, 1], 0.0*p[:, 2])
    v = ecef(true[:, 0], true[:, 1], 0.0*true[:, 2])
    ang = np.arctan2(np.linalg.norm(np.cross(u, v), axis=1), np.sum(u*v, 1))
    return (R - true[:, 2])*ang


iters = sorted(int(d.split('_')[-1]) for d in glob.glob('tloc/iteration_*')
               if os.path.exists(d + '/model/updated_slon.bin'))
misfit = np.atleast_2d(np.loadtxt('tloc/data_misfit.txt'))
prev = init.copy()
print('iter  misfit(rel)  max move lon(m)  max move lat(m)  max move depth(m)  horiz err median/max (m)  depth err median/max (m)')
he = horizontal_error(init)
de = np.abs(init[:, 2] - true[:, 2])
print('%4d  %10.4f %16s %16s %18s %12.0f %8.0f %13.0f %8.0f' % (0, 1.0, '', '', '', np.median(he), he.max(), np.median(de), de.max()))
mmax = 0.0
for it in iters:
    p = np.column_stack([np.fromfile('tloc/iteration_%d/model/updated_%s.bin' % (it, n), np.float32)
                         for n in ['slon', 'slat', 'sdepth']]).astype(float)
    r = R - prev[:, 2]
    mlon = np.max(r*np.cos(geoc(prev[:, 1]))*np.abs(np.radians(p[:, 0] - prev[:, 0])))
    mlat = np.max(r*np.abs(geoc(p[:, 1]) - geoc(prev[:, 1])))
    mdep = np.max(np.abs(p[:, 2] - prev[:, 2]))
    mmax = max(mmax, mlon, mlat)
    he = horizontal_error(p)
    de = np.abs(p[:, 2] - true[:, 2])
    rel = misfit[misfit[:, 0] == it, 2][0] if np.any(misfit[:, 0] == it) else np.nan
    print('%4d  %10.4f %16.1f %16.1f %18.1f %12.0f %8.0f %13.0f %8.0f'
          % (it, rel, mlon, mlat, mdep, np.median(he), he.max(), np.median(de), de.max()))
    prev = p
print('largest horizontal move of any event in any iteration: %.1f m (limit 2000 m)' % mmax)
