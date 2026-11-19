"""Velocity depending only on depth, against exact spherical traveltimes.

v = 5800 + 0.025 d (m/s) for depth d, on a 12 x 5 degree x 200 km grid of
0.05 x 0.05 degree x 5 km cells. Source at (99, 20.5, 10 km); surface
receivers every 0.5 degree along its latitude.

Exact spherical first arrivals come from a ray-parameter integral in the
Earth-flattened model, which is exact for a velocity that depends only on
depth. The same integral without flattening gives the flat-Earth times that
a Cartesian code would compute. Rows whose exact ray bottoms below the
200 km floor of the grid test the grid's depth, not the solver, and are
marked.

Usage: python test_radial.py --bin <spherical bindir> [--workdir DIR]
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

ap = argparse.ArgumentParser()
ap.add_argument('--bin', required=True)
ap.add_argument('--workdir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work_radial'))
a = ap.parse_args()
sph = os.path.join(os.path.abspath(a.bin), 'x_eikonal3')
os.makedirs(a.workdir, exist_ok=True)
os.chdir(a.workdir)

hs = 10000.0
floor = 200000.0


def v_of_d(d):
    return 5800.0 + 0.025*d


def ray_branches(spherical):
    """The two ray families from the source, each as rows of (distance, time,
    bottoming depth) ordered by ray parameter: the direct rays going up to
    the surface, and the turning rays going down, turning, and coming up."""
    dz = 20.0
    edges = np.arange(0, 400000 + dz, dz)
    dmid = 0.5*(edges[1:] + edges[:-1])
    if spherical:
        r = R - dmid
        h = np.diff(R*np.log(R/(R - edges)))
        u = r/(R*v_of_d(dmid))
    else:
        h = np.diff(edges)
        u = 1/v_of_d(dmid)
    above = dmid < hs
    us = u[above].min()
    direct = []
    turning = []
    for p in np.linspace(1e-7, us*0.999999, 20000):
        ok = u > p
        q = np.sqrt(np.where(ok, u**2 - p**2, 1))
        dx = np.where(ok, h*p/q, 0)
        dt = np.where(ok, h*u**2/q, 0)
        xu = dx[above].sum()
        tu = dt[above].sum()
        direct.append((xu, tu, 0.0))
        idx = np.where((~above) & (u <= p))[0]
        if idx.size:
            seg = (~above).copy()
            seg[idx[0]:] = False
            turning.append((xu + 2*dx[seg].sum(), tu + 2*dt[seg].sum(), dmid[idx[0]]))
    return [np.array(direct), np.array(turning)]


def first_arrival(branches, x):
    """Earliest time among the ray segments of either family that cross
    distance x, and the bottoming depth of that ray."""
    best = (np.inf, 0.0)
    for b in branches:
        X, T, B = b[:, 0], b[:, 1], b[:, 2]
        for i in range(len(X) - 1):
            if (X[i] - x)*(X[i + 1] - x) <= 0 and X[i] != X[i + 1]:
                t = T[i] + (T[i + 1] - T[i])*(x - X[i])/(X[i + 1] - X[i])
                if t < best[0]:
                    best = (t, max(B[i], B[i + 1]))
    return best


g = (98.0, 18.0, 0.0, 0.05, 0.05, 5000.0, 241, 101, 41)
L, A, D = grid(*g)
write_model('model/vp.bin', v_of_d(D))
src = (99.0, 20.5, hs)
recs = np.array([(99.0 + x, 20.5, 0.0) for x in np.arange(0.5, 9.51, 0.5)])
write_geometry('geometry', src, recs)
write_param('param.rb', g, 'out', 'snap')
run(sph, 'param.rb', 'log.txt')
tr = read_receivers('out')

Xs0 = ecef(src[0], src[1], 0.0)
Xr = ecef(recs[:, 0], recs[:, 1], recs[:, 2])
arc = R*np.arccos(np.clip(Xr @ Xs0/(np.linalg.norm(Xr, axis=1)*np.linalg.norm(Xs0)), -1, 1))
bs = ray_branches(True)
bf = ray_branches(False)
print(' distance   solver     exact      solver-exact   flat Earth   flat-exact   ray bottom')
for x, t in zip(arc, tr):
    te, bottom = first_arrival(bs, x)
    tf, _ = first_arrival(bf, x)
    note = '' if bottom <= floor else '  below grid floor, not a solver test'
    print('%6.0f km  %8.3f s  %8.3f s  %8.3f s     %8.3f s   %8.3f s   %5.0f km%s'
          % (x/1000, t, te, t - te, tf, tf - te, bottom/1000, note))
