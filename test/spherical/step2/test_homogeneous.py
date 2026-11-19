"""Homogeneous Earth (6000 m/s), 10 x 8 degrees, off-node source at 15.3 km.

The exact traveltime is the straight-line distance between Earth-centred
positions over the velocity. Prints the error of the spherical solver on two
grids: 0.1 x 0.1 degree x 2 km cells and 0.1 x 0.1 degree x 10 km cells.

With --main-bin, also runs the Cartesian x_eikonal3 of the main branch on
the matching flat grids (10 x 10 x 2 km and 10 km cubes) for comparison.
Both codes are inexact below the source, where the local solver falls back
to its 1D and 2D updates; docs/spherical_02_forward_eikonal.md explains why.

Usage: python test_homogeneous.py --bin <spherical bindir> [--main-bin <main bindir>] [--workdir DIR]
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

ap = argparse.ArgumentParser()
ap.add_argument('--bin', required=True)
ap.add_argument('--main-bin', default=None)
ap.add_argument('--workdir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work_homogeneous'))
a = ap.parse_args()
sph = os.path.join(os.path.abspath(a.bin), 'x_eikonal3')
os.makedirs(a.workdir, exist_ok=True)
os.chdir(a.workdir)

src = (104.337, 23.713, 15321.0)
v = 6000.0


def stats(e):
    e = np.abs(e)
    return 'max |e| = %.4f s, mean |e| = %.2e s, 99th percentile = %.4f s' % (e.max(), e.mean(), np.percentile(e, 99))


print('Spherical solver, homogeneous Earth, exact T = |X - X_s|/v')
for label, ddep, ndep in [('0.1 x 0.1 deg x 2 km ', 2000.0, 41), ('0.1 x 0.1 deg x 10 km', 10000.0, 11)]:
    g = (100.0, 20.0, 0.0, 0.1, 0.1, ddep, 101, 81, ndep)
    L, A, D = grid(*g)
    write_model('model/vp.bin', np.full(L.shape, v))
    write_geometry('geometry', src, [(105.0, 24.0, 0.0)])
    write_param('param.rb', g, 'out', 'snap')
    run(sph, 'param.rb', 'log_sph_%d.txt' % ddep)
    Te = np.linalg.norm(ecef(L, A, D) - ecef(*src), axis=-1)/v
    e = read_field('snap', L.shape) - Te
    i, j, k = np.unravel_index(np.abs(e).argmax(), e.shape)
    print('  %s: %s; worst node at (%.2f, %.2f, %.0f m)' % (label, stats(e), L[i, j, k], A[i, j, k], D[i, j, k]))

if a.main_bin:
    cart = os.path.join(os.path.abspath(a.main_bin), 'x_eikonal3')
    print('Cartesian solver of the main branch, matching flat grids, exact T = |x - x_s|/v')
    xs = (433700.0, 371300.0, 15321.0)
    for label, dz, nz in [('10 x 10 x 2 km ', 2000.0, 41), ('10 x 10 x 10 km', 10000.0, 11)]:
        nx, ny = 101, 81
        os.makedirs('cart/model', exist_ok=True)
        os.makedirs('cart/geometry', exist_ok=True)
        np.full((nx, ny, nz), v, np.float32).tofile('cart/model/vp.bin')
        with open('cart/geometry/geometry.txt', 'w') as f:
            f.write('shot_1.txt\n')
        with open('cart/geometry/shot_1.txt', 'w') as f:
            f.write('1\n1\n%.1f %.1f %.1f 0.0\n1\n500000.0 400000.0 0.0 1.0\n' % xs)
        with open('cart/param.rb', 'w') as f:
            f.write(f"nx = {nx}\nny = {ny}\nnz = {nz}\ndx = 10000\ndy = 10000\ndz = {dz}\nns = 1\n"
                    "file_geometry = geometry/geometry.txt\nmodel_name = vp\nfile_vp = model/vp.bin\n"
                    "snaps = 1\ndir_synthetic = out\ndir_snapshot = snap\nsweep_niter_max = 20\n"
                    "sweep_stop_threshold = 1e-6\n")
        cwd = os.getcwd()
        os.chdir('cart')
        run(cart, 'param.rb', 'log_cart_%d.txt' % dz)
        os.chdir(cwd)
        X, Y, Z = np.meshgrid(np.arange(nx)*1e4, np.arange(ny)*1e4, np.arange(nz)*dz, indexing='ij')
        Tc = np.sqrt((X - xs[0])**2 + (Y - xs[1])**2 + (Z - xs[2])**2)/v
        e = read_field('cart/snap', X.shape) - Tc
        print('  %s: %s' % (label, stats(e)))
else:
    print('Cartesian comparison skipped: pass --main-bin with the bin directory of a main-branch build to run it.')
