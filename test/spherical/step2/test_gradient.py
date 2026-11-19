"""Linear velocity gradient, grid refinement.

Velocity grows linearly along the downward direction n at the centre
X0 = (102, 22, 0) of a 4 x 4 degree x 60 km grid:

    v(X) = 5000 + 0.02 (X - X0) . n    (m/s)

This medium has the exact traveltime

    T = (1/G) arccosh( 1 + G^2 |X - X_s|^2 / (2 v(X) v(X_s)) )

Rays are circular arcs. For distant pairs of points the exact ray dips
below the 60 km floor of the grid, which the solver cannot follow, so the
error statistics count only nodes whose exact ray stays above 50 km. The
error there must fall with each refinement.

Usage: python test_gradient.py --bin <spherical bindir> [--levels N] [--workdir DIR]

--levels 4 (default) includes a 5-million-node grid; use 3 for a quick run.
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

ap = argparse.ArgumentParser()
ap.add_argument('--bin', required=True)
ap.add_argument('--levels', type=int, default=4)
ap.add_argument('--workdir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work_gradient'))
a = ap.parse_args()
sph = os.path.join(os.path.abspath(a.bin), 'x_eikonal3')
os.makedirs(a.workdir, exist_ok=True)
os.chdir(a.workdir)

G = 0.02
v0 = 5000.0
X0 = ecef(102.0, 22.0, 0.0)
n = -X0/np.linalg.norm(X0)
src = (101.7313, 22.4127, 12345.0)
Xs = ecef(*src)


def vel(X):
    return v0 + G*np.tensordot(X - X0, n, axes=([-1], [0]))


def texact(X):
    d2 = np.sum((X - Xs)**2, -1)
    return np.arccosh(1 + G**2*d2/(2*vel(X)*vel(Xs)))/G


def ray_bottom(X):
    """Largest distance below X0 along n reached by the exact ray (a circular
    arc centred where v = 0) between the source and X."""
    d = X - Xs
    sp = np.dot(Xs - X0, n)
    sq = np.tensordot(X - X0, n, axes=([-1], [0]))
    w = d - np.multiply.outer(np.tensordot(d, n, axes=([-1], [0])), n)
    aa = np.maximum(np.linalg.norm(w, axis=-1), 1e-6)
    sc = -v0/G
    xc = (aa**2 + (sq - sc)**2 - (sp - sc)**2)/(2*aa)
    rho = np.sqrt(xc**2 + (sp - sc)**2)
    return np.where((xc >= 0) & (xc <= aa), sc + rho, np.maximum(sp, sq))


rng = np.random.default_rng(3)
recs = np.column_stack([rng.uniform(100.1, 103.9, 40), rng.uniform(20.1, 23.9, 40),
                        np.r_[np.zeros(20), rng.uniform(0, 55000, 20)]])
levels = [(0.1, 10000.0), (0.05, 5000.0), (0.025, 2500.0), (0.0125, 1250.0)][:a.levels]
prev = None
print('grid                  nodes      inside: mean|e|   max|e|    order(max)   ray-below-floor: max|e|   receivers: max|e|')
for dl, dd in levels:
    nl = int(round(4/dl)) + 1
    nd = int(round(60000/dd)) + 1
    g = (100.0, 20.0, 0.0, dl, dl, dd, nl, nl, nd)
    L, A, D = grid(*g)
    X = ecef(L, A, D)
    write_model('model/vp.bin', vel(X))
    write_geometry('geometry', src, recs)
    write_param('param.rb', g, 'out', 'snap')
    run(sph, 'param.rb', 'log_%g.txt' % dl)
    e = np.abs(read_field('snap', L.shape) - texact(X))
    inside = ray_bottom(X) < 50000
    er = np.abs(read_receivers('out') - texact(ecef(recs[:, 0], recs[:, 1], recs[:, 2])))
    order = '' if prev is None else '%.2f' % np.log2(prev/e[inside].max())
    print('%.4f deg x %5.0f m %9d   %.3e s  %.3e s  %6s       %.3f s                   %.3e s'
          % (dl, dd, L.size, e[inside].mean(), e[inside].max(), order, e[~inside].max(), er.max()))
    prev = e[inside].max()
