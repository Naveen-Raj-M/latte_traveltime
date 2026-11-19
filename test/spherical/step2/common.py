"""Shared helpers of the step-2 tests: the spherical geometry of
src/module_spherical.f90 reimplemented in numpy, and writers of models,
geometry and parameter files for x_eikonal3."""
import os
import subprocess
import numpy as np

R = 6371000.0
F = 1/298.257223563
E2 = F*(2 - F)


def geoc(latg):
    """Geocentric latitude (radians) of a WGS84 geodetic latitude (degrees)."""
    p = np.radians(latg)
    return np.arctan2((1 - E2)*np.sin(p), np.cos(p))


def ecef(lon, lat, dep):
    """Earth-centred position of (lon, geodetic lat, depth), last axis xyz."""
    r = R - np.asarray(dep, dtype=float)
    t = geoc(lat)
    l = np.radians(lon)
    return np.stack([r*np.cos(t)*np.cos(l), r*np.cos(t)*np.sin(l), r*np.sin(t)], -1)


def grid(olon, olat, odep, dlon, dlat, ddep, nlon, nlat, ndep):
    """Node coordinates with shape (nlon, nlat, ndepth); in C order this is
    the Fortran order (ndepth, nlat, nlon) of LATTE's binary files."""
    lon = olon + np.arange(nlon)*dlon
    lat = olat + np.arange(nlat)*dlat
    dep = odep + np.arange(ndep)*ddep
    return np.meshgrid(lon, lat, dep, indexing='ij')


def write_model(path, v):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.asarray(v, dtype=np.float32).tofile(path)


def write_geometry(directory, src, recs):
    """One shot with one source point (lon, lat, depth) and receivers."""
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, 'geometry.txt'), 'w') as f:
        f.write('shot_1.txt\n')
    with open(os.path.join(directory, 'shot_1.txt'), 'w') as f:
        f.write('1\n1\n%.6f %.6f %.3f 0.0\n' % tuple(src))
        f.write('%d\n' % len(recs))
        for r in recs:
            f.write('%.6f %.6f %.3f 1.0\n' % tuple(r))


def write_param(name, g, outdir, snapdir):
    """Spherical x_eikonal3 parameter file for grid g =
    (olon, olat, odepth, dlon, dlat, ddepth, nlon, nlat, ndepth)."""
    olon, olat, odep, dlon, dlat, ddep, nlon, nlat, ndep = g
    with open(name, 'w') as f:
        f.write(f"""nlon = {nlon}
nlat = {nlat}
ndepth = {ndep}
dlon = {dlon}
dlat = {dlat}
ddepth = {ddep}
olon = {olon}
olat = {olat}
odepth = {odep}
ns = 1
file_geometry = geometry/geometry.txt
model_name = vp
file_vp = model/vp.bin
snaps = 1
sweep_niter_max = 20
sweep_stop_threshold = 1e-6
dir_synthetic = {outdir}
dir_snapshot = {snapdir}
""")


def run(executable, param, log):
    """Run an executable on a parameter file in the current directory."""
    with open(log, 'w') as f:
        subprocess.run([executable, param], stdout=f, stderr=subprocess.STDOUT, check=True)


def read_field(snapdir, shape):
    return np.fromfile(os.path.join(snapdir, 'shot_1_traveltime_p.bin'), np.float32).reshape(shape)


def read_receivers(outdir):
    return np.fromfile(os.path.join(outdir, 'shot_1_traveltime_p.bin'), np.float32)
