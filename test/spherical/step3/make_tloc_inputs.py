"""Inputs of the small TLOC run of step 3, made after make_fatt_inputs.py
and the x_eikonal3 run of run_fatt.sh (which give the model and the observed
data, with origin times 0). geometry_tloc holds the same events and stations
with origin times of 0.3 s, which TLOC must bring back toward 0.
"""
import os
import shutil

os.makedirs('geometry_tloc', exist_ok=True)
shutil.copy('geometry/geometry.txt', 'geometry_tloc/geometry.txt')
for s in range(1, 10):
    lines = open('geometry/shot_%d.txt' % s).read().splitlines()
    parts = lines[2].split()
    parts[3] = '0.3'
    lines[2] = ' '.join(parts)
    open('geometry_tloc/shot_%d.txt' % s, 'w').write('\n'.join(lines) + '\n')
