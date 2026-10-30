
nx = 301
nz = 301
dx = 10
dz = 10

ns = 40
file_geometry = geometry/geometry_no_t0.txt

step_max_sx = 1:5, 10:5
step_max_sz = 1:5, 10:5
step_max_st0 = 1:0.5, 10:0.1

min_sx = 400
max_sx = 2600
min_sz = 400
max_sz = 2600
min_st0 = 0
max_st0 = 3

model_update = sx, sz, st0
file_vp = model/vp.bin
file_vs = model/vs.bin
file_sx = model/sx_init.bin
file_sz = model/sz_init.bin
#file_st0 = model/st0_init.bin

misfit_type = ad
niter_max = 200

# data masked - p only
dir_record = data_masked_no_t0
which_medium = acoustic-iso
model_aux = vp
dir_working = reloc_ad_p

yn_continue = y
