
nx = 301
nz = 301
dx = 10
dz = 10

ns = 40
file_geometry = geometry/geometry_no_t0.txt


step_max_sx = 1:500, 10:100
step_max_sz = 1:500, 10:100

min_sx = 400
max_sx = 2600
min_sz = 400
max_sz = 2600

model_update = sx, sz
file_vp = model/vp.bin
file_vs = model/vs.bin
file_sx = model/sx_init.bin
file_sz = model/sz_init.bin

misfit_type = dd

niter_max = 30


## data masked - p only
#dir_record = data_masked
#which_medium = acoustic-iso
#model_aux = vp
#dir_working = reloc_dd_p


## data masked - p + s
#dir_record = data_masked
#which_medium = elastic-iso
#model_aux = vp, vs
#dir_working = reloc_dd_ps


## data masked more - p only
#dir_record = data_masked_more
#which_medium = acoustic-iso
#model_aux = vp
#dir_working = reloc_dd_p_more


# data masked - p + s
dir_record = data_masked_more
which_medium = elastic-iso
model_aux = vp, vs
dir_working = reloc_dd_ps_more

