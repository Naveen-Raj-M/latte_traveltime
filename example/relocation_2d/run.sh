
bindir=$HOME/src/latte/bin

export OMP_NUM_THREADS=4

## forward modeling
#mpirun -np 10 $bindir/x_eikonal2 param_eikonal.rb

## dd tloc
#mpirun -np 10 $bindir/x_tloc2 param_tloc_dd.rb

# ad tloc
mpirun -np 10 $bindir/x_tloc2 param_tloc_ad.rb
