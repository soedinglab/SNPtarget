#EXTRA_GRID_ARGS="-R cbscratch"
#EXTRA_GRID_ARGS="$EXTRA_GRID_ARGS -m hh"

EXTRA_GRID_ARGS="$EXTRA_GRID_ARGS -q mpi -m 'mpi mpi2 mpi3_all hh sa' -R scratch"
HIC_RESOLUTION=25000

for chrom in {1..22}; do
	CHROMS="${CHROMS} chr$chrom"
done
