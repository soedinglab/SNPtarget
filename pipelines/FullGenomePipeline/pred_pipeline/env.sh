#EXTRA_GRID_ARGS="-R cbscratch -R span[hosts=1]"

# only run on soeding cluster
#EXTRA_GRID_ARGS="$EXTRA_GRID_ARGS -m hh"
EXTRA_GRID_ARGS="$EXTRA_GRID_ARGS -q mpi -m 'mpi mpi2 mpi3_all hh sa' -R scratch"

TSS_MIN_DIST=2000
HIC_RESOLUTION=5000

export TMPDIR=/local
for chrom in {1..22}; do
	CHROMS="${CHROMS} chr$chrom"
done
