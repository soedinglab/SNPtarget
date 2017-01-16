EXTRA_GRID_ARGS="-R cbscratch -R span[hosts=1]"

# only run on soeding cluster
#EXTRA_GRID_ARGS="$EXTRA_GRID_ARGS -m hh"

TSS_MIN_DIST=2000

export TMPDIR=/local
for chrom in {1..22}; do
	CHROMS="${CHROMS} chr$chrom"
done
