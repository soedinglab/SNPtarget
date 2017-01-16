EXTRA_GRID_ARGS="-R cbscratch"
# EXTRA_GRID_ARGS="$EXTRA_GRID_ARGS -m hh"
for chrom in {1..22}; do
	CHROMS="${CHROMS} chr$chrom"
done
