source ../env.sh

grid_params="-W 1:00 -n 2 -R span[hosts=1] ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/mask_tss_$(date +%s)_XXX)
mkdir -p $outdir/{logs,access_db}
for chrom in $CHROMS; do
	cmd="bash job.sh $outdir $chrom"
	bsub $grid_params -J mask_tss -o $outdir/logs/job.log $cmd
	ln -snrf $outdir latest_run
done
