source ../env.sh

grid_params="-W 2:00 -n 16 -R span[hosts=1] ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/access_db_$(date +%s)_XXX)
mkdir -p $outdir/logs $outdir/access_db
for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir"
	bsub $grid_params -J norm_all -o $outdir/logs/${chrom}_job.log $cmd
done
ln -snrf $outdir latest_run
