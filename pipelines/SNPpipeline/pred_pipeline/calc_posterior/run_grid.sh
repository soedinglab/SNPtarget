source ../env.sh

grid_params="-W 24:00 -n 8 -m hh ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/post_pred_$(date +%s)_XXX)
mkdir -p $outdir/{logs,post_pred}

for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir"
	bsub $grid_params -J calc_posterior_$chrom -o $outdir/logs/${chrom}_job.log $cmd
done

ln -snrf $outdir latest_run
