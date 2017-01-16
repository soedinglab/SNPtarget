source ../env.sh

grid_params="-W 4:00 -n 8 -R span[hosts=1] ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/post_pred_$(date +%s)_XXX)
mkdir -p $outdir/{logs,post_pred}

for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir"
	eval "bsub $grid_params -J calc_posterior_$chrom -o $outdir/logs/${chrom}_job.log \"$cmd\""
	sleep 1
done

ln -snrf $outdir latest_run
