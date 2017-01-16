source ../env.sh

grid_params="-W 2:00 -n 2 -R span[hosts=1] ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/pred_bf_$(date +%s)_XXX)
mkdir -p $outdir/logs/hic
mkdir -p $outdir/{pval,dist}/plots

for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir"
	eval "bsub $grid_params -J calc_bf_$chrom -o $outdir/logs/${chrom}_job.log \"$cmd\""
	sleep 2
done

ln -snrf $outdir latest_run
