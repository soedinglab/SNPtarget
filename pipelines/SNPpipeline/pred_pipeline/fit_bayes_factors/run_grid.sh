source ../env.sh

grid_params="-W 2:00 -n 2 ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/pred_bf_$(date +%s)_XXX)
mkdir -p $outdir/logs/hic 
mkdir -p $outdir/{pval,dist,hic}/plots

for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir"
	bsub $grid_params -J calc_bf_$chrom -o $outdir/logs/${chrom}_job.log $cmd
done

cmd="bash job_single.sh $chrom $outdir"
bsub $grid_params -J calc_hicbf -o $outdir/logs/hic/job.log $cmd

ln -snrf $outdir latest_run
