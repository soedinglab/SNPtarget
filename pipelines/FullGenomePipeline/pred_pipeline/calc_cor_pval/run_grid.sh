source ../env.sh

grid_params="-W 24:00 -n 16 -R span[hosts=1] -m hh ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/corr_pval_$(date +%s)_XXX)
mkdir -p $outdir/{logs,cor_pval}

for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir 15"
	bsub $grid_params -J corr_pval_$chrom -o $outdir/logs/${chrom}_job.log $cmd
	sleep 1
done
ln -snrf $outdir latest_run
