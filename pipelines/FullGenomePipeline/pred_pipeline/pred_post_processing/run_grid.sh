source ../env.sh

grid_params=" -W 4:00 -n 8 -R span[hosts=1] $EXTRA_GRID_ARGS"

outdir=$(mktemp -d $(pwd)/pred_pp_$(date +%s)_XXX)
mkdir -p $outdir/pred_db
mkdir -p $outdir/logs/{hic_pval,extr_hic}

for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir 15"
	eval "bsub $grid_params -J pred_pp_$chrom -o $outdir/logs/${chrom}_job.log \"$cmd\""
	sleep 2
done
ln -snrf $outdir latest_run
