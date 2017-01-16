source ../env.sh

grid_params="-W 24:00 -n 16 -R span[hosts=1] ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/hic_mat_$(date +%s)_XXX)
mkdir -p $outdir/logs $outdir/hic_db

for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir"
	eval "bsub $grid_params -J extr_hic_$chrom -o $outdir/logs/${chrom}_job.log \"$cmd\""
	sleep 2
done

ln -snrf $outdir latest_run
