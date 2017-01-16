source ../env.sh

n_proc=16
grid_params="-W 8:00 -n $n_proc -R span[hosts=1] ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/extract_cuts_$(date +%s)_XXX)
mkdir -p $outdir/logs $outdir/cut_db

for chrom in ${CHROMS}; do
	cmd="bash job.sh $chrom $outdir $n_proc"
	bsub $grid_params -J extr_$chrom -o $outdir/logs/${chrom}_job.log $cmd
done

ln -snrf $outdir latest_run
