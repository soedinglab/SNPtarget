source ../env.sh

grid_params="-W 6:00 -n 2 -R span[hosts=1] ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/hic_bf_$(date +%s)_XXX)
mkdir -p $outdir/logs
mkdir -p $outdir/plots

cmd="bash job.sh $outdir"
eval "bsub $grid_params -J calc_hicbf -o $outdir/logs/job.log \"$cmd\""

ln -snrf $outdir latest_run
