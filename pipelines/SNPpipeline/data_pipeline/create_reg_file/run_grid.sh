source ../env.sh

grid_params="-W 1:00 -n 2 ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/create_regions_$(date +%s)_XXX)
mkdir -p $outdir/logs

cmd="bash job.sh $outdir"
bsub $grid_params -J create_regions -o $outdir/logs/job.log $cmd

ln -snrf $outdir latest_run
