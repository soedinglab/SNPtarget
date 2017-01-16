source ../env.sh

grid_params="-W 24:00 -n 16 -R span[hosts=1] ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/norm_access_$(date +%s)_XXX)
mkdir -p $outdir/logs
cmd="bash job.sh $outdir"
bsub $grid_params -J norm_all -o $outdir/logs/job.log $cmd
ln -snrf $outdir latest_run
