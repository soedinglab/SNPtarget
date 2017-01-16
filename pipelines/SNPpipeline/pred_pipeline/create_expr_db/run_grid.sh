source ../env.sh

grid_params="-W 0:30 -n 2 ${EXTRA_GRID_ARGS}"

outdir=$(mktemp -d $(pwd)/expr_db_$(date +%s)_XXX)
mkdir -p $outdir/logs -p $outdir/expr_db

cmd="bash job.sh $outdir"
bsub $grid_params -J generate_expr_db -o $outdir/logs/job.log $cmd

ln -snrf $outdir latest_run
