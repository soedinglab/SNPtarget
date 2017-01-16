source ../env.sh

OUTDIR=$1

source activate sci34

# using rpy2 here
module load R/3.2.3

PRED_INPUT=pred_pp

HIC_PLOT_DIR=$OUTDIR/plots

fit_hic_bf\
	--plot_dir $HIC_PLOT_DIR\
	--log_file $OUTDIR/logs/hic_bf.log\
	pred_pp\
	dist_bf_db\
	dist_bins.txt\
	$OUTDIR/all_chroms.bfpkl
