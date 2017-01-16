source ../env.sh

OUTDIR=$2

source activate sci34

# using rpy2 here
module load R/3.2.3

PRED_INPUT=pred_pp

HIC_PLOT_DIR=$OUTDIR/hic/plots/

fit_hic_bf\
	--plot_dir $HIC_PLOT_DIR\
	--log_file $OUTDIR/logs/hic/hic_bf.log\
	$PRED_INPUT\
	dist_bins.txt\
	$OUTDIR/hic/all_chroms.bfpkl
