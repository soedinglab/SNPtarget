source ../env.sh

CHROM=$1
OUTDIR=$2

source activate sci34

calc_posterior\
	--log_file $OUTDIR/logs/${CHROM}.log\
	pred_pp/${CHROM}.predpp\
	$OUTDIR/post_pred/${CHROM}.post_pred\
	pval_bf/${CHROM}.bfpkl\
	dist_bf/${CHROM}.bfpkl\
	hic_bf/all_chroms.bfpkl
