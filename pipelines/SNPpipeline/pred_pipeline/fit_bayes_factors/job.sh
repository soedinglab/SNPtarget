source ../env.sh

CHROM=$1
OUTDIR=$2

source activate sci34

# using rpy2 here
module load R/3.2.3

PREDFILE=pred_pp/${CHROM}.predpp

fit_pval_bf\
	--plot_file $OUTDIR/pval/plots/${CHROM}_pval.png\
	$PREDFILE\
	$OUTDIR/pval/${CHROM}.bfpkl
	
fit_dist_bf\
	--min_dist ${TSS_MIN_DIST}\
	--max_dist 2000000\
	--plot_file $OUTDIR/dist/plots/${CHROM}_dist.png\
	--pval_enrich_thresh 0.1\
	--pval_bg_thresh 0.6\
	--lsq_pval_estim \
	$PREDFILE\
	$OUTDIR/dist/${CHROM}.bfpkl
