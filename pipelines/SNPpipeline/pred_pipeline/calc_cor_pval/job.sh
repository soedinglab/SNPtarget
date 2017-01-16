source ../env.sh

CHROM=$1
OUTDIR=$2
THREADS=$3
source activate sci34

PRED_FILE=$OUTDIR/cor_pval/${CHROM}.pred

corr_pval\
	--seed 20151130\
	--ks_split 0.7\
	--n_threads $THREADS\
	--log_file $OUTDIR/logs/${CHROM}.log\
	--ks_probe_threshold 1e-3\
	--ks_bg_threshold 1e-1\
	--dist_ks_thresh 0\
	--min_bg_samples 2000\
	--min_n_high_dist 5\
	--min_bg_samples 10\
	--tss_min_dist ${TSS_MIN_DIST}\
	regions.reg\
	$CHROM\
	$OUTDIR/cor_pval/${CHROM}.bgstat\
	${PRED_FILE}\
	access_db\
	expr_db
