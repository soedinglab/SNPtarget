source ../env.sh

CHROM=$1
OUTDIR=$2
source activate sci34

INFILE=pred_db/${CHROM}.pred
OUTFILE=$OUTDIR/pred_db/${CHROM}.predpp

DIST_PRED=$TMPDIR/${CHROM}_dist.tmp
HIC_PRED=$TMPDIR/${CHROM}_hic.tmp
HICPVAL_PRED=$TMPDIR/${CHROM}_hicpval.tmp

calc_dist\
	$INFILE\
	$DIST_PRED

extract_hic\
	--bg_bandwidth 200000\
	--full_bg\
	hic_db/${CHROM}.hicbmat\
	$DIST_PRED\
	$HIC_PRED

hic_pval\
	--log_file $OUTDIR/logs/hic_pval/${chrom}.log\
	$HIC_PRED\
	1000\
	dist_bins.txt\
	$HICPVAL_PRED

mv $HICPVAL_PRED $OUTFILE
rm $DIST_PRED $HIC_PRED
