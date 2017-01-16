CHROM=$1
OUTDIR=$2

source activate sci34
noise_control\
	--type dhs\
	--value 0\
	--mode min_value\
	access_db/$CHROM.norm_access\
	> $OUTDIR/access_db/$CHROM.norm_access	
