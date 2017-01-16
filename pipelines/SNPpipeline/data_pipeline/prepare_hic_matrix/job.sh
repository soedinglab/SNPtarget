source ../env.sh

CHROM=$1
OUTPUT_DIR=$2
RESOLUTION=1000

source activate sci34
rao_hic_to_matrix\
	--type band\
	--max_distance 2000000\
	hic_db/${CHROM}.raomat\
	$RESOLUTION\
	$OUTPUT_DIR/hic_db/${CHROM}.hicbmat
