source ../env.sh

CHROM=$1
OUTPUT_DIR=$2

source activate sci34
rao_hic_to_matrix\
	--type band\
	--max_distance 2050000\
	hic_db/${CHROM}.raomat\
	$HIC_RESOLUTION\
	$OUTPUT_DIR/hic_db/${CHROM}.hicbmat
