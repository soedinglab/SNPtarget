source ../env.sh

OUTPUT_DIR=$1
CHROM=$2

source activate sci34
filter_promoters\
	access_db/${CHROM}.norm_access\
	trscr.db\
	${CHROM}\
	--updown_bp 500\
	> ${OUTPUT_DIR}/access_db/${CHROM}.norm_access
