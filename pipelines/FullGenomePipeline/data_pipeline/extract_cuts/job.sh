source ../env.sh

CHROM=$1
OUTPUT_DIR=$2
N_PROC=$3
WIDTH=100

source activate sci34
extract_cuts_intervals\
	--bam_db bam_db\
	--n_proc $N_PROC\
	--log_file $OUTPUT_DIR/logs/${CHROM}.log\
	${CHROM}\
	${WIDTH}\
	bam_files.txt\
	$OUTPUT_DIR/cut_db/${CHROM}_rawcuts.pkl
