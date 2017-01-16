source ../env.sh

CHROM=$1
OUTPUT_DIR=$2
N_PROC=$3

source activate sci34
extract_cuts\
	--bam_db bam_db\
	--n_proc $N_PROC\
	--log_file $OUTPUT_DIR/logs/${CHROM}.log\
	snp_db/${CHROM}.snpdb\
	snp_lists/${CHROM}_snp.list\
	bam_files.txt\
	$OUTPUT_DIR/cut_db/${CHROM}_rawcuts.pkl
