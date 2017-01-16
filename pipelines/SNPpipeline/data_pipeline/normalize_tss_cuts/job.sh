source ../env.sh

OUTPUT_DIR=$1
PERC_CUTOFF=99.85

source activate sci34
normalize_cuts\
	--log_file ${OUTPUT_DIR}/logs/norm.log\
	--perc_cutoff ${PERC_CUTOFF}\
	cut_matrices\
	${OUTPUT_DIR}/norm_access\
	${OUTPUT_DIR}/error\
	--max_d_dump 100000\
	--dump_raw_data ${OUTPUT_DIR}/raw_dump.pkl\
	--dump_proc_data ${OUTPUT_DIR}/norm_data.pkl
