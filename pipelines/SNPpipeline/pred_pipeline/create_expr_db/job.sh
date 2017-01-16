OUTPUT_DIR=$1

SMOOTH_EXPR=${OUTPUT_DIR}/smooth_expr.expr
EXPR_DB=${OUTPUT_DIR}/expr_db

source activate sci34
noise_control\
	--type expr\
	--mode min_value\
	--value 4\
	expr.expr\
	> $SMOOTH_EXPR
sync

split_expr\
	--output_dir ${EXPR_DB}\
	$SMOOTH_EXPR
