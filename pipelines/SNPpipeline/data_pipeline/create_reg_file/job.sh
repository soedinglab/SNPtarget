source ../env.sh

OUTPUT_DIR=$1

source activate sci34
select_tss\
	gc.trscr\
	tss_access_db\
	> $OUTPUT_DIR/gc.filtered

trscr2reg\
	--range 2000000\
	$OUTPUT_DIR/gc.filtered\
	$OUTPUT_DIR/gc.reg


