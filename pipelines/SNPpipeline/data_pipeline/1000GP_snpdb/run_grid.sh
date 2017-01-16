source ../env.sh

grid_params="-W 2:00 -n 2 ${EXTRA_GRID_ARGS}"

MIN_MAF="0.00"
outdir=$(mktemp -d $(pwd)/snpdb_$(date +%s)_XXX)
mkdir -p $outdir/{logs,snp_db}

for chrom in $CHROMS; do
	infile=$(ls 1000gp_vcfdb/ALL.${chrom}.*.vcf.gz)
	cmd="bash job.sh $infile $outdir/snp_db/${chrom}.snpdb $MIN_MAF"
	bsub $grid_params -J vcf2db_${chrom} -o $outdir/logs/${chrom}.log $cmd
done
ln -snrf $outdir latest_run
