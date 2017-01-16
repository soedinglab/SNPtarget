source ../env.sh

source activate sci34
zcat $1 | tgvcf2snpdb - --minimum_maf $3 > $2
