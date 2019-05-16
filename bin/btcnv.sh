#!/bin/bash

sample_name=$1
orig_bam=$2
out_dir=`readlink -f $3`

bindir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
datadir=${bin_dir%bin}"DataFile"

panel=${datadir}"/ct_panel_152_auto_withg_uniq.bed"
exclude_region=${datadir}"/exclude_regions.bed"
script=${bindir}"/normalize.R"
filter_script=${bindir}"/filter_gam.pl"
insize=${bindir}"/insertsize.R"

SAMBAMBA=$(which sambamba)
RSCRIPT=$(which Rscript)
DATAMASH=$(which datamash)
SAMTOOLS=$(which samtools)
BEDTOOLS=$(which bedtools)

mkdir -p ${out_dir}/${sample_name}_btcnv
cd ${out_dir}/${sample_name}_btcnv

$SAMBAMBA depth region -t 3 -m -L $panel $orig_bam | sort -k1,1 -k2,2n | uniq > $sample_name.interval &
$SAMTOOLS stats -c 100,10000,100 --trim-quality 1 $orig_bam | grep ^IS | cut -f 2- > $sample_name.insertstat &
wait
$RSCRIPT $script -i $sample_name.interval -p $sample_name -d $datadir
$RSCRIPT $insize -i $sample_name.insertstat -p $sample_name

#grep -v "BCL2L11\|CYP2D6" $sample_name.tsv > $sample_name.nogerm.tsv
tail -n +2 $sample_name.tsv | $BEDTOOLS intersect -v -a stdin -b $exclude_region > $sample_name.care.tsv
perl $filter_script $sample_name.care.tsv $sample_name.insertstat $bindir > $sample_name.cnv.xls
