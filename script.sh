#!/bin/bash

ARRAY=(/xchip/gistic/Jeremiah/tracks/repeat_masker_hg19_SINE.filtered.bed \
/xchip/gistic/Jeremiah/tracks/UCSC_predicted_genes.trim.bed \
/xchip/gistic/Jeremiah/Projects/ICGC/RefAndAnnot/Roadmap_regulatory_elements/2f9970ee-f311-4e4d-aa50-4b1ac986c031/roadmap_DNase_all.nochr.bed \
/xchip/gistic/Jeremiah/tracks/repeat_masker_hg19_LINE.filtered.bed \
/xchip/gistic/Jeremiah/tracks/repeat_masker_hg19_LTR.filtered.bed)

INPUTS=""
for i in "${ARRAY[@]}"
do
   INPUTS="${INPUTS} -I $i"
done

EVENTS=/xchip/gistic/Jeremiah/Projects/ICGC/1d_breaks_161014.bed
MASK=/xchip/gistic/Jeremiah/Projects/ICGC/covered_1d.bed
WIDTH=100000
SLOP=1000
##echo "src/ginseng fishhook $EVENTS $INPUTS -m $MASK -v -w $WIDTH -s $SLOP -b $na1 -G $REFHG19"
src/ginseng fishhook $EVENTS $INPUTS -m $MASK -v -w $WIDTH -s $SLOP -b $na1 -G $REFHG19 -p 10
