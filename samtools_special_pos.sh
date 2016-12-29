#!/bin/bash

CONFIG_FILE=$1
source $CONFIG_FILE

FILENAME=$2

OUT_PATH=$3

IFS='/' read -ra FILE_COMP <<< "$FILENAME"

PANEL_POS=${#FILE_COMP[@]}

PANEL=${FILE_COMP[$(($PANEL_POS - 3))]}

OUTNAME=$(basename $FILENAME)
OUTNAME=${OUTNAME%_MERGED.bam}
OUTNAME=${OUTNAME}_$PANEL

OUT=${OUT_PATH}${OUTNAME}


REFERENCE=$REFERENCE_PATH/bwa06_1KGRef/hs37d5.fa
echo $OUT

samtools mpileup -D -S -Q 8 -r 5:1295228-1295228 -uf $REFERENCE $FILENAME | bcftools view -bcg - > ${OUT}_TERT1.bcf
samtools mpileup -D -S -Q 8 -r 5:1295250-1295250 -uf $REFERENCE $FILENAME | bcftools view -bcg - > ${OUT}_TERT2.bcf

samtools mpileup -D -S -Q 8 -r 2:209113111-209113113 -uf $REFERENCE $FILENAME | bcftools view -bcg - > ${OUT}_IDH1_R132.bcf
samtools mpileup -D -S -Q 8 -r 15:90631837-90631839 -uf $REFERENCE $FILENAME | bcftools view -bcg - > ${OUT}_IDH2_R172.bcf

samtools mpileup -D -S -Q 8 -r 7:140453134-140453136 -uf $REFERENCE $FILENAME | bcftools view -bcg - > ${OUT}_BRAF_V600.bcf

samtools mpileup -D -S -Q 8 -r 1:226252134-226252136 -uf $REFERENCE $FILENAME | bcftools view -bcg - > ${OUT}_H3F3A_K27.bcf
samtools mpileup -D -S -Q 8 -r 1:226252155-226252157 -uf $REFERENCE $FILENAME | bcftools view -bcg - > ${OUT}_H3F3A_G34.bcf

bcftools view  ${OUT}_TERT1.bcf > ${OUT}_TERT1.vcf
bcftools view  ${OUT}_TERT2.bcf > ${OUT}_TERT2.vcf
bcftools view  ${OUT}_IDH1_R132.bcf > ${OUT}_IDH1_R132.vcf
bcftools view  ${OUT}_IDH2_R172.bcf > ${OUT}_IDH2_R172.vcf
bcftools view  ${OUT}_BRAF_V600.bcf > ${OUT}_BRAF_V600.vcf
bcftools view  ${OUT}_H3F3A_K27.bcf > ${OUT}_H3F3A_K27.vcf
bcftools view  ${OUT}_H3F3A_G34.bcf > ${OUT}_H3F3A_G34.vcf

#CALL R Summary
Rscript $SCRIPT_PATH/dna/Summary_special.R $OUT_PATH $OUT_PATH.. $PANEL $OUT

echo "$OUT finished"
