#!/bin/bash
CONFIG_FILE=$1
source $CONFIG_FILE

OUT_PATH=$2
BED_FILE=$3
FILENAME=$4

echo "run samtools"
OUT=$OUT_PATH/$(basename $FILENAME)
OUT=${OUT%_MERGED.bam}
OUT=${OUT}_MPILEUP

echo $FILENAME $OUT
REFERENCE=$REFERENCE_PATH/bwa06_1KGRef/hs37d5.fa

samtools mpileup -C50 -D -S -uf $REFERENCE $FILENAME | bcftools view -bvcg - > ${OUT}.bcf 
bcftools view  ${OUT}.bcf > ${OUT}.vcf

vcftools --vcf ${OUT}.vcf --recode --recode-INFO-all --bed $BED_FILE --out  ${OUT}_comp
mv ${OUT}_comp.recode.vcf  ${OUT}_comp.vcf
 
vcftools --vcf  ${OUT}_comp.vcf --remove-indels --recode --recode-INFO-all --out ${OUT}_comp_SNP
bgzip -c ${OUT}_comp_SNP.recode.vcf > ${OUT}_comp_SNP.recode.vcf.gz
tabix ${OUT}_comp_SNP.recode.vcf.gz

vcftools --vcf  ${OUT}_comp.vcf --keep-only-indels --recode --recode-INFO-all --out ${OUT}_comp_INDEL
bgzip -c ${OUT}_comp_INDEL.recode.vcf > ${OUT}_comp_INDEL.recode.vcf.gz
tabix ${OUT}_comp_INDEL.recode.vcf.gz

#Filter
echo "FILTER $OUT"

Rscript $SCRIPT_PATH/dna/VariantFilter_mpileup_SNP.R 40 0.1 0.1 ${OUT}_comp_SNP.recode.vcf.gz
Rscript $SCRIPT_PATH/dna/VariantFilter_mpileup_INDEL.R 40 ${OUT}_comp_INDEL.recode.vcf.gz

#variant annotation
echo "VARIANT ANNOTATION"
bash $SCRIPT_PATH/dna/anovar_annotate.sh $CONFIG_FILE  $(ls -d -1 ${OUT}_comp_SNP.recode_filtered*.vcf)
bash $SCRIPT_PATH/dna/anovar_annotate.sh $CONFIG_FILE  $(ls -d -1 ${OUT}_comp_INDEL.recode_filtered*.vcf)

Rscript $SCRIPT_PATH/dna/annovar_csv2xlsx_mpileup.R  $(ls -d -1 ${OUT}_comp_SNP.recode_filtered*.hg19_multianno.csv)
Rscript $SCRIPT_PATH/dna/annovar_csv2xlsx_mpileup.R  $(ls -d -1 ${OUT}_comp_INDEL.recode_filtered*.hg19_multianno.csv)

echo "$OUT finished"

