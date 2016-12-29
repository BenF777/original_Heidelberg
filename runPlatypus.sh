#!/bin/bash
CONFIG_FILE=$1
source $CONFIG_FILE

OUT_PATH=$2
BED_FILE=$3
BAM_FILE=$4

echo $OUT_PATH
echo $COMPLETE_BED_FILE
echo $BAM_FILE

OUT=$OUT_PATH/$(basename $BAM_FILE)
OUT=${OUT%_MERGED.bam}
OUT=${OUT}_PLATYPUS
REFERENCE=$REFERENCE_PATH/bwa06_1KGRef/hs37d5.fa

python $PLATYPUS_BINARY callVariants --bamFiles=$BAM_FILE --refFile=$REFERENCE --output=$OUT.vcf --filterDuplicates=0 --nCPU=$DNA_PARALLEL_ALIGNMENT --logFileName=$OUT.log

vcftools --vcf ${OUT}.vcf --recode --recode-INFO-all --bed $BED_FILE --out  ${OUT}_comp
mv ${OUT}_comp.recode.vcf  ${OUT}_comp.vcf
 
vcftools --vcf  ${OUT}_comp.vcf --remove-indels --recode --recode-INFO-all --out ${OUT}_comp_SNP
bgzip -c ${OUT}_comp_SNP.recode.vcf > ${OUT}_comp_SNP.recode.vcf.gz
tabix ${OUT}_comp_SNP.recode.vcf.gz

vcftools --vcf  ${OUT}_comp.vcf --keep-only-indels --recode --recode-INFO-all --out ${OUT}_comp_INDEL
bgzip -c ${OUT}_comp_INDEL.recode.vcf > ${OUT}_comp_INDEL.recode.vcf.gz
tabix ${OUT}_comp_INDEL.recode.vcf.gz

echo "FILTER Platypus"
Rscript $SCRIPT_PATH/dna/VariantFilter_platypus_SNP.R 40 0.1 0.1 ${OUT}_comp_SNP.recode.vcf.gz
Rscript $SCRIPT_PATH/dna/VariantFilter_platypus_INDEL.R 40 -1 -1 ${OUT}_comp_INDEL.recode.vcf.gz

#variant annotation
echo "VARIANT ANNOTATION"
bash $SCRIPT_PATH/dna/anovar_annotate.sh $CONFIG_FILE $(ls -d -1 ${OUT}_comp_SNP.recode_filtered*.vcf)
bash $SCRIPT_PATH/dna/anovar_annotate.sh $CONFIG_FILE $(ls -d -1 ${OUT}_comp_INDEL.recode_filtered*.vcf)

Rscript $SCRIPT_PATH/dna/annovar_csv2xlsx_platypus.R  $(ls -d -1 ${OUT}_comp_SNP.recode_filtered*.hg19_multianno.csv)
Rscript $SCRIPT_PATH/dna/annovar_csv2xlsx_platypus.R  $(ls -d -1 ${OUT}_comp_INDEL.recode_filtered*.hg19_multianno.csv)

echo "$OUT finished"
