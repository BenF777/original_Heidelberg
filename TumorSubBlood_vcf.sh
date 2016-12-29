#!/bin/bash
CONFIG_FILE=$1
source $CONFIG_FILE

VAR_PATH=$2
CALLER_TYPE=$3
SAMPLE_A=$4
SAMPLE_B=$5


FILE_A=$(ls -d -1 $VAR_PATH/${SAMPLE_A}*_comp_SNP.recode_filtered*.vcf)
FILE_B=$(ls -d -1 $VAR_PATH/${SAMPLE_B}*_comp_SNP.recode_filtered*.vcf)
FILE_OUT_SNP=$VAR_PATH/${SAMPLE_A}_sub_${SAMPLE_B}_comp_SNP.recode_filtered.vcf

bedtools intersect -a $FILE_A -b $FILE_B -v -header > $FILE_OUT_SNP
bgzip -c $FILE_OUT_SNP > $FILE_OUT_SNP.gz
tabix -f $FILE_OUT_SNP.gz

FILE_A=$(ls -d -1 $VAR_PATH/${SAMPLE_A}*_comp_INDEL.recode_filtered*.vcf)
FILE_B=$(ls -d -1 $VAR_PATH/${SAMPLE_B}*_comp_INDEL.recode_filtered*.vcf)
FILE_OUT_INDEL=$VAR_PATH/${SAMPLE_A}_sub_${SAMPLE_B}_comp_INDEL.recode_filtered.vcf

bedtools intersect -a $FILE_A -b $FILE_B -v -header > $FILE_OUT_INDEL
bgzip -c $FILE_OUT_INDEL > $FILE_OUT_INDEL.gz
tabix -f $FILE_OUT_INDEL.gz

echo "VARIANT ANNOTATION"
bash $SCRIPT_PATH/dna/anovar_annotate.sh $CONFIG_FILE $FILE_OUT_SNP
bash $SCRIPT_PATH/dna/anovar_annotate.sh $CONFIG_FILE $FILE_OUT_INDEL


if [ "$CALLER_TYPE" == "MPILEUP" ]; then   
    Rscript $SCRIPT_PATH/dna/annovar_csv2xlsx_mpileup.R $VAR_PATH/${SAMPLE_A}_sub_${SAMPLE_B}_comp_SNP.recode_filtered.hg19_multianno.csv
    Rscript $SCRIPT_PATH/dna/annovar_csv2xlsx_mpileup.R $VAR_PATH/${SAMPLE_A}_sub_${SAMPLE_B}_comp_INDEL.recode_filtered*hg19_multianno.csv
elif  [ "$CALLER_TYPE" == "PLATYPUS" ]; then
     Rscript $SCRIPT_PATH/dna/annovar_csv2xlsx_platypus.R $VAR_PATH/${SAMPLE_A}_sub_${SAMPLE_B}_comp_SNP.recode_filtered.hg19_multianno.csv
     Rscript $SCRIPT_PATH/dna/annovar_csv2xlsx_platypus.R $VAR_PATH/${SAMPLE_A}_sub_${SAMPLE_B}_comp_INDEL.recode_filtered.hg19_multianno.csv
fi
