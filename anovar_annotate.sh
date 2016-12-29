#!/bin/bash
CONFIG_FILE=$1
source $CONFIG_FILE

FILENAME=$2
AVINPUT=${FILENAME%.vcf}.avinput
CSV_OUT=${FILENAME%.vcf}
$ANNOVAR_HOME/convert2annovar.pl -format vcf4 $FILENAME -outfile $AVINPUT -include -withzyg -includeinfo
$ANNOVAR_HOME/table_annovar.pl $AVINPUT $ANNOVAR_DB -buildver hg19 -out $CSV_OUT -remove -otherinfo -protocol refGene,cytoBand,snp138,ljb_all,1000g2012apr_all,cosmic68 -operation g,r,f,f,f,f -nastring NA -csvout

