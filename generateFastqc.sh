#!/bin/sh

CONFIG_FILE=$1
source $CONFIG_FILE

RUN_PATH=$2
FASTQ_PATH=$RUN_PATH/out/fastq
FASTQC_PATH=$RUN_PATH/out/fastqc
LOG_PATH=$RUN_PATH/out/log

LOG_FILE=$LOG_PATH/fastqc.log

echo "Fastqc-Analysis" > $LOG_FILE
date >> $LOG_FILE
echo "" >> $LOG_FILE
echo "FASTQC-VERSION:" >> $LOG_FILE

$FASTQC_BINARY --version >> $LOG_FILE

#create the folder for specific run
mkdir -p $FASTQC_PATH

echo "" >> $LOG_FILE
echo "RUN-Details:" >> $LOG_FILE

find $FASTQ_PATH -name "*.fastq.gz" | grep -v I1 | sort | grep -v NuGen | parallel -X $FASTQC_BINARY -t 2 --nogroup -o $FASTQC_PATH 


