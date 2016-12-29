#!/bin/bash
CONFIG_FILE=$1
source $CONFIG_FILE

BAM_PATH=$2
OUT_PATH=$3/calls

mkdir -p $OUT_PATH

FILENAMES=($(ls -d -1 $BAM_PATH/*MERGED.bam))

i=0

for FILENAME in "${FILENAMES[@]}" 
do 
    :
    IFS='/' read -ra FILE_COMP <<< "$FILENAME"
    bash $SCRIPT_PATH/dna/samtools_special_pos.sh $CONFIG_FILE "$FILENAME" "$OUT_PATH/" &

    if (( $(($i % $DNA_PARALLEL_ALIGNMENT)) == 0 )) 
    then
        wait
    fi    
    i=$(($i+1))
done
