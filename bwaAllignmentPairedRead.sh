#!/bin/bash
CONFIG_FILE=$1
source $CONFIG_FILE

BAM_DIR=$2
LOG_DIR=$3
RUN_ID=$4
RAW_SEQ_1=$5
RAW_SEQ_2=$6

#TODO check RAW_SEQ R1/R2

FILE_NAME=$(basename "${RAW_SEQ_1}")
echo $FILE_NAME
IFS='_' read -ra FILE_NAME_COMP <<< "$FILE_NAME"



SAMPLE_ID=${FILE_NAME_COMP[0]}
LANE_ID=${FILE_NAME_COMP[2]}

echo $SAMPLE_ID
echo $LANE_ID

RUN_GROUP_HEADER="@RG\tID:${RUN_ID}_${LANE_ID}\tLB:Lib01\tSM:${SAMPLE_ID}\tPL:ILLUMINA"
echo $RUN_GROUP_HEADER

RESULT_FILENAME=${SAMPLE_ID}_${LANE_ID}

SCRATCH_DIR=${BAM_DIR}/tmp/${SAMPLE_ID}_${LANE_ID}

echo "SEQ1: ${RAW_SEQ_1}"
echo "SEQ2: ${RAW_SEQ_2}"
echo "SCRATCH_DIR: ${SCRATCH_DIR}"
echo "BAM_DIR: ${BAM_DIR}"
echo "RESULT_PATH: ${RESULT_PATH}"
echo "FILE_NAME: ${RESULT_FILENAME}"

#create the folder for specific run
mkdir -p ${SCRATCH_DIR}
mkdir -p ${BAM_DIR}

INDEX_PREFIX=$REFERENCE_PATH/bwa06_1KGRef
REFERENCE=${INDEX_PREFIX}/hs37d5.fa

alnThreadOptions=$DNA_BWA_MEM_PARALLEL

baseBWACall="${BWA_BINARY} mem -t ${alnThreadOptions} -R ${RUN_GROUP_HEADER} ${REFERENCE}" 

echo "bwacall $baseBWACall $RAW_SEQ_1 $RAW_SEQ_2 $INDEX_SEQ"

UNZIPTOOL_OPTIONS="-c"

FNPIPE1=$SCRATCH_DIR/NAMED_PIPE1
FNPIPE2=$SCRATCH_DIR/NAMED_PIPE2

nice ${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1}  > $FNPIPE1 &
nice ${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2}  > $FNPIPE2 &
wait

BWA_LOG=${LOG_DIR}/bwa_${RESULT_FILENAME}.log
date > ${BWA_LOG}
echo "" >> ${BWA_LOG}
SAM_LOG=${LOG_DIR}/samtools_${RESULT_FILENAME}.log
date > ${SAM_LOG}
echo "" >> ${SAM_LOG}


BAM_FILE=${SCRATCH_DIR}/${RESULT_FILENAME}.sorted
SAMPESORT_MEMSIZE=2000000000

MBUFFER_2="mbuffer -q -m 500M -l /dev/null"

#alignment and creation of sorted bam
${baseBWACall} ${FNPIPE1} ${FNPIPE2} | ${MBUFFER_2} | ${SAMTOOLS_SORT_BINARY} view -uSbh - | ${MBUFFER_2} |  ${SAMTOOLS_SORT_BINARY} sort -@ 8 -m ${SAMPESORT_MEMSIZE} - ${BAM_FILE}

${SAMTOOLS_SORT_BINARY} index ${BAM_FILE}.bam 

BAM_DUPREM=${SCRATCH_DIR}/${RESULT_FILENAME}.dupsMarked.bam
PICARD_LOG=${LOG_DIR}/picard_${RESULT_FILENAME}.log
picard-tools MarkDuplicates I=${BAM_FILE}.bam O=${BAM_DUPREM} M=$PICARD_LOG REMOVE_DUPLICATES=false

${SAMTOOLS_SORT_BINARY} index ${BAM_DUPREM}


#move and remove files
rm $FNPIPE1
rm $FNPIPE2
mv ${SCRATCH_DIR}/* ${BAM_DIR}

eval "nice $cmd" # execute the command
