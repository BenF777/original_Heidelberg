#!/bin/bash


CONFIG_FILE=$1
source $CONFIG_FILE

OUTPUT_PATH=$2
LOG_PATH=$3

LANE_1=$(basename $4)

LANE_2=$(basename $5)
LANE_3=$(basename $6)
LANE_4=$(basename $7)
echo $LANE_1
echo $LANE_2
echo $LANE_3
echo $LANE_4
IFS='_' read -ra LANE1_COMP <<< "$LANE_1"
IFS='_' read -ra LANE2_COMP <<< "$LANE_2"
IFS='_' read -ra LANE3_COMP <<< "$LANE_3"
IFS='_' read -ra LANE4_COMP <<< "$LANE_4"

if [ ${LANE1_COMP[0]} !=  ${LANE2_COMP[0]} ] || \
       [ ${LANE1_COMP[0]} !=  ${LANE3_COMP[0]} ] || \
       [ ${LANE1_COMP[0]} !=  ${LANE4_COMP[0]} ]
then
    echo "EXIT: IDs not equal!"
    exit
fi

if [ ${LANE1_COMP[1]} ==  ${LANE2_COMP[1]} ] || \
       [ ${LANE1_COMP[1]} ==  ${LANE3_COMP[1]} ] || \
       [ ${LANE1_COMP[1]} ==  ${LANE4_COMP[1]} ] || \
       [ ${LANE2_COMP[1]} ==  ${LANE3_COMP[1]} ] || \
       [ ${LANE2_COMP[1]} ==  ${LANE4_COMP[1]} ] || \
       [ ${LANE3_COMP[1]} ==  ${LANE4_COMP[1]} ]
then
    echo "EXIT: Equal lanes!"
    exit
fi

OUTPUT_FILE=$OUTPUT_PATH/${LANE1_COMP[0]}_MERGED.bam
LOG_FILE=$LOG_PATH/${LANE1_COMP[0]}_MERGED.log
echo "Merge BAM Files" > $LOG_FILE
date >> $LOG_FILE
echo "FILE 1: $LANE_1" >> $LOG_FILE
echo "FILE 2: $LANE_2" >> $LOG_FILE
echo "FILE 3: $LANE_3" >> $LOG_FILE
echo "FILE 4: $LANE_4" >> $LOG_FILE
echo "OUTPUT FILE: $OUTPUT_FILE" >> $LOG_FILE
echo "" >> $LOG_FILE
echo "Run samtools merge:" >> $LOG_FILE

$SAMTOOLS_SORT_BINARY  merge -f $OUTPUT_FILE $OUTPUT_PATH/$LANE_1 $OUTPUT_PATH/$LANE_2 $OUTPUT_PATH/$LANE_3 $OUTPUT_PATH/$LANE_4 2>>$LOG_FILE

echo "" >> $LOG_FILE
echo "Run samtools index:" >> $LOG_FILE

$SAMTOOLS_SORT_BINARY index $OUTPUT_FILE 2>> $LOG_FILE 

echo "" >> $LOG_FILE
echo "Run igv-tools count:" >> $LOG_FILE

time sh ${IGV_TOOL_BINARY} count ${OUTPUT_FILE} ${OUTPUT_FILE}.tdf $IGV_TOOL/genomes/1kg_v37.chrom.sizes >> $LOG_FILE

