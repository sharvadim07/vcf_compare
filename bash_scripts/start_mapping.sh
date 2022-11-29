#!/bin/bash

set -e

INPUT_DIR=/gpfs/AfterDemultiplex/cleaning_witch_broom/clean_v2_2/psyl_reads/

#Tree 1
./start_hisat_PE.sh 10smple_may ${INPUT_DIR}/start/10/10/sortmerna/sortmerna_out_R1.fastq ${INPUT_DIR}/start/10/10/sortmerna/sortmerna_out_R2.fastq 
./start_hisat_PE.sh 11smple_may ${INPUT_DIR}/start/11/11/sortmerna/sortmerna_out_R1.fastq ${INPUT_DIR}/start/11/11/sortmerna/sortmerna_out_R2.fastq

./start_hisat_PE.sh 10smple_aug ${INPUT_DIR}/end/10/10/sortmerna/sortmerna_out_R1.fastq ${INPUT_DIR}/end/10/10/sortmerna/sortmerna_out_R2.fastq
./start_hisat_PE.sh 11smple_aug ${INPUT_DIR}/end/11/11/sortmerna/sortmerna_out_R1.fastq ${INPUT_DIR}/end/11/11/sortmerna/sortmerna_out_R2.fastq


#Tree 2
./start_hisat_PE.sh 14smple_may ${INPUT_DIR}/start/14/14/sortmerna/sortmerna_out_R1.fastq ${INPUT_DIR}/start/14/14/sortmerna/sortmerna_out_R2.fastq
./start_hisat_PE.sh 15smple_may ${INPUT_DIR}/start/15/15/sortmerna/sortmerna_out_R1.fastq ${INPUT_DIR}/start/15/15/sortmerna/sortmerna_out_R2.fastq

./start_hisat_PE.sh 14smple_aug ${INPUT_DIR}/end/14/14/sortmerna/sortmerna_out_R1.fastq ${INPUT_DIR}/end/14/14/sortmerna/sortmerna_out_R2.fastq
./start_hisat_PE.sh 15smple_aug ${INPUT_DIR}/end/15/15/sortmerna/sortmerna_out_R1.fastq ${INPUT_DIR}/end/15/15/sortmerna/sortmerna_out_R2.fastq
