#!/bin/bash
#PBS -N RefAsmKedrPila
#PBS -j oe
#PBS -o /gpfs/AssembleData/kedr_reference_pila_asm/asm.log
#PBS -l nodes=genom:ppn=24
#PBS -l pvmem=1000gb
#PBS -d /gpfs/AssembleData/kedr_reference_pila_asm

SCRIPT_PATH=`pwd -P`
THREADS_NUM=40
SAMTOOLS_0=/gpfs/SOFT/bio/samtools-0.1.19/samtools
BCFTOOLS_0=/gpfs/SOFT/bio/samtools-0.1.19/bcftools/bcftools
VCFUTILS_0=/gpfs/SOFT/bio/samtools-0.1.19/bcftools/vcfutils.pl
HISAT=hisat2

set -e

cd $SCRIPT_PATH

# Aligning reads on assemblies
#INPUT_READS_DIR=/gpfs/AfterDemultiplex/cleaning_larch/NEW/arboretum_IL_needles/
#REF_PATH="/gpfs/AssembleData/larch_CLC_BIO/NEW/sum_1_2_3_4_all_data_v2/cov_all_data_16x/result_rmRed_rmCont_NCBI_Larch.fasta"
ALN_PREFIX=$1

#Paired
PE_R1=$2
PE_R2=$3

#Unpaired
#PE_U=`cat list_U`

echo -e "\n### Aligning PE reads on reference"

#build index
#echo "bowtie2-build $REF_PATH"
#bowtie2-build --large-index  $REF_PATH ${ALN_PREFIX} 2>>./log.stderr

#start bowtie2

echo $PE_R1
echo $PE_R2
#echo $PE_U

INDEX=/gpfs/AssembleDATA_2/Pinus_sylvestr_WB_Trinity/hisat2-index/pinsylv_wb_tr

echo "hisat2_align"
$HISAT -x $INDEX --threads $THREADS_NUM -t \
-1 $PE_R1 \
-2 $PE_R2 \
-S ${ALN_PREFIX}.sam  2>>./log.stderr

#The next steps involve creating a bam files which is a binary sam file (takes less space). 
#Samtools is a critical tool in the arsenal for dealing with the alignment file (sam file). 
echo "samtools view"
samtools view -S ${ALN_PREFIX}.sam -@ $THREADS_NUM -b -o ${ALN_PREFIX}.bam
samtools flagstat ${ALN_PREFIX}.bam > flagstat_${ALN_PREFIX}.stat

#You can open up the bam file mapped onto the reference in tablet and visualize. 
#the alignment but first you need to sort it and create and index of the bam file 
#echo "sambamba sort"
#sambamba sort --memory-limit=1000GB --tmpdir=./tmp_sort -o ${ALN_PREFIX}_sort.bam -p -t $THREADS_NUM ${ALN_PREFIX}.bam
echo "samtools sort"
samtools sort -@ $THREADS_NUM -m 25GB ${ALN_PREFIX}.bam -o ${ALN_PREFIX}_sort.bam
rm ${ALN_PREFIX}.bam
echo "samtools index"
samtools index ${ALN_PREFIX}_sort.bam