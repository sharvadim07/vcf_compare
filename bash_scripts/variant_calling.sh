#!/bin/bash

set -e

for f in `ls *_sort.bam`;
do
#smpl=${f%.*}
smpl=${f%%.*}
#samtools mpileup -uf ../Trinity.fasta ${smpl}.bam | bcftools call -Am | bcftools filter -s LowQual -e '%QUAL<20 || DP>2' > ${smpl}_raw.vcf
samtools mpileup -uf ../Trinity.fasta ${smpl}.bam | bcftools call -c | bcftools view  -i 'MIN(INFO/DP)>1' > ${smpl}_raw.vcf
done