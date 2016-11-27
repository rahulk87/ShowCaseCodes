######################################################
## This code generates logr and baf files for ASCAT ##
######################################################

#!/bin/bash
#BSUB -u rahul.kumar@icr.ac.uk
#BSUB -J seqgenes_BX078
#BSUB -e /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky/analysis/BX078/ascat/seq.err
#BSUB -o /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky/analysis/BX078/ascat/seq.out
#BSUB -n 16
#BSUB -P DBCDOBZAK
#BSUB -q normal
#

module load samtools/1.3
module load bedtools

path1="/scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky/analysis/BX078/BX078_germline"
path2="/scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky/analysis/BX078/ascat"
path3="/scratch/DBC/GENFUNC/NGS_Projects/scripts"
path4="/scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky/analysis/BX078/gatk_joint_calling"

samtools view $path1/BX078_germline.recaliberated.bam > $path2/BX078_germline.recaliberated.sam

$path3/sam2wig.py $path2/BX078_germline.recaliberated.sam $path2/BX078_germline.recaliberated.wig

$path3/rpkm.py -o human -m cnv BX078_germline.recaliberated.wig $path2/BX078_germline.recaliberated.rpkm

awk '$1 ~/^[0-9]/ || ($1 ~/^chro/) || ($1 ~/^X/) || ($1 ~/^Y/)' $path2/BX078_germline.recaliberated.rpkm \
 > $path2/BX078_germline.recaliberated.rpkm.clean

cut -f 1,3,4,7,10 -s $path2/BX078_germline.recaliberated.rpkm.clean | tail -n +2 \
 > $path2/BX078_germline.recaliberated.rpkm.clean.bed

intersectBed -a $path4/BX078.Recalibrator.vcf \
 -b $path2/BX078_germline.recaliberated.rpkm.clean.bed -wo > $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm

sort -k1n -k2n $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm \
 > $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm.sorted

echo BX078_germline
COL=`grep CHROM /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky/analysis/BX078/gatk_joint_calling/BX078.Recalibrator.vcf | tr '\t' '\n' | grep BX078_germline -n | cut -f 1 -d ':'`
echo $COL

COL=$(expr $COL - 1)
echo $COL

$path3/get_allele_freqs.pl --vcf $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm.sorted --columns $COL --out $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm.sorted.vaf

FRE=`head -1 /scratch/DBC/GENFUNC/NGS_Projects/PDXs_Becky/analysis/BX078/ascat/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm.sorted.vaf | tr '\t' '\n' | wc -l`
GENE=`expr $FRE - 5`
LOG=`expr $FRE - 4`

cut -f 1,2,$GENE,$LOG,$FRE $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm.sorted.vaf \
 > $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm.sorted.vaf.format

awk '{print $3 "_" $1 "_" $2 "\t" $1 "\t" $2 "\t" $5}' $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm.sorted.vaf.format | uniq \
 > $path2/BX078_germline.baf

awk '{print $3 "_" $1 "_" $2 "\t" $1 "\t" $2 "\t" $4}' \
 $path2/BX078_germline.recaliberated.rpkm.clean.bed.vcfWithRpkm.sorted.vaf.format | uniq > $path2/BX078_germline.logr

