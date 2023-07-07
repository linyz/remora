#! /usr/bin/env bash

env

# STAR index loc.
genomeDir=~/STAR_index/GRCh38_deaminases
# parent folder dir
pdir=~/myscratch/Novaseq202301/
# each pair of reads are in a same subfolder, folder name were inferred from fastq.gz names
subdir=${PWD##*/}

echo $PWD

# load mudules
module load CBI
module load star/2.7.9a
module load samtools/1.13

# remove tmp dir if already there
rm -r ${pdir}${subdir}/_STARtmp

# unzip
gunzip *.fastq.gz

# run STAR
STAR --genomeDir $genomeDir \
        --outTmpDir ${pdir}${subdir}/_STARtmp \
        --readFilesIn ${pdir}${subdir}/${subdir}_L002_R1_001.fastq ${pdir}${subdir}/${subdir}_L002_R2_001.fastq \
        --outFileNamePrefix ${pdir}${subdir}/${subdir}_ \
        --outReadsUnmapped Fastx \
        --clip3pAdapterSeq "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
        --clip3pAdapterMMp 0.1 0.1 \
        --runThreadN 4
        
samtools view -@ 4 -bS -o ${pdir}${subdir}/${subdir}_Aligned.out.bam ${pdir}${subdir}/${subdir}_Aligned.out.sam

samtools sort -@ 4 -o ${pdir}${subdir}/${subdir}_Aligned.out.sorted.bam ${pdir}${subdir}/${subdir}_Aligned.out.bam

samtools index -@ 4 ${pdir}${subdir}/${subdir}_Aligned.out.sorted.bam

# remove dups and calibrate
module load picard/2.27.1
module load gatk/4.2.6.1

picard MarkDuplicates -I  ${pdir}${subdir}/${subdir}_Aligned.out.sorted.bam \
     -M  ${pdir}${subdir}/${subdir}_picardMarkDup.txt \
     -O  ${pdir}${subdir}/${subdir}_Aligned.out.sortedMD.bam

gatk SplitNCigarReads \
      -R ${genomeFa} \
      -I ${pdir}${subdir}/${subdir}_Aligned.out.sortedMD.bam \
      -O ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit.bam

# seperate bam files by reads strand
samtools view -@ 4 -b -f 99 ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit.bam > ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_R1F.bam
samtools view -@ 4 -b -f 147 ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit.bam > ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_R2R.bam

samtools view -@ 4 -b -f 83 ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit.bam > ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_R1R.bam
samtools view -@ 4 -b -f 163 ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit.bam > ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_R2F.bam

# merge bam files by reads strand
samtools merge -@ 4 -o ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_FWD.bam \
 ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_R1F.bam ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_R2R.bam

samtools merge -@ 4 -o ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_REV.bam \
 ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_R1R.bam ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_R2F.bam

# samtools index
samtools index -@ 4 ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_FWD.bam
samtools index -@ 4 ${pdir}${subdir}/${subdir}_Aligned.out.sortedMDSplit_REV.bam

# end job
qstat -j $JOB_ID

