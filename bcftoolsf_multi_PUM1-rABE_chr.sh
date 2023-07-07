#!/bin/bash

module load CBI
module load bcftools

bcftools mpileup -f~/genomes/human/GRCh38.primary_assembly_mod.genome.fa -R ~/chrBed/chr$1.bed -d 10000000 -I -a DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR \
    S01-1-rABE-Flag_S74_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S01-2-rABE-Flag_S75_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S01-3-rABE-Flag_S76_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S01-4-rABE-Flag_S77_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S01-5-rABE-Flag_S78_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S01-6-rABE-Flag_S79_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S01-7-rABE-Flag_S80_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S01-8-rABE-Flag_S81_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S01-9-rABE-Flag_S82_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-1-PUM1-rABE_S128_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-2-PUM1-rABE_S129_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-3-PUM1-rABE_S130_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-4-PUM1-rABE_S131_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-5-PUM1-rABE_S132_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-6-PUM1-rABE_S133_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-7-PUM1-rABE_S134_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-8-PUM1-rABE_S135_v3-Aligned.out.sortedMDSplit_FWD.bam \
    S07-9-PUM1-rABE_S136_v3-Aligned.out.sortedMDSplit_FWD.bam \
    | bcftools filter -i 'INFO/AD[1-]>2 & MAX(FORMAT/DP)>20' -O v - > mpileupfilter_neg_PUM1-rABE_chr$1.vcf


bcftools mpileup -f~/genomes/human/GRCh38.primary_assembly_mod.genome.fa -R ~/chrBed/chr$1.bed -d 10000000 -I -a DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR \
    S01-1-rABE-Flag_S74_v3-Aligned.out.sortedMDSplit_REV.bam \
    S01-2-rABE-Flag_S75_v3-Aligned.out.sortedMDSplit_REV.bam \
    S01-3-rABE-Flag_S76_v3-Aligned.out.sortedMDSplit_REV.bam \
    S01-4-rABE-Flag_S77_v3-Aligned.out.sortedMDSplit_REV.bam \
    S01-5-rABE-Flag_S78_v3-Aligned.out.sortedMDSplit_REV.bam \
    S01-6-rABE-Flag_S79_v3-Aligned.out.sortedMDSplit_REV.bam \
    S01-7-rABE-Flag_S80_v3-Aligned.out.sortedMDSplit_REV.bam \
    S01-8-rABE-Flag_S81_v3-Aligned.out.sortedMDSplit_REV.bam \
    S01-9-rABE-Flag_S82_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-1-PUM1-rABE_S128_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-2-PUM1-rABE_S129_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-3-PUM1-rABE_S130_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-4-PUM1-rABE_S131_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-5-PUM1-rABE_S132_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-6-PUM1-rABE_S133_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-7-PUM1-rABE_S134_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-8-PUM1-rABE_S135_v3-Aligned.out.sortedMDSplit_REV.bam \
    S07-9-PUM1-rABE_S136_v3-Aligned.out.sortedMDSplit_REV.bam \
    | bcftools filter -i 'INFO/AD[1-]>2 & MAX(FORMAT/DP)>20' -O v - > mpileupfilter_pos_PUM1-rABE_chr$1.vcf