1. Making STAR index and other refs:
    - concat deaminases sequence with GRCh38
    - concat deaminases gtf annotation with gencode v40 basic gtf
    - samtools faidx
    - picard CreateSequenceDictionary
    - make chromosome bed files (for run mpileup in parallel)
2. run run_STAR2.sh
    - STAR alignment: need to trim 3' illumina universal adapter
    - samtools sorting etc
    - picard MarkDuplicate
    - gatk SplitNCigarReads
3. Use Rsubreads for gene expression count and DE analysis
4. run bcftools mpileup for variants calling
5. Custom R script (vcf_preprocess.R) to reformat vcf files and annotate genes
    - Remove ##contig lines to only keep chr1-23, X, Y, M; 