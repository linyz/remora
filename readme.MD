1. Making STAR index and other refs:
    - (optional for DE analysis) concat deaminases sequence with GRCh38
    - (optional for DE analysis) concat deaminases gtf annotation with gencode v40 basic gtf
    - samtools faidx
    - picard CreateSequenceDictionary
    - make chromosome bed files (for run mpileup in parallel)
        ./chrBed, chromosome length from GRCh38
2. run run_STAR.sh
    - STAR alignment: need to trim 3' illumina universal adapter
    - samtools sorting etc
    - picard MarkDuplicate
    - gatk SplitNCigarReads
    - split by read orientation

3. run bcftools mpileup for variants calling. I split up files by chromosome and strand for parallel computing.
4. Custom R script (vcf_preprocess.R) to reformat vcf files, convert DNA alt to RNA alt, and basic filtering
    - If refs have extra chromosomes, remove ##contig lines to only keep chr1-23, X, Y, M;
5. concat all [strand_chr].vcf files for downstream editing site calling. (so multiple test correction is done tx-wide)
6. Edit_call_cmh.R
