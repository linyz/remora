remove(list = ls())
gc()
setwd('D://recolor/NovaseqAll/vcf/')
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
library(tidyr)
library(dplyr)
# library(ggplot2)
# library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
s=as.numeric(args[1]) ## strand, pos or neg
k=as.numeric(args[2])  # sample group
j=as.numeric(args[3])  # chr

s=1
k=1
j=25

if (s == 1){
  st = 'pos'
} else if (s == 2){
  st = 'neg'
}


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

chrList = seqlevels0(txdb)
seqlevels(txdb) <- chrList[j]

sList = c('C.N.1','C.N.2','C.N.3','C.L.1','C.L.2','C.L.3','C.H.1','C.H.2','C.H.3',
          'T.N.1','T.N.2','T.N.3','T.L.1','T.L.2','T.L.3','T.H.1','T.H.2','T.H.3')

# hard coded file names for easier qsub
if (k==1){
  setName = 'RBFox2-hyperTRIBE'
} else if (k==2){
  setName = 'RBFox2-rABE'
} else if (k==3){
  setName = 'RBFox2-BE4'
} else if (k==4){
  setName = 'PUM1-rABE'
} else if (k==5){
  setName = 'PUM1-BE4'
} else if (k==6){
  setName = 'PUM2-rABE'
} else if (k==7){
  setName = 'PUM2-BE4'
} else if (k==8){
  setName = 'dual'
  sList = c('C.N.A.1','C.N.A.2','C.N.A.3','C.L.A.1','C.L.A.2','C.L.A.3','C.H.A.1','C.H.A.2','C.H.A.3',
            'C.N.B.1','C.N.B.2','C.N.B.3','C.L.B.1','C.L.B.2','C.L.B.3','C.H.B.1','C.H.B.2','C.H.B.3',
            'T.N.A.1','T.N.A.2','T.N.A.3','T.L.A.1','T.L.A.2','T.L.A.3','T.H.A.1','T.H.A.2','T.H.A.3',
            'T.N.B.1','T.N.B.2','T.N.B.3','T.L.B.1','T.L.B.2','T.L.B.3','T.H.B.1','T.H.B.2','T.H.B.3')
} else if (k==9){
  setName = 'RBFox2-rABE-12h'
} else if (k==10){
  setName = 'RBFox2-rABE-24h'
} else if (k==11){
  setName = 'RBFox2-rABE-48h'
}



# read vcf file from mpileup, annotate strand, type (cds/utr etc.)
# Fix mutation by strand
# extract read depth (DP) and allele depth (AD)
vcf <- readVcf( sprintf("mpileup/v3_mpileupfilter_%s_%s_chr%s.vcf", st, setName, j))
rd <- rowRanges(vcf)
# read in as table for easy cleaning up
vcfraw <- read.table(sprintf("mpileup/v3_mpileupfilter_%s_%s_chr%s.vcf", st, setName, j))

# read gene anntations on current chr
alltx <- transcripts(txdb,columns=c('TXNAME', "GENEID", "TXSTRAND", "TXTYPE"), filter=NULL)
# TODO: may need to check strand orientation depends on library prep, assign strand accordingly
if (st=='pos'){
  alltx = alltx[strand(alltx) == "-"]
  strand = 'pos'
} else if (st=='neg') {
  alltx = alltx[strand(alltx) == "+"]
  strand = 'neg'
}
vcfraw$strand <- strand

# keytypes(txdb)

# intersct vcf with txs
# hits = findOverlaps(rd, alltx)
# vcfraw$TXNAME <- sapply(split(mcols(alltx)$TXNAME, queryHits(hits)), paste, collapse = ",")
# tmp = gsub('\\..*', "", vcfraw$TXNAME)
# vcfraw$TXNAME <- tmp

# vcfraw$GENENAME <- mapIds(edb,keys =tmp,keytype="TXNAME", column=c("GENENAME")) 
# vcfraw$ensGENEID <- mapIds(edb,keys =tmp,keytype="TXNAME", column=c("GENEID")) 
# vcfraw$SYMBOL <- mapIds(edb,keys =tmp,keytype="TXNAME", column=c("SYMBOL")) 

# only keep first ALT 
vcfraw$V5 <- substr(vcfraw$V5, 0, 1)

loccds <- locateVariants(rd, txdb, CodingVariants())
locfiveutr <- locateVariants(rd, txdb, FiveUTRVariants())
locthreeutr <- locateVariants(rd, txdb, ThreeUTRVariants())
# MT has no intron
if (j!=25){
locintron <- locateVariants(rd, txdb, IntronVariants())
}


vcfraw$cds <- seq(length(rd)) %in% loccds$QUERYID
vcfraw$fiveutr <- seq(length(rd)) %in% locfiveutr$QUERYID
vcfraw$threeutr <- seq(length(rd)) %in% locthreeutr$QUERYID
if (j!=25){
vcfraw$intron <- seq(length(rd)) %in% locintron$QUERYID
} else {vcfraw$intron = FALSE}
remove(list=c("loccds", "locfiveutr", "locthreeutr", "locintron"))
gc()

ncol = ncol(vcfraw)
out <- vcfraw %>% drop_na()

# filter out ambi. strand site
# out <- out[(out$posstrand + out$negstrand) == 1, ]
# replace ATCG to rv on neg strand
out$RNAREF <- out$V4
out$RNAALT <- out$V5
if (strand=='neg'){
  out[out$V4=='A', 'RNAREF'] <- 'T'
  out[out$V4=='C', 'RNAREF'] <- 'G'
  out[out$V4=='G', 'RNAREF'] <- 'C'
  out[out$V4=='T', 'RNAREF'] <- 'A'
  out[out$V5=='A', 'RNAALT'] <- 'T'
  out[out$V5=='C', 'RNAALT'] <- 'G'
  out[out$V5=='G', 'RNAALT'] <- 'C'
  out[out$V5=='T', 'RNAALT'] <- 'A'
}
# colnames(out)

# set up treated vs control, need to match file order in vcf, which is bcftools mpileup input order
vList = (colnames(out))[10:ncol]
print(vList)
for (i in seq(length(sList))){
  Vn = vList[i]
  Sp = sList[i]
  out <- out %>% 
    separate(Vn,sep=':', c(NA, paste(Sp,"DP",sep='.'), NA, NA, NA, paste(Sp,"AD",sep='.')), remove=TRUE,convert=TRUE) %>%
    separate(paste(Sp,"AD",sep='.'),sep=',', c(NA, paste(Sp,"AD",sep='.')), remove=TRUE, convert = TRUE)
}
# 

# calc mutation rate 
vars <- seq(10, (9+length(sList)*2))
varsDP = vars[seq(1,length(vars),2)]
varsAD = vars[seq(2,length(vars),2)]
for (i in seq(length(varsDP))){
  Sp = sList[i]
  out$ratio <- as.numeric(out[,varsAD[i]])/as.numeric(out[,varsDP[i]])
  names(out)[names(out) == 'ratio'] <- paste(Sp,'ratio', sep = ".")
}

# remove muts exists in no dox controls (germline / endo editing)
cutoff = 0.01
if (k != 8) {
  out <- out[(out$C.N.1.ratio<cutoff & out$C.N.2.ratio <cutoff & out$C.N.3.ratio<cutoff), ]
} else if (k == 2) {
  # RBFox2-rABE-nodox rep1 is actually low dox
  out <- out[(out$C.N.2.ratio <cutoff & out$C.N.3.ratio<cutoff), ]
} else {
  out <- out[(out$C.N.A.1.ratio<cutoff & out$C.N.A.2.ratio <cutoff & out$C.N.A.3.ratio<cutoff & out$C.N.B.1.ratio<cutoff & out$C.N.B.2.ratio <cutoff & out$C.N.B.3.ratio<cutoff), ]
}
out$avgDP <- rowMeans(as.matrix(out[,varsDP]))

write.table(out,file=sprintf('%s_%s_anno_chr%s.vcf',setName, strand, j), quote = FALSE,sep='\t', row.names=FALSE,col.names = TRUE)

