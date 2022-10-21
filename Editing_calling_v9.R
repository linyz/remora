## Yizhu Lin, 10/21/2022 
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(ggplot2)
library(reshape2)
library(GGally)
library("BSgenome.Hsapiens.UCSC.hg38")
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
hgenome <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(tidyr)
library(dplyr)

library(magrittr)
library(tidyverse)
library(gtools)

# cmh test
rowcmh <- function(x, myorder=c(10,7,4,1, 11,8,5,2, 12,9,6,3), dim=c(2,2,3)){
  # DP-AD
  # c(7,1,8,2,9,3,  10,4,11,5,12,6)c(10,7,4,1, 11,8,5,2, 12,9,6,3)
  ratios <- x[7:12]/x[1:6]
  Tratios <- mean(ratios[4:6]) 
  Cratios <- mean(ratios[1:3])
  diffratios = Tratios - Cratios
  x[1:6 ] = x[1:6]-x[7:12]
  # reformat to 3d array
  te <- array(t(as.numeric(x[myorder])), dim=dim)
  # print(te)
  cmh <- mantelhaen.test(te)
  return(c(cmh$statistic, cmh$p.value, cmh$estimate, Tratios, Cratios, diffratios))
}

##================ cmh test for all A-to-I editings ========================##
for (setName in c('RBFox2-rABE', 'PUM2-rABE','PUM1-rABE', 'dual-S15', 'dual-S16')){
  # for (setName in c( 'dual-S15', 'dual-S16')){
  refnt = 'A'
  mutnt = 'G'
  for (dox in c('H', 'L')){
    
    setNameInput = setName
    if (setName == 'dual-S15' | setName == 'dual-S16'){setNameInput = 'dual'}
    df=read.csv(sprintf("%s_anno_withexon_all.vcf",setNameInput))

    # only keeps most abundat exon annotation
    df <- df %>%
      group_by(exonidx, seqnames) %>%
      mutate(avgExonDP = mean(avgDP)) %>%
      group_by(vcfidx, seqnames) %>%
      dplyr::slice(which.max(avgExonDP)) %>%
      ungroup()
    
    # filtering for A-to-I / C-to-U
    outallg <- df[df$RNAREF==refnt & df$RNAALT==mutnt, ]
    
    dfcount <- outallg %>%
      dplyr::select(contains('.DP') | contains('.AD')) %>%
      dplyr::select(contains(sprintf('.%s.', dox)))
    cols <- colnames(dfcount)
    
    dfinfos <- outallg %>%
      dplyr::select(!contains('.DP')) %>%
      dplyr::select(!contains('.AD')) %>%
      dplyr::select(!contains('.ratio'))
    cmht <- data.frame(t(apply(as.matrix(dfcount), 1, FUN=rowcmh)))
    if (setName == 'dual-S15'){
      cmht <- data.frame(t(apply(as.matrix(dfcount[,c(1,2,3,7,8,9,13,14,15,19,20,21)]), 1, FUN=rowcmh))) # sample 15
    } else if (setName == 'dual-S16'){
      cmht <- data.frame(t(apply(as.matrix(dfcount[,c(1,2,3,10,11,12,13,14,15,22,23,24)]), 1, FUN=rowcmh))) #sample 16
    }
    colnames(cmht) <-c("chi-squared", "pval", "OR","Tratios", "Cratios", "diffratios")
    cmht <- cmht %>% 
      mutate(padj = p.adjust(pval, method='BH'))
    cmht$de <- "ns"
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > 1, "de"] <- "up"
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1, "de"] <- "down"
    table(cmht$de)
    rownames(cmht) <- rownames(outallg)
    dfcmh <- cbind(dfinfos, dfcount, cmht) %>%
      drop_na()

    tmp=gsub("\\..*","",dfcmh$TXNAME)
    dfcmh$GENENAME <- mapIds(edb,keys = tmp,keytype="TXNAME", column=c("GENENAME"))
    write.csv(dfcmh,sprintf("%s_cmh_%s-%s-dox%s.csv",setName, refnt, mutnt, dox))

  }
}

##================ cmh test for all C-to-U editings ========================##
for (setName in c('RBFox2-BE4', 'PUM2-BE4','PUM1-BE4', 'dual-S15', 'dual-S16')){
  refnt = 'C'
  mutnt = 'T'
  for (dox in c('H', 'L')){
    
    setNameInput = setName
    if (setName == 'dual-S15' | setName == 'dual-S16'){setNameInput = 'dual'}
    df=read.csv(sprintf("%s_anno_withexon_all.vcf",setNameInput))
    
    # only keeps most abundant exon annotation
    df <- df %>%
      group_by(exonidx, seqnames) %>%
      mutate(avgExonDP = mean(avgDP)) %>%
      group_by(vcfidx, seqnames) %>%
      dplyr::slice(which.max(avgExonDP)) %>%
      ungroup()
    
    # filtering for A-to-I / C-to-U
    outallg <- df[df$RNAREF==refnt & df$RNAALT==mutnt, ]
    
    dfcount <- outallg %>%
      dplyr::select(contains('.DP') | contains('.AD')) %>%
      dplyr::select(contains(sprintf('.%s.', dox)))
    cols <- colnames(dfcount)
    
    dfinfos <- outallg %>%
      dplyr::select(!contains('.DP')) %>%
      dplyr::select(!contains('.AD')) %>%
      dplyr::select(!contains('.ratio'))
    cmht <- data.frame(t(apply(as.matrix(dfcount), 1, FUN=rowcmh)))
    if (setName == 'dual-S15'){
      cmht <- data.frame(t(apply(as.matrix(dfcount[,c(4,5,6,7,8,9,16,17,18,19,20,21)]), 1, FUN=rowcmh))) # sample 15
    } else if (setName == 'dual-S16'){
      cmht <- data.frame(t(apply(as.matrix(dfcount[,c(4,5,6,10,11,12,16,17,18,22,23,24)]), 1, FUN=rowcmh))) #sample 16
    }
    colnames(cmht) <-c("chi-squared", "pval", "OR","Tratios", "Cratios", "diffratios")
    cmht <- cmht %>% 
      mutate(padj = p.adjust(pval, method='BH'))
    cmht$de <- "ns"
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > 1, "de"] <- "up"
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1, "de"] <- "down"
    table(cmht$de)
    rownames(cmht) <- rownames(outallg)
    dfcmh <- cbind(dfinfos, dfcount, cmht) %>%
      drop_na()
    
    tmp=gsub("\\..*","",dfcmh$TXNAME)
    dfcmh$GENENAME <- mapIds(edb,keys = tmp,keytype="TXNAME", column=c("GENENAME"))
    write.csv(dfcmh,sprintf("%s_cmh_%s-%s-dox%s.csv",setName, refnt, mutnt, dox))
  }
}
