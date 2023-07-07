remove(list = ls())
gc()
setwd('D://recolor/NovaseqAll/vcf/')
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
library(ggseqlogo)
library(tidyr)
library(dplyr)
library(ggpointdensity)

library(magrittr)
library(edgeR)
library(tidyverse)
library(VennDiagram)
library(gtools)
library(epitools)
library(ggpubr)

### -------------mut type summary ----- #####
for (setName in c('RBFox2-hyperTRIBE','RBFox2-rABE', 'RBFox2-BE4','PUM2-rABE','PUM2-BE4','PUM1-rABE','PUM1-BE4','RBFox2-rABE-12h', 'RBFox2-rABE-24h','RBFox2-rABE-48h', 'dual')){
  setName = 'allPUM'
  df <- read.table(sprintf("%s_all_v3.vcf",setName), header = TRUE)

  ## mut type summary
  for (cond in c('T.H', 'T.L', 'T.N', 'C.H', 'C.L')){
  # for (cond in c('T.H.A', 'T.L.A', 'T.N.A', 'T.H.B', 'T.L.B', 'T.N.B')){
    # for (cond in c('T.H.D', 'T.L.D', 'T.N.D', 'C.H.B', 'C.L.B', 'C.N.B')){
    cutoffAR=0.005
    cutoffDP= 100
    out2 <- df %>% 
      dplyr::select(contains(cond)) %>%
      dplyr::select(contains("ratio")) %>%
      mutate_if(is.character, as.numeric) %>%
      mutate(avgRatio = rowMeans(.,na.rm=TRUE))
    out2$RNAREF <- df$RNAREF
    out2$RNAALT <- df$RNAALT
    out2$avgDP <- df$avgDP
    ggplot(out2[out2$avgDP>cutoffDP & out2$avgRatio>cutoffAR,], aes(x=RNAALT, y=avgRatio)) +
      # ggplot(tp, aes(group, YL0001)) +
      # geom_bar(stat='identity') +
      ggtitle(paste(setName, cond)) +
      geom_boxplot() +
      facet_grid(.~RNAREF) +
      ylab("Eiditing Rate")+
      # theme_classic()+
      scale_y_continuous(trans='log10', limits=c(0.001,1))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(sprintf("boxplot_avgEditRate_%sDP_%sAR_%s_%s.pdf",cutoffDP, cutoffAR,setName, cond), width = 15, height = 6, units = "cm")
    out3 <- out2[out2$avgRatio>cutoffAR & out2$avgDP > cutoffDP,]
    print(paste(setName, cond))
    t <- table(out3[,c('RNAREF', 'RNAALT')]); sum(t); t/sum(t)
    print(t)
  }
}


rowcmh <- function(x){
  nrep = length(x)/4
  # c(7,1,8,2,9,3,  10,4,11,5,12,6)
  ratios <- x[(2*nrep+1):(4*nrep)]/x[1:(2*nrep)]
  Tratios <- mean(ratios[(nrep+1):(2*nrep)]) 
  Cratios <- mean(ratios[1:nrep])
  diffratios = Tratios - Cratios
  
  if (nrep == 1){# chi-squared
    re <- chisq.test(x) # not necessary to reorder
    if (x[3] == 0 & x[4] == 0){return(c(NA, 1, NA, 0,0,0))}
    if (x[2]*x[3] == 0){OR = Inf}else{
      OR <- (x[4]*x[1]) / (x[2]*x[3])
    }
    return(c(re$statistic, re$p.value, OR, Tratios, Cratios, diffratios))
  } 
  else {# if n >= 2, Cochran-Mantel-Haenszel
    myorder = c(1+3*nrep, 1+2*nrep, 1+nrep, 1)
    for (i in seq(2,nrep)){
      myorder=c(myorder, i+3*nrep, i+2*nrep, i+nrep, i)
    }
    dim = c(2,2,nrep)
    # DP-AD
    x[1:(2*nrep) ] = x[1:(2*nrep)]-x[(2*nrep+1):(4*nrep)]
    # reformat to 3d array
    te <- array(t(as.numeric(x[myorder])), dim=dim)
    # print(te)
    re <- mantelhaen.test(te)
    return(c(re$statistic, re$p.value, re$estimate, Tratios, Cratios, diffratios))
  }
  
}


DPcutoff=50
ORcutoff=1.2
dEcutoff = 0.005 # diff editing ratio cutoff

##================ cmh test for all A-to-I editings ========================##

for (setName in c('RBFox2-hyperTRIBE','RBFox2-rABE','PUM2-rABE','PUM1-rABE','RBFox2-rABE-12h', 'RBFox2-rABE-24h','RBFox2-rABE-48h','dual-S15', 'dual-S16')){

  refnt = 'A'
  mutnt = 'G'
  for (dox in c('H', 'L')){
    dox = 'H'
    setNameInput = setName
    if (setName == 'dual-S15' | setName == 'dual-S16'){setNameInput = 'dual'}
    df=read.csv(sprintf("%s_all_v3.vcf",setNameInput), header=TRUE)
    
    # filtering for A-to-I / C-to-U
    outallg <- df[df$RNAREF==refnt & df$RNAALT==mutnt, ] 
    
    DPs <- outallg %>% 
      dplyr::select(contains('.DP')) %>%
      dplyr::select(contains(sprintf('.%s.', dox))) 
    outallg$DPmin <- apply(DPs, MARGIN =  1, FUN = min, na.rm = T)
    outallg <- outallg[outallg$DPmin >0 & outallg$avgDP>DPcutoff, ]
    
    dfcount <- outallg %>%
      dplyr::select(contains('.DP') | contains('.AD')) %>%
      dplyr::select(contains(sprintf('.%s.', dox)))
    
    cols <- colnames(dfcount)
    
    dfinfos <- outallg %>%
      dplyr::select(!contains('.DP')) %>%
      dplyr::select(!contains('.AD')) %>%
      dplyr::select(!contains('.ratio'))

    if (setName == 'dual-S15'){
      cmht <- data.frame(t(apply(as.matrix(dfcount[,c(1,2,3,7,8,9,13,14,15,19,20,21)]), 1, FUN=rowcmh))) # sample 15
    } else if (setName == 'dual-S16'){
      cmht <- data.frame(t(apply(as.matrix(dfcount[,c(1,2,3,10,11,12,13,14,15,22,23,24)]), 1, FUN=rowcmh))) #sample 16
    } else {
      ### see bcftools step for column order
      cmht <- data.frame(t(apply(as.matrix(dfcount), 1, FUN=rowcmh)))
    }
    colnames(cmht) <-c("chi-squared", "pval", "OR","Tratios", "Cratios", "diffratios")
    cmht <- cmht %>% 
      mutate(padj = p.adjust(pval, method='BH'))
    cmht$de <- "ns"
    # cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > ORcutoff, "de"] <- "up"
    # cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1.0/ORcutoff, "de"] <- "down"
    
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > ORcutoff & cmht$diffratios > dEcutoff, "de"] <- "up"
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1/ORcutoff & cmht$diffratios < -dEcutoff, "de"] <- "down"
    
    print(setName)
    print(table(cmht$de))
    rownames(cmht) <- rownames(outallg)
    dfcmh <- cbind(dfinfos, dfcount, cmht)

    write.table(dfcmh,sprintf("%s_cmh_%s-%s-dox%s_DP%s_OR%s_dE%s.vcf",setName, refnt, mutnt, dox, DPcutoff, ORcutoff, dEcutoff))
    
    dfcmh$color <- 'grey'
    if (nrow(dfcmh[dfcmh$de=='up', ])>0){
      dfcmh[dfcmh$de=='up', ]$color = 'red'}
    if (nrow(dfcmh[dfcmh$de=='down', ])>0){
      dfcmh[dfcmh$de=='down', ]$color = 'blue'}
    ggplot(data = dfcmh, aes(x=avgDP, y=diffratios)) +
      geom_point(alpha=0.5,size=0.3, color=dfcmh$color) +    
      labs( y = "Editing rate", x = "Read depth") + 
      scale_x_continuous(trans='log10') + 
      ylim(-0.8, 0.8)+
      #scale_fill_manual(values=c(FALSE="black",TRUE="red")) + 
      theme_bw()
    ggsave(sprintf("diffEditing_scatter_cmh_padj005_%s_%s-%s-dox%s_DP%s_OR%s_dE%s.png", setName, refnt, mutnt, dox, DPcutoff, ORcutoff, dEcutoff), width = 6, height = 5, units = "cm")
  }
}


##================ cmh test for all C-to-U editings ========================##
for (setName in c('RBFox2-BE4', 'PUM2-BE4','PUM1-BE4', 'dual-S15', 'dual-S16')){
  refnt = 'C'
  mutnt = 'T'
  for (dox in c('H', 'L')){
    setNameInput = setName
    if (setName == 'dual-S15' | setName == 'dual-S16'){setNameInput = 'dual'}
    df=read.csv(sprintf("%s_all_v3.vcf",setNameInput), header=TRUE)
    
    # filtering for A-to-I / C-to-U
    outallg <- df[df$RNAREF==refnt & df$RNAALT==mutnt, ] 
    
    DPs <- outallg %>% 
      dplyr::select(contains('.DP')) %>%
      dplyr::select(contains(sprintf('.%s.', dox))) 
    outallg$DPmin <- apply(DPs, MARGIN =  1, FUN = min, na.rm = T)
    outallg <- outallg[outallg$DPmin >0 & outallg$avgDP>DPcutoff, ]
    
    dfcount <- outallg %>%
      dplyr::select(contains('.DP') | contains('.AD')) %>%
      dplyr::select(contains(sprintf('.%s.', dox)))
    
    cols <- colnames(dfcount)
    
    dfinfos <- outallg %>%
      dplyr::select(!contains('.DP')) %>%
      dplyr::select(!contains('.AD')) %>%
      dplyr::select(!contains('.ratio'))
    
    if (setName == 'dual-S15'){
      cmht <- data.frame(t(apply(as.matrix(dfcount[,c(4,5,6,7,8,9,16,17,18,19,20,21)]), 1, FUN=rowcmh))) # sample 15
    } else if (setName == 'dual-S16'){
      cmht <- data.frame(t(apply(as.matrix(dfcount[,c(4,5,6,10,11,12,16,17,18,22,23,24)]), 1, FUN=rowcmh))) #sample 16
    } else {
      cmht <- data.frame(t(apply(as.matrix(dfcount), 1, FUN=rowcmh)))
    }
    colnames(cmht) <-c("chi-squared", "pval", "OR","Tratios", "Cratios", "diffratios")
    cmht <- cmht %>% 
      mutate(padj = p.adjust(pval, method='BH'))
    cmht$de <- "ns"
    # cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > ORcutoff, "de"] <- "up"
    # cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1.0/ORcutoff, "de"] <- "down"
    
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > ORcutoff & cmht$diffratios > dEcutoff, "de"] <- "up"
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1/ORcutoff & cmht$diffratios < -dEcutoff, "de"] <- "down"
    
    print(setName)
    print(table(cmht$de))
    rownames(cmht) <- rownames(outallg)
    dfcmh <- cbind(dfinfos, dfcount, cmht)
    
    write.table(dfcmh,sprintf("%s_cmh_%s-%s-dox%s_DP%s_OR%s_dE%s.vcf",setName, refnt, mutnt, dox, DPcutoff, ORcutoff, dEcutoff))

    dfcmh$color <- 'grey'
    if (nrow(dfcmh[dfcmh$de=='up', ])>0){
      dfcmh[dfcmh$de=='up', ]$color = 'red'}
    if (nrow(dfcmh[dfcmh$de=='down', ])>0){
      dfcmh[dfcmh$de=='down', ]$color = 'blue'}
    ggplot(data = dfcmh, aes(x=avgDP, y=diffratios)) +
      geom_point(alpha=0.5,size=0.3, color=dfcmh$color) +    
      labs( y = "Editing rate", x = "Read depth") + 
      scale_x_continuous(trans='log10') + 
      ylim(-0.8, 0.8)+
      #scale_fill_manual(values=c(FALSE="black",TRUE="red")) + 
      theme_bw()
    ggsave(sprintf("diffEditing_scatter_cmh_padj005_%s_%s-%s-dox%s_DP%s_OR%s_dE%s.png", setName, refnt, mutnt, dox, DPcutoff, ORcutoff, dEcutoff), width = 6, height = 5, units = "cm")
  }
}


###===========================deaminase only (background edits) A>G ======================###
for (setName in c('PUM1-rABE', 'RBFox2-hyperTRIBE')){
  #setName = 'RBFox2-hyperTRIBE'
  #setName = 'PUM1-rABE'
  
  refnt = 'A'
  mutnt = 'G'
  for (dox in c('H', 'L')){
    #dox='L'
    df <- read.csv(sprintf("%s_all_v3.vcf",setName), header=TRUE)
    # filtering for A-to-I / C-to-U
    outallg <- df[df$RNAREF==refnt & df$RNAALT==mutnt, ] 
    
    DPs <- outallg %>% 
      dplyr::select(contains('.DP')) %>%
      dplyr::select(contains('C.')) %>%
      dplyr::select(contains('.N.') | contains(sprintf('.%s.', dox))) 
    outallg$DPmin <- apply(DPs, MARGIN =  1, FUN = min, na.rm = T)
    outallg <- outallg[outallg$DPmin >0, ]
    
    dfcount <- outallg %>%
      dplyr::select(contains('C.')) %>%
      dplyr::select(contains('.N.') | contains(sprintf('.%s.', dox))) %>%
      dplyr::select(contains('.DP') | contains('.AD'))
    
    cols <- colnames(dfcount)
    
    dfinfos <- outallg %>%
      dplyr::select(!contains('.DP')) %>%
      dplyr::select(!contains('.AD')) %>%
      dplyr::select(!contains('.ratio'))
    
    cmht <- data.frame(t(apply(as.matrix(dfcount), 1, FUN=rowcmh)))

    colnames(cmht) <-c("chi-squared", "pval", "OR","Tratios", "Cratios", "diffratios")
    cmht <- cmht %>% 
      mutate(padj = p.adjust(pval, method='BH'))
    cmht$de <- "ns"
    #cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > 1, "de"] <- "up"
    #cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1, "de"] <- "down"
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > ORcutoff & cmht$diffratios > dEcutoff, "de"] <- "up"
    cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1/ORcutoff & cmht$diffratios < -dEcutoff, "de"] <- "down"
    print(setName)
    print(table(cmht$de))
    rownames(cmht) <- rownames(outallg)
    dfcmh <- cbind(dfinfos, dfcount, cmht)
    
    # write.table(dfcmh,sprintf("bg_%s_cmh_%s-%s-dox%s.vcf",setName, refnt, mutnt, dox))
    write.table(dfcmh,sprintf("bg_%s_cmh_%s-%s-dox%s_DP%s_OR%s_dE%s.vcf",setName, refnt, mutnt, dox, DPcutoff, ORcutoff, dEcutoff))
    #volcano
    ggplot(data=dfcmh, aes(x=log2(OR), y=-log10(pval)),size=2, col=factor(de))+
      geom_point()+
      geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="red") +
      geom_hline(yintercept=-log10(0.05), col="red") + 
      xlim(-10, 10)+
      labs(title = sprintf("%s_%s-%s_dox%s",setName, refnt, mutnt, dox))
    ggsave(sprintf("bg_volcanoplot_cmh_padj005_%s_%s-%s-dox%s.png", setName, refnt, mutnt, dox), width = 12, height = 6, units = "cm")
    
    dfcmh$color <- 'grey'
    if (nrow(dfcmh[dfcmh$de=='up', ])>0){
      dfcmh[dfcmh$de=='up', ]$color = 'red'}
    if (nrow(dfcmh[dfcmh$de=='down', ])>0){
      dfcmh[dfcmh$de=='down', ]$color = 'blue'}
    ggplot(data = dfcmh, aes(x=avgDP, y=diffratios)) +
      geom_point(alpha=0.5,size=0.5, color=dfcmh$color) +    
      # labs(title = sprintf("rABE", refnt, mutnt), y = "Editing rate", x = "Read depth") + 
      labs( y = "Editing rate", x = "Read depth") + 
      # scale_y_continuous(trans='log10') +
      scale_x_continuous(trans='log10') + 
      ylim(-0.8, 0.8)+
      #scale_fill_manual(values=c(FALSE="black",TRUE="red")) + 
      theme_bw()
    ggsave(sprintf("diffEditing_scatter_bg_cmh_padj005_%s_%s-%s-dox%s.png", setName, refnt, mutnt, dox), width = 6, height = 5, units = "cm")
  }}


#####------------------deaminase only (background edits) , C > U----------------#############
setName = 'PUM1-BE4'

refnt = 'C'
mutnt = 'T'
dox='L'
df <- read.csv(sprintf("%s_all_v3.vcf",setName), header=TRUE)
# filtering for A-to-I / C-to-U
outallg <- df[df$RNAREF==refnt & df$RNAALT==mutnt, ] 

DPs <- outallg %>% 
  dplyr::select(contains('.DP')) %>%
  dplyr::select(contains('C.')) %>%
  dplyr::select(contains('.N.') | contains(sprintf('.%s.', dox))) 
outallg$DPmin <- apply(DPs, MARGIN =  1, FUN = min, na.rm = T)
outallg <- outallg[outallg$DPmin >0, ]

dfcount <- outallg %>%
  dplyr::select(contains('C.')) %>%
  dplyr::select(contains('.N.') | contains(sprintf('.%s.', dox))) %>%
  dplyr::select(contains('.DP') | contains('.AD'))

cols <- colnames(dfcount)

dfinfos <- outallg %>%
  dplyr::select(!contains('.DP')) %>%
  dplyr::select(!contains('.AD')) %>%
  dplyr::select(!contains('.ratio'))

cmht <- data.frame(t(apply(as.matrix(dfcount), 1, FUN=rowcmh)))

colnames(cmht) <-c("chi-squared", "pval", "OR","Tratios", "Cratios", "diffratios")
cmht <- cmht %>% 
  mutate(padj = p.adjust(pval, method='BH'))
cmht$de <- "ns"
#cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > 1, "de"] <- "up"
#cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1, "de"] <- "down"
cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR > ORcutoff & cmht$diffratios > dEcutoff, "de"] <- "up"
cmht[(!is.na(cmht$padj)) & cmht$padj<0.05 & cmht$OR < 1/ORcutoff & cmht$diffratios < -dEcutoff, "de"] <- "down"
print(setName)
print(table(cmht$de))
rownames(cmht) <- rownames(outallg)
dfcmh <- cbind(dfinfos, dfcount, cmht)

# write.table(dfcmh,sprintf("bg_%s_cmh_%s-%s-dox%s.vcf",setName, refnt, mutnt, dox))
write.table(dfcmh,sprintf("bg_%s_cmh_%s-%s-dox%s_DP%s_OR%s_dE%s.vcf",setName, refnt, mutnt, dox, DPcutoff, ORcutoff, dEcutoff))
#volcano
ggplot(data=dfcmh, aes(x=log2(OR), y=-log10(pval)),size=2, col=factor(de))+
  geom_point()+
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  xlim(-10, 10)+
  labs(title = sprintf("%s_%s-%s_dox%s",setName, refnt, mutnt, dox))
ggsave(sprintf("bg_volcanoplot_cmh_padj005_%s_%s-%s-dox%s.png", setName, refnt, mutnt, dox), width = 12, height = 6, units = "cm")

dfcmh$color <- 'grey'
if (nrow(dfcmh[dfcmh$de=='up', ])>0){
  dfcmh[dfcmh$de=='up', ]$color = 'red'}
if (nrow(dfcmh[dfcmh$de=='down', ])>0){
  dfcmh[dfcmh$de=='down', ]$color = 'blue'}
ggplot(data = dfcmh, aes(x=avgDP, y=diffratios)) +
  geom_point(alpha=0.5,size=0.5, color=dfcmh$color) +    
  # labs(title = sprintf("rABE", refnt, mutnt), y = "Editing rate", x = "Read depth") + 
  labs( y = "Editing rate", x = "Read depth") + 
  # scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') + 
  ylim(-0.8, 0.8)+
  #scale_fill_manual(values=c(FALSE="black",TRUE="red")) + 
  theme_bw()
ggsave(sprintf("diffEditing_scatter_bg_cmh_padj005_%s_%s-%s-dox%s.png", setName, refnt, mutnt, dox), width = 6, height = 5, units = "cm")

