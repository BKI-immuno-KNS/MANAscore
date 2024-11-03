## MANAscorehi TIL define for NSCLC
## zzeng
rm(list=ls())
library(magrittr)
library(Seurat)
library(dplyr)
library(data.table)
library(scales)
library(ggbeeswarm)
library(tidyverse)
library(gplots)
library(ggpubr)
library(ggsci)
library(reshape2)
library(heatmap3)
library(cowplot)
library(matrixStats)
library(ggpointdensity)
library(stringr)
library(ggsci)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

setwd("D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/")
umap = readRDS('D:/Project/ManaScore/01.Result/00.TCR_info/00.umap_new.rds')

Cutoff = NULL
pa = umap[umap$cancer=='NSCLC' & umap$tissue=='tumor',]$patient_id %>% unique()
m = unique(umap[,c('patient_id','imid')])

for(ps in pa){
  cutoff= NULL
  imid = m[m$patient_id==ps,]$imid
  d1 = read.csv(paste0('voting_i_',ps,'_predict_score.csv'),row.names = 1)
  colnames(d1)[2] = 'score_i'
  d2 = read.csv(paste0('voting_ni_',ps,'_predict_score.csv'),row.names = 1)
  colnames(d2)[2] = 'score_ni'
  
  d = d1 %>% left_join(d2)
  
  di <- density(d$score_i)
  trough1 <- NULL
  for (i in 2:(length(di$y)-1)) {
    if (di$y[i-1] >= di$y[i] & di$y[i] <= di$y[i+1]) {
      trough1 <- cbind(trough1, c(di$x[i], di$y[i]))
    }
  }
  cutoff1 =  trough1[,dim(trough1)[2]][1]
  d$it = ifelse(d$score_i>cutoff1,'hi','lo')
  
  dni <- density(d$score_ni)
  trough2 <- NULL
  for (i in 2:(length(dni$y)-1)) {
    if (dni$y[i-1] >= dni$y[i] & dni$y[i] <= dni$y[i+1]) {
      trough2 <- cbind(trough2, c(dni$x[i], dni$y[i]))
    }
  }
  cutoff2 =  trough2[,dim(trough2)[2]][1]
  d$nt = ifelse(d$score_ni>cutoff2,'hi','lo')
  
  write.csv(d,paste0('../02.MANAscore_hi/00.NSCLC/',ps,'_predict_scores_type.csv'),quote=FALSE)
  cutoff = data.frame(patient_ID = ps,
                      cutoff_i = cutoff1,
                      cutoff_ni = cutoff2)
  Cutoff = rbind(Cutoff,cutoff)
  
  ##
  plot1 <- ggscatter(d, x = "score_i", y = "score_ni",size = 0.5,title = ps,
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "#008B45FF", fill = "lightgray"), # Customize reg. line
                     conf.int = TRUE # Add confidence interval
  )+
    geom_hline(yintercept=cutoff2,color='red',linetype="dashed",size=1.5)+
    geom_vline(xintercept=cutoff1,color='red',linetype="dashed",size=1.5)+
    scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    scale_x_continuous(expand = c(0,0),limits = c(0,1))+
    # geom_rug()+
    stat_cor(method = "pearson",label.x = 0.05, label.y = 0.95,color="#008B45FF",size=4) +
    labs(x='MANAscore_i',y='MANAscore_ni') +
    theme(axis.title = element_text(size=14,face='bold'),
          axis.text = element_text(size=12),
          plot.title = element_blank())
  
  dens1 <- ggplot(d, aes(x = score_i)) + 
    geom_density() + 
    theme_classic()+
    geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="red",size=3)+
    # geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="blue",size=2)+
    labs(title=imid)+
    scale_x_continuous(expand=c(0,0),limits = c(0,1))+
    theme(legend.position = "none",axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=16,hjust=0.5,face='bold'))
  
  
  dens2 <- ggplot(d, aes(x = score_ni)) + 
    geom_density(alpha = 0.4) + 
    theme_classic() + 
    geom_point(aes(x=trough2[,dim(trough2)[2]][1],y=trough2[,dim(trough2)[2]][2]),colour="red",size=3)+
    scale_x_continuous(expand=c(0,0),limits = c(0,1))+
    theme(legend.position = "none",axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    coord_flip()
  dens2
  
  dens1 + plot_spacer() + plot1 + dens2 + 
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  ggsave(paste0('../02.MANAscore_hi/00.NSCLC/',ps,'_correlation.tiff'),height=6,width=6)
  
  
}

###
write.csv(Cutoff,"../02.MANAscore_hi/00.NSCLC/specific_patient_cutoff.csv",quote=FALSE)

## MANAscore UMAP
ser <- readRDS("D:/Project/ManaScore/01.Result/00.Lung_melanoma_integration/ser_lung_me_integrated_fi.rds")
ser$response = ifelse(ser$patient_id=='JS10','NR',ser$response)
ser$imid = ifelse(ser$patient_id=='JS10','MD01-019',ser$imid)

Score = NULL
for(p in pa){
  dat = read.csv(paste0('../02.MANAscore_hi/00.NSCLC/',p,'_predict_scores_type.csv'),row.names=1)
  Score = rbind(Score,dat)
}

NSCLC_predict_ser = ser[,Score$barcode]
rownames(Score) = Score$barcode
NSCLC_predict_ser$score_i = Score[rownames(NSCLC_predict_ser[[]]),]$score_i
NSCLC_predict_ser$score_ni = Score[rownames(NSCLC_predict_ser[[]]),]$score_ni

p1 = FeaturePlot(NSCLC_predict_ser,features = c('score_i'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  labs(title='Imputation MANAscore') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes() & NoLegend()
p2 = FeaturePlot(NSCLC_predict_ser,features = c('score_ni'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  labs(title='Non-imputation MANAscore') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes()

p1|p2
ggsave("../02.MANAscore_hi/00.NSCLC/01.lung_tumor_manascore.tiff",height=3.5,width=8)

##
## filter by TCR 
TCR_not_shared = umap %>% filter(is.na(TRB_aa_1) == F) %>%
  group_by(TRB_aa_1,patient_id) %>%summarise(n=n()) %>%
  group_by(TRB_aa_1) %>% summarise(n=n()) %>% filter(n==1)

MANAhi <- umap %>% filter(barcode %in% Score$barcode) %>% filter(is.na(TRB_aa_1) == F ) %>% 
  left_join(Score) %>% filter(it=='hi' & nt == 'hi') %>% filter(TRB_aa_1 %in% TCR_not_shared$TRB_aa_1) ## 19753 cell

# saveRDS(MANAhi,'../02.MANAscore_hi/00.NSCLC/MANAhiTcell_umap.rds')
# MANAhi <- readRDS('../02.MANAscore_hi/00.NSCLC/MANAhiTcell_umap.rds')

NSCLC_predict_ser$it = Score[rownames(NSCLC_predict_ser@meta.data),]$it
NSCLC_predict_ser$nt = Score[rownames(NSCLC_predict_ser@meta.data),]$nt

NSCLC_predict_ser$MANAscoreType = ifelse(NSCLC_predict_ser$barcode %in% MANAhi$barcode,'MANAscorehi','MANAscorelo')
NSCLC_predict_ser$MANAscoreType = ifelse(NSCLC_predict_ser$it == 'hi' &  NSCLC_predict_ser$nt == 'hi' & 
                                           NSCLC_predict_ser$MANAscoreType == 'MANAscorelo', 'Null',NSCLC_predict_ser$MANAscoreType)

## 
sub_ser <- subset(NSCLC_predict_ser, subset=MANAscoreType!='Null')
## visualization
sub_ser$MANAscore = ifelse(sub_ser$MANAscoreType == 'MANAscorehi',1,0)
sub_ser_R = subset(sub_ser,subset = response=='R')
sub_ser_NR = subset(sub_ser,subset = response=='NR')
p1 = FeaturePlot(sub_ser_R,features = c('MANAscore'),cols = c('grey','red'),pt.size = 0.01) & 
  NoLegend() & NoAxes() &
  labs(title = 'R')&
  theme(plot.margin = margin(0,0,0,0))

p2 = FeaturePlot(sub_ser_NR,features = c('MANAscore'),cols = c('grey','red'),pt.size = 0.01) & 
  NoLegend() & NoAxes() &
  labs(title = 'NR') &
  theme(plot.margin = margin(0,0,0,0))

p1/p2
ggsave('../02.MANAscore_hi/00.NSCLC/03.MANAscorehi_T_cell_umap_by_response.jpeg',height=5.5,width=3)


## celltype contribution
st = MANAhi %>% group_by(celltype) %>% summarise(n=n()) 
st$celltype <- ifelse(st$celltype=='Effector(III)', 'NK-like',st$celltype)
st$per = st$n/19753 

a = data.frame(celltype = c("Effector(I)","Effector(II)","NK-like",
                            "TRM(I)","TRM(II)","TRM(III)","TRM(IV)","TRM(V)","TRM(VI)",
                            "NK-like TRM","Naïve/stem-like","CD4+CD8+","Proliferating","MAIT"),
               number = c(1,3,9,
                          5,2,10,13,4,8,7,11,6,14,12))
st$celltype <- factor(st$celltype,levels = a$celltype)
rownames(a) <- a$celltype
order = st%>% arrange(-per) %>% as.data.frame() %>%.[,'celltype']
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
col = getPalette(14) 

ggplot(st,aes(x=celltype,y=per,fill=celltype))+
  geom_col() +
  scale_x_discrete(limits=order)+
  scale_fill_manual(values = getPalette(14)[c(1,3,9,
                                              5,2,10,13,4,8,7,11,6,14,12)]) +
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  # scale_y_continuous(labels = scales::percent,position="right",expand=c(0,0))+
  theme_classic()+
  ylab(bquote('Proportion of MANAscore' ^ 'hi' ~ 'T cell'))+
  xlab('Cell clusters')+
  # coord_flip() +
  theme(axis.text.x=element_text(angle=30,size=10,hjust=1,vjust=1,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        axis.title.y = element_text(size=12,face='bold'),
        axis.title.x = element_text(size=14,face='bold'),
        legend.position = 'none')
ggsave("../02.MANAscore_hi/00.NSCLC/02.MANAhi_cell_proportion.tiff",height=3.5,width=5)

MANAhiTCR = MANAhi$TRB_aa_1 %>% unique() ## 2144 TCR
NSCLC_predict_ser$TRB_aa_1 = umap[NSCLC_predict_ser$barcode,]$TRB_aa_1
NSCLC_predict_ser$TCRtype = ifelse(NSCLC_predict_ser$TRB_aa_1 %in% MANAhi$TRB_aa_1,'MANAhi','MANAlo')

##
## visualization
s = NSCLC_predict_ser
s$TCR = ifelse(s$TCRtype == 'MANAhi',1,0)

s_R = subset(s,subset = response=='R')
s_NR = subset(s,subset = response=='NR')
p3 = FeaturePlot(s_R,features = c('TCR'),cols = c('grey','red'),pt.size = 0.01) & 
  NoLegend() & NoAxes() &
  labs(title = 'R')&
  theme(plot.margin = margin(0,0,0,0))

p4 = FeaturePlot(s_NR,features = c('TCR'),cols = c('grey','red'),pt.size = 0.01) & 
  NoLegend() & NoAxes() &
  labs(title = 'NR') &
  theme(plot.margin = margin(0,0,0,0))

p3/p4
ggsave('../02.MANAscore_hi/00.NSCLC/04.MANAhi_TCR_cell_umap_by_response.jpeg',height=5.5,width=3)

saveRDS(NSCLC_predict_ser,'../02.MANAscore_hi/00.NSCLC/NSCLC_predict_ser_with_MANAscoreType_MANATCRType.rds')

##
NSCLC_umap = umap %>% filter(cancer=='NSCLC') 
MANAscoreType = NSCLC_predict_ser[[]] %>% .[,c('barcode','MANAscoreType')]
NSCLC_umap <- NSCLC_umap %>% left_join(MANAscoreType) 
NSCLC_umap$TCRType = ifelse(NSCLC_umap$TRB_aa_1 %in% MANAhiTCR,'MANAhi','MANAlo')

NSCLC_umap$TCRType = ifelse(NSCLC_umap$antigen %in% c('EBV','MANA','InfluenzaA','Viral (CMV, EBV, Influenza A)'),
                            NSCLC_umap$antigen,
                            NSCLC_umap$TCRType)

saveRDS(NSCLC_umap,'../02.MANAscore_hi/00.NSCLC/NSCLC_umap_MANAscoreType_TCRType.rds')

## cutoff 5
NSCLC_ser <- NSCLC_predict_ser
new_pTRC <- NSCLC_ser[[]] %>% filter(TCRType=='MANAhi' & tissue=='tumor') %>% 
  group_by(imid,TRB_aa_1,response) %>% mutate(total=n()) %>% 
  filter(MANAscoreType=='MANAscorehi') %>%
  group_by(imid,TRB_aa_1,response,total) %>% summarise(n=n()) %>%
  filter(n>=5) ## 443

NSCLC_ser$TCRType_new <- NSCLC_ser$TCRType

NSCLC_ser$TCRType_new <- ifelse(NSCLC_ser$TRB_aa_1 %in% new_pTRC$TRB_aa_1,'pTRC',NSCLC_ser$TCRType_new)
NSCLC_ser$TCRType_new <- ifelse(!NSCLC_ser$TRB_aa_1 %in% new_pTRC$TRB_aa_1 & NSCLC_ser$TCRType == 'MANAhi','non-pTRC' ,NSCLC_ser$TCRType_new)
NSCLC_ser$TCRType_new <- ifelse(NSCLC_ser$TCRType=='MANAlo','non-pTRC',NSCLC_ser$TCRType_new)
NSCLC_ser$TCRType_new = ifelse(NSCLC_ser$TCRType_new=='non-pTRC' & NSCLC_ser$imid =='MD043-003' & NSCLC_ser$TRB_aa_1=='CAILRQGLNEQYF', 'EBV',NSCLC_ser$TCRType_new)

saveRDS(NSCLC_ser,'NSCLC_ser_new_cutoff5.rds')

d <- NSCLC_ser[[]] %>% dplyr::filter(tissue=='tumor') %>% dplyr::filter(TCRType_new %in% c('pTRC','non-pTRC')) %>%
  mutate(MANAscoreType = ifelse(MANAscoreType=='MANAscorehi','MANAscorehi','MANAscorelo')) %>%
  mutate(MANAscoreType = ifelse(is.na(MANAscoreType),'MANAscorelo',MANAscoreType))%>%
  group_by(imid,response,TRB_aa_1,TCRType_new,MANAscoreType) %>%
  summarise(n=n()) %>% 
  reshape2::dcast(imid+response+TRB_aa_1+TCRType_new~MANAscoreType, value.var = 'n') %>% 
  replace(is.na(.),0) 

d$total = d$MANAscorehi + d$MANAscorelo
write.csv(d,'NSCLC_TCRType_new.csv',quote = FALSE,row.names = FALSE)

NSCLC_ser$pTRC <- ifelse(NSCLC_ser$TCRType_new=='pTRC' &NSCLC_ser$tissue=='tumor' ,1,0)
NSCLC_ser$MANAscorehi <- ifelse(NSCLC_ser$TCRType_new == 'pTRC' & NSCLC_ser$MANAscoreType=='MANAscorehi'&NSCLC_ser$tissue=='tumor',1,0)
FeaturePlot(NSCLC_ser, features = c('MANAscorehi','pTRC'),cols = c('grey','#B2182B'),order = T) & NoLegend() & NoAxes()
ggsave('UMAP_MANAscorehipTRC_pTRC_cutoff5.tiff',height=3.5,width=7)

## celltype contribution to pTRC
st0 = NSCLC_ser[[]] %>% filter(MANAscoreType=='MANAscorehi') %>% group_by(celltype) %>% summarise(n=n()) 
st0$celltype <- ifelse(st0$celltype == 'Effector(III)','NK-like', st0$celltype)
order = st0%>% arrange(-n) %>% as.data.frame() %>%.[,'celltype']


st = NSCLC_ser[[]] %>% filter(tissue=='tumor' & TCRType_new=='pTRC') %>% group_by(celltype) %>% summarise(n=n()) 
st$per = st$n/45190 
st$celltype <- ifelse(st$celltype == 'Effector(III)','NK-like', st$celltype)

st$celltype <- factor(st$celltype,
                      levels=c("Effector(I)","Effector(II)","NK-like",
                               "TRM(I)","TRM(II)","TRM(III)","TRM(IV)","TRM(V)","TRM(VI)",
                               "NK-like TRM","Naïve/stem-like","CD4+CD8+","Proliferating","MAIT"))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
ggplot(st,aes(x=celltype,y=per,fill=celltype))+
  geom_col()+
  scale_x_discrete(limits=order)+
  scale_fill_manual(values = getPalette(14)[c(1,3,9,
                                              5,2,10,13,4,8,7,11,6,14,12)]) +
  scale_y_continuous(labels = scales::percent,expand=c(0,0))+
  # scale_y_continuous(labels = scales::percent,position="right",expand=c(0,0))+
  theme_classic()+
  ylab(bquote(bold('Proportion of MANA' ^ 'hi' ~ 'clone T cells')))+
  xlab('Cell clusters')+
  # coord_flip() +
  theme(axis.text.x=element_text(angle=30,size=10,hjust=1,vjust=1,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        axis.title.y = element_text(size=12,face='bold'),
        axis.title.x = element_text(size=14,face='bold'),
        legend.position = 'none')
ggsave("pTRC_cell_proportion_cutoff5.tiff",height=3.5,width=5)

##
TIL = d$type %>% table() %>% as.data.frame()
TIL$TIL = 'TIL'
colnames(TIL)[1] <- 'type'

TIL$Type = factor(TIL$type, 
                  levels = c('MANAscore_hi pTRC','MANAscore_lo pTRC','MANAscore_hi non-pTRC','non-pTRC','unknown'))

TIL$Type[TIL$type=='Null'] <- 'unknown'
ggplot(TIL,aes(x=TIL,y=Freq,fill=Type))+
  geom_col(position='fill') +
  theme_classic()+
  scale_y_continuous(expand= c(0,0),labels = scales::percent)+
  scale_fill_manual(values = c("#B2182B",'skyblue','red',"#2166AC",'grey'))+
  labs(fill='',y='',x='')+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.text.y = element_text(size=11,color='black'))
ggsave('TIL_barplot_cutoff5.tiff',height=4,width=3.5)

## all TIL by patient
d$type <- ifelse(d$type=='Null', 'unknown',d$type)

patient_TIL = d %>% group_by(response,imid,type) %>%
  summarise(n=n())
p <- patient_TIL[,c('response','imid')] %>% unique() %>% as.data.frame () %>% .[,'imid']
patient_TIL$imid <- factor(patient_TIL$imid, levels = p)
patient_TIL$type <- factor(patient_TIL$type, 
                           levels = c('MANAscore_hi pTRC','MANAscore_lo pTRC','MANAscore_hi non-pTRC','non-pTRC','unknown'))

patient_TIL %>% na.omit() %>%
  ggplot(aes(x=imid,y=n,fill=type)) +
  geom_col(position = 'fill',width=0.7)+
  theme_classic()+
  scale_y_continuous(expand= c(0,0),labels = scales::percent)+
  scale_fill_manual(values = c("#B2182B",'skyblue','red',"#2166AC",'grey'))+
  labs(fill='',y='',x='')+
  theme(axis.text.x = element_text(size=12,face='bold',color='black', angle=30,hjust=1,vjust=1),
        axis.text.y = element_text(size=11,color='black'))
NoLegend()
ggsave('TIL_barplot_by_patient_cutoff5.tiff',height=4.3,width=10)

######### 
meta <- readRDS('D:/Project/ManaScore/01.Result/36.12_single_moldes/02.MANAscore_hi/00.NSCLC/NSCLC_umap_MANAscoreType_TCRType.rds')
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
col = getPalette(14)
mycol = data.frame(cluster=factor(c(0:13)),
                   color = col)
a <- meta[,c('integrated_snn_res.0.5','celltype')] %>% unique()
colnames(a)[1] <- 'cluster'
mycol = a %>% left_join(mycol)
rownames(mycol) <- mycol$celltype

d <- meta %>% filter(tissue=='tumor') %>% filter(TCRType %in% c('MANAhi','MANAlo')) %>% 
  group_by(celltype,MANAscoreType) %>%
  summarise(n=n()) %>% 
  na.omit() %>%
  reshape2::dcast(celltype~MANAscoreType,value.var = 'n') %>%
  mutate(MANAhiper=MANAscorehi/(MANAscorehi+MANAscorelo+Null)) %>%
  arrange(-MANAhiper) %>%
  mutate(celltype = factor(celltype,level=celltype))
mycol = mycol[levels(d$celltype),]

ggplot(d,aes(x=celltype,y=MANAhiper,fill=celltype))+
  geom_col() +
  scale_fill_manual(values=mycol$color)+
  xlab('Cell clusters')+
  theme_classic()+
  ylab(bquote('Proportion of MANAscore' ^ 'hi' ~ ' in each cluster'))+
  scale_y_continuous(expand = c(0,0),labels = scales::percent) +
  theme(axis.text.x=element_text(angle=30,size=10,hjust=1,vjust=1,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        axis.title.y = element_text(size=12,face='bold'),
        axis.title.x = element_text(size=14,face='bold'),
        legend.position = 'none')
ggsave("../02.MANAscore_hi/00.NSCLC/02.MANAhi_cell_proportion_in_each_cluster.tiff",height=4.5,width=7)
