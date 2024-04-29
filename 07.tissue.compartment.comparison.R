## tissue compartment comparison
## zzeng
rm(list=ls())
library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)

setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/39.tumor_vs_normal_cutoff5/')
NSCLC_ser <- readRDS('../14.specific_genes_1/NSCLC_ser_add_scores.rds')

## tumor vs normal
##
rm(list=ls())
library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(ashr)
library(tidyverse)
library(pheatmap)
library(scales)
library(ggridges)

setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/39.tumor_vs_normal_cutoff5/')
NSCLC_ser <- readRDS('../14.specific_genes_1/NSCLC_ser_add_scores.rds')

shared_clone <- NSCLC_ser[[]] %>% filter(tissue=='normal' & TCRType_new =='pTRC') %>%  ## 3620 cells in normal
  group_by(imid,TRB_aa_1) %>%
  summarise(normal=n()) ## 156 MANAhi clone shared with tumor

shared_clone_in_tumor_normal <- NSCLC_ser[[]] %>% filter(tissue=='tumor' & TCRType_new=='pTRC') %>%  ## 30436 cells in tumor
  group_by(imid,TRB_aa_1) %>%
  summarise(tumor=n()) %>%
  right_join(shared_clone)

write.csv(shared_clone_in_tumor_normal,'Shared_clone_in_tumor_and_normal_cutoff5.csv',quote=FALSE)

shared_NSCLC_ser.t.n <- subset(NSCLC_ser,
                               subset = TRB_aa_1 %in% shared_clone_in_tumor_normal$TRB_aa_1 & tissue %in% c('tumor','normal'))
## 34,056 cells
FetchData(shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response')) %>%
  group_by(imid,TRB_aa_1,tissue,response) %>%
  summarise(mean.CD39=mean(ENTPD1)) %>%
  ggplot(aes(x=tissue,y=mean.CD39,fill=tissue))+
  geom_boxplot(outlier.size = 0.2)+
  scale_fill_manual(values = c('#008B45FF','#EE0000FF'))+
  stat_compare_means(comparisons = list(c('tumor','normal')))+
  theme_classic() +
  # scale_y_continuous(expand = c(0,0))+
  labs(x='',y='Mean CD39 expression within tissue clone',fill='')+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold',color='black')) &
  NoLegend()
ggsave('CD39_exp_in_shared_tumor_vs_normal_cutoff5.tiff',heigh=4,width=3.2)

##
FetchData(shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response')) %>%
  group_by(imid,TRB_aa_1,tissue,response) %>%
  summarise(mean.CD39=mean(ENTPD1)) %>%
  ggplot(aes(x=tissue,y=mean.CD39,fill=tissue))+
  geom_boxplot(outlier.size = 0.2)+
  facet_wrap(.~response)+
  scale_fill_manual(values = c('#008B45FF','#EE0000FF')) +
  stat_compare_means(comparisons = list(c('tumor','normal')))+
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='Mean CD39 expression within tissue clone',fill='')+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold',color='black'),
        strip.text = element_text(size=13,face='bold',color='black')) &
  NoLegend()
ggsave('CD39_exp_in_shared_tumor_vs_normal_by_response_cutoff5.tiff',heigh=4,width=4)


FetchData(shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response')) %>%
  group_by(imid,TRB_aa_1,tissue,response) %>%
  summarise(mean.CD39=mean(ENTPD1)) %>%
  ggplot(aes(x=response,y=mean.CD39,fill=response))+
  geom_boxplot(outlier.size = 0.2)+
  facet_wrap(.~tissue)+
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('R','NR')))+
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='Average CD39 expression within tissue clone',fill='')+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold',color='black'),
        strip.text = element_text(size=13,face='bold',color='black')) &
  NoLegend()
ggsave('CD39_exp_in_shared_R_vs_NR_in_tumor_and_normal_cutoff5.tiff',height=4,width=4)

##
## signature score
NSCLC_ser_score <- readRDS('D:/Project/ManaScore/01.Result/36.12_single_moldes/14.specific_genes_1/NSCLC_meta_add_scores.rds')
NSCLC_ser_score %>% filter(barcode %in% shared_NSCLC_ser.t.n$barcode) %>%
  group_by(imid,TRB_aa_1,response,tissue) %>%
  summarise(mean=mean(checkpoint1)) %>%
  ggplot(aes(x=tissue,y=mean,fill=tissue)) +
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c('#008B45FF','#EE0000FF')) +
  stat_compare_means(comparisons = list(c('tumor','normal')))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='Mean checkpoint score of clones')+
  NoLegend()+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_checkpoint_score_in_shared_pTRC_in_tumor_vs_normal_cutoff5.tiff',height=4,width=3)

NSCLC_ser_score %>% filter(barcode %in% shared_NSCLC_ser.t.n$barcode) %>%
  group_by(imid,TRB_aa_1,response,tissue) %>%
  summarise(mean=mean(TRM_sig1)) %>%
  ggplot(aes(x=tissue,y=mean,fill=tissue)) +
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c('#008B45FF','#EE0000FF')) +
  stat_compare_means(comparisons = list(c('tumor','normal')))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='Mean TRM signature score of clones')+
  NoLegend()+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_TRM_sig_score_in_shared_pTRC_in_tumor_vs_normal_cutoff5.tiff',height=4,width=3)

#
NSCLC_ser_score %>% filter(barcode %in% shared_NSCLC_ser.t.n$barcode) %>%
  group_by(imid,TRB_aa_1,response,tissue) %>%
  summarise(mean=mean(cytotoxicity1)) %>%
  ggplot(aes(x=tissue,y=mean,fill=tissue)) +
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c('#008B45FF','#EE0000FF')) +
  stat_compare_means(comparisons = list(c('tumor','normal')))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='Mean cytotoxicity score of clones')+
  NoLegend()+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_cytotoxicity_score_in_shared_pTRC_in_tumor_vs_normal_cutoff5.tiff',height=4,width=3)

##
NSCLC_ser_score %>% filter(barcode %in% shared_NSCLC_ser.t.n$barcode) %>%
  group_by(imid,TRB_aa_1,response,tissue) %>%
  summarise(mean=mean(TCR_signaling1)) %>%
  ggplot(aes(x=tissue,y=mean,fill=tissue)) +
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c('#008B45FF','#EE0000FF')) +
  stat_compare_means(comparisons = list(c('tumor','normal')))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='Mean TCR signaling score of clones')+
  NoLegend()+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_TCR_signaling_score_in_shared_pTRC_in_tumor_vs_normal_cutoff5.tiff',height=4,width=3)

##
NSCLC_ser_score %>% filter(barcode %in% shared_NSCLC_ser.t.n$barcode) %>%
  group_by(imid,TRB_aa_1,response,tissue) %>%
  summarise(mean=mean(stemness1)) %>%
  ggplot(aes(x=tissue,y=mean,fill=tissue)) +
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c('#008B45FF','#EE0000FF')) +
  stat_compare_means(comparisons = list(c('tumor','normal')))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='Mean stemness score of clones')+
  NoLegend()+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_stemness_score_in_shared_pTRC_in_tumor_vs_normal_cutoff5.tiff',height=4,width=3)


## % CD39-
FetchData(shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response')) %>%
  group_by(imid,TRB_aa_1,tissue,response) %>%
  mutate(tissue.clone = n()) %>%
  filter(ENTPD1 == 0) %>%
  group_by(imid,TRB_aa_1,tissue,response,tissue.clone) %>%
  summarise(CD39neg = n()) %>%
  mutate(per = CD39neg/tissue.clone) %>%
  ggplot(aes(x=response,y=per,fill=response))+
  geom_boxplot(outlier.size = 0.2)+
  facet_wrap(.~tissue)+
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('R','NR')))+
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)),labels = scales::percent)+
  labs(x='',y='% CD39neg within tissue clone',fill='')+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold',color='black'),
        strip.text = element_text(size=13,face='bold',color='black')) &
  NoLegend()
ggsave('CD39neg_proportion_within_clone_cutoff5.tiff',height=4,width=4)  

## 
patient_tisue_number <- NSCLC_ser[[]] %>% group_by(imid,tissue) %>% summarise(tissue.n = n())
##
FetchData(shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response','TCRType')) %>%
  group_by(imid,TRB_aa_1,tissue,response) %>%
  summarise(tissue.clone = n()) %>%
  reshape2::dcast(imid+TRB_aa_1+response ~ tissue, value.var = 'tissue.clone') %>%
  mutate(ratio=tumor/normal) %>%
  # filter(tumor+normal>10)%>%
  ggplot(aes(x=response,y=log(ratio),fill=response))+
  geom_violin()+
  geom_boxplot(outlier.size = 0.2,width=0.4) +
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('R','NR'))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='tumor:normal (log)')+
  theme_classic()+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold',color='black'))+
  NoLegend()
ggsave('ratio_of_cell_in_normal_tumor_each_clone.tiff',height=4,width=3.2)  

## ratio for MANAlo
shared_clone_MANAlo <- NSCLC_ser[[]] %>% filter(tissue=='normal' & TCRType_new=='non-pTRC') %>%  # 19676 cells
  group_by(imid,TRB_aa_1) %>%
  summarise(normal=n()) ## 22845 MANAlo clone shared with tumor

shared_clone_MANAlo_in_tumor_normal <- NSCLC_ser[[]] %>% filter(tissue=='tumor' & TCRType_new=='non-pTRC') %>%  
  group_by(imid,TRB_aa_1) %>%
  summarise(tumor=n()) %>%
  right_join(shared_clone_MANAlo)

write.csv(shared_clone_MANAlo_in_tumor_normal,'Shared_clone_non_pTRC_in_tumor_and_normal_cutoff5.csv',quote=FALSE)

cell <- NSCLC_ser[[]] %>% right_join(shared_clone_MANAlo_in_tumor_normal) %>% filter(tissue %in% c('tumor','normal'))
MANAlo_shared_NSCLC_ser.t.n <- NSCLC_ser[,cell$barcode]

FetchData(MANAlo_shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response')) %>%
  group_by(imid,TRB_aa_1,tissue,response) %>%
  summarise(tissue.clone = n()) %>%
  reshape2::dcast(imid+TRB_aa_1+response ~ tissue, value.var = 'tissue.clone') %>%
  mutate(ratio=tumor/normal) %>%
  ggplot(aes(x=response,y=log(ratio),fill=response))+
  geom_violin()+
  geom_boxplot(outlier.size = 0.2,width=0.3) +
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('R','NR'))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='tumor:normal (log)')+
  theme_classic()+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold',color='black'))+
  NoLegend()
ggsave('MANAlo_ratio_of_cell_in_normal_tumor_each_clone_cutoff5.tiff',height=4,width=3.2)  

##
FetchData(MANAlo_shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response')) %>%
  group_by(imid,TRB_aa_1,tissue,response) %>%
  summarise(tissue.clone = n()) %>%
  reshape2::dcast(imid+TRB_aa_1+response ~ tissue, value.var = 'tissue.clone') %>%
  mutate(ratio=tumor/normal) %>%
  # filter(tumor+normal>10) %>%
  ggplot(aes(x=response,y=log(ratio),fill=response))+
  geom_violin()+
  geom_boxplot(outlier.size = 0.2,width=0.3) +
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('R','NR'))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='tumor:normal (log)')+
  theme_classic()+
  theme(axis.text.x = element_text(size=12,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold',color='black'))+
  NoLegend()

patient_tissue_total <- NSCLC_ser[[]] %>% filter(tissue %in% c('tumor','normal')) %>%
  group_by(imid,tissue) %>%
  summarise(tissue.n = n())

df_hi <- FetchData(shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response','TCRType_new')) %>%
  group_by(imid,TRB_aa_1,tissue,response,TCRType_new) %>%
  summarise(tissue.clone = n()) %>% left_join(patient_tissue_total) %>%
  mutate(per = tissue.clone/tissue.n) %>%
  reshape2::dcast(imid+TRB_aa_1+response + TCRType_new ~ tissue, value.var = 'per') %>%
  mutate(ratio=tumor/normal)

df_lo <- FetchData(MANAlo_shared_NSCLC_ser.t.n,vars = c('tissue','ENTPD1','TRB_aa_1','imid','response','TCRType_new')) %>%
  group_by(imid,TRB_aa_1,tissue,response,TCRType_new) %>%
  summarise(tissue.clone = n()) %>% left_join(patient_tissue_total) %>%
  mutate(per = tissue.clone/tissue.n) %>%
  reshape2::dcast(imid+TRB_aa_1+response + TCRType_new~ tissue, value.var = 'per') %>%
  mutate(ratio=tumor/normal)

rbind(df_hi, df_lo) %>%
  filter(TCRType_new %in% c('non-pTRC','pTRC')) %>%
  ggplot(aes(x=TCRType_new,y=log(ratio),fill = TCRType_new))+
  geom_violin()+
  geom_boxplot(outlier.size = 0.2,width=0.2) +
  scale_fill_manual(values = c("#2166AC","#B2182B")) +
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC'))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  labs(x='',y='tumor:normal (log)')+
  theme_classic()+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.title.y = element_text(size=12,face='bold',color='black'))+
  NoLegend()
ggsave('ratio_tumor_normal_pTRC_vs_non-pTRC_clones_cutoff5.tiff',height=3.5,width=2.5)


## shared clone in JM6 and PP3
JM6 <- subset(NSCLC_ser, subset = tissue != 'mettumor' & imid == 'MD043-011')
JM6_clone <- JM6[[]] %>% group_by(tissue) %>% mutate(tissue.total=n()) %>%
  group_by(tissue,tissue.total,TRB_aa_1,TCRType_new)  %>%
  summarise(n=n()) %>%
  filter(TCRType_new %in% c('pTRC','non-pTRC')) %>%
  reshape2::dcast(TRB_aa_1+TCRType_new~tissue, value.var = 'n')
JM6_clone[is.na(JM6_clone)] <- 0
JM6_shared_clone <- JM6_clone %>% filter(LN>0 & normal>0 & tumor>0)
JM6_shared_clone$TCRType_new %>% table() ## 38 pTRC 222 non-pTRC

JM6[[]] %>% group_by(tissue) %>% mutate(tissue.total=n()) %>%
  filter(TRB_aa_1 %in% JM6_shared_clone$TRB_aa_1) %>%
  group_by(TCRType_new, tissue, tissue.total, TRB_aa_1) %>%
  summarise(n=n()) %>%
  mutate(per=n/tissue.total) %>%
  mutate(TCRType_new= factor(TCRType_new, levels=c('pTRC','non-pTRC'))) %>%
  mutate(tissue = factor(tissue, levels = c('normal','LN','tumor'))) %>%
  ggplot(aes(x=tissue, y=log(per), fill = tissue))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(size=0.8)+
  facet_wrap(.~TCRType_new, scales='free_y')+
  scale_fill_manual(values = c('#008B45FF','#3B4992FF','#EE0000FF'))+
  theme_classic() + NoLegend()+
  labs(x='',y='Frequency of clone in tissue (log)') +
  theme(axis.text.x = element_text(size=13,color='black',face='bold'),
        axis.text.y = element_text(size=13,color='black'),
        axis.title.y = element_text(size=14,color='black',face='bold'),
        strip.background = element_blank(),
        strip.text = element_text(size=13,face='bold',color='black')) +
  stat_compare_means(comparisons = list(c('normal','LN'),c('tumor','LN'),c('tumor','normal')),paired=TRUE)

ggsave('JM6_shared_clone_freq.tiff',height=4,width=6)

## PP3
PP3 <- subset(NSCLC_ser, subset = tissue != 'mettumor' & imid == 'MD01-005')
PP3_clone <- PP3[[]] %>% group_by(tissue) %>% mutate(tissue.total=n()) %>%
  group_by(tissue,tissue.total,TRB_aa_1,TCRType_new)  %>%
  summarise(n=n()) %>%
  filter(TCRType_new %in% c('pTRC','non-pTRC')) %>%
  reshape2::dcast(TRB_aa_1+TCRType_new~tissue, value.var = 'n')
PP3_clone[is.na(PP3_clone)] <- 0
PP3_shared_clone <- PP3_clone %>% filter(LN>0 & normal>0 & tumor>0)
PP3_shared_clone$TCRType_new %>% table() ## 12 pTRC 260 non-pTRC

PP3[[]] %>% group_by(tissue) %>% mutate(tissue.total=n()) %>%
  filter(TRB_aa_1 %in% PP3_shared_clone$TRB_aa_1) %>%
  group_by(TCRType_new, tissue, tissue.total, TRB_aa_1) %>%
  summarise(n=n()) %>%
  mutate(per=n/tissue.total) %>%
  mutate(TCRType_new= factor(TCRType_new, levels=c('pTRC','non-pTRC'))) %>%
  mutate(tissue = factor(tissue, levels = c('normal','LN','tumor'))) %>%
  ggplot(aes(x=tissue, y=log(per), fill = tissue))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(size=0.8)+
  facet_wrap(.~TCRType_new, scales='free_y')+
  scale_fill_manual(values = c('#008B45FF','#3B4992FF','#EE0000FF'))+
  theme_classic() + NoLegend()+
  labs(x='',y='Frequency of clone in tissue (log)') +
  theme(axis.text.x = element_text(size=13,color='black',face='bold'),
        axis.text.y = element_text(size=13,color='black'),
        axis.title.y = element_text(size=14,color='black',face='bold'),
        strip.background = element_blank(),
        strip.text = element_text(size=13,face='bold',color='black')) +
  stat_compare_means(comparisons = list(c('normal','LN'),c('tumor','LN'),c('tumor','normal')),paired=TRUE)

ggsave('PP3_shared_clone_freq.tiff',height=4,width=6)