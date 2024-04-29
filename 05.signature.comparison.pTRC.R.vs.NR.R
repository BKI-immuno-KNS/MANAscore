## pTRC R vs NR
## zzeng
rm(list=ls())
library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)

getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
col = getPalette(15) 
setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/43.cutoff5/')
NSCLC_ser <- readRDS('../14.specific_genes_1/NSCLC_ser_add_scores.rds')
# No.1  IL7R 
FetchData(NSCLC_ser,vars = c('rna_IL7R','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(rna_IL7R)) %>%
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.1,outlier.size = 0.3) +
  # geom_point()+
  stat_compare_means(comparisons = list(c('NR','R')))+
  scale_fill_aaas()+
  theme_classic()+
  NoLegend()+
  labs(x='',y='Mean IL7R expression of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))
ggsave('mean_IL7R_expression_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)

##
## No.2 checkpoint

FetchData(NSCLC_ser,vars = c('checkpoint1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(checkpoint1)) %>%
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.2,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('NR','R')))+
  labs(x='',y='Mean checkpoint score of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))+
  NoLegend()
ggsave('mean_checkpoint_score_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)


##
## No. 3 exhaustion 

FetchData(NSCLC_ser,vars = c('Exhaustion1','TCRType_new','response','TRB_aa_1','imid', 'tissue')) %>%
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(Exhaustion1)) %>%
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('NR','R')))+
  labs(x='',y='Mean exhaustion score of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))+
  NoLegend()
ggsave('mean_exhaustion_score_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)

##
## No. 4 cytotoxicity 

FetchData(NSCLC_ser,vars = c('cytotoxicity1','TCRType_new','response','TRB_aa_1','imid', 'tissue')) %>%
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(cytotoxicity1)) %>%
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('NR','R')))+
  labs(x='',y='Mean cytotoxicity score of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))+
  NoLegend()
ggsave('mean_cytotoxicity_score_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)


##  No. 5 TCR signaling
FetchData(NSCLC_ser,vars = c('TCR_signaling1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(TCR_signaling1)) %>%
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  scale_fill_aaas()+
  stat_compare_means(comparisons = list(c('NR','R')))+
  labs(x='',y='Mean TCR signaling score of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))+
  NoLegend()
ggsave('mean_TCR_signaling_score_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)


## No. 6 TRM sig
FetchData(NSCLC_ser,vars = c('TRM_sig1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>% 
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(TRM_sig1)) %>% 
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  theme_classic()+
  scale_fill_aaas() +
  stat_compare_means(comparisons = list(c('NR','R')))+
  labs(x='',y='Mean TRM sig score of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))+
  NoLegend()
ggsave('mean_TRM_sig_score_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)

## No. 8 IFNG 
FetchData(NSCLC_ser,vars = c('IFNG','TCRType_new','response','TRB_aa_1','imid','tissue')) %>% 
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(IFNG)) %>% 
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point() +
  theme_classic()+
  scale_fill_aaas() +
  stat_compare_means(comparisons = list(c('NR','R')))+
  labs(x='',y='Mean IFNG expression of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))+
  NoLegend()
ggsave('mean_IFNG_expression_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)

## No. 9 IFNG sig
FetchData(NSCLC_ser,vars = c('IFNG_sig1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>% 
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(IFNG_sig1)) %>% 
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point() +
  theme_classic()+
  scale_fill_aaas() +
  stat_compare_means(comparisons = list(c('NR','R')))+
  labs(x='',y='Mean IFNG signature score of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))+
  NoLegend()
ggsave('mean_IFNG_signature_score_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)

##
## No. 10 ITGB1 
FetchData(NSCLC_ser,vars = c('ITGB1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>% 
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(ITGB1)) %>% 
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point() +
  theme_classic()+
  scale_fill_aaas() +
  stat_compare_means(comparisons = list(c('NR','R')))+
  labs(x='',y='Mean ITGB expression of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))+
  NoLegend()
ggsave('mean_ITGB1_expression_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)

# No.11  stemness 
FetchData(NSCLC_ser,vars = c('stemness1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(TCRType_new=='pTRC' & tissue=='tumor') %>%
  group_by(imid,TRB_aa_1,response) %>%
  summarise(mean=mean(stemness1)) %>%
  ggplot(aes(x=response,y=mean,fill=response))+
  geom_violin()+
  geom_boxplot(width=0.1,outlier.size = 0.3) +
  # geom_point()+
  stat_compare_means(comparisons = list(c('NR','R')))+
  scale_fill_aaas()+
  theme_classic()+
  NoLegend()+
  labs(x='',y='Mean stemness score of pTRCs')+
  theme(axis.text.x = element_text(size=14,face='bold',color='black'),
        axis.title.y = element_text(size=12,face='bold'))
ggsave('mean_stemness_expression_in_MANAhi_clones_cutoff_5.tiff',height=3.5,width=2.5)

## pTRC in R and NR
NSCLC_ser[[]] %>% filter(tissue=='tumor') %>% group_by(imid,tissue) %>% mutate(total = n()) %>% filter(TCRType_new =='pTRC') %>%
  group_by(patient_id,imid,response,resi_tumor,total,tissue) %>%
  summarise(n = n()) %>%
  mutate(prop = n/total) %>% as.data.frame() %>%
  ggplot(aes(x=response,y=prop)) +
  geom_boxplot(aes(fill=response),width=0.7,outlier.colour=NA)+
  geom_point(aes(color=imid),size=3,position='jitter') +
  scale_color_manual(values = col) +
  scale_fill_aaas(guide="none")+
  labs(color='',y='Proportion of pTRC in tumor',x='') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  stat_compare_means(comparisons = list(c("R","NR")), method = 't.test') +
  # stat_pvalue_manual(prop.stat.test,size=5,hide.ns = F,label = 'p.adj') +
  theme(axis.title.y = element_text(size=14,color='black',face='bold'),
        axis.text.x = element_text(face='bold',size=14,color='black'),
        axis.text.y = element_text(size=12,color='black'),
        legend.key.size = unit(0.5,'cm'))+
  NoLegend()
ggsave('pTRC_proportion_R_NR_cutoff5.tiff',height=4.1,width=3.5)
##
NSCLC_ser[[]] %>% filter(tissue=='tumor') %>% group_by(imid,tissue) %>% mutate(total = n()) %>% filter(TCRType_new =='pTRC') %>%
  group_by(patient_id,imid,response,resi_tumor,total,tissue) %>%
  summarise(n = n()) %>%
  mutate(prop = n/total) %>% as.data.frame() %>% 
  ggscatter(x = "resi_tumor", y = "prop",size = 3,color="imid",
            add = "reg.line",  # Add regressin line
            add.params = list(color = "black"), # Customize reg. line
            conf.int = FALSE # Add confidence interval
  )+
  stat_cor(method = "spearman",label.x = 0.4, label.y = 0.75,color="red")+
  scale_color_manual(values = col)+
  scale_y_continuous(labels=scales::percent,limits = c(0,0.8))+
  scale_x_continuous(labels=scales::percent)+
  labs(x='resi_tumor',y='Proportion of pTRC in tumor',color='')+
  theme(legend.position = "right",legend.key.size = unit(0.5, 'cm'),
        axis.title = element_text(size=14,color='black',face='bold'),
        axis.text = element_text(size=13,color='black')) +
  NoLegend()

ggsave('pTRC_proportion_vs_resi_tumor_cutoff5.tiff',height=4.1,width=4)
