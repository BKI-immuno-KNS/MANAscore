## signature compare
## pTRC vs non-pTRC
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

NSCLC_ser = readRDS("D:/Project/ManaScore/01.Result/36.12_single_moldes/02.MANAscore_hi/00.NSCLC/NSCLC_ser.rds")
score_ser = readRDS('D:/Project/ManaScore/01.Result/36.12_single_moldes/02.MANAscore_hi/00.NSCLC/NSCLC_predict_ser_with_MANAscoreType_MANATCRType.rds')
scores = score_ser[[]] %>% select(barcode,score_i,score_ni)
scores$score = sqrt(scores$score_i^2+scores$score_ni^2)

checkpoint = c('CTLA4', 'PDCD1', 'LAG3', 'HAVCR2', 'TIGIT','ENTPD1')
Exhaustion = c('ENTPD1','LAYN','ITGAE','BATF')
sigs <- readRDS('D:/Project/ManaScore/01.Result/35.downstream_analysis/06.signature/sigs.rds')
TCR_signaling = sigs[['TCR Signaling']] %>% intersect(rownames(NSCLC_ser))
IFNG_sig = c('IDO1','CXCL10','CXCL9','HLA-DRA','STAT1','IFNG') %>% intersect(rownames(NSCLC_ser))
TRM_sig = c('CA10', 'ITGA1','ITGAE','IL2','IL10','CXCR6','CXCL13','KCNK5','RGS1','CRTAM','DUSP6','PDCD1','IL23R','ZNF683') %>%
  intersect(rownames(NSCLC_ser))
cytotoxicity = c('PRF1','GZMB','GZMA','GZMH','NKG7','GNLY')  %>% intersect(rownames(NSCLC_ser))
NK_cell_receptor = c('KLRD1','FGFBP2','FCGR3A','S1PR5','KLRC1','KLRC3','KLRB1','KLRC2') %>% intersect(rownames(NSCLC_ser))


NSCLC_ser = AddModuleScore(NSCLC_ser,features = list(checkpoint),name='checkpoint')
NSCLC_ser = AddModuleScore(NSCLC_ser,features = list(Exhaustion), name ='Exhaustion' )
NSCLC_ser = AddModuleScore(NSCLC_ser,features = list(TCR_signaling), name ='TCR_signaling')
NSCLC_ser = AddModuleScore(NSCLC_ser,features = list(IFNG_sig), name ='IFNG_sig')
NSCLC_ser = AddModuleScore(NSCLC_ser,features = list(TRM_sig), name ='TRM_sig')
NSCLC_ser = AddModuleScore(NSCLC_ser,features = list(cytotoxicity), name ='cytotoxicity')
NSCLC_ser = AddModuleScore(NSCLC_ser,features = list(NK_cell_receptor), name ='NK_cell_receptor')


## No.1 checkpoint
##
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(checkpoint1)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.2,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')))+
  labs(x='',y='Mean checkpoint score of clones')+
  NoLegend()+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_checkpoint_score_in_MANAhi_vs_MANAlo_clones_all_cutoff_5.tiff',height=4.5,width=3)

## No. 2 Exhaustion
##
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(Exhaustion1)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  labs(x='',y='Mean exhaustion score of clones')+
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')),label = 'p.format')+
  NoLegend()+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_exhaustion_score_in_MANAhi_vs_MANAlo_clones_all_cutoff_5.tiff',height=4.5,width=3)

## No. 3 TCR signaling
##
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(TCR_signaling1)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  labs(x='',y='Mean TCR signaling score of clones')+
  NoLegend() +
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')),label = 'p.format')+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))

ggsave('mean_TCR_signaling_score_in_MANAhi_vs_MANAlo_clones_all_cutoff_5.tiff',height=4.5,width=3)

## No. 4 IFNG
##
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','TCRType_new','response','TRB_aa_1','imid','tissue','IFNG')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(IFNG)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.1,outlier.size = 0.3,alpha=0.5) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  labs(x='',y='Mean IFNG expression of clones')+
  NoLegend() +
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')),label = 'p.format')+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))

ggsave('mean_IFNG_in_MANAhi_vs_MANAlo_clones_all_cutoff_5.tiff',height=4.5,width=3)

## No. 5 IFNG_sig
##
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','TCRType_new','response','TRB_aa_1','imid','tissue','IFNG_sig1')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(IFNG_sig1)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  labs(x='',y='Mean IFNG signature expression of clones')+
  NoLegend() +
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')),label = 'p.format')+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))

ggsave('mean_IFNG_sig_in_MANAhi_vs_MANAlo_clones_all_cutff_5.tiff',height=4.5,width=3)


## No. 6 ITGB1
##
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','TCRType_new','response','TRB_aa_1','imid','tissue','ITGB1')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(ITGB1)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.1,outlier.size = 0.3,alpha=0.5) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  labs(x='',y='Mean ITGB1 expression of clones')+
  NoLegend() +
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')),label = 'p.format')+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))

ggsave('mean_ITGB1_in_MANAhi_vs_MANAlo_clones_all_cutoff_5.tiff',height=4.5,width=3)

## No. 7 Cytotoxicity
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','NK_cell_receptor1',
                             'cytotoxicity1','TRM_sig1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(cytotoxicity1)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  labs(x='',y='Mean cytotoxicity score of clones')+
  NoLegend() +
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')),label = 'p.format')+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_cytotoxicity_score_in_MANAhi_vs_MANAlo_clones_all_cutoff_5.tiff',height=4.5,width=3)

## No. 8 TRM
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','NK_cell_receptor1',
                             'cytotoxicity1','TRM_sig1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(TRM_sig1)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  labs(x='',y='Mean TRM signature score of clones')+
  NoLegend() +
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')),label = 'p.format')+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_TRM_signature_score_in_MANAhi_vs_MANAlo_clones_all_cutoff_5.tiff',height=4.5,width=3)


##
## No.9 stemness
FetchData(NSCLC_ser,vars = c('Exhaustion1','checkpoint1','TCR_signaling1','NK_cell_receptor1','stemness1',
                             'cytotoxicity1','TRM_sig1','TCRType_new','response','TRB_aa_1','imid','tissue')) %>%
  filter(tissue=='tumor' & TCRType_new %in% c('pTRC','non-pTRC')) %>%
  # mutate(TCRType = factor(TCRType,levels = c('MANAlo','MANAhi')))%>%
  # mutate(TCRType = ifelse(TCRType=='MANAhi','pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,response,TCRType_new) %>%
  summarise(mean=mean(stemness1)) %>%
  ggplot(aes(x=TCRType_new,y=mean,fill=TCRType_new))+
  geom_violin()+
  geom_boxplot(width=0.1,outlier.size = 0.3) +
  # geom_point()+
  theme_classic()+
  # facet_wrap(.~response)+
  scale_fill_manual(values = c("#2166AC","#B2182B"))+
  labs(x='',y='Mean stemness score of clones')+
  NoLegend() +
  stat_compare_means(comparisons = list(c('pTRC','non-pTRC')),label = 'p.format')+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('mean_stemness_score_in_MANAhi_vs_MANAlo_clones_all_cutoff_5.tiff',height=4.5,width=3)

## UMAP for MANAscorelo pTRC 
NSCLC_ser$MANAscorelo_pTRC <- ifelse(NSCLC_ser$TCRType_new=='pTRC' & NSCLC_ser$MANAscoreType == 'MANAscorelo',1,0)
FeaturePlot(NSCLC_ser, features = c('MANAscorelo_pTRC'),cols = c('grey','#B2182B'),order = T) & NoLegend() & NoAxes()
ggsave('UMAP_MANAscorelopTRC_pTRC_cutoff5.tiff',height=3.5,width=3.5)

## clonal size
##
NSCLC_ser[[]] %>% filter(tissue=='tumor') %>% group_by(imid) %>% mutate(patient.total=n()) %>%
  filter(TCRType_new %in%c('pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,TCRType_new,patient.total) %>% 
  summarise(n=n()) %>%
  mutate(per=n/patient.total) %>%
  ggplot(aes(x=TCRType_new,y=log(per),fill=TCRType_new))+
  geom_boxplot(outlier.size = 0.4)+
  theme_classic()+
  scale_fill_manual(values = c("#2166AC","#B2182B")) +
  NoLegend()+
  labs(x='',y='Log clone frequency in tumor')+
  stat_compare_means(comparisons = list(c('non-pTRC','pTRC')))+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('Clone_frequency_nonpTRC_vs_pTRC_cutoff5.tiff',height=4.5, width=3)

## clonal size >5
NSCLC_ser[[]] %>% filter(tissue=='tumor') %>% group_by(imid) %>% mutate(patient.total=n()) %>%
  filter(TCRType_new %in%c('pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,TCRType_new,patient.total) %>% 
  summarise(n=n()) %>%
  mutate(per=n/patient.total) %>%
  filter(n>=5) %>% 
  ggplot(aes(x=TCRType_new,y=log(per),fill=TCRType_new))+
  geom_boxplot(outlier.size = 0.4)+
  theme_classic()+
  scale_fill_manual(values = c("#2166AC","#B2182B")) +
  NoLegend()+
  labs(x='',y='Log clone frequency in tumor')+
  stat_compare_means(comparisons = list(c('non-pTRC','pTRC')))+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('Clone_frequency_nonpTRC_vs_pTRC_cutoff5_over5.tiff',height=4.5, width=3)

## clonal size by patient
NSCLC_ser[[]] %>% filter(tissue=='tumor') %>% group_by(imid) %>% mutate(patient.total=n()) %>%
  filter(TCRType_new %in%c('pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,TCRType_new,patient.total) %>% 
  summarise(n=n()) %>%
  mutate(per=n/patient.total) %>%
  filter(n>=5) %>% 
  ggplot(aes(x=TCRType_new,y=log(per),fill=TCRType_new))+
  geom_boxplot(outlier.size = 0.4)+
  theme_classic()+
  scale_fill_manual(values = c("#2166AC","#B2182B")) +
  facet_wrap(.~imid, ncol=8)+
  scale_y_continuous(expand = c(0.1,0.15))+
  NoLegend()+
  labs(x='',y='Log clone frequency in tumor')+
  stat_compare_means(comparisons = list(c('non-pTRC','pTRC')))+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('Clone_frequency_nonpTRC_vs_pTRC_cutoff5_over5_by_patient.tiff',height=7, width=11)

## by cluster
NSCLC_ser$celltype <- ifelse(NSCLC_ser$celltype=='Effector(III)','NK-like',NSCLC_ser$celltype)
NSCLC_ser[[]] %>% filter(tissue=='tumor') %>% group_by(imid) %>% mutate(patient.total=n()) %>%
  filter(TCRType_new %in%c('pTRC','non-pTRC')) %>%
  group_by(imid,TRB_aa_1,TCRType_new,patient.total) %>% 
  mutate(n=n()) %>%
  mutate(per=n/patient.total) %>%
  filter(n>=5) %>% 
  ggplot(aes(x=TCRType_new,y=log(per),fill=TCRType_new))+
  geom_boxplot(outlier.size = 0.4)+
  theme_classic()+
  scale_fill_manual(values = c("#2166AC","#B2182B")) +
  facet_wrap(.~celltype, ncol=7)+
  scale_y_continuous(expand = c(0.1,0.15))+
  NoLegend()+
  labs(x='',y='Log clone frequency of cell in tumor')+
  stat_compare_means(comparisons = list(c('non-pTRC','pTRC')))+
  theme(axis.text.x = element_text(size=12,face='bold',angle=30,hjust=1,vjust=1,color='black'),
        axis.text.y = element_text(size=12,face='bold'),
        axis.title = element_text(size=12,face='bold'),
        strip.text = element_text(size=12,face='bold'))
ggsave('Clone_frequency_nonpTRC_vs_pTRC_cutoff5_over5_by_cluster.tiff',height=7, width=11)

saveRDS('../14.specific_genes_1/NSCLC_ser_add_scores.rds')