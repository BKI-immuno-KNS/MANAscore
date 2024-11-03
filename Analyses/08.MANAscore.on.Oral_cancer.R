## Oral cancer
## zzeng
rm(list=ls())
library(Seurat)
library(magrittr)
library(dplyr)
library(reshape2)
library(data.table)
library(patchwork)
library("RColorBrewer")


setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/17.Oral_cancer/')
ser <- readRDS('D:/Project/ManaScore/data/Oral_cancer/CD8.tumor.TCRs.RDS')
meta <- ser[[]]
meta$barcode = rownames(meta)

TCR <- readRDS('Oral_cancer_TCR.rds')
meta <- meta %>% left_join(TCR)

### cutoff for post
sample = meta %>% 
  group_by(Patient.ID,pre_post) %>%
  summarise(n=n())

Cutoff = NULL
sample <- sample[sample$Patient.ID!='P30',]
for(k in 1:24){
  
  cutoff= NULL
  ps = sample[k,]$Patient.ID
  tr = sample[k,]$pre_post
  d1 = read.csv(paste0('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/voting_i_Oral_',ps,'_',tr,'_predict_score.csv'),row.names = 1)
  colnames(d1)[2] = 'score_i'
  d2 = read.csv(paste0('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/voting_ni_Oral_',ps,'_',tr,'_predict_score.csv'),row.names = 1)
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
  
  write.csv(d,paste0('../02.MANAscore_hi/02.Oral_cancer/',ps,'_',tr,'_predict_scores_type.csv'),quote=FALSE)
  cutoff = data.frame(patient_ID = ps,
                      pre_post = tr,
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
    labs(x='Imputation MANAscore',y='Non-imputation MANAscore') +
    theme(axis.title = element_text(size=14,face='bold'),
          axis.text = element_text(size=12),
          plot.title = element_blank())
  
  dens1 <- ggplot(d, aes(x = score_i)) + 
    geom_density() + 
    theme_classic()+
    geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="red",size=3)+
    # geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="blue",size=2)+
    labs(title=paste0(ps, '_', tr))+
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
  ggsave(paste0('../02.MANAscore_hi/02.Oral_cancer/',ps,'_',tr,'_correlation.jpeg'),height=6,width=6)
  
  
}
Cutoff
write.csv(Cutoff,'MANAscore_cutoff.csv',quote = FALSE)

##
Cutoff = read.csv('MANAscore_cutoff.csv',row.names = 1)
tcr_color <- meta %>% filter(Tscan.TCR== 'Tumor.reactive' & pre_post=='post') %>% group_by(Patient.ID,TRB_aa_1) %>% summarise(n=n())
tcr_color$col <-  c('red','blue','purple','green') 

for(ps in c('P20','P21','P32')){
  cutoff= NULL
  tr = 'post'
  d1 = read.csv(paste0('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/voting_i_Oral_',ps,'_',tr,'_predict_score.csv'),row.names = 1)
  colnames(d1)[2] = 'score_i'
  d2 = read.csv(paste0('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/voting_ni_Oral_',ps,'_',tr,'_predict_score.csv'),row.names = 1)
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
  
  d <- d %>% left_join(meta)
  
  ##
  plot1 <- ggscatter(d, x = "score_i", y = "score_ni",size = 0.5,title = ps,color='grey16',
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "#008B45FF", fill = "lightgray"), # Customize reg. line
                     conf.int = TRUE # Add confidence interval
  )+
    # geom_point(data= d %>% filter(Tscan.TCR=='Tumor.reactive') , aes(x=score_i,y=score_ni,color=TRB_aa_1),shape=15,size=2)+
    geom_hline(yintercept=cutoff2,color='red',linetype="dashed",size=1)+
    geom_vline(xintercept=cutoff1,color='red',linetype="dashed",size=1)+
    scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    scale_x_continuous(expand = c(0,0),limits = c(0,1))+
    scale_color_manual(values = d %>% filter(Tscan.TCR=='Tumor.reactive')  %>% 
                         group_by(Patient.ID,TRB_aa_1) %>% 
                         summarise(n=n()) %>% 
                         left_join(tcr_color) %>% 
                         as.data.frame() %>% .[,'col']) +
    # geom_rug()+
    stat_cor(method = "pearson",label.x = 0.05, label.y = 0.95,color="#008B45FF",size=4) +
    labs(x='Imputation MANAscore',y='Non-imputation MANAscore',color='') +
    theme(axis.title = element_text(size=14,face='bold'),
          axis.text = element_text(size=12),
          plot.title = element_blank(),
          legend.position = 'bottom',
          legend.key.size = unit(0.25, 'cm'))
  
  dens1 <- ggplot(d, aes(x = score_i)) + 
    geom_density() + 
    theme_classic()+
    geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="red",size=3)+
    # geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="blue",size=2)+
    labs(title=paste0(ps, '_', tr))+
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
  ggsave(paste0('../02.MANAscore_hi/02.Oral_cancer/',ps,'_',tr,'_correlation.tiff'),height=7.5,width=7.5)
  
}

##
Cutoff_post <- Cutoff %>% filter(pre_post == 'post')
SCORE = NULL
for(i in 1:18){
  ps = Cutoff_post[i,]$patient_ID
  tr = Cutoff_post[i,]$pre_post
  score <- read.csv(paste0('../02.MANAscore_hi/02.Oral_cancer/',ps,'_',tr,'_predict_scores_type.csv'),row.names = 1)
  SCORE <- rbind(SCORE,score)
}
SCORE$MANAscoreType[SCORE$it=='hi' & SCORE$nt == 'hi'] = 'MANAscorehi'
SCORE$MANAscoreType[!(SCORE$it=='hi' & SCORE$nt == 'hi')] = 'MANAscorelo'
saveRDS(SCORE,'SCORE.rds')
SCORE <- readRDS('SCORE.rds')

meta_post = meta %>% filter(pre_post=='post')
meta_post = meta_post %>% left_join(SCORE) %>% filter(Patient.ID!='P30') 

##
TCR_shared = meta_post %>% filter(is.na(TRB_aa_1)==FALSE) %>% group_by(TRB_aa_1,Patient.ID) %>%
  summarise(n=n()) %>%
  group_by(TRB_aa_1) %>%
  summarise(n=n()) %>% filter(n>1) %>% as.data.frame() %>% .[,'TRB_aa_1']

meta_post$MANAscoreType[meta_post$TRB_aa_1 %in% TCR_shared] <- 'Null'

MANAhi_clones = meta_post[meta_post$MANAscoreType=='MANAscorehi',] %>% .[,'TRB_aa_1'] %>%  na.omit() %>%
  unique()

MANAlo_clones = meta_post %>% filter(is.na(TRB_aa_1)==FALSE) %>% .[,'TRB_aa_1'] %>% unique() %>%
  setdiff(MANAhi_clones)

meta_post$MANAtype[is.na(meta_post$TRB_aa_1)== FALSE & meta_post$TRB_aa_1 %in% MANAhi_clones] <- 'MANAhi'
meta_post$MANAtype[is.na(meta_post$TRB_aa_1)== FALSE & meta_post$TRB_aa_1 %in% MANAlo_clones] <- 'MANAlo'
meta_post$MANAtype[is.na(meta_post$TRB_aa_1)==TRUE] <- NA


meta_pre = meta %>% filter(pre_post=='pre')
meta_pre$MANAscoreType = 'Null'
meta_pre$MANAtype[is.na(meta_pre$TRB_aa_1)== FALSE  & meta_pre$TRB_aa_1 %in% MANAhi_clones] <- 'MANAhi'
meta_pre$MANAtype[is.na(meta_pre$TRB_aa_1)== FALSE & meta_pre$TRB_aa_1 %in% MANAlo_clones] <- 'MANAlo'
meta_pre %>% filter(MANAtype=='MANAhi') %>% .[,'TRB_aa_1'] %>% unique() ## 65

meta_new = rbind(meta_post[,c(1:27,32,33)],meta_pre)


ser_new = ser[,meta_new$barcode]
ser_new@meta.data = meta_new
saveRDS(ser_new,'Oral_cancer_MANAtype_MANAscoretype_ser.rds')

##
ser_post <- subset(ser,pre_post=='post')
rownames(SCORE) <- SCORE$barcode

ser_post.f <- ser_post[,SCORE$barcode]

ser_post.f$score_i <- SCORE[rownames(ser_post.f@meta.data),]$score_i
ser_post.f$score_ni <- SCORE[rownames(ser_post.f@meta.data),]$score_ni


ser_post.f$MANAscore <- sqrt(ser_post.f$score_i^2+ser_post.f$score_ni^2)
FeaturePlot(ser_post.f,features = c('MANAscore'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  labs(title='MANAscore') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes()
ggsave('Oral_MANAscore_combined.tiff',height=3.2,width=4)

##

ser_post.f$MANA <- ifelse(ser_post.f$Tscan.TCR == 'Tumor.reactive',1,0)
MANA = ser_post.f[[]] %>% filter(Tscan.TCR == 'Tumor.reactive')
umap = ser_post.f@reductions$umap@cell.embeddings %>% as.data.frame()
umap$barcode <- rownames(umap)

ggplot(umap,aes(x = UMAP_1,y=UMAP_2))+
  geom_point(col='grey',size=0.1)+
  geom_point(data = umap %>% filter(barcode %in% rownames(MANA)), aes(x = UMAP_1,y=UMAP_2),color='red',shape=15,size=2)+
  NoLegend()+
  theme_classic()+
  NoLegend()+
  NoAxes()
ggsave('00.MANA_1.tiff',height=3.2,width=3.5)  
saveRDS(ser_post.f,'ser_post_f.rds')
