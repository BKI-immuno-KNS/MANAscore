## MCPyV+ Merkel
## zzeng
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(bslib)
library(dplyr)
library(Seurat)
library("ggsci")
library(patchwork)
library(magrittr)
library(dplyr)
library(ggpubr)
library(scales)
library(RColorBrewer)

setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/44.Merkel/')
ser <- readRDS('D:/MANAscore/share/object/Merkel/ser_antigenSpeCllInfo_outCellRm.rds')
score_i <- read.csv('D:/Project/tmp/Merkel/MANAscore/voting_i_Merkel_predict_score.csv')
colnames(score_i)[2] <- 'score_i'
score_ni <- read.csv('D:/Project/tmp/Merkel/MANAscore/voting_ni_Merkel_predict_score.csv')
colnames(score_ni)[2] <- 'score_ni'

score <- score_i %>% left_join(score_ni) %>%
  mutate(score=sqrt(score_i^2 + score_ni^2))
rownames(score) <- score$barcode

ser$score_i <- score[rownames(ser[[]]),]$score_i
ser$score_ni <- score[rownames(ser[[]]),]$score_ni
ser$score <- score[rownames(ser[[]]),]$score

ser_pos <- subset(ser, subset= patient %in% c('Z1513','Z1504','Z1525'))

ser_pos_noC5 <- subset(ser_pos, subset=seurat_clusters!=5)

FeaturePlot(ser_pos_noC5,features = c('score'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  # labs(title='Imputation MANAscore') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes() 
ggsave('Merkel_pos_MANAscore_combined_noC5.tiff',height=3.2, width=4)

##
p1 = FeaturePlot(ser_pos_noC5,features = c('score_i'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  labs(title='MANAscore_i') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes() & NoLegend()
p2 = FeaturePlot(ser_pos_noC5,features = c('score_ni'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  labs(title='MANAscore_ni') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes()
p1|p2
ggsave('Merkel_pos_MANAscore_noC5.tiff', height=3.2, width=7)


## MCPyV
ser_pos_noC5$MCPyV_clone <- ifelse(ser_pos_noC5$`MCPyV_specific Clone T cells c1`=='1', 1, 0)
FeaturePlot(ser_pos_noC5, features = c('MCPyV_clone'), cols = c('grey','red'), order=T) + NoLegend() + NoAxes()
ggsave('Merkel_pos_MCPyV.tiff', height=3.2, width=3.5)
##
ser_pos_noC5$MCPyV <- ifelse(ser_pos_noC5$`MCPyV_specific T cells`=='1', 1, 0)

##
umap <- ser_pos_noC5@reductions$umap@cell.embeddings %>% as.data.frame()
umap$barcode <- rownames(umap) 
umap$seurat_clusters <- ser_pos_noC5$seurat_clusters

umap$MCPyV_specific <- ser_pos_noC5[[]]$`MCPyV_specific T cells`
umap$MCPyV_specific_clone <- ser_pos_noC5[[]]$`MCPyV_specific Clone T cells c1`
umap %>% 
  # filter(seurat_clusters!=5) %>%
  filter(MCPyV_specific=='1' & MCPyV_specific_clone=='1') %>% dim() ## 823
umap %>% 
  # filter(seurat_clusters!=5) %>% 
  filter(MCPyV_specific_clone=='1') %>% dim() ## 4473

umap %>% 
  # filter(seurat_clusters!=5) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(col='grey',size=0.1) +
  geom_point(data=umap%>% filter(MCPyV_specific=='1' & MCPyV_specific_clone=='1'),aes(x=UMAP_1, y=UMAP_2), color='red', size=0.5)+
  NoLegend()+
  theme_classic() +
  NoAxes()
ggsave('Merkel_pos_MCPyV_by_CITEseq_no_C5.tiff', height=3.2, width=3.2)

##
umap %>% 
  # filter(seurat_clusters!=5) %>% 
  ggplot(aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(col='grey',size=0.1) +
  geom_point(data=umap%>% filter(MCPyV_specific_clone=='1'),aes(x=UMAP_1, y=UMAP_2), color='red', size=0.5)+
  NoLegend()+
  theme_classic() +
  NoAxes()
ggsave('Merkel_pos_MCPyV_clone_by_CITEseq.tiff', height=3.2, width=3.2)

TCR <- readRDS('trab.wide.rds')
# ser_pos$TRB_aa_1 <- ser_pos[[]] %>% left_join(TCR) %>% .[,'TRB_aa_1']
ser_pos_noC5$TRB_aa_1 <- ser_pos_noC5[[]] %>% left_join(TCR) %>% .[,'TRB_aa_1']
## scatter plot

Cutoff = NULL
D=NULL
# MANAhi_putative_overlapped = c()
# pTCR_cutoff = 1 #<---adjust the cutoff of MANA clone size here

sample <- c("Z1504", "Z1513", "Z1525")

for(k in 1:length(sample)){
  cutoff= NULL
  ps = sample[k]
  print(ps)
  
  d = ser_pos_noC5[[]] %>%  filter(patient==ps) 
  
  di <- density(d$score_i)
  trough1 <- NULL
  for (i in 2:(length(di$y)-1)) {
    if (di$y[i-1] >= di$y[i] & di$y[i] <= di$y[i+1]) {
      trough1 <- cbind(trough1, c(di$x[i], di$y[i]))
    }
  }
  cutoff1 =  trough1[,dim(trough1)[2]][1]
  cutoff1_y =  trough1[,dim(trough1)[2]][2]
  d$it = ifelse(d$score_i>cutoff1,'hi','lo')
  dni <- density(d$score_ni)
  trough2 <- NULL
  for (i in 2:(length(dni$y)-1)) {
    if (dni$y[i-1] >= dni$y[i] & dni$y[i] <= dni$y[i+1]) {
      trough2 <- cbind(trough2, c(dni$x[i], dni$y[i]))
    }
  }
  cutoff2 =  trough2[,dim(trough2)[2]][1]
  cutoff2_y = trough2[,dim(trough2)[2]][2]
  d$nt = ifelse(d$score_ni>cutoff2,'hi','lo')
  
  cutoff = data.frame(patient_ID = ps,
                      cutoff_i = cutoff1,
                      cutoff_ni = cutoff2,
                      cutoff_i_y = cutoff1_y,
                      cutoff_ni_y = cutoff2_y)
  Cutoff = rbind(Cutoff,cutoff)
  D = rbind(D,d)
  
  MANAscorehi_TCR = d %>% filter(MCPyV_clone==1 & is.na(TRB_aa_1)==FALSE & it=='hi' & nt=='hi') %>% .[,'TRB_aa_1'] %>% unique() 
  
  hi_light = d %>% filter(MCPyV_clone==1) %>%
    mutate(shared = ifelse(TRB_aa_1 %in% MANAscorehi_TCR, 'Yes','No'))
  print(table(hi_light$shared))
  
  ##
  plot1 <- ggscatter(d %>% filter(!barcode %in% hi_light$barcode), x = "score_i", y = "score_ni",size = 0.5,title = ps,
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "#008B45FF", fill = "lightgray"), # Customize reg. line
                     conf.int = TRUE  # Add confidence interval
  )+ 
    geom_point(data=hi_light %>% filter(shared=='Yes'), aes(x=score_i,y=score_ni), color='red',size=1.5) +
    geom_point(data=hi_light %>% filter(shared=='No'), aes(x=score_i,y=score_ni), color='green',size=1.5)+ 
    theme(legend.position = 'bottom') + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    geom_hline(yintercept=cutoff2,color='blue',linetype="dashed",size=1.5)+
    geom_vline(xintercept=cutoff1,color='blue',linetype="dashed",size=1.5)+
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
    geom_point(aes(x=cutoff1,y=cutoff1_y),colour="blue",size=3)+
    # geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="blue",size=2)+
    labs(title=paste0(ps))+
    scale_x_continuous(expand=c(0,0),limits = c(0,1))+
    theme(legend.position = "none",axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=16,hjust=0.5,face='bold'))
  
  dens2 <- ggplot(d, aes(x = score_ni)) +
    geom_density(alpha = 0.4) +
    theme_classic() +
    geom_point(aes(x=cutoff2,y=cutoff2_y),colour="blue",size=3)+
    scale_x_continuous(expand=c(0,0),limits = c(0,1))+
    theme(legend.position = "none",axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    coord_flip()
  
  dens1 + plot_spacer() + plot1 + dens2 +
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  ggsave(paste0(ps,'_correlation.tiff'),height=6,width=7.5)
}

##
D$MANAscoreType <- ifelse(D$it=='hi' & D$nt=='hi', 'MANAscorehi','MANAscorelo')
D$clonal.size = D %>% group_by(patient, TRB_aa_1) %>% mutate(clonal.size=n()) %>% as.data.frame() %>% .[,'clonal.size']

D$clonal.size = ifelse(is.na(D$TRB_aa_1)==TRUE, D$TRB_aa_1, D$clonal.size)

##
MANAscorehi_TCR <- D %>% filter(MANAscoreType=='MANAscorehi' & is.na(TRB_aa_1)==FALSE) %>% 
  group_by(patient,TRB_aa_1) %>%
  summarise(n=n()) %>%
  filter(n>4) ## 74

patient_TCR_MANAscore_num <- D %>% filter(is.na(TRB_aa_1)==FALSE) %>% group_by(patient, TRB_aa_1, MANAscoreType) %>%
  summarise(n=n()) %>%
  reshape2::dcast(patient+TRB_aa_1~MANAscoreType,value.var = 'n') %>%
  replace(is.na(.),0)

D <- D %>% left_join(patient_TCR_MANAscore_num[,c('patient','TRB_aa_1','MANAscorehi')])
rownames(D) <- D$barcode
##
D$TCRType = ifelse(D$TRB_aa_1 %in% MANAscorehi_TCR$TRB_aa_1 & D$MANAscorehi>4,'pTRC','non-pTRC')

D$TCRType = ifelse(is.na(D$TRB_aa_1)==TRUE, D$TRB_aa_1, D$TCRType)

TCRType_count <- D %>% filter(is.na(TRB_aa_1)==FALSE) %>% group_by(patient, TRB_aa_1, TCRType, MANAscoreType) %>% summarise(n=n()) %>% 
  reshape2::dcast(patient+TRB_aa_1+TCRType~MANAscoreType, value.var = 'n') %>%
  replace(is.na(.),0) %>% arrange(-MANAscorehi) %>%
  mutate(total=MANAscorehi+MANAscorelo)

write.csv(TCRType_count,'Merkel_pTRC_nonPTRC_count.csv',quote=FALSE, row.names = FALSE)

##

ser_pos_noC5$MANAscoreType <- D[rownames(ser_pos_noC5[[]]),]$MANAscoreType
ser_pos_noC5$clonal.size <- D[rownames(ser_pos_noC5[[]]),]$clonal.size
ser_pos_noC5$clonal.size <- as.numeric(ser_pos_noC5$clonal.size)
ser_pos_noC5$TCRType <- D[rownames(ser_pos_noC5[[]]),]$TCRType

ser_pos_noC5$MANAscorehi_pTRC <- ifelse(ser_pos_noC5$TCRType=='pTRC' & ser_pos_noC5$MANAscoreType=='MANAscorehi',1,0)
ser_pos_noC5$pTRC <- ifelse(ser_pos_noC5$TCRType=='pTRC',1,0)

ser_pos_noC5$pTRC  <- ifelse(is.na(ser_pos_noC5$pTRC)==TRUE,0, ser_pos_noC5$pTRC)

FeaturePlot(ser_pos_noC5, features = c('MANAscorehi_pTRC','pTRC'), cols = c('grey','red'),order=T) &
  NoLegend() & NoAxes() 
ggsave('UMAP_MANAscoreehi_pTRC_pTRC_noC5.tiff',height=4,width = 7)

## 
ser_pos_noC5[[]] %>% filter(clonal.size > 4 & MCPyV_clone==1) %>% dim() ## 4443
ser_pos_noC5[[]] %>% filter(clonal.size > 4 &pTRC==1 & MCPyV_clone==1) %>% dim() ## 4401
ser_pos_noC5[[]] %>% filter(pTRC==1 & MCPyV_clone==0) %>% dim() ## 4239

MCPyV_clone_TCR <- ser_pos_noC5[[]] %>% filter(MCPyV_clone==1) %>% group_by(patient,TRB_aa_1) %>% summarise(n=n()) # 33 clones
MCPyV_clone_TCR_over5 <- ser_pos_noC5[[]] %>% filter(MCPyV_clone==1 & clonal.size>=5) %>% group_by(patient,TRB_aa_1) %>% summarise(n=n()) ## 22 clones
pTRC_TCR <- ser_pos_noC5[[]] %>% filter(pTRC==1) %>% group_by(patient,TRB_aa_1) %>% summarise(n=n()) ## 74
ser_pos_noC5[[]] %>% filter(pTRC==1 & MCPyV_clone==1) %>% group_by(patient,TRB_aa_1) %>% summarise(n=n()) ## 16

##
install.packages('eulerr')
library(eulerr)

fit <- euler(c(pTRC = 4239,  `MCPyV-specific T cell`= 42,
               "pTRC&MCPyV-specific T cell" = 4401))
tiff('Merkel_pos_pTRC_MCPyV_venn_noC5.tiff', height=400,width = 450, units = 'px')
plot(fit,
     fills = c("skyblue", "red"),
     edges = TRUE,
     labels = NULL,
     # fontsize = 10,
     quantities = list(fontsize = 0))
dev.off()

##
ser_pos_downsample <- subset(ser_pos_noC5, downsample=1000)
FeaturePlot(ser_pos_downsample, features = c('CXCL13','ENTPD1'), 
            cols = c('grey','red'), order=T) & NoLegend() & NoAxes()
ggsave('Merkel_pos_CXCL13_CD39.tiff', height=3, width=6)

genes <- c('CD8A','CD4','ITGAE','ZNF683','GZMK','GZMB','PRF1',
           'IL7R','TCF7','CCR7','SELL','GNLY','IFNG','NKG7',
           'PDCD1','CTLA4','HAVCR2','TIGIT','LAG3','CXCL13','ENTPD1')
FeaturePlot(ser_pos_downsample, features = genes, cols = c('grey','red'),order=T, ncol = 7) & 
  NoLegend() &NoAxes()
ggsave('Merkel_pos_marker.tiff', height=7.5, width=13.75)
saveRDS(ser_pos_noC5[[]], 'ser_pos_noC5_meta.rds')

# install.packages('wesanderson')
library(wesanderson)
table(ser$integrated_snn_res.0.5) 
DimPlot(ser_pos_noC5, group.by = 'seurat_clusters') &
  scale_color_manual(values =brewer.pal(10, 'Paired')[-6]) & NoLegend() & NoAxes() 
ggsave('Merkel_cluster_umap_without_C5.tiff', height=5, width=4.8)
Idents(ser_pos) <- 'seurat_clusters'
marker <- FindAllMarkers(ser_pos_noC5, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)

marker_sig <- marker %>% filter(p_val_adj<0.05) 
write.csv(marker_sig, 'Merkle_cluster_marker.csv', quote=FALSE, row.names = FALSE)
