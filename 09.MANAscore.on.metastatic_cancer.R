## rosenberg's data
rm(list=ls())
library(magrittr)
library(Seurat)
library(preprocessCore)
library(dplyr)
library(data.table)
library(scales)
library(ggbeeswarm)
library(tidyverse)
library(gplots)
library(dplyr)
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


setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/37.rosenberg_data/')
fig2 <- read.csv('D:/TCR/102018-Neountx/Data/public/Steve Rosenberg/fig2.csv')
fig3 <- read.csv('D:/TCR/102018-Neountx/Data/public/Steve Rosenberg/fig3.csv')

d1<-readRDS("D:/TCR/102018-Neountx/Data/public/Steve Rosenberg/12gmacstil_0.85.rds")

d2 <- readRDS('D:/TCR/102018-Neountx/Data/public/Steve Rosenberg/4393til 10x.rds')
d3 <- readRDS('D:/TCR/102018-Neountx/Data/public/Steve Rosenberg/4394til10x.rds')
d4 <- readRDS('D:/TCR/102018-Neountx/Data/public/Steve Rosenberg/4400til10x.rds')
d5 <- readRDS('D:/TCR/102018-Neountx/Data/public/Steve Rosenberg/4421til10x.rds')

fig2$Previously.known.or.newly.identified %>% table()

MANA <- c(fig2[fig2$Previously.known.or.newly.identified %in% c('Known','Newly identified'),]$Barcode,
          fig3[fig3$Previously.known.or.newly.identified %in% c('Known','Newly identified'),]$Barcode) %>% unique()

DimPlot(til_4421)

ser <- readRDS('D:/TCR/102018-Neountx/Result/rosenberg.integrate/integrate_20PCs_2KHVGs.rds')
DimPlot(ser,split.by = 'orig.ident',ncol = 5)

## only CD8 T cell
df <- FetchData(ser,vars=c('CD8A','barcode','sample.n')) 

df %>%
  filter(!sample.n %in% c('WL44211','WL44212')) %>% ## WL4421 introduce strong batch effect
  ggplot(aes(x=CD8A))+
  geom_density()

cell = df %>%
  filter(!sample.n %in% c('WL44211','WL44212')) %>% .[,'barcode']
ser1 <- subset(ser,subset = barcode %in% cell)

df1 <- FetchData(ser1,vars=c('CD8A','barcode','sample.n')) 
d1 <- density(df1$CD8A)
trough <- NULL
for (i in 2:(length(d1$y)-1)) {
  if (d1$y[i-1] >= d1$y[i] & d1$y[i] <= d1$y[i+1]) {
    trough <- cbind(trough, c(d1$x[i], d1$y[i]))
  }
}
cutoff =  trough[1,1] ## 0.6526353

ser_sub <- subset(ser1, CD8A>0.6526353)
DimPlot(ser_sub,split.by = 'orig.ident',ncol = 5)

##
ser_sub$sample <- ser_sub$sample.n
ser_sub$sample <- sub('FrTu|WL|.h5','',ser_sub$sample)

ser_sub$sample <- ifelse(ser_sub$sample %in% c('42981','42982'), '4298',ser_sub$sample)
ser_sub$sample <- ifelse(ser_sub$sample %in% c('43421','43422'), '4342',ser_sub$sample)

table(ser_sub$sample)

## integration again  ###
#########################
library(cowplot)
library(matrixStats)
library(ggpointdensity)

hypervar <- function(data, span = 0.5, showplot = TRUE, font_size = 14){
  
  gene_mean_all <- rowMeans(data)
  gene_var_all <- rowVars(data)
  
  data_filter <- data[gene_mean_all > 0 & gene_var_all > 0,]
  
  gene_mean <- rowMeans(data_filter)
  gene_var <- rowVars(data_filter)
  
  data_fit <- data.frame(X=gene_mean,Y=gene_var)
  fit_model <- loess(formula = log2(x=Y) ~ log2(x=X),
                     data = data_fit,
                     span = span)
  
  gene_var_expect <- fit_model$fitted
  gene_hyper_var <- log2(gene_var) - gene_var_expect
  
  result <- data.frame(feature=row.names(data_filter), mean=gene_mean, var=gene_var,
                       var_expect_log2=gene_var_expect,hypervar_log2=gene_hyper_var)
  
  p1 <- ggplot(result, aes(log2(mean), log2(var))) +
    geom_pointdensity() +
    scale_color_viridis_c() +
    geom_point(data=result,aes(log2(mean),var_expect_log2),color="red") +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
  
  p2 <- ggplot(result, aes(log2(mean), hypervar_log2)) +
    geom_pointdensity() +
    scale_color_viridis_c() +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
  
  combined_plot <- plot_grid(p1, p2, labels = c('A', 'B'))
  
  if(showplot){
    print(combined_plot)
  }
  
  return(list(data=result,plot=combined_plot))
}

theme_set(theme_classic())


### 

ser.list<-SplitObject(ser_sub,split.by = 'sample')

dissociation <- read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/dissociation.csv")  # dissociation genes from biorxiv paper 
ifn<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/ifn.csv")
ig<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/ig.csv")
mito<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/mito.csv")
stress<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/stress.rds")
ncrna<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/ncrna.csv")

exclude<-c(as.character(dissociation$gene_symbol),
           as.character(ncrna$ncRNA_module),
           as.character(stress),
           as.character(ig$Ig_module),
           as.character(mito$mito),
           as.character(stress))


for (i in 1:length(ser.list)) {
  ser.list[[i]] <- NormalizeData(ser.list[[i]] , verbose = FALSE)
  
  ###### WZ 04/16/2020 #######
  #ser.list[[i]]  <- FindVariableFeatures(ser.list[[i]] , selection.method="vst",nfeatures=3000,verbose = FALSE)
  #########get variable feature##########
  gene_count <- as.matrix(ser.list[[i]]@assays$RNA@counts)
  gene_count_norm <- sweep(gene_count,2,colSums(gene_count),FUN="/")*1e4
  gene_hypervar <- hypervar(gene_count_norm,showplot=FALSE)
  gene_hypervar_sort <- gene_hypervar$data %>% arrange(.,desc(hypervar_log2))
  VariableFeatures(ser.list[[i]]) <- setdiff(gene_hypervar_sort$feature,exclude) %>% .[1:3000]
  #######################################
}



###### integrate data #######
# use sample with top N as reference group
temp<-ser_sub[[]]
sample <-temp %>% group_by(sample) %>% summarise(n=n()) %>% filter(n>100)
remain_sample<-which(names(ser.list) %in% sample$sample)


ser.list[remain_sample]

anchors <- FindIntegrationAnchors(object.list = ser.list[remain_sample],
                                  dims = 1:30)  # specify anchoring dataset here

ser.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)


ser.integrated <- ScaleData(ser.integrated, verbose = FALSE) #vars.to.regress = c("cc_score1","S.Score","G2M.Score")
# DefaultAssay(ser.integrated)<-"integrated"

ser.integrated <- RunPCA(object = ser.integrated,npcs =30,verbose = FALSE,
                         features = setdiff(VariableFeatures(object = ser.integrated),exclude))

ElbowPlot(ser.integrated, ndims = 30)
print(ser.integrated[["pca"]], dims = 1:30, nfeatures = 5)
# saveRDS(ser.integrated,'MANAhi_reclustering_mid.rds')
# ser.integrated <- readRDS('MANAhi_reclustering_mid.rds')

ser.integrated <- FindNeighbors(ser.integrated, reduction = "pca", dims = 1:15, nn.eps = 0.5)

ser.integrated <- FindClusters(ser.integrated, resolution = c(0.3,0.4,0.45,0.5,0.55,0.6,0.7), n.start = 10)

ser.integrated <- RunUMAP(object = ser.integrated,reduction='pca',dims = 1:15)
ser.integrated <- RunTSNE(object = ser.integrated,reduction='pca',dims = 1:15)
saveRDS(ser.integrated,'rosenber_CD8_integrated.rds')

ser.integrated <- readRDS('D:/Project/ManaScore/01.Result/36.12_single_moldes/37.rosenberg_data/rosenber_CD8_integrated.rds')

##
DimPlot(ser.integrated,group.by = 'integrated_snn_res.0.5')
Idents(ser.integrated) <- 'integrated_snn_res.0.5'

DimPlot(ser.integrated,group.by = 'integrated_snn_res.0.5',label = T,label.size = 8) &
  NoLegend()& NoAxes()
ggsave('00.UMAP_res0.5.jpeg',height=4,width=4.5)

##
ser.integrated$MANA <- ifelse(ser.integrated$barcode %in% MANA,1,0)

ggplot(umap,aes(x = UMAP_1,y=UMAP_2))+
  geom_point(col='grey',size=0.1)+
  geom_point(data = umap %>% filter(barcode %in% MANA), aes(x = UMAP_1,y=UMAP_2),color='red',shape=17,size=2)+
  NoLegend()+
  theme_classic()+
  NoLegend()+
  NoAxes()
ggsave('00.MANA_1.tiff',height=3.2,width=3.5)  

meta <- ser.integrated@meta.data
meta[meta$barcode %in% MANA,] %>% group_by(integrated_snn_res.0.5) %>% summarise(n=n()) %>%
  arrange(-n) %>%
  mutate(integrated_snn_res.0.5 = factor(integrated_snn_res.0.5,levels = c(3,0,1,2,6,5,7))) %>% 
  ggplot(aes(x=integrated_snn_res.0.5,y=n))+
  geom_col()+
  scale_y_continuous(expand = c(0,0))+
  labs(x='cluster',y='MANA number')+
  theme_classic()+
  theme(axis.title = element_text(face='bold',size=12,color='black'))
ggsave('MANA_cluster_contribution.tiff',height=3.5,width = 5)

## MANAscore prediction
## data preparation
##
saver <- readRDS('rosenber_CD8_integrated_saver.rds')
## imputation
d <- ser.integrated[[]] %>%
  group_by(sample) %>% summarise(n=n()) %>% as.data.frame()

###
gene = c('CXCL13','ENTPD1','IL7R')
RES = NULL
for(i in 1:length(d[,1])){
  print(i)
  p = d[i,]$sample
  
  patient_tissue_barcode <- ser.integrated[[]] %>% filter(sample ==p) %>% 
    rownames() 
  
  if(length(patient_tissue_barcode)>0){
    saver_sub <- saver[,patient_tissue_barcode]
    
    ## quantile normalization
    saver_sub_norm <- normalize.quantiles(saver_sub, copy = TRUE)
    colnames(saver_sub_norm) <- colnames(saver_sub)
    rownames(saver_sub_norm) <- rownames(saver_sub)
    ## scale
    saver_sub_norm_scale <- scale(saver_sub_norm)
    
    subdat = saver_sub_norm_scale[gene,intersect(colnames(saver_sub_norm_scale),patient_tissue_barcode)] %>% t() %>% as.data.frame()
    RES = rbind(RES,subdat)
    
  }
}
write.csv(RES,paste0("../../24.normalize_for_each_patient/rosenberg/rosenberg_saver.csv"))

## non-imputation

df = ser.integrated@assays$RNA@data[gene,] %>% as.matrix() %>% t() %>% as.data.frame() 

write.csv(df,"../../24.normalize_for_each_patient/rosenberg/rosenberg_RNA.csv")

##
score_i <- read.csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/rosenberg/voting_i_rosenberg_predict_score.csv',row.names = 1)
colnames(score_i)[2] <- 'score_i'
score_ni <- read.csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/rosenberg/voting_ni_rosenberg_predict_score.csv',row.names = 1)
colnames(score_ni)[2] <- 'score_ni'

score <- score_i %>% left_join(score_ni)
rownames(score) <- score$barcode

ser.integrated$score_i <- score[rownames(ser.integrated@meta.data),]$score_i
ser.integrated$score_ni <- score[rownames(ser.integrated@meta.data),]$score_ni

p1 = FeaturePlot(ser.integrated,features = c('score_i'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  labs(title='MANAscore_i') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes() & NoLegend()
p2 = FeaturePlot(ser.integrated,features = c('score_ni'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  labs(title='MANAscore_ni') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes()

p1|p2
ggsave('rosenberg_MANAscore.tiff',height=3.2,width=7)

ser.integrated$MANAscore <- sqrt(ser.integrated$score_i^2+ser.integrated$score_ni^2)
FeaturePlot(ser.integrated,features = c('MANAscore'),order=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  labs(title='MANAscore') &
  theme(plot.title=element_text(size=12),
        plot.margin = margin(c(0,0,0,0))) &
  NoAxes()
ggsave('rosenberg_MANAscore_combined.tiff',height=3.2,width=4)

##################
DefaultAssay(ser.integrated) <- 'RNA'
gene = c('CD3E','CD8A','CD8B','ITGAE','ZNF683','GZMK',
         'CD4','CXCL13','ENTPD1',
         'GNLY','KLRD1','NKG7',
         'TCF7','SELL','LEF1','CCR7','IL7R',
         'LAG3','TIGIT','PDCD1','HAVCR2','CTLA4',
         'IL2','GZMA','PRF1','GZMB',
         'CD28','TNFRSF14','ICOS','TNFRSF9',
         'EOMES','HOPX','TBX21','ZEB2','HIF1A','TOX')

FeaturePlot(ser.integrated,features = gene,cols=c('grey','red'),order = T,ncol = 6) & 
  NoLegend() &
  NoAxes()
ggsave('00.feature_gene.jpg',height=12,width=13)

FeaturePlot(ser.integrated,features = c('CXCL13','ENTPD1'),cols=c('grey','red'),order = T,ncol = 2) & 
  NoLegend() &
  NoAxes()
ggsave('00.feature_CXCL13_CD39.tiff',height=3,width=6)


##
markers <- FindAllMarkers(ser.integrated, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DefaultAssay(ser.integrated) <- 'integrated'
DoHeatmap(ser.integrated,features = top10$gene)+NoLegend()

saveRDS(score,'predicted_score.rds')
top10 %>% as.data.frame()

## pTRC
meta <- ser.integrated@meta.data
TCR = rbind(fig2[,c('Barcode','CDR3B','Tumor.ID','Clonotype','Previously.known.or.newly.identified')],
            fig3[,c('Barcode','CDR3B','Tumor.ID','Clonotype','Previously.known.or.newly.identified')])
colnames(TCR)[1] <- 'barcode'

meta <- meta %>% left_join(TCR) 
head(meta)

##
Cutoff = NULL
pa = meta %>% filter(MANA==1) %>% group_by(sample) %>% summarise(n=n()) %>% as.data.frame() %>% .[,'sample'] 

META = NULL
for(ps in pa){
  cutoff= NULL
  d = meta %>% filter(sample==ps) 
  
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
  
  META <- rbind(META, d)
  # cutoff = data.frame(patient_ID = ps,
  #                     cutoff_i = cutoff1,
  #                     cutoff_ni = cutoff2)
  # Cutoff = rbind(Cutoff,cutoff)
  
  ##
  plot1 <- ggscatter(d, x = "score_i", y = "score_ni",size = 0.5,title = ps,
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "#008B45FF", fill = "lightgray"), # Customize reg. line
                     conf.int = TRUE # Add confidence interval
  )+
    geom_point(data= d %>% filter(MANA==1) %>% left_join(TCR), aes(x=score_i,y=score_ni,color=CDR3B),shape=17,size=3)+
    geom_hline(yintercept=cutoff2,color='red',linetype="dashed",size=1.5)+
    geom_vline(xintercept=cutoff1,color='red',linetype="dashed",size=1.5)+
    scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    scale_x_continuous(expand = c(0,0),limits = c(0,1))+
    # geom_rug()+
    stat_cor(method = "pearson",label.x = 0.05, label.y = 0.95,color="#008B45FF",size=4) +
    labs(x='MANAscore_i',y='MANAscore_ni',color='') +
    theme(axis.title = element_text(size=14,face='bold'),
          axis.text = element_text(size=12),
          plot.title = element_blank(),
          legend.position = 'bottom',
          legend.key.size = unit(0.25, 'cm'))+
    guides(color = guide_legend(ncol = 3))
  
  dens1 <- ggplot(d, aes(x = score_i)) + 
    geom_density() + 
    theme_classic()+
    geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="red",size=3)+
    # geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="blue",size=2)+
    labs(title=ps)+
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
  ggsave(paste0('./cutoff/',ps,'_correlation.tiff'),height=7.5,width=7.5)
  
}

for(ps in c('4283','4324','4385')){
  cutoff= NULL
  d = meta %>% filter(sample==ps) 
  META = META %>% filter(sample != ps)
  # Cutoff = Cutoff %>% filter(patient_ID != ps)
  
  di <- density(d$score_i)
  trough1 <- NULL
  for (i in 2:(length(di$y)-1)) {
    if (di$y[i-1] >= di$y[i] & di$y[i] <= di$y[i+1]) {
      trough1 <- cbind(trough1, c(di$x[i], di$y[i]))
    }
  }
  cutoff1 =  trough1[,1][1]
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
  
  META <- rbind(META, d)
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
    geom_point(data= d %>% filter(MANA==1) %>% left_join(TCR), aes(x=score_i,y=score_ni,color=CDR3B),shape=17,size=3)+
    geom_hline(yintercept=cutoff2,color='red',linetype="dashed",size=1.5)+
    geom_vline(xintercept=cutoff1,color='red',linetype="dashed",size=1.5)+
    scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    scale_x_continuous(expand = c(0,0),limits = c(0,1))+
    # geom_rug()+
    stat_cor(method = "pearson",label.x = 0.05, label.y = 0.95,color="#008B45FF",size=4) +
    labs(x='MANAscore_i',y='Non-imputation MANAscore_ni',color='') +
    theme(axis.title = element_text(size=14,face='bold'),
          axis.text = element_text(size=12),
          plot.title = element_blank(),
          legend.position = 'bottom',
          legend.key.size = unit(0.25, 'cm'))+
    guides(color = guide_legend(ncol = 3))
  
  dens1 <- ggplot(d, aes(x = score_i)) + 
    geom_density() + 
    theme_classic()+
    geom_point(aes(x=trough1[,1][1],y=trough1[,1][2]),colour="red",size=3)+
    # geom_point(aes(x=trough1[,dim(trough1)[2]][1],y=trough1[,dim(trough1)[2]][2]),colour="blue",size=2)+
    labs(title=ps)+
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
  ggsave(paste0('./cutoff/',ps,'_correlation.tiff'),height=7.5,width=7.5)
  
}
saveRDS(META,'META.rds')
META <- readRDS('META.rds')
