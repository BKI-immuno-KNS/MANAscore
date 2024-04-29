## integration NSCLC and melanoma
## zzeng
rm(list = ls()) 
library(Seurat)
library(data.table)
library(scales)
library(ggbeeswarm)
library(tidyverse)
library(gplots)
library(dplyr)
library(ggpubr)
library(ggsci)
library(reshape2)
library(ggplot2)
library(heatmap3)

## integrated data of NSCLC and melanoma
# data from melanoma 
me<-readRDS("D:/Dropbox/Scdata/data/antigen specific t cells/Melanoma/tils.CD8.R0.6.harmonized.20200422.rds")
me.count<-me@assays$RNA@counts
dissociation <- read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/dissociation.csv")  # dissociation genes from biorxiv paper 
#ifn<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/ifn.csv")
#ig<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/ig.csv")
#mito<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/mito.csv")
stress<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/stress.rds")
ncrna<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/ncrna.csv")
genes<-rownames(me.count)
exclude<-genes[grepl('TRBV|TRBJ|TRBD|TRBC|TRAV|TRAJ|TRAD|TRDV|TRDJ|TRDC|TRGV|TRGJ|TRGC|MT-|RP[SL]',substr(genes,1,4))]  # remove genes with TRB/TRA/TRD/TRG
exclude<-c(as.character(dissociation$gene_symbol),
           as.character(ncrna$ncRNA_module),
           as.character(stress),
           exclude) # genes to exclude

me.count<-me.count[!rownames(me.count) %in% exclude,]


ser.me<-CreateSeuratObject(me.count,min.cells =5,min.features = 250,meta.data = me.meta)
ser.me$orig.ident<-ser.me$patient

####get variable feature####
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

dir.create("D:/TCR/102018-Neountx/Result/lung.melanoma.integrate/")
path<-"D:/TCR/102018-Neountx/Result/lung.melanoma.integrate/"
annotate<-read.csv("D:/TCR/102018-Neountx/Data/Processed/Scseq/annotate.new.dge.csv")
annotate$sample.n<-paste0(annotate$patient_id,':',annotate$tissue,'-',annotate$batch)
trab.wide<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/trab.wide.public.rds")
ser.count<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/ser.count.rds")
ser.count<-ser.count@assays$RNA
cells<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/cells.pass.qc.rds")
qc.meta<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/qc.rds")
qc.meta<-qc.meta %>% left_join(unique(annotate[,c('sample','sample.n')]))
rownames(qc.meta)<-qc.meta$barcode
file<-annotate %>% filter(mana.include==1)
cd8<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/DGE/saver/cutoff/CD8A.cutoff.rds")
cd4<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/DGE/saver/cutoff/CD4.sd.cutoff.rds")

check<-qc.meta %>% left_join(trab.wide.neotx) %>% filter(tissue=='tumor') %>%
  group_by(antigen, patient_id) %>% summarise(n=n())
# genes to be excluded from clustering
ifn<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/ifn.csv")
ig<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/ig.csv")
mito<-read.csv("D:/TCR/102018-Neountx/Data/Raw/Published gene list/mito.csv")
cluster.exclude<-c(ifn$IFN_module,ig$Ig_module,mito$mito)

ser.count<-ser.count[,sub('_.*','',colnames(ser.count)) %in% file$sample&colnames(ser.count) %in% cells&colnames(ser.count) %in% cd8]
ser<-CreateSeuratObject(ser.count,min.cells =5,min.features = 250,meta.data = qc.meta)
ser$orig.ident<-ser$sample.n

ser.combined <- merge(ser, y = ser.me)

test<-ser.combined[[]] %>% group_by(orig.ident) %>% summarise(n=n())
ser.combined<-subset(ser.combined,subset=orig.ident!='p6') ## remove p6
ser.list<-SplitObject(ser.combined,split.by = 'orig.ident')

for (i in 1:length(ser.list)) {
  ser.list[[i]] <- NormalizeData(ser.list[[i]] , verbose = FALSE)
  
  #########get variable feature##########
  gene_count <- as.matrix(ser.list[[i]]@assays$RNA@counts)
  gene_count_norm <- sweep(gene_count,2,colSums(gene_count),FUN="/")*1e4
  gene_hypervar <- hypervar(gene_count_norm,showplot=FALSE)
  gene_hypervar_sort <- gene_hypervar$data %>% arrange(.,desc(hypervar_log2))
  VariableFeatures(ser.list[[i]]) <- gene_hypervar_sort$feature[1:3000]
  VariableFeatures(ser.list[[i]])<-setdiff(VariableFeatures(ser.list[[i]]),cluster.exclude)
  #######################################
}

###### integrate data #######
# use sample with top N as reference group
temp<-ser.combined[[]]
top.sample<-temp %>% group_by(sample.n,tissue,response) %>% summarise(n=n()) %>%
  group_by(tissue,response) %>% top_n(1,wt=n)
ref<-which(names(ser.list) %in% c(top.sample$sample.n,ser.me$orig.ident))
ref<-as.numeric(ref)

ser.list[ref]

anchors <- FindIntegrationAnchors(object.list = ser.list,reference=ref,
                                  dims = 1:30)  # specify anchoring dataset here

ser.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)


ser.integrated <- ScaleData(ser.integrated, verbose = FALSE) #vars.to.regress = c("cc_score1","S.Score","G2M.Score")
DefaultAssay(ser.integrated)<-"integrated"

ser.integrated <- RunPCA(object = ser.integrated,npcs =30,verbose = FALSE,
                         features =setdiff(VariableFeatures(object = ser.integrated),cluster.exclude))

ElbowPlot(ser.integrated, ndims = 100)
print(ser.integrated[["pca"]], dims = 1:30, nfeatures = 5)
ser.integrated <- FindNeighbors(ser.integrated, reduction = "pca", dims = 1:30, nn.eps = 0.5)

ser.integrated <- FindClusters(ser.integrated, resolution = c(0.7,0.3,0.6,0.5), n.start = 10)
ser.integrated <- RunUMAP(object = ser.integrated,reduction='pca',dims = 1:30)

ser.integrated<-subset(ser.integrated,subset=integrated_snn_res.0.5 %in% c(0:13))

top3000 <- head(VariableFeatures(ser.integrated), 3000) 
dir.create(paste0(path,'/'))
saveRDS(top3000,paste0(path,"/top3000.hvg.rds"))

ser.integrated[[]] %>% group_by(integrated_snn_res.0.5) %>% summarise(n=n())
# find diffrential markers
ser.integrated.small<-subset(ser.integrated,downsample=5000)
markers <- FindAllMarkers(ser.integrated.small, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
tail20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = -avg_logFC)
write.csv(markers,paste0(path,"/markers.rawcount.csv"),row.names = F)
write.csv(top20,paste0(path,"/top20.rawcount.csv"),row.names = F)
write.csv(tail20,paste0(path,"/tail20.rawcount.csv"),row.names = F)
write.csv(top10,paste0(path,"/top10.rawcount.csv"),row.names = F)



DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.3' ,label = F)
ggsave(paste0(path, "/cluster.0.3.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.5' ,label = F)
ggsave(paste0(path, "/cluster.0.5.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.7' ,label = F)
ggsave(paste0(path, "/cluster.0.7.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.3' ,label = T)
ggsave(paste0(path, "/cluster.0.3.label.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.5' ,label = T)
ggsave(paste0(path, "/cluster.0.5.label.tiff"),height = 6,width =7,dpi = 300)
ggsave('D:/Dropbox/Scdata/result/manascore/cluster.0.5.label.tiff',width = 5,height = 3.5)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.6' ,label = T)
ggsave(paste0(path, "/cluster.0.6.label.tiff"),height = 6,width =7,dpi = 300)

DimPlot(ser.integrated,reduction = 'umap',group.by = 'integrated_snn_res.0.7' ,label = T)
ggsave(paste0(path, "/cluster.0.7.label.tiff"),height = 6,width =7,dpi = 300)


cluster<-Idents(ser.integrated)
umap<-ser.integrated@reductions$umap@cell.embeddings

saveRDS(top10, paste0(path,"/top10.rds"))
saveRDS(ser.integrated,paste0(path,"/ser.integrate.rds"))
saveRDS(cluster,paste0(path,"/cluster.rds"))
saveRDS(umap, paste0(path,"/umap.rds"))

######### plot ############
ser <- readRDS("D:/TCR/102018-Neountx/Result/lung.melanoma.integrate/ser.integrate.rds")
ser <- subset(ser,subset=(integrated_snn_res.0.5 %in% c(0:13)))
dir.create("D:/Project/ManaScore/01.Result/00.Lung_melanoma_integration/")
setwd("D:/Project/ManaScore/01.Result/00.Lung_melanoma_integration/")
# saveRDS(ser,"ser_lung_me_integrated_fi.rds")
ser <- readRDS("ser_lung_me_integrated_fi.rds")

ser@meta.data
DimPlot(ser,label = T) + NoLegend()+
  theme(plot.title=element_blank())
ggsave("00.cluster.0.5.lable.tiff",height=5,width=5.5)


FeaturePlot(ser,features = c("CD8A","CD4","FOXP3","GZMK","ITGAE",
                             "ZNF683","TCF7","CXCL13","SLC4A10",
                             "PDCD1","CTLA4","HAVCR2","TIGIT","ENTPD1","LAG3"),
            cols = c("grey","red"),ncol = 3) & NoLegend() & NoAxes()

ggsave("01.specific_gene_plot.tiff",height=10,width=6)

##heatmap plot
celltype<-read.csv('D:/TCR/102018-Neountx/Result/lung.melanoma.integrate/celltype.csv')
celltype[11,2] <- 'Naive/stem-like'
celltype$CellType <- sub(" \\(","\\(",celltype$CellType)
ser.small<-subset(ser,downsample=5000)
meta <- ser.small@meta.data
colnames(celltype)[1] <- "integrated_snn_res.0.5"
meta <- merge(meta,celltype,by='integrated_snn_res.0.5')
rownames(meta) <- meta$barcode
ser.small$CellType <- meta[rownames(ser.small@meta.data),]$CellType

markers<-read.csv("D:/TCR/102018-Neountx/Result/lung.melanoma.integrate/markers.rawcount.csv")
top5<-markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)  


DefaultAssay(ser.small)<-'integrated'
ser.small$CellType <- factor(ser.small$CellType,
                             levels = c("Effector(I)","Effector(II)",
                                        "TRM(I)","TRM(II)","TRM(III)","TRM(IV)","TRM(V)","TRM(VI)",
                                        "NK-like","NK-like TRM","Naive/stem-like","CD4CD8","Proliferating","MAIT"))
clu <- c(0,2,
         4,1,9,12,3,7,8,6,10,5,13,11)
genes = NULL
for(i in clu){
  g = top5 %>% filter(cluster==i) %>%as.data.frame() %>% .[,'gene'] 
  genes = c(genes,g)
  
}

DoHeatmap(ser.small,features = genes,angle=90,group.by = 'CellType')+NoLegend()
ggsave("02.top5_marker_gene.tiff",height=9,width=7)

##
trab.wide<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/trab.wide.public.rds")

umap <- Embeddings(ser,reduction = "umap")
umap<-as.data.frame(umap)
umap$barcode<-rownames(umap)
newmeta<-ser[[]]

newmeta<-newmeta[,c(1:180,183:189)]
umap<-newmeta %>% left_join(umap) %>% left_join(trab.wide,by=c('barcode','patient_id'))

umap$antigen[umap$antigen=='Viral (CMV, EBV, Influenza A)'&umap$Type=='EBV']<-'EBV'
umap$antigen[umap$antigen=='VIRAL'&umap$Type=='EBV']<-'EBV'
umap$antigen[umap$antigen=='Viral (Influenza)']<-'InfluenzaA'
umap$antigen[umap$antigen=='VIRAL'&umap$Type=='InfluenzaA']<-'InfluenzaA'
umap$antigen[umap$antigen=='Flu']<-'InfluenzaA'

umap$antigen[umap$activity=='nontumor'&is.na(umap$Cathegory)==T]<-'Non.Tumor.antigen.NA'
umap$antigen[umap$activity=='tumor'&is.na(umap$Cathegory)==T]<-'Tumor.antigen.NA'

umap$antigen[umap$Cathegory=='NeoAg'|umap$Cathegory=='Multi']<-'MANA'
umap$antigen[umap$Cathegory=='Viral'&umap$NeoAg=='EBV']<-'EBV'
umap$antigen[umap$Cathegory=='Viral'&umap$NeoAg=='Flu']<-'InfluenzaA'
umap$antigen[umap$antigen=='Viral (CMV, EBV, Influenza A)'&umap$Type=='InfluenzaA']<-'InfluenzaA'
umap$antigen[umap$Cathegory=='MAA']<-'MAA'
umap$antigen[is.na(umap$antigen)]<-'antigen.NA'
umap$cancer[is.na(umap$tissue)==T]<-'melanoma'
umap$cancer[is.na(umap$tissue)==F]<-'NSCLC'
umap$tissue[is.na(umap$tissue)==T]<-'melanoma tumor bed'

rownames(umap) <- umap$barcode
table(umap$Type)
##
celltype<-read.csv('D:/TCR/102018-Neountx/Result/lung.melanoma.integrate/celltype.csv')
names(celltype)[1]<-'seurat_clusters'
celltype$CellType <- sub("Na\xefve/stem-like","Naive/stem-like",celltype$CellType)
celltype$CellType <- sub("Effector ","Effector",celltype$CellType)

celltype$seurat_clusters<-as.factor(celltype$seurat_clusters)
umap<-umap %>% left_join(celltype)

umap$patient_id[is.na(umap$patient_id)==T]<- umap$patient[is.na(umap$patient_id)==T]
rownames(umap) <- umap$barcode
saveRDS(umap,'umap.rds')
# umap <- readRDS("umap.rds")
umap <- readRDS("D:/Project/ManaScore/01.Result/00.TCR_info/00.umap_new.rds")
rownames(umap) <- umap$barcode
ser <- AddMetaData(ser,metadata = umap[rownames(ser@meta.data),c("cancer","tissue","celltype","antigen")],col.name = c("cancer","tissue","CellType","antigen"))

umap %>% group_by(antigen,Type,Cathegory) %>% dplyr::summarise(n=n()) 

##
ser$cancer <- factor(x = ser$cancer, levels = c("NSCLC","melanoma"))
DimPlot(ser,reduction="umap",split.by = "cancer") + NoLegend() + NoAxes()
ggsave("03.cluster_split_by_cancer.tiff",height=5,width=11)

ser$tissue <- factor(x = ser$tissue, levels = c("normal","tumor","mettumor","LN","melanoma tumor bed"))
DimPlot(ser,reduction="umap",split.by = "tissue",ncol = 2) + NoLegend() + NoAxes()
ggsave("04.cluster_split_by_tissue.tiff",height=8,width=5.3)
##
head(umap)
umap$cancer <- factor(umap$cancer,levels = c("NSCLC","melanoma"))

umap$dataset <- ifelse(umap$cancer=='NSCLC','Caushi et al., NSCLC','Oliveira et al., melanoma')
ggplot(umap,aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(aes(color=seurat_clusters),alpha=0.2,size=0.8,color='grey')+
  geom_point(data=umap %>% filter(antigen=='InfluenzaA'),col='blue',size=2,shape=16)+
  geom_point(data=umap %>% filter(antigen=='EBV'),col='purple',size=2,shape=16)+
  geom_point(data=umap %>% filter(antigen=='MANA'),col='red',size=2,shape=17)+
  geom_point(data=umap %>% filter(Cathegory=='MAA'),col='green',size=2,shape=16)+
  theme_classic() + 
  facet_grid(.~dataset)+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size=20,face='bold'))
ggsave("05.lung.mana.cef.cluster.colored.tiff",width =10, height = 5)

head(celltype)
library(RColorBrewer)

# col<-colorRampPalette(brewer.pal(name="RdYlBu", n = 8))(14)
col<-c(pal_npg("nrc",alpha = 0.5)(8),pal_aaas(alpha = 0.8)(6))

umap$imid <- ifelse(umap$cancer=="melanoma",umap$patient,umap$imid)

umap %>% filter(antigen %in% c('MANA','mana_score')) %>% group_by(imid,cancer,tissue,seurat_clusters) %>%
  summarise(n=n()) %>% group_by(imid) %>% mutate(total=sum(n)) %>%
  filter(total>=30 & tissue %in% c('Mixed','tumor')) %>%
  left_join(celltype) %>%
  ggplot(aes(x=imid,y=n,fill=CellType))+
  geom_col(position="fill") +
  theme_classic()+
  scale_fill_manual(values = col)+
  scale_y_continuous(expand = c(0,0.01),labels = scales::percent)+
  labs(y='',x='')+
  facet_wrap(~cancer,scales='free_x',nrow=1)+
  theme(axis.text.x = element_text(angle =40, hjust = 1,size = 14, face = "bold"),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14,face = "bold"),
        axis.title.x =element_text(size = 14,face = "bold"),
        axis.title.y=element_text(size = 12,face = "bold")
  )
ggsave("06.mana.prop.tiff",width = 6,height =5.5)

total = umap %>%  group_by(imid,cancer,tissue) %>%
  summarise(total=n()) %>% 
  filter(tissue %in% c('Mixed','tumor'))

tumor_MANA_proportion = umap %>% filter(antigen %in% c('MANA')) %>% group_by(imid,cancer,tissue) %>%
  summarise(n=n()) %>% group_by(imid)  %>%
  filter(tissue %in% c('Mixed','tumor')) %>% left_join(total[,c("imid","total")])
tumor_MANA_proportion$prop = tumor_MANA_proportion$n/tumor_MANA_proportion$total


umap %>% filter(antigen %in% c('MANA'))  %>%
  filter(tissue %in% c('Mixed','tumor'))%>% group_by(cancer,seurat_clusters) %>%
  summarise(n=n()) %>% 
  left_join(celltype) %>% as.data.frame()

# saveRDS(ser,"ser_lung_me_integrated_fi.rds")