## dge analysis
## zzeng
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
library(ggrepel)

setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/20.pseudobulk_dge/')
## dge analysis between pTRC vs non-pTRC
NSCLC_ser <- readRDS('D:/Project/ManaScore/01.Result/36.12_single_moldes/43.cutoff5/NSCLC_ser_new_cutoff5.rds')
MANA_tumor_ser <- subset(NSCLC_ser,subset = TCRType %in% c('MANAhi','MANAlo') & tissue == 'tumor') 
## filter clonalsize >1
tcr_f <- MANA_tumor_ser[[]] %>% group_by(TRB_aa_1,imid,response) %>% summarise(n=n()) %>% filter(n>1)
MANA_tumor_ser_f <- subset(MANA_tumor_ser,subset = TRB_aa_1 %in% tcr_f$TRB_aa_1) ## 6087

meta <- MANA_tumor_ser_f[[]]

rownames(meta) <- meta$barcode
meta$group <- paste0(meta$batch,meta$center)
meta$comb.lane <- paste0(sub('-[0-9]*_.*','',rownames(meta)),'_',meta$group)
aa <- meta %>% group_by(sample.n) %>% tally
meta.new = left_join(meta,aa)
TCR = meta$TRB_aa_1
names(TCR) = meta1$TRB_aa_1

meta <- meta.new %>% dplyr::select(barcode,orig.ident,imid,patient_id,tissue,comb.lane,sample.n,resi_tumor,response,center) %>%
  mutate(patient.tissue = paste0(patient_id,':',tissue))

##
count <- MANA_tumor_ser_f@assays$RNA@counts
TCR <- MANA_tumor_ser_f$TRB_aa_1
uTCR <- unique(TCR) %>% sort()

s <- count

res = NULL
for(sTCR in uTCR){
  print(sTCR)
  tmp <- s[,colnames(s) %in% names(TCR)[which(TCR==sTCR)]]
  new <- inner_join(data.frame(barcode = colnames(tmp)),meta)
  p <- new$imid
  t = sapply(unique(p),function(sp) rowSums(tmp[,p==sp,drop=F]))
  colnames(t) = paste0(colnames(t),":",sTCR)
  res = cbind(res,t)
}

####
tcr_f$tcr = paste0(tcr_f$imid,":",tcr_f$TRB_aa_1)
tcr_f <- as.data.frame(tcr_f)
rownames(tcr_f) <- tcr_f$tcr
tcr_f <- tcr_f %>% na.omit()
res <- res[,tcr_f$tcr]

tcr_f$TCRType <- ifelse(tcr_f$TRB_aa_1 %in% tcr_f$TRB_aa_1,'MANAhi','MANAlo')

TCR_ser <- CreateSeuratObject(counts = res,meta.data = tcr_f)
TCR_ser<-NormalizeData(TCR_ser, verbose = FALSE)
TCR_ser <- ScaleData(TCR_ser, verbose = FALSE) 
saveRDS(TCR_ser,'MANA_6072_clones_pseudobulk.rds')
##
TCRser <- readRDS('D:/Project/ManaScore/01.Result/36.12_single_moldes/20.pseudobulk_dge/MANA_6072_clones_pseudobulk.rds')
tcr <- NSCLC_ser[[]] %>% filter(is.na(TRB_aa_1)==FALSE) %>%
  filter(TCRType_new %in% c('pTRC','non-pTRC') & tissue=='tumor') %>%
  group_by(imid,TCRType_new,TRB_aa_1) %>%
  summarise(n=n()) %>%
  mutate(tcr = paste0(imid,':',TRB_aa_1))  %>% as.data.frame() ## 18818
rownames(tcr) <- tcr$tcr
##
TCRser$TCRType_new <- tcr[rownames(TCRser[[]]),]$TCRType_new


TCRser_cutoff5 <- subset(TCRser, subset=n>=5)

Idents(TCRser_cutoff5) <- 'TCRType_new'
dge1 <- FindMarkers(TCRser_cutoff5, ident.1 = "pTRC"  , 
                    ident.2 =   "non-pTRC",logfc.threshold=0.1)

dge1_select <- dge1 %>% arrange(-avg_log2FC) %>% filter(p_val_adj<0.01) %>% head(20) %>%
  rbind(dge1 %>% arrange(-avg_log2FC) %>% filter(p_val_adj<0.01) %>% tail(20))
saveRDS(dge1,'MANAhi_vs_MANAlo_dge_cutoff5.rds')
##
dge1$Regulation = 'not. sig' 
dge1 <- dge1 %>% mutate(Regulation = ifelse(p_val_adj< 0.05 & avg_log2FC> 0, 'upregulated in pTRC', Regulation)) %>%
  mutate(Regulation = ifelse(p_val_adj< 0.05 & avg_log2FC< 0, 'downregulated in pTRC', Regulation)) 
write.csv(dge1,'dge_pseudobulk_pTRC_vs_non-pTRC.csv',quote=FALSE, row.names = TRUE)

##
side_MANA<- TCRser_cutoff5@meta.data[,c('TCRType_new','imid')] %>% arrange(TCRType_new)
colnames(side_MANA) <- c('TCRType','patient')
annot_colors1=list(TCRType = c(pTRC='#B2182B',
                               `non-pTRC`='#2166AC'),
                   patient=c(`MD01-005`="#D95F02",
                             `MD01-010`="#A6A810",
                             `MD043-003` = '#1B9E77',
                             `MD043-008` = '#E7298A',
                             `NY016-022` = '#A7675A',
                             `NY016-025` = '#66A61E',
                             `MD01-004`="#C5900F",
                             `MD01-019`="#7570B3",
                             `MD01-024` = '#666666',
                             `MD043-006` = '#A6761D',
                             `MD043-011` = '#A66753',
                             `NY016-007` = '#AD4C9E',
                             `NY016-014` ='#7A7E3C',
                             `NY016-015` = '#E6AB02',
                             `NY016-021` = '#866E41')
)
plot_m1 <- TCRser_cutoff5@assays$RNA@scale.data[rownames(dge1_select),rownames(side_MANA)] %>% as.matrix()
library(pheatmap)
tiff("MANAhi_vs_MANAlo_top20_cutoff5.tiff", width =6.5, height = 7, units = 'in', res = 300)
pheatmap(MinMax(plot_m1,min=-2.5,max = 2.5),cluster_rows = F,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols = F,
         # scale = 'row',
         show_colnames =F,
         annotation_col = side_MANA[,c('patient','TCRType')],
         # cutree_rows = 2,
         annotation_colors=annot_colors1
)
dev.off()


## pTRC R vs NR
## dge between pTRC in R and NR
pTRC_TCR <- subset(TCRser_cutoff5, subset = TCRType_new=='pTRC')
Idents(pTRC_TCR) <- 'response'
dge2 <- FindMarkers(pTRC_TCR, ident.1 = "R"  , 
                    ident.2 =   "NR",logfc.threshold=0.1)
dge2 %>% arrange(-avg_log2FC) %>% filter(p_val_adj<0.01) %>% head(10)

dge2_select <- dge2 %>% arrange(-avg_log2FC) %>% filter(p_val_adj<0.01) %>% head(20) %>%
  rbind(dge2 %>% arrange(-avg_log2FC) %>% filter(p_val_adj<0.01) %>% tail(20))

library(viridis)

DoHeatmap(pTRC_TCR,features = rownames(dge2_select),group.colors = c('#00468BFF','#ED0000FF'))+
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  guides(color='none')
ggsave('MANAhi_R_vs_NR_clones_top20_cutoff5.jpeg',height=6,width=5.5)  
saveRDS(dge2,"MANAhi_R_vs_NR_dge_cutoff5.rds")

dge2$Regulation = 'not. sig'
dge2 <- dge2 %>% mutate(Regulation = ifelse(p_val_adj< 0.05 & avg_log2FC> 0, 'upregulated in R', Regulation)) %>%
  mutate(Regulation = ifelse(p_val_adj< 0.05 & avg_log2FC< 0, 'downregulated in R', Regulation)) 
write.csv(dge2,'dge_pseudobulk_pTRC_R_vs_NR.csv',quote=FALSE, row.names = TRUE)


#dge2 <- readRDS('MANAhi_R_vs_NR_dge_cutoff5.rds')
##
library(pheatmap)
side_MANAhi<- pTRC_TCR[[]][,c('response','imid')] %>% arrange(response)
colnames(side_MANAhi) <- c('response','patient')
annot_colors=list(response = c(R='#EE0000FF',
                               NR='#3B4992FF'),
                  patient=c(`MD01-005`="#D95F02",
                            `MD01-010`="#A6A810",
                            `MD043-003` = '#1B9E77',
                            `MD043-008` = '#E7298A',
                            `NY016-022` = '#A7675A',
                            `NY016-025` = '#66A61E',
                            `MD01-004`="#C5900F",
                            `MD01-019`="#7570B3",
                            `MD01-024` = '#666666',
                            `MD043-006` = '#A6761D',
                            `MD043-011` = '#A66753',
                            `NY016-007` = '#AD4C9E',
                            `NY016-014` ='#7A7E3C',
                            `NY016-015` = '#E6AB02',
                            `NY016-021` = '#866E41')
)
plot_m <- pTRC_TCR@assays$RNA@scale.data[rownames(dge2_select),rownames(side_MANAhi)] %>% as.matrix()

tiff("pTRC_R_vs_NR_top20.tiff", width =6.5, height = 7, units = 'in', res = 300)
pheatmap(MinMax(plot_m,min=-2.5,max = 2.5),cluster_rows = F,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols = F,
         # scale = 'row',
         show_colnames =F,
         annotation_col = side_MANAhi[,c('patient','response')],
         # cutree_rows = 2,
         annotation_colors=annot_colors
)
dev.off()

###
dge2$Gene = rownames(dge2)
dd <- dge2 %>% arrange(-avg_log2FC) %>% mutate(rank=1:length(dge2[,1])) 

gene_interest = c('GZMK','GIMAP4','GIMAP5','GIMPA7','CXCR4','DUSP2','TRBC1','NFATC2','PLEK','GPR183','EOMES',
                  'GIMAP1','KDSR','TAGAP','FOS','TNFSF9','FOSB','KLRD1','LAYN','CXCL13','ENTPD1','RBPJ','KLRC2','SNX9',
                  'ID2','GALNT2','ZNF683','GZMB','HLA-DRA','IL7R','KLF2','TXNIP','JUNB','ZFP36','PIK3R1',
                  'ANXA1','FOS','FOSB','CTLA4','ITGAE','GPR25','GNLY','CCL3','LAG3','HAVCR2','GZMA','PRF1','PDCD1','TIGIT','RUNX2',
                  'CXCR6','RGS1','IFNG','EOMES','ZFP36','KLF2','CCR7','TOX') %>% unique()

highlight = dd %>% filter(Gene %in% gene_interest)
geneType <- data.frame(Gene = c('CXCL13','IL7R','CCR7','GZMK','PDCD1','PIK3R1','SNX9' ,'EOMES','CXCR6','TXNIP',
                                'GZMB','CCL3','LAYN','GIMAP5','ITGAE','KLRC2','GIMAP1','TIGIT','GIMAP4',
                                'LAG3','ZNF683','IFNG','GPR25','GZMA','GNLY','TOX'),
                       type = c('tumor reactive','Renewal/memory','Renewal/memory','Effector','chekpoint','TCR signaling',
                                'TCR signaling','Renewal/memory','TRM','Renewal/memory','Effector','TRM',
                                'Effector','TCR signaling','TRM','Effector','TCR signaling','chekpoint','TCR signaling',
                                'chekpoint','TRM','Effector','TRM','Effector','Effector','Exhaustion'))

highlight <- highlight %>% left_join(geneType) %>% na.omit()

ggplot() +
  geom_point(data = dd,
             aes(x=rank,y=avg_log2FC,color=p_val_adj)) +
  theme_classic() +
  scale_color_gradient(low="red", high="yellow") +
  geom_label_repel(data= highlight,
                   aes(x=rank,y=avg_log2FC,label = Gene,fill=type),
                   segment.colour = "black",
                   fontface = 'bold',
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"))+
  geom_hline(yintercept=0,linetype='dashed')+
  scale_fill_aaas(alpha=0.6)+
  theme(axis.title = element_text(size=14,color='black'),
        axis.text = element_text(size=13,color='black'))

ggsave('pTRC_R_vs_NR_highlight_cutoff5.tiff',height=5.5,width=5.5)

##


