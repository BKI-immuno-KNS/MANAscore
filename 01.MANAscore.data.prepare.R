## MANAscore data prepare
## zzeng
library(magrittr)
library(Seurat)
library(preprocessCore)
library(dplyr)

setwd('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/')
saver<-readRDS("D:/TCR/102018-Neountx/Data/Processed/Scseq/saver.lung.melanoma.rds") ## impute data
umap<-readRDS("D:/Project/ManaScore/01.Result/00.TCR_info/00.umap_new.rds")
ser <- readRDS("D:/Project/ManaScore/01.Result/00.Lung_melanoma_integration/ser_lung_me_integrated_fi.rds")

# training data and validation data 

pa = c('p2','p11','p15', 'IS2','PP3')
me_known <- read.csv("D:/Project/ManaScore/01.Result/01.MANAscore_model_data_preparation/01.me_train_me_trainging_cell_label_gene_f_expression.csv",
                     check.names = FALSE,row.names = 1)

NSCLC_known = read.csv("D:/Project/ManaScore/01.Result/01.MANAscore_model_data_preparation/02.NSCLC_train_NSCLC_training_cell_label_gene_f_expression.csv",
                       check.names = FALSE,row.names = 1)

known <- rbind(me_known,NSCLC_known)
p_info <- known[,c("patient","Label")]
p_info$barcode = rownames(p_info)
###
gene = c('CXCL13','ENTPD1','IL7R')
## imputation
for(p in pa){
  info = p_info %>% filter(patient==p)
  patient_tumor_barcode <- umap %>% filter(patient_id ==p) %>% 
    rownames() 
  
  if(length(patient_tumor_barcode)>0){
    saver_sub <- saver[,patient_tumor_barcode]
    
    ## quantile normalization
    saver_sub_norm <- normalize.quantiles(saver_sub, copy = TRUE)
    colnames(saver_sub_norm) <- colnames(saver_sub)
    rownames(saver_sub_norm) <- rownames(saver_sub)
    ## scale
    saver_sub_norm_scale <- scale(saver_sub_norm)

    subdat_known = saver_sub_norm_scale[gene,
                                        intersect(patient_tumor_barcode,
                                                  rownames(known))] %>% 
      t() %>% as.data.frame() %>%
      merge(info %>% select(Label),by='row.names')
    
    write.csv(subdat_known,paste0(p,"_known.csv"))

    antigen_NA_barcode = umap %>% filter(tissue=="tumor" &
                                           patient_id==p & 
                                           antigen=="antigen.NA") %>% .[,"barcode"] 
    
    subdat_unknown = saver_sub_norm_scale[gene,
                                          intersect(colnames(saver_sub_norm_scale),
                                                    antigen_NA_barcode)] %>% 
      t() %>% as.data.frame()
    
    write.csv(subdat_known,paste0(p,"_unknown.csv"))
    
  }
}

## non-imputation
for(p in pa){
  info = p_info %>% filter(patient==p)
  patient_tumor_barcode <-  ser[[]] %>% filter(patient==p) %>% rownames()
  
  if(length(patient_tumor_barcode)>0){
    subdat <- ser@assays$RNA[,patient_tumor_barcode]
    
    subdat_known <- subdat[gene,intersect(patient_tumor_barcode,
                                          rownames(me_known))] %>% 
      t() %>% as.data.frame() %>%
      merge(info %>% select(Label),by='row.names')
  
    write.csv(sub_known,paste0(p,"_known_RNA.cvs"))
    
    antigen_NA_barcode = umap %>% filter(tissue=="tumor" &
                                           patient_id==p & 
                                           antigen=="antigen.NA") %>% .[,"barcode"] 
    
    
    subdat_unknown <- subdat[gene,intersect(patient_tumor_barcode,
                                            antigen_NA_barcode)] %>% 
      t() %>% as.data.frame()
    write.csv(sub_known,paste0(p,"_unknown_RNA.cvs"))
  }
}

## unknown for other NSCLC patient #######
NSCLC_other_pa <- umap %>% filter(tissue=='tumor' & cancer=='NSCLC') %>%
  .[,'patient_id'] %>% unique() %>% setdiff(c('IS2','PP3'))

## imputation
for(p in NSCLC_other_pa){
  patient_tumor_barcode <- umap %>% filter(patient_id ==p) %>% 
    rownames() 
  
  if(length(patient_tumor_barcode)>0){
    saver_sub <- saver[,patient_tumor_barcode]
    
    ## quantile normalization
    saver_sub_norm <- normalize.quantiles(saver_sub, copy = TRUE)
    colnames(saver_sub_norm) <- colnames(saver_sub)
    rownames(saver_sub_norm) <- rownames(saver_sub)
    ## scale
    saver_sub_norm_scale <- scale(saver_sub_norm)
    
    antigen_NA_barcode = umap %>% filter(tissue=="tumor" &
                                           patient_id==p & 
                                           antigen=="antigen.NA") %>% .[,"barcode"] 
    
    subdat_unknown = saver_sub_norm_scale[gene,
                                          intersect(patient_tumor_barcode,
                                                    antigen_NA_barcode)] %>% 
      t() %>% as.data.frame()
    
    write.csv(subdat_known,paste0(p,"_unknown.csv"))
  }
}
## non-imputation
for(p in NSCLC_other_pa){
  info = p_info %>% filter(patient==p)
  patient_tumor_barcode <-  ser[[]] %>% filter(patient==p) %>% rownames()
  
  if(length(patient_tumor_barcode)>0){
    subdat <- ser@assays$RNA[,patient_tumor_barcode]
    
    antigen_NA_barcode = umap %>% filter(tissue=="tumor" &
                                           patient_id==p & 
                                           antigen=="antigen.NA") %>% .[,"barcode"] 
    
    subdat_unknown <- subdat[gene,intersect(patient_tumor_barcode,
                                            antigen_NA_barcode)] %>% 
      t() %>% as.data.frame()
    write.csv(sub_known,paste0(p,"_unknown_RNA.cvs"))
  }
}

#####

## NSCLC data for normal and LN samples #######
setwd('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/NSCLC_tissue/')

d <- umap %>% filter(cancer=='NSCLC' & tissue %in% c('normal','LN')) %>%
  group_by(imid,tissue) %>% summarise(n=n()) %>% as.data.frame()

###

for(i in 1:length(d[,1])){
  print(i)
  p = d[i,]$imid
  t = d[i,]$tissue
  
  patient_tissue_barcode <- umap %>% filter(imid ==p & tissue==t) %>% 
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
    write.csv(subdat,paste0(p,"_",t,"_unknown.csv"))
  }
}

## non-imputation
cell = umap %>% filter(cancer=='NSCLC' & tissue %in% c('normal','LN')) 
df = ser@assays$RNA@data[gene,cell$barcode] %>% as.matrix() %>% t() %>% as.data.frame() 

write.csv(df,"NSCLC_normal_LN_RNA.csv")


