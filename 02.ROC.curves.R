##ROC curves - MANAscore using 12 models
## zzeng
library(ROCR)
library(pROC)
library(ggplot2)
library(magrittr)
library(tidyverse)

ssgsea <- readRDS('D:/Project/ManaScore/01.Result/35.downstream_analysis/04.signature_comparison/ssgsea.rds')
res = ssgsea %>% t() %>% as.data.frame()
res$barcode = rownames(res)
##
## NeoTCR8
setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/01.model_statistics/')
d1 = read.csv('D:/Project/ManaScore/01.Result/34.combine_model_trained_on_me/voting_6_models_i_test_score.csv',row.name=1)
d11 = d1 %>% as.data.frame() %>% left_join(res)  %>% dplyr::select(barcode,label,NeoTCR8)

d11_pred <- prediction(d11$NeoTCR8,d11$label)

d11_perf <- performance(d11_pred, "tpr", "fpr" )

d11_cost_perf = performance(d11_pred, "cost") 
d11_pred@cutoffs[[1]][which.min(d11_cost_perf@y.values[[1]])]

## CXCL13_i
d12 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CXCL13_voting_6_models_i_test_score.csv',row.names = 1)
colnames(d12)[3] <- 'CXCL13_i'
d12_pred <- prediction(d12$CXCL13_i,d12$label)
d12_perf <- performance(d12_pred, "tpr", "fpr" )
d12_cost_perf = performance(d12_pred, "cost") 
d12_pred@cutoffs[[1]][which.min(d12_cost_perf@y.values[[1]])]

d13 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CXCL13_voting_6_models_ni_test_score.csv',row.names = 1)
colnames(d13)[3] <- 'CXCL13_ni'
d13_pred <- prediction(d13$CXCL13_ni,d12$label)
d13_perf <- performance(d13_pred, "tpr", "fpr" )
d13_cost_perf = performance(d13_pred, "cost") 
d13_pred@cutoffs[[1]][which.min(d13_cost_perf@y.values[[1]])]

d14 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/voting_models_i_test_score.csv',row.names = 1)
colnames(d14)[3] <- 'MANAscore_i'
d14_pred <- prediction(d14$MANAscore_i,d14$label)
d14_perf <- performance(d14_pred, "tpr", "fpr" )
d14_cost_perf = performance(d14_pred, "cost") 
d14_pred@cutoffs[[1]][which.min(d14_cost_perf@y.values[[1]])] # 0.4257571
#
d15 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/voting_models_ni_test_score.csv',row.names = 1)
colnames(d15)[3] <- 'MANAscore_ni'
d15_pred <- prediction(d15$MANAscore_ni,d15$label)
d15_perf <- performance(d15_pred, "tpr", "fpr" )
d15_cost_perf = performance(d15_pred, "cost") 
d15_pred@cutoffs[[1]][which.min(d15_cost_perf@y.values[[1]])] # 0.4728286

auc(d11$label,d11$NeoTCR8) ## 0.5826
auc(d12$label,d12$CXCL13_i) # 0.8278
auc(d13$label,d13$CXCL13_ni) # 0.7266
auc(d14$label,d14$MANAscore_i) ## 0.891
auc(d15$label,d15$MANAscore_ni) ## 0.803
## plot
tiff("00.all_5_model_20per_MANA_nonMANA_in_me.tiff",width=5,height=4,units = "in",res=300)
par(mar = (c(4,4.5,1,1)))
plot(d14_perf, colorize = FALSE,col="firebrick2",lwd=2.5,cex.lab = 1.5)
plot(d15_perf, add = TRUE,colorize = FALSE,col="royalblue3",lwd=2)

plot(d11_perf, add = TRUE,colorize = FALSE,col="seagreen",lwd=2)
plot(d12_perf, add = TRUE,colorize = FALSE,col="darkorange",lwd=2)
plot(d13_perf,add = TRUE, colorize = FALSE,col="darkorange3",lwd=2)

legend(x=0.48,y=0.4, c('MANAscore_i (0.89)','MANAscore_ni (0.80)',
                       "NeoTCR8 (0.58)",'CXCL13_i (0.83)',
                       'CXCL13_ni (0.73)'), lty=1,lwd = 2, 
       col = c("firebrick2","royalblue3","seagreen",'darkorange','darkorange3'), 
       bty="n", inset=c(0,-0.15),cex=0.9,text.font=1,title = 'AUC score')
abline(0,1,lty='dashed')
dev.off()

##
d2 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/MAA_voting_i_test_score.csv',row.name=1)
d21 = d2 %>% as.data.frame() %>% left_join(res)  %>% dplyr::select(barcode,label,NeoTCR8)
d21_pred <- prediction(d21$NeoTCR8,d21$label)
d21_perf <- performance(d21_pred, "tpr", "fpr" )
d21_cost_perf = performance(d21_pred, "cost") 
d21_pred@cutoffs[[1]][which.min(d21_cost_perf@y.values[[1]])]

d22 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CXCL13_MAA_voting_i_test_score.csv',row.names = 1)
colnames(d22)[3] <- 'CXCL13_i'
d22_pred <- prediction(d22$CXCL13_i,d22$label)
d22_perf <- performance(d22_pred, "tpr", "fpr" )
d22_cost_perf = performance(d22_pred, "cost") 
d22_pred@cutoffs[[1]][which.min(d22_cost_perf@y.values[[1]])]

d23 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CXCL13_MAA_voting_ni_test_score.csv',row.names = 1)
colnames(d23)[3] <- 'CXCL13_ni'
d23_pred <- prediction(d23$CXCL13_ni,d23$label)
d23_perf <- performance(d23_pred, "tpr", "fpr" )
d23_cost_perf = performance(d23_pred, "cost") 
d23_pred@cutoffs[[1]][which.min(d23_cost_perf@y.values[[1]])]

d24 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/MAA_voting_i_test_score.csv',row.names = 1)
colnames(d24)[3] <- 'MANAscore_i'
d24_pred <- prediction(d24$MANAscore_i,d24$label)
d24_perf <- performance(d24_pred, "tpr", "fpr" )
d24_cost_perf = performance(d24_pred, "cost") 
d24_pred@cutoffs[[1]][which.min(d24_cost_perf@y.values[[1]])] #0.1314989

d25 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/MAA_voting_ni_test_score.csv',row.names = 1)
colnames(d25)[3] <- 'MANAscore_ni'
d25_pred <- prediction(d25$MANAscore_ni,d25$label)
d25_perf <- performance(d25_pred, "tpr", "fpr" )
d25_cost_perf = performance(d25_pred, "cost") 
d25_pred@cutoffs[[1]][which.min(d25_cost_perf@y.values[[1]])] # 0.3024408


auc(d21$label,d21$NeoTCR8) ## 0.7758
auc(d22$label,d22$CXCL13_i) # 0.8659
auc(d23$label,d23$CXCL13_ni) # 0.6987
auc(d24$label,d24$MANAscore_i) # 0.8531
auc(d25$label,d25$MANAscore_ni) # 0.8044


tiff("00.all_5_model_MAA_nonMANA_in_p2.tiff",width=5,height=4,units = "in",res=300)
par(mar=c(4,4.5,1,1))
plot(d24_perf, colorize = FALSE,col="firebrick2",lwd=2.5,cex.lab = 1.5)
plot(d25_perf, add = TRUE,colorize = FALSE,col="royalblue3",lwd=2)

plot(d21_perf, add = TRUE,colorize = FALSE,col="seagreen",lwd=2)
plot(d22_perf, add = TRUE,colorize = FALSE,col="darkorange",lwd=2)
plot(d23_perf,add = TRUE, colorize = FALSE,col="darkorange3",lwd=2)

legend(x=0.48,y=0.4, c('MANAscore_i (0.85)','MANAscore_ni (0.80)',
                       "NeoTCR8 (0.78)",'CXCL13_i (0.87)',
                       'CXCL13_ni (0.70)'), lty=1,lwd = 2, 
       col = c("firebrick2","royalblue3","seagreen",'darkorange','darkorange3'), 
       bty="n", inset=c(0,-0.15),cex=0.9,text.font=1,title = 'AUC score')
abline(0,1,lty='dashed')
dev.off()

## IS2
d3 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/IS2_voting_i_test_score.csv',row.name=1)
d31 = d3 %>% as.data.frame() %>% left_join(res)  %>% dplyr::select(barcode,label,NeoTCR8)
d31_pred <- prediction(d31$NeoTCR8,d31$label)
d31_perf <- performance(d31_pred, "tpr", "fpr" )
d31_cost_perf = performance(d31_pred, "cost") 
d31_pred@cutoffs[[1]][which.min(d31_cost_perf@y.values[[1]])]

d32 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CXCL13_IS2_voting_i_test_score.csv',row.names = 1)
colnames(d32)[3] <- 'CXCL13_i'
d32_pred <- prediction(d32$CXCL13_i,d32$label)
d32_perf <- performance(d32_pred, "tpr", "fpr" )
d32_cost_perf = performance(d32_pred, "cost") 
d32_pred@cutoffs[[1]][which.min(d32_cost_perf@y.values[[1]])]

d33 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CXCL13_IS2_voting_ni_test_score.csv',row.names = 1)
colnames(d33)[3] <- 'CXCL13_ni'
d33_pred <- prediction(d33$CXCL13_ni,d33$label)
d33_perf <- performance(d33_pred, "tpr", "fpr" )
d33_cost_perf = performance(d33_pred, "cost") 
d33_pred@cutoffs[[1]][which.min(d33_cost_perf@y.values[[1]])]

d34 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/IS2_voting_i_test_score.csv',row.names = 1)
colnames(d34)[3] <- 'MANAscore_i'
d34_pred <- prediction(d34$MANAscore_i,d34$label)
d34_perf <- performance(d34_pred, "tpr", "fpr" )
d34_cost_perf = performance(d34_pred, "cost") 
d34_pred@cutoffs[[1]][which.min(d34_cost_perf@y.values[[1]])] ## 0.04055744

d35 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/IS2_voting_ni_test_score.csv',row.names = 1)
colnames(d35)[3] <- 'MANAscore_ni'
d35_pred <- prediction(d35$MANAscore_ni,d35$label)
d35_perf <- performance(d35_pred, "tpr", "fpr" )
d35_cost_perf = performance(d35_pred, "cost") 
d35_pred@cutoffs[[1]][which.min(d35_cost_perf@y.values[[1]])] ## 0.442155


auc(d31$label,d31$NeoTCR8) ## 0.8611
auc(d32$label,d32$CXCL13_i) ## 0.9475
auc(d33$label,d33$CXCL13_ni) ## 0.8819
auc(d34$label,d34$MANAscore_i) # 0.966
auc(d35$label,d35$MANAscore_ni) # 0.9703

tiff("00.all_5_model_IS2.tiff",width=5,height=4,units = "in",res=300)
par(mar=c(4,4.5,1,1))
plot(d34_perf, colorize = FALSE,col="firebrick2",lwd=2.5,cex.lab = 1.5)
plot(d35_perf, add = TRUE,colorize = FALSE,col="royalblue3",lwd=2)

plot(d31_perf, add = TRUE,colorize = FALSE,col="seagreen",lwd=2)
plot(d32_perf, add = TRUE,colorize = FALSE,col="darkorange",lwd=2)
plot(d33_perf,add = TRUE, colorize = FALSE,col="darkorange3",lwd=2)

legend(x=0.4,y=0.45, c('MANAscore_i (0.97)','MANAscore_ni (0.97)',
                       "NeoTCR8 (0.86)",'CXCL13_i (0.95)',
                       'CXCL13_ni (0.88)'), lty=1,lwd = 2, 
       col = c("firebrick2","royalblue3","seagreen",'darkorange','darkorange3'), 
       bty="n", inset=c(0,-0.15),cex=1,text.font=1,title = 'AUC score')
abline(0,1,lty='dashed')
dev.off()

##
d4 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/PP3_voting_i_test_score.csv',row.name=1)
d41 = d4 %>% as.data.frame() %>% left_join(res)  %>% dplyr::select(barcode,label,NeoTCR8)
d41_pred <- prediction(d41$NeoTCR8,d41$label)
d41_perf <- performance(d41_pred, "tpr", "fpr" )
d41_cost_perf = performance(d41_pred, "cost") 
d41_pred@cutoffs[[1]][which.min(d41_cost_perf@y.values[[1]])]

d42 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CXCL13_PP3_voting_i_test_score.csv',row.names = 1)
colnames(d42)[3] <- 'CXCL13_i'
d42_pred <- prediction(d42$CXCL13_i,d42$label)
d42_perf <- performance(d42_pred, "tpr", "fpr" )
d42_cost_perf = performance(d42_pred, "cost") 
d42_pred@cutoffs[[1]][which.min(d42_cost_perf@y.values[[1]])]

d43 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CXCL13_PP3_voting_ni_test_score.csv',row.names = 1)
colnames(d43)[3] <- 'CXCL13_ni'
d43_pred <- prediction(d43$CXCL13_ni,d43$label)
d43_perf <- performance(d43_pred, "tpr", "fpr" )
d43_cost_perf = performance(d43_pred, "cost") 
d43_pred@cutoffs[[1]][which.min(d43_cost_perf@y.values[[1]])]


d44 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/PP3_voting_i_test_score.csv',row.names = 1)
colnames(d44)[3] <- 'MANAscore_i'
d44_pred <- prediction(d44$MANAscore_i,d44$label)
d44_perf <- performance(d44_pred, "tpr", "fpr" )
d44_cost_perf = performance(d44_pred, "cost") 
d44_pred@cutoffs[[1]][which.min(d44_cost_perf@y.values[[1]])] ## 0.04170888

d45 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/PP3_voting_ni_test_score.csv',row.names = 1)
colnames(d45)[3] <- 'MANAscore_ni'
d45_pred <- prediction(d45$MANAscore_ni,d45$label)
d45_perf <- performance(d45_pred, "tpr", "fpr" )
d45_cost_perf = performance(d45_pred, "cost") 
d45_pred@cutoffs[[1]][which.min(d45_cost_perf@y.values[[1]])] # 0.5225799

auc(d44$label,d44$MANAscore_i) ## 0.9564
auc(d45$label,d45$MANAscore_ni) ## 0.9429
auc(d41$label,d41$NeoTCR8) ## 0.8293
auc(d42$label,d42$CXCL13_i) ## 0.77
auc(d43$label,d43$CXCL13_ni) ## 0.614


tiff("00.all_5_model_PP3.tiff",width=5,height=4,units = "in",res=300)
par(mar=c(4,4.5,1,1))
plot(d44_perf, colorize = FALSE,col="firebrick2",lwd=2.5,cex.lab = 1.5)
plot(d45_perf, add = TRUE,colorize = FALSE,col="royalblue3",lwd=2)

plot(d41_perf, add = TRUE,colorize = FALSE,col="seagreen",lwd=2)
plot(d42_perf, add = TRUE,colorize = FALSE,col="darkorange",lwd=2)
plot(d43_perf,add = TRUE, colorize = FALSE,col="darkorange3",lwd=2)

legend(x=0.4,y=0.45, c('MANAscore_i (0.96)','MANAscore_ni (0.94)',
                       "NeoTCR8 (0.83)",'CXCL13_i (0.77)',
                       'CXCL13_ni (0.61)'), lty=1,lwd = 2, 
       col = c("firebrick2","royalblue3","seagreen",'darkorange','darkorange3'), 
       bty="n", inset=c(0,-0.15),cex=1,text.font=1,title = 'AUC score')
abline(0,1,lty='dashed')
dev.off()

#####

#####################
#####
### ###### add single gene model 
# 1: NeoTCR8
# 2: CXCL13_i
# 3: CXCL13_ni
# 4: MANAscore_i
# 5: MANAscore_ni
# 6: ENTPD1_i
# 7: ENTPD1_ni
# 8: IL7R_i
# 9: IL7R_ni


setwd('D:/Project/ManaScore/01.Result/36.12_single_moldes/01.model_statistics/')


#
d16 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CD39_voting_6_models_i_test_score.csv',row.names = 1)
colnames(d16)[3] <- 'ENTPD1_i'
d16_pred <- prediction(d16$ENTPD1_i,d16$label)
d16_perf <- performance(d16_pred, "tpr", "fpr" )

d17 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CD39_voting_6_models_ni_test_score.csv',row.names = 1)
colnames(d17)[3] <- 'ENTPD1_ni'
d17_pred <- prediction(d17$ENTPD1_ni,d17$label)
d17_perf <- performance(d17_pred, "tpr", "fpr" )

d18 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/IL7R_voting_6_models_i_test_score.csv',row.names = 1)
colnames(d18)[3] <- 'IL7R_i'
d18_pred <- prediction(d18$IL7R_i,d18$label)
d18_perf <- performance(d18_pred, "tpr", "fpr" )

d19 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/IL7R_voting_6_models_ni_test_score.csv',row.names = 1)
colnames(d19)[3] <- 'IL7R_ni'
d19_pred <- prediction(d19$IL7R_ni,d19$label)
d19_perf <- performance(d19_pred, "tpr", "fpr" )

auc(d11$label,d11$NeoTCR8) ## 0.5826
auc(d12$label,d12$CXCL13_i) # 0.8278
auc(d13$label,d13$CXCL13_ni) # 0.7266
auc(d14$label,d14$MANAscore_i) ## 0.891
auc(d15$label,d15$MANAscore_ni) ## 0.803
auc(d16$label,d16$ENTPD1_i) ## 0.8584
auc(d17$label,d17$ENTPD1_ni) # 0.5857
auc(d18$label,d18$IL7R_i) # 0.799
auc(d19$label,d19$IL7R_ni) ## 0.6854

## plot
tiff("02.all_8_model_20per_MANA_nonMANA_in_me.tiff",width=6.5,height=5.5,units = "in",res=300)
par(mar = (c(4,4.5,1,1)))
plot(d14_perf, colorize = FALSE,col="firebrick2",lwd=2.5,cex.lab = 1.5)
plot(d15_perf, add = TRUE,colorize = FALSE,col="royalblue3",lwd=2,lty='dashed')

plot(d12_perf, add = TRUE,colorize = FALSE,col="darkorange",lwd=2)
plot(d13_perf,add = TRUE, colorize = FALSE,col="darkorange3",lwd=2,lty='dashed')

plot(d16_perf, add = TRUE,colorize = FALSE,col="cyan3",lwd=2)
plot(d17_perf,add = TRUE, colorize = FALSE,col="darkcyan",lwd=2,lty='dashed')

plot(d18_perf, add = TRUE,colorize = FALSE,col="darkorchid",lwd=2)
plot(d19_perf,add = TRUE, colorize = FALSE,col="darkorchid4",lwd=2,lty='dashed')

legend(x=0.55,y=0.4, c('MANAscore_i (0.89)','MANAscore_ni (0.80)',
                       'CXCL13_i (0.83)','CXCL13_ni (0.73)',
                       'ENTPD1_i (0.86)','ENTPD1_ni (0.58)',
                       'IL7R_i (0.80)','IL7R_ni (0.69)'), 
       lty = c(1,2,1,2,1,2,1,2),lwd = 2, 
       col = c("firebrick2","royalblue3","darkorange",'darkorange3','cyan3','darkcyan',
               'darkorchid','darkorchid4'), 
       bty="n", inset=c(0,-0.15),cex=0.9,text.font=1,title = 'AUC score')
abline(0,1,lty='dashed')
dev.off()

##
d26 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CD39_MAA_voting_i_test_score.csv',row.names = 1)
colnames(d26)[3] <- 'ENTPD1_i'
d26_pred <- prediction(d26$ENTPD1_i,d26$label)
d26_perf <- performance(d26_pred, "tpr", "fpr" )

d27 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CD39_MAA_voting_ni_test_score.csv',row.names = 1)
colnames(d27)[3] <- 'ENTPD1_ni'
d27_pred <- prediction(d27$ENTPD1_ni,d27$label)
d27_perf <- performance(d27_pred, "tpr", "fpr" )

d28 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/IL7R_MAA_voting_i_test_score.csv',row.names = 1)
colnames(d28)[3] <- 'IL7R_i'
d28_pred <- prediction(d28$IL7R_i,d28$label)
d28_perf <- performance(d28_pred, "tpr", "fpr" )

d29 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/IL7R_MAA_voting_ni_test_score.csv',row.names = 1)
colnames(d29)[3] <- 'IL7R_ni'
d29_pred <- prediction(d29$IL7R_ni,d29$label)
d29_perf <- performance(d29_pred, "tpr", "fpr" )


auc(d21$label,d21$NeoTCR8) ## 0.7758
auc(d22$label,d22$CXCL13_i) # 0.8659
auc(d23$label,d23$CXCL13_ni) # 0.6987
auc(d24$label,d24$MANAscore_i) ## 0.8531
auc(d25$label,d25$MANAscore_ni) ## 0.8044
auc(d26$label,d26$ENTPD1_i) ## 0.8193
auc(d27$label,d27$ENTPD1_ni) # 0.6176
auc(d28$label,d28$IL7R_i) # 0.7882
auc(d29$label,d29$IL7R_ni) ## 0.6909

## plot
tiff("02.all_8_model_MAA_20per_nonMANA_in_p2.tiff",width=6.5,height=5.5,units = "in",res=300)
par(mar = (c(4,4.5,1,1)))
plot(d24_perf, colorize = FALSE,col="firebrick2",lwd=2.5,cex.lab = 1.5)
plot(d25_perf, add = TRUE,colorize = FALSE,col="royalblue3",lwd=2,lty='dashed')

plot(d22_perf, add = TRUE,colorize = FALSE,col="darkorange",lwd=2)
plot(d23_perf,add = TRUE, colorize = FALSE,col="darkorange3",lwd=2,lty='dashed')

plot(d26_perf, add = TRUE,colorize = FALSE,col="cyan3",lwd=2)
plot(d27_perf,add = TRUE, colorize = FALSE,col="darkcyan",lwd=2,lty='dashed')

plot(d28_perf, add = TRUE,colorize = FALSE,col="darkorchid",lwd=2)
plot(d29_perf,add = TRUE, colorize = FALSE,col="darkorchid4",lwd=2,lty='dashed')


legend(x=0.55,y=0.4, c('MANAscore_i (0.85)','MANAscore_ni (0.80)',
                       'CXCL13_i (0.87)','CXCL13_ni (0.70)',
                       'ENTPD1_i (0.86)','ENTPD1_ni (0.58)',
                       'IL7R_i (0.80)','IL7R_ni (0.69)'), 
       lty = c(1,2,1,2,1,2,1,2),lwd = 2, 
       col = c("firebrick2","royalblue3","darkorange",'darkorange3','cyan3','darkcyan',
               'darkorchid','darkorchid4'), 
       bty="n", inset=c(0,-0.15),cex=0.9,text.font=1,title = 'AUC score')
abline(0,1,lty='dashed')
dev.off()

## IS2
##
d36 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CD39_IS2_voting_i_test_score.csv',row.names = 1)
colnames(d36)[3] <- 'ENTPD1_i'
d36_pred <- prediction(d36$ENTPD1_i,d36$label)
d36_perf <- performance(d36_pred, "tpr", "fpr" )

d37 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CD39_IS2_voting_ni_test_score.csv',row.names = 1)
colnames(d37)[3] <- 'ENTPD1_ni'
d37_pred <- prediction(d37$ENTPD1_ni,d37$label)
d37_perf <- performance(d37_pred, "tpr", "fpr" )

d38 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/IL7R_IS2_voting_i_test_score.csv',row.names = 1)
colnames(d38)[3] <- 'IL7R_i'
d38_pred <- prediction(d38$IL7R_i,d38$label)
d38_perf <- performance(d38_pred, "tpr", "fpr" )

d39 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/IL7R_IS2_voting_ni_test_score.csv',row.names = 1)
colnames(d39)[3] <- 'IL7R_ni'
d39_pred <- prediction(d39$IL7R_ni,d39$label)
d39_perf <- performance(d39_pred, "tpr", "fpr" )

auc(d31$label,d31$NeoTCR8) ## 0.8611
auc(d32$label,d32$CXCL13_i) # 0.9475
auc(d33$label,d33$CXCL13_ni) # 0.8819
auc(d34$label,d34$MANAscore_i) ## 0.966
auc(d35$label,d35$MANAscore_ni) ## 0.9705
auc(d36$label,d36$ENTPD1_i) ## 0.8889
auc(d37$label,d37$ENTPD1_ni) # 0.6763
auc(d38$label,d38$IL7R_i) # 0.9599
auc(d39$label,d39$IL7R_ni) ## 0.9313

## plot
tiff("02.all_8_model_IS2.tiff",width=6.5,height=5.5,units = "in",res=300)
par(mar = (c(4,4.5,1,1)))
plot(d34_perf, colorize = FALSE,col="firebrick2",lwd=2.5,cex.lab = 1.5)
plot(d35_perf, add = TRUE,colorize = FALSE,col="royalblue3",lwd=2,lty='dashed')

plot(d32_perf, add = TRUE,colorize = FALSE,col="darkorange",lwd=2)
plot(d33_perf,add = TRUE, colorize = FALSE,col="darkorange3",lwd=2,lty='dashed')

plot(d36_perf, add = TRUE,colorize = FALSE,col="cyan3",lwd=2)
plot(d37_perf,add = TRUE, colorize = FALSE,col="darkcyan",lwd=2,lty='dashed')

plot(d38_perf, add = TRUE,colorize = FALSE,col="darkorchid",lwd=2)
plot(d39_perf,add = TRUE, colorize = FALSE,col="darkorchid4",lwd=2,lty='dashed')


legend(x=0.55,y=0.4, c('MANAscore_i (0.97)','MANAscore_ni (0.97)',
                       'CXCL13_i (0.95)','CXCL13_ni (0.88)',
                       'ENTPD1_i (0.89)','ENTPD1_ni (0.68)',
                       'IL7R_i (0.96)','IL7R_ni (0.93)'), 
       lty = c(1,2,1,2,1,2,1,2),lwd = 2, 
       col = c("firebrick2","royalblue3","darkorange",'darkorange3','cyan3','darkcyan',
               'darkorchid','darkorchid4'), 
       bty="n", inset=c(0,-0.15),cex=0.9,text.font=1,title = 'AUC score')
abline(0,1,lty='dashed')
dev.off()

## PP3
##
d46 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CD39_PP3_voting_i_test_score.csv',row.names = 1)
colnames(d46)[3] <- 'ENTPD1_i'
d46_pred <- prediction(d46$ENTPD1_i,d46$label)
d46_perf <- performance(d46_pred, "tpr", "fpr" )

d47 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/CD39_PP3_voting_ni_test_score.csv',row.names = 1)
colnames(d47)[3] <- 'ENTPD1_ni'
d47_pred <- prediction(d47$ENTPD1_ni,d47$label)
d47_perf <- performance(d47_pred, "tpr", "fpr" )

d48 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/IL7R_PP3_voting_i_test_score.csv',row.names = 1)
colnames(d48)[3] <- 'IL7R_i'
d48_pred <- prediction(d48$IL7R_i,d48$label)
d48_perf <- performance(d48_pred, "tpr", "fpr" )

d49 = read.csv('D:/Project/ManaScore/01.Result/36.12_single_moldes/00.models/single_gene_model/IL7R_PP3_voting_ni_test_score.csv',row.names = 1)
colnames(d49)[3] <- 'IL7R_ni'
d49_pred <- prediction(d49$IL7R_ni,d49$label)
d49_perf <- performance(d49_pred, "tpr", "fpr" )


auc(d41$label,d41$NeoTCR8) ## 0.8293
auc(d42$label,d42$CXCL13_i) # 0.77
auc(d43$label,d43$CXCL13_ni) # 0.614
auc(d44$label,d44$MANAscore_i) ## 0.9546
auc(d45$label,d45$MANAscore_ni) ## 0.9429
auc(d46$label,d46$ENTPD1_i) ## 0.7017
auc(d47$label,d47$ENTPD1_ni) # 0.6018
auc(d48$label,d48$IL7R_i) # 0.956
auc(d49$label,d49$IL7R_ni) ## 0.9216
## plot
tiff("02.all_8_model_PP3.tiff",width=6.5,height=5.5,units = "in",res=300)
par(mar = (c(4,4.5,1,1)))
plot(d44_perf, colorize = FALSE,col="firebrick2",lwd=2.5,cex.lab = 1.5)
plot(d45_perf, add = TRUE,colorize = FALSE,col="royalblue3",lwd=2,lty='dashed')

plot(d42_perf, add = TRUE,colorize = FALSE,col="darkorange",lwd=2)
plot(d43_perf,add = TRUE, colorize = FALSE,col="darkorange3",lwd=2,lty='dashed')

plot(d46_perf, add = TRUE,colorize = FALSE,col="cyan3",lwd=2)
plot(d47_perf,add = TRUE, colorize = FALSE,col="darkcyan",lwd=2,lty='dashed')

plot(d48_perf, add = TRUE,colorize = FALSE,col="darkorchid",lwd=2)
plot(d49_perf,add = TRUE, colorize = FALSE,col="darkorchid4",lwd=2,lty='dashed')


legend(x=0.55,y=0.4, c('MANAscore_i (0.95)','MANAscore_ni (0.94)',
                       'CXCL13_i (0.77)','CXCL13_ni (0.61)',
                       'ENTPD1_i (0.70)','ENTPD1_ni (0.60)',
                       'IL7R_i (0.96)','IL7R_ni (0.92)'), 
       lty = c(1,2,1,2,1,2,1,2),lwd = 2, 
       col = c("firebrick2","royalblue3","darkorange",'darkorange3','cyan3','darkcyan',
               'darkorchid','darkorchid4'), 
       bty="n", inset=c(0,-0.15),cex=0.9,text.font=1,title = 'AUC score')
abline(0,1,lty='dashed')
dev.off()
