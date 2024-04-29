cd D:/Project/ManaScore/01.Result/34.combine_model_trained_on_me/

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_score

patient = ['p2','p11','p15']

dat1 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/p2_known.csv',index_col=0)
dat2 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/p2_known_RNA.csv',index_col=0)
dat3 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/p11_known.csv',index_col=0)
dat4 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/p11_known_RNA.csv',index_col=0)
dat5 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/p15_known.csv',index_col=0)
dat6 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/p15_known_RNA.csv',index_col=0)

dats = [dat1,dat2,dat3,dat4,dat5,dat6]

for i in range(6):
    dat = dats[i]
    X_train, X_test, y_train, y_test = train_test_split(dat.iloc[:,dat.columns !="Label"][['CXCL13','ENTPD1','IL7R']],
                                                        dat['Label'],
                                                        test_size=0.2,
                                                        random_state=42)
    myStr1 = 'X_train' + str(i+1)
    myVars1 = vars()
    myVars1.__setitem__(myStr1,X_train)
    
    myStr2 = 'X_test' + str(i+1)
    myVars2 = vars()
    myVars2.__setitem__(myStr2,X_test)
    
    myStr3 = 'y_train' + str(i+1)
    myVars3 = vars()
    myVars3.__setitem__(myStr3,y_train)
    
    myStr4 = 'y_test' + str(i+1)
    myVars4 = vars()
    myVars4.__setitem__(myStr4,y_test)

X_TRAIN = [X_train1,X_train2,X_train3,X_train4,X_train5,X_train6]
X_TEST = [X_test1,X_test2,X_test3,X_test4,X_test5,X_test6]
y_TRAIN = [y_train1,y_train2,y_train3,y_train4,y_train5,y_train6]
y_TEST = [y_test1,y_test2,y_test3,y_test4,y_test5,y_test6]

## ground truth
## LM
## LMi_p2, LMni_p2, LMi_p11,LMni_p11,LMi_p15,LMni_p15, LMi_IS2, LMni_IS2,LMi_pp, LMni_pp3

## RF
##  RFi_p11,RFni_p11, RFi_IS2, RFi_IS2, RFni_pp3
LM = ['LMi_p2', 'LMni_p2', 'LMi_p11','LMni_p11','LMi_p15','LMni_p15']
#RF = [RFi_p11,RFni_p11, RFi_IS2, RFi_IS2, RFni_pp3]
LM_models = []
dats = [dat1,dat2,dat3,dat4,dat5,dat6]
Res = []
for i in range(6):
    res= {}
    res['model']=LM[i]
    X_train = X_TRAIN[i]
    X_test = X_TEST[i]
    y_train = y_TRAIN[i]
    y_test = y_TEST[i]
    
    model = LogisticRegression()
    scores = cross_val_score(model, X_train, y_train, cv=5, scoring='roc_auc')
    res['score'] = scores.mean()
    
    # Make predictions on the testing data
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    # Calculate the probabilities of the positive class
    y_prob = model.predict_proba(X_test)[:,1]

    # Calculate the AUC score
    auc_score = roc_auc_score(y_test, y_prob)
    res['auc'] = auc_score
    Res.append(res)

    myStr = LM[i]
    myVars = vars()
    myVars.__setitem__(myStr,model)



LM_models = [LMi_p2, LMni_p2, LMi_p11,LMni_p11,LMi_p15,LMni_p15]

### RF = ['RFi_p11','RFni_p11', 'RFi_IS2', 'RFni_IS2', 'RFni_pp3'] ## selected RF models
RF = ['RFi_p2', 'RFni_p2', 'RFi_p11','RFni_p11','RFi_p15','RFni_p15']

Res1=[]
for i in range(6):
    res= {}
    res['model']=RF[i]
    
    X_train = X_TRAIN[i]
    X_test = X_TEST[i]
    y_train = y_TRAIN[i]
    y_test = y_TEST[i]

    ### Using GridSearchCV
    max_features = [1,2,3]
    n_estimators = [100,500,1000,2000,2500]
    param_distributions = dict(max_features = max_features,
                               n_estimators = n_estimators)
    ##
    dt = RandomForestClassifier(random_state=42)
            
    # Build grid search
    grid = GridSearchCV(estimator=dt,
                        param_grid = param_distributions,
                        cv=10)
    grid_result = grid.fit(X_train,y_train)
    
    res['best_pa'] = grid_result.best_params_
    
    best_model = grid_result.best_estimator_
            
    # Make predictions on the testing data
    best_model.fit(X_train,y_train)
    y_pred = best_model.predict(X_test)

    # Calculate the probabilities of the positive class
    y_prob = best_model.predict_proba(X_test)[:,1]

    # Calculate the AUC score
    auc_score = roc_auc_score(y_test, y_prob)
    res['auc'] = auc_score
    Res1.append(res)
    
    myStr = RF[i]
    myVars = vars()
    myVars.__setitem__(myStr,best_model)

RF_models = [RFi_p2, RFni_p2, RFi_p11,RFni_p11,RFi_p15,RFni_p15]       

#import voting classifier 6+6 models
from sklearn.ensemble import VotingClassifier

## imputation
# create a voting classifier with soft voting
voting_i_classifier = VotingClassifier(
    estimators = [('LMi_p2',LMi_p2),
                  #('LMni_p2', LMni_p2 ),
                  ('LMi_p11',LMi_p11),
                  #('LMni_p11',LMni_p11),
                  ('LMi_p15',LMi_p15),
                  #('LMni_p15',LMni_p15),
                  ('RFi_p2',RFi_p2),
                  #('RFni_p2',RFni_p2),
                  ('RFi_p11',RFi_p11),
                  #('RFni_p11',RFni_p11),
                  ('RFi_p15',RFi_p15),
                  #('RFni_p15',RFni_p15)
                 ], 
    voting='soft')


### non imputation

# create a voting classifier with soft voting
voting_ni_classifier = VotingClassifier(
    estimators = [('LMni_p2',LMni_p2),
                  #('LMi_p2', LMi_p2 ),
                  #('LMi_p11',LMi_p11),
                  ('LMni_p11',LMni_p11),
                  #('LMi_p15',LMi_p15),
                  ('LMni_p15',LMni_p15),
                  #('RFi_p2',RFi_p2),
                  ('RFni_p2',RFni_p2),
                  #('RFi_p11',RFi_p11),
                  ('RFni_p11',RFni_p11),
                  #('RFi_p15',RFi_p15),
                  ('RFni_p15',RFni_p15)
                 ], 
    voting='soft')

## imputation
# make predictions with the hard soft model
mydats = [dat1,dat3,dat5]

voting_i_classifier.fit(pd.concat([X_train1, X_train3, X_train5]), 
               np.concatenate((y_train1, y_train3, y_train5)))

Y_PRED_i = voting_i_classifier.predict(pd.concat([X_test1, X_test3, X_test5]))
Y_PROB_i = voting_i_classifier.predict_proba(pd.concat([X_test1, X_test3, X_test5]))[:, 1]

AUC_SCORE_i = roc_auc_score(pd.concat([y_test1, y_test3, y_test5]), Y_PROB_i)    

di = pd.DataFrame({'barcode':pd.concat([X_test1, X_test3, X_test5]).index,
                   'label':pd.concat([y_test1, y_test3, y_test5]),
                 'score':Y_PROB_i})
di.to_csv('voting_6_models_i_test_score.csv')


AUC_SCORE_i

## non-imputation
# make predictions with the soft voting model
mydats = [dat2,dat4,dat6]
voting_ni_classifier.fit(pd.concat([X_train2, X_train4, X_train6]), 
               np.concatenate((y_train2, y_train4, y_train6)))

Y_PRED_ni = voting_ni_classifier.predict(pd.concat([X_test2, X_test4, X_test6]))
Y_PROB_ni = voting_ni_classifier.predict_proba(pd.concat([X_test2, X_test4, X_test6]))[:, 1]

AUC_SCORE_ni = roc_auc_score(pd.concat([y_test2, y_test4, y_test6]), Y_PROB_ni)

dni = pd.DataFrame({'barcode':pd.concat([X_test2, X_test4, X_test6]).index,
                    'label':pd.concat([y_test2, y_test4, y_test6]),
                    'score':Y_PROB_ni})

dni.to_csv('voting_6_models_ni_test_score.csv')


AUC_SCORE_ni

## test on IS2
## imputation
IS2i = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/IS2_known.csv',index_col=0)
X_TEST_i = IS2i.iloc[:,IS2i.columns !="Label"][['CXCL13','ENTPD1','IL7R']]
Y_TEST_i = IS2i['Label']

Y_PRED_IS2i = voting_i_classifier.predict(X_TEST_i)

# Calculate the probabilities of the positive class
Y_PROB_IS2i = voting_ni_classifier.predict_proba(X_TEST_i)[:, 1]

AUC_SCORE_IS2i = roc_auc_score(Y_TEST_i, Y_PROB_IS2i)

IS2_di = pd.DataFrame({'barcode':X_TEST_i.index,
                    'label':Y_TEST_i,
                    'score':Y_PROB_IS2i})
IS2_di.to_csv('IS2_voting_i_test_score.csv')

AUC_SCORE_IS2i

## test on IS2
## non-imputation
IS2ni = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/IS2_known_RNA.csv',index_col=0)
X_TEST_ni = IS2ni.iloc[:,IS2ni.columns !="Label"][['CXCL13','ENTPD1','IL7R']]
Y_TEST_ni = IS2ni['Label']

Y_PRED_IS2ni = voting_ni_classifier.predict(X_TEST_i)

# Calculate the probabilities of the positive class
Y_PROB_IS2ni = voting_ni_classifier.predict_proba(X_TEST_ni)[:, 1]

AUC_SCORE_IS2ni = roc_auc_score(Y_TEST_ni, Y_PROB_IS2ni)

IS2_dni = pd.DataFrame({'barcode':X_TEST_ni.index,
                    'label':Y_TEST_ni,
                    'score':Y_PROB_IS2ni})
IS2_dni.to_csv('IS2_voting_ni_test_score.csv')

AUC_SCORE_IS2ni

## test on PP3
## imputation
PP3i = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/PP3_known.csv',index_col=0)
X_TEST_i = PP3i.iloc[:,PP3i.columns !="Label"][['CXCL13','ENTPD1','IL7R']]
Y_TEST_i = PP3i['Label']

Y_PRED_PP3i = voting_i_classifier.predict(X_TEST_i)

# Calculate the probabilities of the positive class
Y_PROB_PP3i = voting_ni_classifier.predict_proba(X_TEST_i)[:, 1]

AUC_SCORE_PP3i = roc_auc_score(Y_TEST_i, Y_PROB_PP3i)

PP3_di = pd.DataFrame({'barcode':X_TEST_i.index,
                    'label':Y_TEST_i,
                    'score':Y_PROB_PP3i})
PP3_di.to_csv('PP3_voting_i_test_score.csv')

AUC_SCORE_PP3i

## test on PP3
## non-imputation
PP3ni = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/PP3_known_RNA.csv',index_col=0)
X_TEST_ni = PP3ni.iloc[:,PP3ni.columns !="Label"][['CXCL13','ENTPD1','IL7R']]
Y_TEST_ni = PP3ni['Label']

Y_PRED_PP3ni = voting_ni_classifier.predict(X_TEST_i)

# Calculate the probabilities of the positive class
Y_PROB_PP3ni = voting_ni_classifier.predict_proba(X_TEST_ni)[:, 1]

AUC_SCORE_PP3ni = roc_auc_score(Y_TEST_ni, Y_PROB_PP3ni)

IS2_dni = pd.DataFrame({'barcode':X_TEST_ni.index,
                    'label':Y_TEST_ni,
                    'score':Y_PROB_PP3ni})
IS2_dni.to_csv('PP3_voting_ni_test_score.csv')

AUC_SCORE_PP3ni

### test on Head & Neck data
## imputation
HNi = pd.read_csv('D:/Project/ManaScore/01.Result/10.HN/single_patient_data/Nivo1_known.csv',index_col=0)
X_TEST_i = HNi.iloc[:,HNi.columns !="Label"][['CXCL13','ENTPD1','IL7R']]
Y_TEST_i = HNi['Label']

Y_PRED_HNi = voting_i_classifier.predict(X_TEST_i)

# Calculate the probabilities of the positive class
Y_PROB_HNi = voting_ni_classifier.predict_proba(X_TEST_i)[:, 1]

AUC_SCORE_HNi = roc_auc_score(Y_TEST_i, Y_PROB_HNi)

HN_di = pd.DataFrame({'barcode':X_TEST_i.index,
                    'label':Y_TEST_i,
                    'score':Y_PROB_HNi})
HN_di.to_csv('Nivo1_voting_i_test_score.csv')



AUC_SCORE_HNi

HNni = pd.read_csv('D:/Project/ManaScore/01.Result/10.HN/single_patient_data/Nivo1_known_RNA.csv',index_col=0)
X_TEST_ni = HNni.iloc[:,HNni.columns !="Label"][['CXCL13','ENTPD1','IL7R']]
Y_TEST_ni = HNni['Label']

Y_PRED_HNni = voting_ni_classifier.predict(X_TEST_ni)

# Calculate the probabilities of the positive class
Y_PROB_HNni = voting_ni_classifier.predict_proba(X_TEST_ni)[:, 1]

AUC_SCORE_HNni = roc_auc_score(Y_TEST_ni, Y_PROB_HNni)

HN_dni = pd.DataFrame({'barcode':X_TEST_ni.index,
                    'label':Y_TEST_ni,
                    'score':Y_PROB_HNni})
HN_dni.to_csv('Nivo1_voting_ni_test_score.csv')


AUC_SCORE_HNni

## p2:MAA-1 nonMANA-0 imputation
p2_MAA_nonMANA = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/p2_MAA_nonMANA.csv',index_col=0)
p2_MAA_nonMANA_1 = p2_MAA_nonMANA[p2_MAA_nonMANA['Label']==1]
p2_MAA_nonMANA_0 = p2_MAA_nonMANA[p2_MAA_nonMANA['Label']==0]
p2_MAA_nonMANA_0 = p2_MAA_nonMANA_0.loc[[i for i in p2_MAA_nonMANA_0.index if i in y_test1.index]]
p2_MAA_nonMANA = pd.concat([p2_MAA_nonMANA_1,p2_MAA_nonMANA_0])

X_TEST_i = p2_MAA_nonMANA.iloc[:,p2_MAA_nonMANA.columns !="Label"][['CXCL13','ENTPD1','IL7R']]
Y_TEST_i = p2_MAA_nonMANA['Label']

Y_PRED_MAAi = voting_i_classifier.predict(X_TEST_i)

# Calculate the probabilities of the positive class
Y_PROB_MAAi = voting_ni_classifier.predict_proba(X_TEST_i)[:, 1]

AUC_SCORE_MAAi = roc_auc_score(Y_TEST_i, Y_PROB_MAAi)

MAA_di = pd.DataFrame({'barcode':X_TEST_i.index,
                    'label':Y_TEST_i,
                    'score':Y_PROB_MAAi})
MAA_di.to_csv('MAA_voting_i_test_score.csv')


AUC_SCORE_MAAi

## p2:MAA-1 nonMANA-0 non-imputation
p2_MAA_nonMANAni = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/p2_MAA_nonMANA_RNA.csv',index_col=0)
p2_MAA_nonMANA_1 = p2_MAA_nonMANAni[p2_MAA_nonMANAni['Label']==1]
p2_MAA_nonMANA_0 = p2_MAA_nonMANAni[p2_MAA_nonMANAni['Label']==0]
p2_MAA_nonMANA_0 = p2_MAA_nonMANA_0.loc[[i for i in p2_MAA_nonMANA_0.index if i in y_test2.index]]
p2_MAA_nonMANAni = pd.concat([p2_MAA_nonMANA_1,p2_MAA_nonMANA_0])

X_TEST_ni = p2_MAA_nonMANAni.iloc[:,p2_MAA_nonMANAni.columns !="Label"][['CXCL13','ENTPD1','IL7R']]
Y_TEST_ni = p2_MAA_nonMANAni['Label']

Y_PRED_MAAni = voting_ni_classifier.predict(X_TEST_ni)

# Calculate the probabilities of the positive class
Y_PROB_MAAni = voting_ni_classifier.predict_proba(X_TEST_ni)[:, 1]

AUC_SCORE_MAAni = roc_auc_score(Y_TEST_ni, Y_PROB_MAAni)

MAA_dni = pd.DataFrame({'barcode':X_TEST_ni.index,
                    'label':Y_TEST_ni,
                    'score':Y_PROB_MAAni})
MAA_dni.to_csv('MAA_voting_ni_test_score.csv')


AUC_SCORE_MAAni

## predict for unseen data
prediction_patient = ['PM1','NY016-015','EF7','NY016-022','JS10','PL5','NY016-007','JM6','NY016-025',
  'PP4','NY016-021','AH12','PP3','IS2','NY016-014']
## imputation data
for ps in prediction_patient:
    pdat1 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/'+ ps+'_unknown.csv',index_col=0)
    pdat2 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/'+ ps+'_unknown_RNA.csv',index_col=0)
    prob1 = voting_i_classifier.predict_proba(pdat1[['CXCL13','ENTPD1','IL7R']])[:,1]
    d1 = pd.DataFrame({'barcode':pdat1.index,
                       'score':prob1})
    d1.to_csv('voting_i_'+ps+'_predict_score.csv')
    
    prob2 = voting_ni_classifier.predict_proba(pdat2[['CXCL13','ENTPD1','IL7R']])[:,1]
    d2 = pd.DataFrame({'barcode':pdat2.index,
                       'score':prob2})
    d2.to_csv('voting_ni_'+ps+'_predict_score.csv')
            

## JM6 total
JM6_all_RNA = pd.read_csv("D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/JM6_tumor_all_RNA.csv",index_col=0)
JM6_all = pd.read_csv("D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/JM6_tumor_all.csv",index_col=0)

prob1 = voting_i_classifier.predict_proba(JM6_all[['CXCL13','ENTPD1','IL7R']])[:,1]
d1 = pd.DataFrame({'barcode':JM6_all.index,
                   'score':prob1})
d1.to_csv('voting_i_JM6_all_predict_score.csv')
    
prob2 = voting_ni_classifier.predict_proba(JM6_all_RNA[['CXCL13','ENTPD1','IL7R']])[:,1]
d2 = pd.DataFrame({'barcode':JM6_all_RNA.index,
                   'score':prob2})
d2.to_csv('voting_ni_JM6_all_predict_score.csv')
            


## JM6 mettumor
JM6_met_RNA = pd.read_csv("D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/JM6_mettumor_unknown_RNA.csv",index_col=0)
JM6_met = pd.read_csv("D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/JM6_mettumor_unknown.csv",index_col=0)

prob1 = voting_i_classifier.predict_proba(JM6_met[['CXCL13','ENTPD1','IL7R']])[:,1]
d1 = pd.DataFrame({'barcode':JM6_met.index,
                   'score':prob1})
d1.to_csv('voting_i_JM6_met_predict_score.csv')
    
prob2 = voting_ni_classifier.predict_proba(JM6_met_RNA[['CXCL13','ENTPD1','IL7R']])[:,1]
d2 = pd.DataFrame({'barcode':JM6_met_RNA.index,
                   'score':prob2})
d2.to_csv('voting_ni_JM6_met_predict_score.csv')


#NSCLC normal and LN
prediction_patient = ["MD01-004_LN","MD01-005_LN","MD01-005_normal","MD01-010_normal","MD01-019_normal","MD043-003_normal","MD043-006_normal","MD043-008_normal","MD043-011_LN",
                      "MD043-011_normal","NY016-007_normal","NY016-015_normal","NY016-016_normal","NY016-022_normal","NY016-025_normal"]
## imputation data
for ps in prediction_patient:
    pdat1 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/NSCLC_tissue/'+ ps+'_unknown.csv',index_col=0)
    prob1 = voting_i_classifier.predict_proba(pdat1[['CXCL13','ENTPD1','IL7R']])[:,1]
    d1 = pd.DataFrame({'barcode':pdat1.index,
                       'score':prob1})
    d1.to_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/NSCLC_tissue/voting_i_'+ps+'_predict_score.csv')

## non-imputation
pdat2 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/NSCLC_tissue/NSCLC_normal_LN_RNA.csv',index_col=0)
prob2 = voting_ni_classifier.predict_proba(pdat2[['CXCL13','ENTPD1','IL7R']])[:,1]
d2 = pd.DataFrame({'barcode':pdat2.index,'score':prob2})
d2.to_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/NSCLC_tissue/voting_ni_NSCLC_normal_LN_predict_score.csv')

## Head Neck
prediction_patient = ['Nivo1','Nivo2','Nivo3','Nivo10',"Nivo4","J12","Nivo18","Nivo17","J11","Nivo15","Nivo23","Nivo26",
  "Nivo27","Nivo21","Nivo19"]
## imputation data
for ps in prediction_patient:
    pdat1 = pd.read_csv('D:/Project/ManaScore/01.Result/10.HN/single_patient_data/'+ ps+'_unknown.csv',index_col=0)
    pdat2 = pd.read_csv('D:/Project/ManaScore/01.Result/10.HN/single_patient_data/'+ ps+'_unknown_RNA.csv',index_col=0)
    prob1 = voting_i_classifier.predict_proba(pdat1[['CXCL13','ENTPD1','IL7R']])[:,1]
    d1 = pd.DataFrame({'barcode':pdat1.index,
                       'score':prob1})
    d1.to_csv('voting_i_'+ps+'_predict_score.csv')
    
    prob2 = voting_ni_classifier.predict_proba(pdat2[['CXCL13','ENTPD1','IL7R']])[:,1]
    d2 = pd.DataFrame({'barcode':pdat2.index,
                       'score':prob2})
    d2.to_csv('voting_ni_'+ps+'_predict_score.csv')
   

## Head Neck Pre
prediction_patient = ['Nivo1','Nivo2','Nivo3','Nivo10',"Nivo4","Nivo17","Nivo15","Nivo23","Nivo27","Nivo21"]
## imputation data
for ps in prediction_patient:
    pdat1 = pd.read_csv('D:/Project/ManaScore/01.Result/10.HN/single_patient_sample/'+ ps+'_Pre_unknown.csv',index_col=0)
    pdat2 = pd.read_csv('D:/Project/ManaScore/01.Result/10.HN/single_patient_sample/'+ ps+'_Pre_unknown_RNA.csv',index_col=0)
    prob1 = voting_i_classifier.predict_proba(pdat1[['CXCL13','ENTPD1','IL7R']])[:,1]
    d1 = pd.DataFrame({'barcode':pdat1.index,
                       'score':prob1})
    d1.to_csv('voting_i_'+ps+'_predict_score.csv')
    
    prob2 = voting_ni_classifier.predict_proba(pdat2[['CXCL13','ENTPD1','IL7R']])[:,1]
    d2 = pd.DataFrame({'barcode':pdat2.index,
                       'score':prob2})
    d2.to_csv('voting_ni_'+ps+'_predict_score.csv')

## Head Neck Post
prediction_patient = ['Nivo1','Nivo2','Nivo3','Nivo10',"Nivo4","Nivo17","Nivo15","Nivo23","Nivo27","Nivo21"]
## imputation data
for ps in prediction_patient:
    pdat1 = pd.read_csv('D:/Project/ManaScore/01.Result/10.HN/single_patient_sample/'+ ps+'_Post_unknown.csv',index_col=0)
    pdat2 = pd.read_csv('D:/Project/ManaScore/01.Result/10.HN/single_patient_sample/'+ ps+'_Post_unknown_RNA.csv',index_col=0)
    prob1 = voting_i_classifier.predict_proba(pdat1[['CXCL13','ENTPD1','IL7R']])[:,1]
    d1 = pd.DataFrame({'barcode':pdat1.index,
                       'score':prob1})
    d1.to_csv('voting_i_'+ps+'_predict_score.csv')
    
    prob2 = voting_ni_classifier.predict_proba(pdat2[['CXCL13','ENTPD1','IL7R']])[:,1]
    d2 = pd.DataFrame({'barcode':pdat2.index,
                       'score':prob2})
    d2.to_csv('voting_ni_'+ps+'_predict_score.csv')

## Oral Cancer
prediction_patient = ["P13","P14","P15","P16","P18","P19","P20","P21","P22","P23","P24","P25","P26","P27","P28","P29","P30","P31","P32"]
## imputation data
for ps in prediction_patient:
    pdat1 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/Oral_cancer/'+ ps+'_unknown.csv',index_col=0)
    pdat2 = pd.read_csv('D:/Project/ManaScore/01.Result/24.normalize_for_each_patient/Oral_cancer/'+ ps+'_unknown_RNA.csv',index_col=0)
    prob1 = voting_i_classifier.predict_proba(pdat1[['CXCL13','ENTPD1','IL7R']])[:,1]
    d1 = pd.DataFrame({'barcode':pdat1.index,
                       'score':prob1})
    d1.to_csv('voting_i_Oral_'+ps+'_predict_score.csv')
    
    prob2 = voting_ni_classifier.predict_proba(pdat2[['CXCL13','ENTPD1','IL7R']])[:,1]
    d2 = pd.DataFrame({'barcode':pdat2.index,
                       'score':prob2})
    d2.to_csv('voting_ni_Oral_'+ps+'_predict_score.csv')
  

# Merkel
pdat1 = pd.read_csv('D:/Project/tmp/Merkel/MANAscore/Merkel_MANAthrGn.csv',index_col=0)
prob1 = voting_i_classifier.predict_proba(pdat1[['CXCL13','ENTPD1','IL7R']])[:,1]
d1 = pd.DataFrame({'barcode':pdat1.index,'score':prob1})
d1.to_csv('D:/Project/tmp/Merkel/MANAscore/voting_i_Merkel_predict_score.csv')

pdat2 = pd.read_csv('D:/Project/tmp/Merkel/MANAscore/raw_thrGn.csv',index_col=0)
prob2 = voting_ni_classifier.predict_proba(pdat2[['CXCL13','ENTPD1','IL7R']])[:,1]
d2 = pd.DataFrame({'barcode': pdat2.index,'score':prob2})
d2.to_csv('D:/Project/tmp/Merkel/MANAscore/voting_ni_Merkel_predict_score.csv')

