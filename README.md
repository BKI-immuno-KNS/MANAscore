# A minimal gene set characterizes multiple classes of tumor-specific TIL among different cancer types
## MANAscore
### Package required
- python==3.9
- scikit-learn==1.4.2
- pandas
- numpy
### Process
- step 1 Construct 6 linear regression models
- step 2 Construct 6 random forest models
- step 3 Construct imputation combine voting models (3 linear regression models and 3 random forest models from imputaion data) and non-imputation combine voting model (3 linear regression models and 3 random forest models from non-imputaion data)
### Files
- 3gene/ Three gene imputed and non-imputed matrix of ground truth for training and test data in melanoma, validation data in Lung cancer. MANA-/MAA-specific TIL are with label of 1, EBV-/flu-specific TIL are with label of 0.
- models/ The imputation combine voting model (voting_i_classifier.pkl.gz) and non-imputation combine voting model (voting_ni_classifier.pkl.gz) saved for MANAscore prediction.
- examples:
  - example1.py: Loading exsited models for predition.
  - example2.py: Starting with loading ground truth to build the voting models then predict the MANAscore.
## Analyses

- ./Analyses/00.integration.R Integration CD8 T cells from NSCLC and melanoma
- ./Analyses/01.MANAscore.data.prepare.R
- ./Analyses/02.ROC.curves.R ROC curves generated for different models including published signatures (NeoTCR8 and CXCL13) on test and validation data
- ./Analyses/03.cutoff.R Patient specific cutoffs for defining MANAscorehi T cells and pTRC
- ./Analyses/04.signature.comparison.pTRC.vs.non-pTRC.R Signatures including checkpoint, cytotoxicity, TCR signaling, TRM and clonal size were compared between pTRC and non-pTRC
- ./Analyses/05.signature.comparison.pTRC.R.vs.NR.R Signatures including checkpoint, cytotoxicity, TCR signaling, TRM and clonal size were compared between pTRC in responder and non-responder
- ./Analyses/06.pseudobulk.dge.analysis.R DGE anlyses between pTRC and non-pTRC, pTRC in R and NR
- ./Analyses/07.tissue.compartment.comparison.R Tissue compartment comparison based on MANAscore
- ./Analyses/08.MANAscore.on.Oral_cancer.R
- ./Analyses/09.MANAscore.on.metastatic_cancer.R
- ./Analyses/10.MANAscore.on.MCPyV_pos_MCC.R





