# A minimal gene set characterizes multiple classes of tumor-specific TIL among different cancer types
## Overview
Pinpointing tumor-specific T cell clones responsible for immunotherapy responses remains challenging. Based on previous findings that validated mutation-associated neoantigen (MANA)-specific CD8+ tumor-infiltrating lymphocytes (TIL) from several cancer types express high levels of CXCL13 and CD39 (ENTPD1), and low IL7R levels, we develop a ‘MANAscore’ using weighted expression levels of these three genes from lung cancer and melanoma single-cell RNAseq, aiming at prospectively identify MANA-specific TIL. Our three-gene MANAscore algorithm outperforms other RNAseq-based algorithms in identifying validated neoantigen-specific CD8+ clones and is further demonstrated to accurately identify TIL recognizing other classes of tumor antigens, including cancer-testis antigens, endogenous retroviruses, and viral oncogenes, and most of these TIL expressed a tissue resident memory program. 

## MANAscore
### Description
Model construction and MANAscore prediction.
### Installation
Users shall also install the following python packages before install MANAscore
- git-lfs



  Users can use [Anaconda](https://www.anaconda.com/download) to create a conda environment:
  ```
  source path/of/your/anaconda/bin/activate ## activate your anaconda
  conda create -n manascore python==3.9
  ```
  Activate environment
  ```
  conda activate manascore
  ```
  Then install git-lfs for downloading large files:
  ```
  conda install -c conda-forge git-lfs
  ```
  Initiate git-lfs for downloading large model files saved:
  ```
  git lfs install
  ```
  Then install MANAscore:
  ```
  git clone https://github.com/BKI-immuno-KNS/MANAscore.git
  cd path/of/MANAscore/MANAscore/
  pip install -e .
  ```
### Process
- step 1 Construct 6 linear regression models
- step 2 Construct 6 random forest models
- step 3 Construct imputation combine voting models (3 linear regression models and 3 random forest models from imputaion data) and non-imputation combine voting model (3 linear regression models and 3 random forest models from non-imputaion data)
### Files
- ./MANAscore/3gene/ Three gene imputed and non-imputed matrix of ground truth for training and test data in melanoma, validation data in lung cancer. MANA-/MAA-specific TIL are with label of 1, EBV-/flu-specific TIL are with label of 0.
- ./MANAscore/models/ The imputation combine voting model (voting_i_classifier.pkl.gz) and non-imputation combine voting model (voting_ni_classifier.pkl.gz) saved for MANAscore prediction.


Users can directly predict MANAscore using saved models under /models for their own data
```
manascore predict input_imputation_3_gene.csv input_non_imputation_3_gene.csv output_imputation_score.csv output_non_imputation_score.csv
```
Users can also construct models and then use the models to predict MANAscore
```
manascore train_predict input_imputation_3_gene.csv input_non_imputation_3_gene.csv output_imputation_score.csv output_non_imputation_score.csv
```

## Tutorial
For the user manual of MANAscore, please refer to: https://bki-immuno-kns.github.io/MANAscore/MANAscore-Tutorial.html
## Analyses
### Description
Including the R scripts for data integration, data preprocessing for MANAscore prediction, model evaluation and differential gene/signature analyses.
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





