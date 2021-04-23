# Robust Machine Learning for Treatment Effects in Multilevel Observational Studies Under Cluster-level Unmeasured Confounding

Youmi Suk and Hyunseung Kang

## Overview

Recently, machine learning (ML) methods have been used in causal inference to estimate treatment effects in order to reduce concerns for model mis-specification. However, many ML methods require that all confounders are measured to consistently estimate treatment effects. In this paper, we propose a family of ML methods that estimate treatment effects in the presence of cluster-level unmeasured confounders, a type of unmeasured confounders that are shared within each cluster and are common in multilevel observational studies. We show through simulation studies that our proposed methods are robust from biases from unmeasured cluster-level confounders in a variety of multilevel observational studies. We also examine the effect of taking an algebra course on math achievement scores from the Early Childhood Longitudinal Study, a multilevel observational educational study, using our methods. The proposed methods are available in the CURobustML R package.

Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis using the kindergarten cohort of ECLS (ECLS-K) data. 

## Simulation Study

* `DGP_twolevel_crossclassified.R`  

   This `R` file includes data generating codes for multilevel observational data with cluster-level unmeasured confounding.
 
```R
twolevel.pop
ccrem.pop
```

* `CURobustML_sepFun.R`  

   This `R` file includes functions for the propensity score and outcome predictions based on our proposed methods.

* `Simulation_robustML.R`
 
   This `R` file includes simulation codes with our proposed methods. For more details on simulation condtions, see our paper, https://psyarxiv.com/t7vbz/.


## ECLS-K Data Study

* `ECLSK_Algebra_complete.csv`

  This is our complete data. The original ECLSK 1998-99 dataset is available at https://nces.ed.gov/ecls/dataproducts.asp. 

* `CURobustML_sepFun.R`  

   Again, this `R` file includes functions for the propensity score and outcome predictions based on our proposed methods.

* `ECLSK_Algebra_robustML.R` 
 
   This `R` file can be used to replicate our data analysis.
