
# RHL4S 

RHL4S is a R package to use a spatial score-based prognostic model to predict risk group for classic Hodgkin lymphoma.

I. Algorithm
-------------------------------

RHL4S is a prognostic model developed for patients with relapsed or refractory classic Hodgkin lymphoma (r/r CHL) who have undergone high dose chemotherapy
followed by autologous stem cell transplantation (HDT/ASCT). The model is based on imaging mass cytometry (IMC) and incorporates information about cellular interactions, 
specific expression features of Hodgkin and Reed-Sternberg (HRS) cells, and the spatial architecture of the tumor microenvironment (TME). 
The study found that CXCR5+ HRS cells, PD1+ CD4+ T cells, macrophages, and CXCR5+ B cells were significantly associated with post-ASCT failure-free survival. 
The RHL4S model was found to be more effective than classical protein percentage-based models in predicting outcomes in r/r CHL patients.

The RHL4S prognostic model is based on four spatial score variables that are inflated into a 0-100 scale range: CXCR5 HRS spatial score, PD1 CD4 spatial score, 
Mac spatial score, and CXCR5 B spatial score. The model was built using our new XGpred algorithm, which combines the machine learning method XGBoost with
traditional statistical techniques such as model-based clustering, spline regression, LPS (Linear Prediction Score), and empirical Bayesian approaches. XGPred functions will be added into R package csmpv at https://github.com/ajiangsfu/csmpv.

To make the model applicable to patients in daily clinical practice, the findings from IMC were translated into simplified data. 
A simplified panel was built to obtain information on the four variables from the predictive model and calculate the spatial score from MC-IF data. 
Calibration was performed by applying this panel to a subset of the IMC cohort, confirming the concordance between spatial score from IMC and MC-IF, 
and setting the optimal cut-off on MC-IF for RHL4S assay. The mean difference in RHL4S scores between IMC and MC-IF for the calibration cohort was used 
as a parameter to adjust RHL4S scores for any MC-IF validation cohorts.

The current RHL4S R package calculates RHL4S scores and predicts RHL4S risk group classification based on MC-IF data.

### References:
Tomohiro Aoki, Aixiang Jiang, Alexander Xu, Alicia Gamboa, Yifan Yin, Katy Milne, Celia Strong, Talia Goodyear, Shaocheng Wu, Lauren C. Chong, Katsuyoshi Takata, Elizabeth Chavez, Tomoko Miyata-Takata, Anthony R. Colombo, Monirath Hav, Adele Telenius, Susana Ben-Neriah, Andrew P. Weng, Kerry J. Savage, David W. Scott, Andrew Roth, Pedro Farinha, Brad H Nelson, Akil Merchant, Christian Steidl; Spatial Tumor Microenvironment Characterization and Outcome of Relapsed/Refractory Classic Hodgkin Lymphoma. Blood 2022; 140 (Supplement 1): 170–171. doi: https://doi.org/10.1182/blood-2022-158709


II Package installation
-------------------------------

RHL4S R package is currently avaiable at GitHub (<https://github.com/ajiangsfu/RHL4S>).

There are at least two ways to install RHL4S R package from GitHub.

### 1. Install RHL4S with R package "devtools""
``` r
install.packages("devtools")  ### Run this line only if "devtools"" is not installed, otherwise please ignore this line
devtools :: install_github(repo = "ajiangsfu/RHL4S",force = TRUE)
```

### 2. Install RHL4S with R package "remotes"

``` r
install.packages("remotes")  ### Run this line only if "remotes"" is not installed, otherwise please ignore this line
remotes :: install_github(repo = "ajiangsfu/RHL4S", force = TRUE)
```

III Example code
-------------------------------------

Example code is shown in the following:

library(RHL4S)

allout = RHL4Scalls(datapath = "../data", patientFile = "patient.xlsx", patientSheet = "B")


RHL4S R package general information
----------------------------------

Version: 0.1.1

Author: Aixiang Jiang, Yifan Yin, and Alex Xu

Maintainer: Maintainer: Aixiang Jiang <aijiang@bccrc.ca>

Depends: R (>= 4.0)

Suggests: knitr

VignetteBuilder: knitr

Imports: data.table, dplyr, tidyr, readxl, spatstat.geom, magrittr,
        stringr
        


