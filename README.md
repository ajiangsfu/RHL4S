
# RHL4S 


RHL4S is a R package to use a spatial score-based prognostic model to predict risk group for classic Hodgkin lymphoma 

Example code is shown in the following:

library(RHL4S)

allout = RHL4Scalls(datapath = "../data", patientFile = "patient.xlsx", patientSheet = "B")

# References:
Tomohiro Aoki, Aixiang Jiang, Alexander Xu, Alicia Gamboa, Yifan Yin, Katy Milne, Celia Strong, Talia Goodyear, Shaocheng Wu, Lauren C. Chong, Katsuyoshi Takata, Elizabeth Chavez, Tomoko Miyata-Takata, Anthony R. Colombo, Monirath Hav, Adele Telenius, Susana Ben-Neriah, Andrew P. Weng, Kerry J. Savage, David W. Scott, Andrew Roth, Pedro Farinha, Brad H Nelson, Akil Merchant, Christian Steidl; Spatial Tumor Microenvironment Characterization and Outcome of Relapsed/Refractory Classic Hodgkin Lymphoma. Blood 2022; 140 (Supplement 1): 170–171. doi: https://doi.org/10.1182/blood-2022-158709
