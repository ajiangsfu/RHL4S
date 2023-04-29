---
output:
  pdf_document: default
  html_document: default
---
RHL4S: 
=======================================================

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

RHL4S is a R package to use a spatial score-based prognostic model to predict risk group for classic Hodgkin lymphoma 

Example code is shown in the following

library(RHL4S)
allout = RHL4Scalls(datapath = "../data", patientFile = "patient.xlsx", patientSheet = "B")

