# MTLRRC
Authors: Akira Okazaki and Shuichi Kawano
## Introduction
This is an implementation of the proposed method in the research paper titled "Multi-task learning via robust regularized clustering with non-convex group penalties".
The R codes of MTLRRC, MTLCVX, and MTLK are provided, and the codes for generating artificial datasets used for simulation studies are also provided.

### Code explanation
- R_sourcefiles
  - base_function.R: Simple functions used in other functions.
  - estimator.R: Main estimation functions.
  - generate_data.R: Functions for generation of artificial data and processing of real datasets.
  - Initialize_R.R: Functions for constructing a weight matrix R. 
- demo codes
  - realdata_MTLRRC_demo.R: A demonstration of MTLRRC for real data analysis studied in the research paper. The analyzed datasets "school data" and "landmine data" can be obtained from the "MALSAR" package for MATLAB (https://jiayuzhou.github.io/MALSAR/) and the landmine data from https://www.cs.columbia.edu/~jebara/code/multisparse/, respectively. In order to read the dataset, please prepare csv files divided by each task and further divided into the response and explanatory variables. Specifically, for the k-th task in the landmine data, the response file should be "landmine_Yk.csv" and the explanatory variables should be "landmine_Xk.csv". For more information, check the ”make_real_dataset” function in "generate_data.R".
  - realdata_MTLCVX_demo.R: A demonstration of MTLCVX for real data analysis studied in the research paper.
  - simulation_MTLRRC_demo.R: A demonstration of MTLRRC for simulation studies in the research paper including Case1 and Case2．

# Import some libraries and our source codes.

  ```
  library(MASS)
  library(glmnet)
  library(tidyverse)
  library(rsample)
  library(ROCR)
  library(truncnorm)

  setwd("~/R_sourcefiles/") #directory of sourcefiles
  source("base_function.R")
  source("estimator.R")
  source("generate_data.R")
  source("Initialize_R.R")

```

