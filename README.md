# MTLRRC
Authors: Akira Okazaki and Shuichi Kawano
## Introduction
This is an implementation of the proposed method in "Multi-task learning via robust regularized clustering with non-convex group penalties
".
The R codes of MTLRRC, MTLCVX, and MTLK are provided, and the codes for generating artificial datasets used for simulation studies are also provided.
We have obtained the school data from the "MALSAR" package for MATLAB (https://jiayuzhou.github.io/MALSAR/) and the landmine data from https://www.cs.columbia.edu/~jebara/code/multisparse/.



Details of 

**1. Import some libraries and our source codes.**

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

