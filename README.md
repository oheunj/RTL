# Reluctant Transfer Learning in Penalized Regressions for Individualized Treatment Rules under Effect Heterogeneity

# Files in this repository
The source code is currently provided in 'rtl_main_code.R'

# Installation
R is a statistical software program, and RStudio is a user interface for R. We recommend that users install both R and R Studio. Both R and RStudio are free and open source.

## Requirements
You may need to install the following dependencies first:
```{r}
library(MASS)
library(survival)
library(glmnet)
library(mice)
library(smcfcs)
library(survcomp)
library(DescTools)
library(Hmisc)
library(knitr)
```
Additionally, you need to save the following R codes from the 'functions' folder to your working directory and import these files:
```{r}
source("run_methods_eval_fun")
source("simul_dat_fun")
```
# License
```{r}
Licensed under the GNU General Public License v3.0 (GPL-3.0);
you may not use this file except in compliance with the License.
```
