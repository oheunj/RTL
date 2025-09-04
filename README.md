# Reluctant Transfer Learning in Penalized Regressions for Individualized Treatment Rules under Effect Heterogeneity

The full manuscript will be available soon:

* __Oh, E. J.,__ and Qian, M. (In Preparation). Reluctant transfer learning in penalized regressions for individualized treatment rules under effect heterogeneity.



# Files in this repository
The source code is currently provided in 'rtl_main_code.R'

# Installation
R is a statistical software program, and RStudio is a user interface for R. We recommend that users install both R and R Studio. Both R and RStudio are free and open source.

## Requirements
You may need to install the following dependencies first:
```{r}
library(MASS)
library(glmnet)
library(ranger)
```
Additionally, you need to save the following R codes from the 'functions' folder to your working directory and import these files:
```{r}
source("adaptive_lasso_with_enet_weights")
source("interaction_method")
source("get_value")
source("get_true_opt_value")
source("get_vsm")
```
In addition, for the additional benchmark model, you need to download the R codes from the Supplementary Materials of Wu & Yang (2023) and import the files:
```{r}
source("all")
source("DRITR")
source("obj_fun")
source("obj")
source("real_data_new")
```

# License
```{r}
Licensed under the GNU General Public License v3.0 (GPL-3.0);
you may not use this file except in compliance with the License.
```
