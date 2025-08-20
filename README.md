# Leveraging Two-Phase Data for Improved Prediction of Survival Outcomes with Application to Nasopharyngeal Cancer
The full manuscript is available below:
* __Oh, E. J.,__ Ahn, S., Tham, T., and Qian, M. (2025). Leveraging two-phase data for improved prediction of survival outcomes with application to nasopharyngeal cancer. _Biometrics_, 81(2), ujaf080. [[link]](https://academic.oup.com/biometrics/article/doi/10.1093/biomtc/ujaf080/8174609?utm_source=authortollfreelink&utm_campaign=biometrics&utm_medium=email&guestAccessKey=3eeab5b5-3196-45c2-a2dc-26cd057f53c8) [[pdf]](https://oheunj.github.io/files/ujaf080.pdf) [[supp]](https://oheunj.github.io/files/BIOM_SuppMaterials_TwoPhaseSurv.pdf)

# Files in this repository
The source code is currently provided in 'tpsurv_main_code.R'

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
