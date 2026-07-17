# Reluctant Transfer Learning in Penalized Regressions for Individualized Treatment Rules under Effect Heterogeneity
<!-- 
The full manuscript will be available soon:
-->

The full manuscript is available below:
* __Oh, E. J.,__ and Qian, M. (2026). Reluctant transfer learning in penalized regressions for individualized treatment rules under effect heterogeneity, _Statistics in Medicine_, 45(15-17), e70671. [[link]](https://onlinelibrary.wiley.com/doi/10.1002/sim.70671) [[pdf]](https://pmc.ncbi.nlm.nih.gov/articles/PMC13349454/pdf/SIM-45-0.pdf) [[supp]](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fsim.70671&file=sim70671-sup-0001-Supinfo.pdf.pdf) [[code]](https://github.com/oheunj/RTL)



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
source("get_value")
source("get_true_opt_value")
source("get_vsm")
```
In addition, for the benchmark model (ITL), you need to download the R codes from the Supplementary Materials of Wu & Yang (2023) to your working directory and import these files:
```{r}
source("all")
source("DRITR")
source("obj_fun")
source("obj")
source("real_data_new")
```
For another benchmark model (TransLasso), you need to download the R codes from https://github.com/saili0103/TransLasso to your working directory and import the following:
```{r}
source("TransLasso-functions")
```

# License
```{r}
Licensed under the GNU General Public License v3.0 (GPL-3.0);
you may not use this file except in compliance with the License.
```
