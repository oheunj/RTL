# the file is licensed under the GNU General Public License v3.0 (GPL-3.0)

# you need to set your own working directory
setwd("Specify-Your-Working-Directory")

# load packages
library(MASS)
library(glmnet)
library(ranger)

# parameters
# - n_sim: # of Monte Carlo simulations
# - n_s: sample size for source data
# - n_t: sample size for target data
# - n_test: sample size for test set
# - q: number of covariates
# - noise_sd: sd of noise
# - nfold: # of folds used for CV

# import several R functions specified in README.md (note: save those to your working directory)
source()
source()
source()
source()
source()
source()
source()
source()
source()
source()
source()

# generate AR(1) covariance matrix
make_ar_cov = function(q, rho = rho) {
  outer(1:q, 1:q, function(i, j) rho^abs(i - j))
}
Sigma_ar = make_ar_cov(q, rho = 0.3)

# weak dense signal
es = 0.5
beta1_dense = c(rep(1.25, q / 2), rep(0, q / 2))
beta2_dense = 4.5 * abs(es - 0.3) * beta1_dense
alpha1 = 1
alpha2_dense = es * sqrt(t(beta1_dense) %*% Sigma_ar %*% beta1_dense + t(beta2_dense) %*% Sigma_ar %*% beta2_dense + noise_sd^2) / 2
beta_s_dense = c(alpha1, alpha2_dense, beta1_dense, beta2_dense)

# sparse strong signal
beta1_sparse = c(seq(.1*q+.5, 1.5, length.out = 0.1*q), rep(0, 0.9*q))
beta2_sparse = 4.5 * abs(es - 0.3) * beta1_sparse
alpha2_sparse = es * sqrt(t(beta1_sparse) %*% Sigma_ar %*% beta1_sparse + t(beta2_sparse) %*% Sigma_ar %*% beta2_sparse + noise_sd^2) / 2
beta_s_sparse = c(alpha1, alpha2_sparse, beta1_sparse, beta2_sparse)

# theta
theta_list = list(
  {tmp = rep(0, q * 2 + 2); tmp[2] <- 2.5; tmp[2 + q + 1] = 2.5; tmp},
  {tmp = rep(0, q * 2 + 2); tmp[(2 + q + 1):(2 + 2 * q)] = c(seq(4.5, 1.5, length.out = q * 0.3), rep(0, q * 0.7)); tmp},
  {tmp = rep(0.5, q * 2 + 2); tmp[(2 + q + 1):(2 + 2 * q)] = c(rep(0.5, q * 0.7), rep(0, q * 0.3)); tmp},
  {tmp = rep(0, q * 2 + 2); tmp[2] <- 2.5; tmp[2 + q + 1] = 2.5; tmp},
  {tmp = rep(0, q * 2 + 2); tmp[(2 + q + 1):(2 + 2 * q)] = c(seq(4.5, 1.5, length.out = q * 0.3), rep(0, q * 0.7)); tmp},
  {tmp = rep(0.5, q * 2 + 2); tmp[(2 + q + 1):(2 + 2 * q)] = c(rep(0.5, q * 0.7), rep(0, q * 0.3)); tmp}
)

beta_s_list = list(beta_s_dense, beta_s_dense, beta_s_dense,
                   beta_s_sparse, beta_s_sparse, beta_s_sparse)

offset = 2 + q
OA_start = offset + 1
OA_end = OA_start + q - 1

get_beta_s_for_case <- function(case_idx) {
  return(beta_s_list[[case_idx]])
}





results = vector("list", length(theta_list))
names(results) = paste0("theta", 1:length(theta_list))

for (theta_case in seq_along(theta_list)){
  theta = theta_list[[theta_case]]
  beta_s = get_beta_s_for_case(theta_case)
  beta_t = beta_s + theta
  true_nonzero = which(beta_t[OA_start:OA_end] != 0)
  
  sim_result = matrix(NA, nrow = n_sim, ncol = 20)
  colnames(sim_result) = c("RTL_Value", "Homo_Value", "Pre_Value", "Post_Value", "Inter_Value", "WuYang_Value",
                           "RTL_RMSE", "Homo_RMSE", "Pre_RMSE", "Post_RMSE", "Inter_RMSE", "WuYang_RMSE",
                           "RTL_RMSEc", "Homo_RMSEc", "Pre_RMSEc", "Post_RMSEc", "Inter_RMSEc", "WuYang_RMSEc",
                           "Obs_Value", "Opt_Value")
  
  FP_rtl = FN_rtl = FP_homo = FN_homo = FP_pre = FN_pre = FP_post = FN_post = FP_inter = FN_inter = FP_wuyang = FN_wuyang = rep(NA, n_sim)
  TP_rtl = TN_rtl = TP_homo = TN_homo = TP_pre = TN_pre = TP_post = TN_post = TP_inter = TN_inter = TP_wuyang = TN_wuyang = rep(NA, n_sim)
  
  set.seed(1)
  O_test = mvrnorm(n_test, rep(0, q), Sigma_ar)
  A_test = rbinom(n_test, 1, 0.5)
  Phi_test = cbind(A_test, O_test, A_test * O_test)
  y_test = model.matrix(~Phi_test) %*% beta_t + rnorm(n_test, sd = noise_sd)
  
  for (sim in 1:n_sim){
    
    set.seed(sim)
    cat(sprintf("\rTheta case %d | Simulation %d / %d (%.0f%%)",
                theta_case, sim, n_sim, 100 * sim / n_sim))
    flush.console()
    
    O_s = mvrnorm(n_s, rep(0, q), Sigma_ar)
    A_s = rbinom(n_s, 1, 0.5)
    Phi_s = cbind(A_s, O_s, A_s * O_s)
    y_s = model.matrix(~Phi_s) %*% beta_s + rnorm(n_s, sd = noise_sd)
    
    O_t = mvrnorm(n_t, rep(0, q), Sigma_ar)
    A_t = rbinom(n_t, 1, 0.5)
    Phi_t = cbind(A_t, O_t, A_t * O_t)
    y_t = model.matrix(~Phi_t) %*% beta_t + rnorm(n_t, sd = noise_sd)
    
    obs_val = mean(model.matrix(~Phi_test) %*% beta_t)
    opt_val = get_true_opt_value(beta_t, O_test)
    
    # RTL
    beta_s_hat = adaptive_lasso(Phi_s, y_s, q = q)
    y_tilde = y_t - model.matrix(~Phi_t) %*% beta_s_hat
    theta_hat = adaptive_lasso(Phi_t, y_tilde, q = q)
    beta_t_hat = beta_s_hat + theta_hat
    RTL_val = get_value(beta_t_hat, O_test, beta_t)
    RTL_rmse = sqrt(mean((model.matrix(~Phi_test) %*% beta_t_hat - y_test)^2))
    
    true_contrast = 1 * beta_t[2] + O_test %*% beta_t[OA_start:OA_end]
    est_contrast_rtl = 1 * beta_t_hat[2] + O_test %*% beta_t_hat[OA_start:OA_end]
    RTL_rmsec = sqrt(mean((est_contrast_rtl - true_contrast)^2))
    
    # homo
    X_homo = rbind(Phi_s, Phi_t)
    y_homo = c(y_s, y_t)
    homo_hat = adaptive_lasso(X_homo, y_homo, q = q)
    Homo_val = get_value(homo_hat, O_test, beta_t)
    Homo_rmse = sqrt(mean((model.matrix(~Phi_test) %*% homo_hat - y_test)^2))
    
    # pre
    pre_hat = adaptive_lasso(Phi_s, y_s, q = q)
    Pre_val = get_value(pre_hat, O_test, beta_t)
    Pre_rmse = sqrt(mean((model.matrix(~Phi_test) %*% pre_hat - y_test)^2))
    
    # post
    post_hat = adaptive_lasso(Phi_t, y_t, q = q)
    Post_val = get_value(post_hat, O_test, beta_t)
    Post_rmse = sqrt(mean((model.matrix(~Phi_test) %*% post_hat - y_test)^2))
    
    # interaction
    inter_hat = interaction_adaptive_lasso(Phi_s, y_s, Phi_t, y_t, q = q)
    Inter_val = get_value(inter_hat, O_test, beta_t)
    Inter_rmse = sqrt(mean((model.matrix(~Phi_test) %*% inter_hat - y_test)^2))
    
    # Wu & Yang (2023)
    y.in.rct = y_s
    x.in.rct = O_s
    a.in.rct = A_s
    x.in.rwe = O_t
    
    w.tilde = learn_weights(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, w.method = 4, misspecify = FALSE, N = NA)
    
    loc.a1 = which(a.in.rct == 1)
    loc.a0 = which(a.in.rct == 0)
    dat = data.frame(y.in.rct = y.in.rct, x.in.rct)
    dat.rwe = data.frame(x.in.rwe)
    fit1 = ranger(y.in.rct ~ ., data = dat[loc.a1,], case.weights = w.tilde[loc.a1])
    reg.mu1.rf = predict(fit1, data = dat.rwe)$predictions
    fit0 = ranger(y.in.rct ~ ., data = dat[loc.a0,], case.weights = w.tilde[loc.a0])
    reg.mu0.rf = predict(fit0, data = dat.rwe)$predictions
    
    x.in.rct = data.matrix(x.in.rct)
    x.in.rwe = data.matrix(x.in.rwe)
    suppressWarnings({
      res = lalonde_learn(nfold = nfold, reg.mu0.rf, reg.mu1.rf, w.tilde,
                          a.in.rct, x.in.rct, y.in.rct, x.in.rwe)
    })
    a.test = A_test
    x.test = O_test
    y.test = y_test
    
    loc.a1 = which(a.test == 1)
    loc.a0 = which(a.test == 0)
    dat = data.frame(y.test = y.test, x.test)
    dat.rwe = data.frame(x.test)
    fit1 = ranger(y.test ~ ., data = dat[loc.a1,])
    reg.mu1.rf = predict(fit1, data = dat.rwe)$predictions
    fit0 = ranger(y.test ~ ., data = dat[loc.a0,])
    reg.mu0.rf = predict(fit0, data = dat.rwe)$predictions
    
    action.weight = as.numeric(as.matrix(x.test) %*% res$beta.hat > 0)
    weighted.value = reg.mu0.rf
    weighted.value[action.weight == 1] = reg.mu1.rf[action.weight == 1]
    WuYang_val = mean(weighted.value)
    WuYang_rmse = sqrt(mean((predict(ranger(y.test ~ ., data = dat), data = dat)$predictions - y_test)^2))
    
    WuYang_rmsec = NA
    
    est_contrast_homo = 1 * homo_hat[2] + O_test %*% homo_hat[OA_start:OA_end]
    Homo_rmsec = sqrt(mean((est_contrast_homo - true_contrast)^2))
    
    est_contrast_pre = 1 * pre_hat[2] + O_test %*% pre_hat[OA_start:OA_end]
    Pre_rmsec = sqrt(mean((est_contrast_pre - true_contrast)^2))
    
    est_contrast_post = 1 * post_hat[2] + O_test %*% post_hat[OA_start:OA_end]
    Post_rmsec = sqrt(mean((est_contrast_post - true_contrast)^2))
    
    est_contrast_inter = 1 * inter_hat[2] + O_test %*% inter_hat[OA_start:OA_end]
    Inter_rmsec = sqrt(mean((est_contrast_inter - true_contrast)^2))
    
    sim_result[sim, ] = c(RTL_val, Homo_val, Pre_val, Post_val, Inter_val, WuYang_val, 
                          RTL_rmse, Homo_rmse, Pre_rmse, Post_rmse, Inter_rmse, WuYang_rmse,
                          RTL_rmsec, Homo_rmsec, Pre_rmsec, Post_rmsec, Inter_rmsec, WuYang_rmsec,
                          obs_val, opt_val)
    
    FP_rtl[sim] = get_vsm(theta_hat)$FP; FN_rtl[sim] = get_vsm(theta_hat)$FN
    FP_homo[sim] = get_vsm(homo_hat)$FP; FN_homo[sim] = get_vsm(homo_hat)$FN
    FP_pre[sim] = get_vsm(pre_hat)$FP; FN_pre[sim] = get_vsm(pre_hat)$FN
    FP_post[sim] = get_vsm(post_hat)$FP; FN_post[sim] = get_vsm(post_hat)$FN
    FP_inter[sim] = get_vsm(inter_hat)$FP; FN_inter[sim] = get_vsm(inter_hat)$FN
    FP_wuyang[sim] = get_vsm(res$beta.hat, is_wuyang = TRUE)$FP; FN_wuyang[sim] = get_vsm(res$beta.hat, is_wuyang = TRUE)$FN
    
    TP_rtl[sim] = get_vsm(theta_hat)$TP; TN_rtl[sim] = get_vsm(theta_hat)$TN
    TP_homo[sim] = get_vsm(homo_hat)$TP; TN_homo[sim] = get_vsm(homo_hat)$TN
    TP_pre[sim] = get_vsm(pre_hat)$TP; TN_pre[sim] = get_vsm(pre_hat)$TN
    TP_post[sim] = get_vsm(post_hat)$TP; TN_post[sim] = get_vsm(post_hat)$TN
    TP_inter[sim] = get_vsm(inter_hat)$TP; TN_inter[sim] = get_vsm(inter_hat)$TN
    TP_wuyang[sim] = get_vsm(res$beta.hat, is_wuyang = TRUE)$TP; TN_wuyang[sim] = get_vsm(res$beta.hat, is_wuyang = TRUE)$TN
  }
  
  results[[theta_case]] = list(
    coeff = data.frame(beta_s = beta_s, beta_t = beta_t, theta = theta),
    summary = sim_result,
    FP_vec = list(RTL = FP_rtl, Homo = FP_homo, Pre = FP_pre, Post = FP_post, Inter = FP_inter, WuYang = FP_wuyang),
    FN_vec = list(RTL = FN_rtl, Homo = FN_homo, Pre = FN_pre, Post = FN_post, Inter = FN_inter, WuYang = FN_wuyang),
    TP_vec = list(RTL = TP_rtl, Homo = TP_homo, Pre = TP_pre, Post = TP_post, Inter = TP_inter, WuYang = TP_wuyang),
    TN_vec = list(RTL = TN_rtl, Homo = TN_homo, Pre = TN_pre, Post = TN_post, Inter = TN_inter, WuYang = TN_wuyang)
  )
  cat("\n")

  wd = getwd() # create 'Output' folder within your working directory to save the output below
  filename = paste0(wd,"/Output/ns",n_s,
                    "_nt",n_t,
                    "_q",q,
                    "_theta",theta_case,
                    ".Rdata")

  output = results[[theta_case]]
  save(output, file = filename)

}






