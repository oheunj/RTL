interaction_adaptive_lasso = function(X_source, y_source, X_target, y_target, q, nfold){
  n_s = nrow(X_source)
  n_t = nrow(X_target)
  X = rbind(X_source, X_target)
  y = c(y_source, y_target)
  
  domain = c(rep(1, n_s), rep(0, n_t))
  X_int = X * domain
  X_full = cbind(X, domain, X_int)
  
  enet.b = cv.glmnet(X_full, y, alpha = 0.5, nfolds = nfold)
  abs.enet = abs(coef(enet.b, s = "lambda.min")[-1]) + 1 / nrow(X_full)
  weights = 1 / abs.enet
  alass = cv.glmnet(X_full, y, alpha = 1, penalty.factor = weights, nfolds = nfold)
  coef_hat = as.vector(coef(alass, s = "lambda.min"))
  beta_target = coef_hat[1:(2*q+2)]
  return(beta_target)
}