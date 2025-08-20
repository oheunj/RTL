adaptive_lasso = function(X, y, method = NULL, q, nfold){
  enet.b = cv.glmnet(X, y, alpha = 0.5, nfolds = nfold)
  abs.enet = abs(coef(enet.b, s = "lambda.min")[-1]) + 1 / nrow(X)
  weights = 1 / abs.enet
  alass = cv.glmnet(X, y, alpha = 1, penalty.factor = weights, nfolds = nfold)
  return(as.vector(coef(alass, s = "lambda.min")))
}
