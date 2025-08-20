get_true_opt_value = function(beta_t, O_test){
  A0 = rep(0, nrow(O_test))
  A1 = rep(1, nrow(O_test))
  q0 = cbind(1, A0, O_test, A0 * O_test) %*% beta_t
  q1 = cbind(1, A1, O_test, A1 * O_test) %*% beta_t
  opt_A = ifelse(q1 > q0, 1, 0)
  OA_opt = O_test * matrix(opt_A, nrow = length(opt_A), ncol = ncol(O_test))
  Phi_opt = cbind(1, opt_A, O_test, OA_opt)
  return(mean(Phi_opt %*% beta_t))
}