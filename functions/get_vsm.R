get_vsm = function(est, is_wuyang = FALSE){
  
  if(is_wuyang == TRUE){
    active = which(est != 0)
    total = length(est)
  }else{
    active = which(est[OA_start:OA_end] != 0)
    total = length(est[OA_start:OA_end])
  }
  
  TP = length(intersect(active, true_nonzero))
  FP = length(setdiff(active, true_nonzero))
  FN = length(setdiff(true_nonzero, active))
  TN = total - TP - FP - FN
  
  return(list(TP = TP, FP = FP, FN = FN, TN = TN))
}
