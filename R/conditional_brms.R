conditional_brms <- function(object, data, centered = FALSE, ...){
  yhat <- fitted(
    object = object, newdata = data,
    re_formula = NA, scale = "linear",
    summary=FALSE)
  
  
  if(centered){
    # projection
    ones <- matrix(rep(1,ncol(yhat)))
    H_matrix <- ones %*% solve(t(ones) %*% ones) %*% t(ones)
    M_matrix <- diag(1,ncol(yhat)) - H_matrix
    predicted <- M_matrix %*% t(yhat)
  }else{
    predicted <- t(yhat)
  }
  
  as.data.table(do.call(rbind, apply(predicted, 1, bsummary, ...)))
  
  
  ## method 2 minus colmean of designmatrix 
  # 
  # 
  # object <- restructure(object)
  # prep <- prepare_predictions(
  #   object, newdata = data, re_formula = NA,
  #   check_response = FALSE)
  # 
  # ## Xs: basis function for smooth term, without penalty (ncol = 1)
  # pred_Xs <- prep$dpars$mu$sm$fe$Xs
  # 
  # ## Zs: basis function for smooth term (ncol = k-2)
  # pred_Zs <- sapply(prep$dpars$mu$sm$re, function(re.)re.$Zs)
  # pred_Zs <- do.call(cbind, pred_Zs)
  # 
  # ## X: linear term, for example intercept + x1 + x2 + x1:x2
  # pred_X <- prep$dpars$mu$fe$X
  # 
  # pred_B <- cbind(pred_X, pred_Xs, pred_Zs)
  # 
  # if(centered) pred_B <- sweep(pred_B,2,colMeans(pred_B),'-') # Return centered smooths.
  # 
  # 
  # 
  # ## get beta from brms
  # ## names of variables
  # post <- as_draws_matrix(object) # posterior samples
  # variables <- variables(object) # variable names
  # 
  # namesX <- variables[grep(pattern = "b_ *", variables)]
  # namesXs <- variables[grep(pattern = paste0("bs_s *"), variables)]
  # namesZs <- variables[grep(pattern = paste0("^s_s *"), variables)]
  # 
  # names <- c(namesX, namesXs, namesZs)
  # 
  # ## beta
  # beta <- post[,names]
  # Predicted <- pred_B %*% t(beta)
  # 
  # ## return
  # tmp2 <- as.data.table(do.call(rbind, apply(Predicted, 1, bsummary, ...)))
}
