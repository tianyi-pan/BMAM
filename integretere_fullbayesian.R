integratere_fullbayesian <- function(d2, sd, L, J, post, usecoef,usegroup, k, yhat, backtrans.fun){
  
    re_list <- seq_along(d2)
    
    sample_i <- function(i){
      Z <- vector("list",length(re_list))
      for(re_id in re_list){
        sd_i <- sd[[re_id]][i,]
        L_i <- L[[re_id]][i,]
        yhat_i <- yhat[1,]
        Z[[re_id]] <- sapply(J[[re_id]], function(J.){
          index <- which(unique(J[[re_id]])==J.) 
          z <- .buildz(data = post, id = J., usecoef,usegroup)
          do.call("cbind",z) %*% t(sd_i %*% matrix(L_i, nrow = 2))
        })
      }
      Z.sum <- do.call("+",Z)
      colMeans(backtrans.fun(yhat_i + Z.sum))
    }
    
    t(sapply(seq_len(nrow(sd[[1]])), sample_i))
}