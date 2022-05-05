integratere_fullbayesian <- function(d2, sd, L, J, post, usecoef,usegroup, k, yhat, backtrans.fun){
  
    re_list <- seq_along(d2)
    nsample <- nrow(post)
    sample_i <- function(i){
      Z <- vector("list",length(re_list))
      r <- vector("list",length(re_list))
      for(re_id in re_list){
        sd_i <- sd[[re_id]][i,]
        L_i <- L[[re_id]][i,]
        yhat_i <- yhat[i,]
        r[[re_id]] <- lapply(J[[re_id]], function(J.){
          index <- which(unique(J[[re_id]])==J.) 
          z <- .buildz(data = post, id = J., usecoef,usegroup)
          z_sam <- sample(seq_len(nsample), min(nsample,k))
          z <- lapply(z, function(z.) z.[z_sam])
          do.call("cbind",z) %*% diag(sd_i) %*% matrix(L_i, nrow = 2)
        })
        Z[[re_id]] <- mapply(r[[re_id]], as.data.frame(t(d2[[re_id]])), 
                             FUN = function(r_tmp, d2_tmp) r_tmp %*% d2_tmp) 
      }
      
      Z.sum <- do.call("+",Z)
      colMeans(backtrans.fun(sweep(Z.sum, 2, yhat_i, "+")))
    }
    
    t(sapply(seq_len(nrow(sd[[1]])), sample_i))
}