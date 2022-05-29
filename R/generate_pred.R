#' @title generate predicted data
#' 
#' @param object Objects of Class 'brms'
#' @param length number of observations in the generated data
#'
#' @return predicted data, used by \code{marginalcoef()}
#'
generate_pred <- function(object, length = 100){
  mf <- model.frame(object) # data in object 
  
  ## smooth term
  smterm <- brmsterms(object$formula)$dpars$mu$sm # smooth term 
  stopifnot(!is.null(smterm)) # check
  smvariable <- lapply(smterm[[2]][-1], function(term.) term.[[2]]) # variable for smooth term
  
  sm_pred <- lapply(smvariable, function(var){
    x <- mf[[var]]
    x <- x[which(x!=0)] # remove 0
    x.max <- max(x)
    x.min <- min(x)
    seq(x.min, x.max, length = length)
  })
  names(sm_pred) <- smvariable
  sm_pred <- do.call("cbind.data.frame", sm_pred)
  
  ## fix effect term
  feterm <- brmsterms(object$formula)$dpars$mu$fe[[2]][-1]
  if(feterm[[1]] == 1) feterm <- feterm[[-1]] # delete the intercept
  if(!is.null(feterm)){
    fevariable <- feterm
    fe_pred <- lapply(fevariable, function(var){
      x <- mf[[var]]
      if(is.factor(x)) x_pred <- factor(levels(x)[1],levels=levels(x))
      if(is.numeric(x)) x_pred <- 0
      rep(x_pred, length = length)
    })
    names(fe_pred) <- list(fevariable)
    fe_pred <- do.call("cbind.data.frame", fe_pred)
  }else{
    fe_pred <- NULL
  }
  
  preddat <- lapply(smvariable, function(var){
    ## pred data for this variable
    sm_pred_design <- sm_pred; sm_pred_design[,] <- 0
    sm_pred_design[[var]] <- sm_pred[[var]]
    preddat_var <- cbind(sm_pred_design, fe_pred)
    preddat_var$varname <- as.character(var) # column to indicate the name of variable
    preddat_var
  })
  
  preddat <- do.call("rbind",preddat)
  
  return(preddat)
}
