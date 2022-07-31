#' @title generate predicted data
#'
#' @param object Objects of Class 'brms'
#' @param length number of observations in the generated data
#'
#' @return predicted data, used by \code{marginalcoef()}
#' @import stringi
generate_pred <- function(object, length = 100){
  mf <- model.frame(object) # data in object

  ## smooth term
  smterm <- brmsterms(object$formula)$dpars$mu$sm # smooth term
  stopifnot(!is.null(smterm)) # check

  smterm <- as.character(smterm)
  smterm <- unlist(stri_split_fixed(smterm, " "))
  smvariable <- as.list(unique(
    smterm[!is.na(stri_extract(smterm, regex = "[a-zA-Z0-9]"))]
  ))
  smvariable <- str_replace_all(smvariable, "s\\(","")
  smvariable <- str_replace_all(smvariable, "\\)","")

  # smvariable <- lapply(smterm[[2]][-1], function(term.) term.[[2]]) # variable for smooth term

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
  feterm <- brmsterms(object$formula)$dpars$mu$fe



  if(length(feterm) > 0){
    # if(feterm[[1]] == 1) feterm <- feterm[[-1]] # delete the intercept

    feterm <- as.character(feterm)
    feterm <- unlist(stri_split_fixed(feterm, " "))

    feterm <- as.list(unique(
      feterm[!is.na(stri_extract(feterm, regex = "[a-zA-Z]"))]
      ))
    feterm <- feterm[is.na(stri_extract(feterm, regex = ":"))]
    if(!is.null(feterm)){
    fevariable <- feterm
    fe_pred <- lapply(fevariable, function(var){
      x <- mf[[var]]
      if(is.factor(x)) x_pred <- factor(levels(x)[1],levels=levels(x))
      if(is.numeric(x)) x_pred <- 0
      rep(x_pred, length = length)
    })
    names(fe_pred) <- fevariable
    fe_pred <- do.call("cbind.data.frame", fe_pred)
    }else{
      fe_pred <- NULL
    }
  }else{
    fe_pred <- NULL
  }

  preddat <- lapply(smvariable, function(var){
    ## pred data for this variable
    sm_pred_design <- sm_pred; sm_pred_design[,] <- 0
    sm_pred_design[[var]] <- sm_pred[[var]]
    preddat_var <- sm_pred_design
    if(!is.null(fe_pred)) preddat_var <- cbind(sm_pred_design, fe_pred)
    preddat_var$varname <- as.character(var) # column to indicate the name of variable
    preddat_var
  })

  preddat <- do.call("rbind",preddat)

  return(preddat)
}
