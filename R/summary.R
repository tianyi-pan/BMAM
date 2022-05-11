#' @title summary() Method for Objects of Class 'bmam'
#'
#' @param x Objects of Class 'bmam'
#'
#' @return a list containing estimates of parameters for smooth terms
#' @import dplyr
summary.bmam <- function(x){
  smooth_est <- lapply(x$Bname, function(name){
    x$Summary[x$Summary$Label %in% name,] %>% select(Parameter = Label,M,Mdn,LL,UL,CI,CIType)
  })
  smterm <- brmsterms(x$Formula)$dpars$mu$sm # smooth term 
  names(smooth_est) <- as.character(smterm[[2]][-1])
  print(x$Family)
  print(x$Formula)
  return(smooth_est)
}


#' @title print() Method for Objects of Class 'bmam'
#'
#' @param x Objects of Class 'bmam'
#'
#' @return invisible() with printing
#' 
print.bmam <- function(x){
  print(x$Family)
  print(x$Formula)
  cat("\n")
  print(summary(x))
  invisible()
}


