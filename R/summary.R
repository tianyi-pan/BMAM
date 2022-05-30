#' @title summary() Method for Objects of Class 'bmam'
#'
#' @param object Objects of Class 'bmam'
#' @param plot_smooth Whether or not to plot the smooth function 
#' @param ... Additional arguments passed to \code{plot.bmam()}
#' @return a list containing estimates of parameters for smooth terms
#' @import dplyr
summary.bmam <- function(object, plot_smooth = FALSE, ...){
  
  print(object$Family)
  print(object$Formula)
  
  
  ### 1. Marginal Model ##################
  ## smooth term
  smooth_est <- lapply(object$Bname, function(name){
    object$Summary[object$Summary$Label %in% name,] %>% select(Parameter = Label,M,Mdn,LL,UL,CI,CIType)
  })
  smterm <- brmsterms(object$Formula)$dpars$mu$sm # smooth term 
  names(smooth_est) <- as.character(smterm[[2]][-1])
  
  ## linear term
  variables <- variables(object$Conditional$Brms) ## get variables of linear term
  names <- sapply(variables[grep(pattern = "b_ *", variables)], function(name) substring(name,3))
  if(!is.null(names)){
    linear <- object$Posterior[names,]
    linear_est <- as.data.table(do.call(rbind,
                                        apply(linear, 1, function(var) do.call("bsummary",c(list(object = var), object$Summary_para)))))
    
    linear_est$Label <- names
    linear_est <- select(linear_est, Parameter = Label,M,Mdn,LL,UL,CI,CIType)
    
  }else{
    linear_est <- NULL
  }
  
  BMAM <- list(Linear = linear_est,
               Smooth = smooth_est)
  
  
  
  ### 2. Conditional Model #############
  post <- as_draws_matrix(object$Conditional$Brms)
  
  ## smooth term
  smooth_est <- vector("list",length = length(object$Bname))
  for(i in seq_along(object$Bname)){
    names1 <- variables[grep(pattern = paste0("bs_s",as.character(smterm[[2]][[i+1]][[2]])," *"), variables)]
    names2<- variables[grep(pattern = paste0("^s_s",as.character(smterm[[2]][[i+1]][[2]]),"_\\d *"), variables)]
    
    variables[grep(pattern = paste0("^s_s *"), variables)]
    
    names <- c(names1,names2)
    smooth <- post[, names]
    smooth_est_i <- as.data.table(do.call(rbind,
                                          apply(smooth, 2, function(var) do.call("bsummary",c(list(object = var), object$Summary_para)))))
    smooth_est_i$Label <- object$Bname[[i]]
    smooth_est[[i]] <- select(smooth_est_i, Parameter = Label,M,Mdn,LL,UL,CI,CIType)
  }
  names(smooth_est) <- as.character(smterm[[2]][-1])
  
  
  ## linear term
  names <- sapply(variables[grep(pattern = "b_ *", variables)], function(name)substring(name,3))
  if(!is.null(names)){
    linear <- post[,variables[grep(pattern = "b_ *", variables)]]
    linear_est <- as.data.table(do.call(rbind,
                                        apply(linear, 2, function(var) do.call("bsummary",c(list(object = var), object$Summary_para)))))
    linear_est$Label <- names
    linear_est <- select(linear_est, Parameter = Label,M,Mdn,LL,UL,CI,CIType)
    
  }else{
    linear_est <- NULL
  }
  
  Conditional_Model <- list(Linear = linear_est,
                            Smooth = smooth_est)
  
  
  
  out <- list(BMAM = BMAM,
              Conditional_Model = Conditional_Model)
  options(digits=3)
  cat("\n Marginal Model \n")
  print(out$BMAM$Linear)
  cat("\n Conditional Model \n")
  print(out$Conditional_Model$Linear)
  
  if(length(list(...)) != 0) plot_smooth <- TRUE
  if(plot_smooth){
    gg <- plot(object, ... )
    out$plot <- gg
  }
  invisible(out)
}


#' @title print() Method for Objects of Class 'bmam'
#'
#' @param object Objects of Class 'bmam'
#'
#' @return invisible() with printing
#' 
print.bmam <- function(object,...){
  summary(object,...)
  invisible()
}


