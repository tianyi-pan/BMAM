#' @title summary() Method for Objects of Class 'bmam'
#'
#' @param x Objects of Class 'bmam'
#' @param plot_smooth Whether or not to plot the smooth function 
#' @param ... Additional arguments passed to \code{plot.bmam()}
#' @return a list containing estimates of parameters for smooth terms
#' @import dplyr
summary.bmam <- function(x, plot_smooth = FALSE, ...){

  print(x$Family)
  print(x$Formula)
  

  ### 1. Marginal Model ##################
  ## smooth term
  smooth_est <- lapply(x$Bname, function(name){
    x$Summary[x$Summary$Label %in% name,] %>% select(Parameter = Label,M,Mdn,LL,UL,CI,CIType)
  })
  smterm <- brmsterms(x$Formula)$dpars$mu$sm # smooth term 
  names(smooth_est) <- as.character(smterm[[2]][-1])

  ## linear term
  variables <- variables(x$Conditional$Brms) ## get variables of linear term
  names <- sapply(variables[grep(pattern = "b_ *", variables)], function(name) substring(name,3))
  if(!is.null(names)){
    linear <- x$Posterior[names,]
    linear_est <- as.data.table(do.call(rbind,
                                        apply(linear, 1, function(var) do.call("bsummary",c(list(x = var), x$Summary_para)))))
    
    linear_est$Label <- names
    linear_est <- select(linear_est, Parameter = Label,M,Mdn,LL,UL,CI,CIType)
    
  }else{
    linear_est <- NULL
  }
  
  BMAM <- list(Linear = linear_est,
               Smooth = smooth_est)
  
  
  
  ### 2. Conditional Model #############
  post <- as_draws_matrix(x$Conditional$Brms)

  ## smooth term
  smooth_est <- vector("list",length = length(x$Bname))
  for(i in seq_along(x$Bname)){
    names1 <- variables[grep(pattern = paste0("bs_s",as.character(smterm[[2]][[i+1]][[2]])," *"), variables)]
    names2<- variables[grep(pattern = paste0("zs_",as.character(i),"_\\d *"), variables)]
    names <- c(names1,names2)
    smooth <- post[, names]
    smooth_est_i <- as.data.table(do.call(rbind,
                                    apply(smooth, 2, function(var) do.call("bsummary",c(list(x = var), x$Summary_para)))))
    smooth_est_i$Label <- x$Bname[[i]]
    smooth_est[[i]] <- select(smooth_est_i, Parameter = Label,M,Mdn,LL,UL,CI,CIType)
  }
  names(smooth_est) <- as.character(smterm[[2]][-1])
  
  
  ## linear term
  names <- sapply(variables[grep(pattern = "b_ *", variables)], function(name)substring(name,3))
  if(!is.null(names)){
    linear <- post[,variables[grep(pattern = "b_ *", variables)]]
    linear_est <- as.data.table(do.call(rbind,
                                        apply(linear, 2, function(var) do.call("bsummary",c(list(x = var), x$Summary_para)))))
    linear_est$Label <- names
    linear_est <- select(linear_est, Parameter = Label,M,Mdn,LL,UL,CI,CIType)
    
  }else{
    linear_est <- NULL
  }
  
  Conditional_Model <- list(Linear = linear_est,
               Smooth = smooth_est)
  
  
  
  out <- list(BMAM = BMAM,
              Conditional_Model = Conditional_Model)
  cat("\n Marginal Model \n")
  print(BMAM$Linear)
  cat("\n Conditional Model \n")
  print(Conditional_Model$Linear)
  
  if(!is.null(list(...))) plot_smooth <- TRUE
  if(plot_smooth){
    gg <- plot(x, ... )
    out$plot <- gg
  }
  invisible(out)
}


#' @title print() Method for Objects of Class 'bmam'
#'
#' @param x Objects of Class 'bmam'
#'
#' @return invisible() with printing
#' 
print.bmam <- function(x,...){
  summary(x,...)
  invisible()
}


