#' @title summary() Method for Objects of Class 'bmamfit'
#'
#' @param object Objects of Class 'bmamfit'
#' @param plot.smooth Whether or not to plot bmam
#' @param digits number of signigicant digits, passing to \code{options()}.
#'   Defaults to 3
#' @param ... Additional arguments passed to \code{plot.bmam()}
#' @return a list containing estimates of parameters for smooth terms
#' @import dplyr
#' @import stringr
#' @export
#'
summary.bmamfit <- function(object, plot.smooth = FALSE, digits = 3, ...){

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
  names <- sapply(variables[grep(pattern = "^b_ *", variables)], function(name) str_replace(substring(name,3), "a_",""))
  if(!is.null(names)){
    # linear <- object$Posterior[names,]
    # if(length(names) == 1) linear <- t(as.data.table(linear))
    # linear_est <- as.data.table(do.call(rbind,
    #                                     apply(linear, 1, function(var) do.call("bsummary",c(list(x = var), object$Summary_para)))))
    linear_est <- object$Summary[Label %in% names,]
    # linear_est$Label <- names
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

    if(length(object$Bname) == 1) {
      names1 <- variables[grep(pattern = paste0("^bs_.*s",as.character(smterm[[2]][2])," *"), variables)]
      names2<- variables[grep(pattern = paste0("^s_.*s",as.character(smterm[[2]][2]),"_\\d *"), variables)]
    }else{
      names1 <- variables[grep(pattern = paste0("^bs_.*s",as.character(smterm[[2]][[i+1]][[2]])," *"), variables)]
      names2<- variables[grep(pattern = paste0("^s_.*s",as.character(smterm[[2]][[i+1]][[2]]),"_\\d *"), variables)]
    }

    variables[grep(pattern = paste0("^s_.*s *"), variables)]

    names <- c(names1,names2)
    smooth <- post[, names]
    smooth_est_i <- as.data.table(do.call(rbind,
                                          apply(smooth, 2, function(var) do.call("bsummary",c(list(x = var), object$Summary_para)))))
    smooth_est_i$Label <- object$Bname[[i]]
    smooth_est[[i]] <- select(smooth_est_i, Parameter = Label,M,Mdn,LL,UL,CI,CIType)
  }
  names(smooth_est) <- as.character(smterm[[2]][-1])


  ## linear term
  names <- sapply(variables[grep(pattern = "b_ *", variables)], function(name)substring(name,3))
  if(!is.null(names)){
    linear <- post[,variables[grep(pattern = "b_ *", variables)]]
    # if(length(names) == 1) linear <- t(as.data.table(linear))
    linear_est <- as.data.table(do.call(rbind,
                                        apply(linear, 2, function(var) do.call("bsummary",c(list(x = var), object$Summary_para)))))
    linear_est$Label <- names
    linear_est <- select(linear_est, Parameter = Label,M,Mdn,LL,UL,CI,CIType)

  }else{
    linear_est <- NULL
  }

  Conditional_Model <- list(Linear = linear_est,
                            Smooth = smooth_est)



  out <- list(BMAM = BMAM,
              Conditional_Model = Conditional_Model)
  options(digits = digits)
  cat("\n Marginal Model \n")
  print(out$BMAM$Linear)
  cat("\n Conditional Model \n")
  print(out$Conditional_Model$Linear)

  if(length(list(...)) != 0) plot.smooth <- TRUE
  if(plot.smooth){
    gg <- plot(object, ... )
    out$plot <- gg
  }
  invisible(out)
}


#' @title print() Method for Objects of Class 'bmamfit'
#'
#' @param object Objects of Class 'bmamfit'
#'
#' @return invisible() with printing
#' @export
print.bmamfit <- function(object,...){
  summary(object,...)
  invisible()
}


#' Bayesian Summary. Function in brmsmargins. Change default.
#' See ?brmsmargins::bsummary
#' @param x The posterior samples of a parameter
#' @param CI Width of the credible interval.
#' @param CIType Type of credible interval, passed on
#'   to the \code{\link[bayestestR]{ci}} function as the method for CIs.
#' @export
#' @seealso  \code{\link{brmsmargins::bsummary}}
bsummary <- function(x, CI = 0.95, CIType = "ETI", ...)
  brmsmargins::bsummary(x, CI = CI, CIType = CIType, ...)
