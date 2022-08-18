#' @title summary() Method for Objects of Class 'bmamfit'
#'
#' @param object Objects of Class 'bmamfit'
#' @param plot.smooth Whether or not to plot bmam
#' @param ... Additional arguments passed to \code{plot.bmam()}
#' @return a list containing estimates of parameters for smooth terms
#' @import dplyr
#' @import stringr
#' @export
#'
summary.bmamfit <- function(object, plot.smooth = FALSE, ...){

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
    linear <- object$Posterior[names,]
    linear_est <- as.data.table(do.call(rbind,
                                        apply(linear, 1, function(var) do.call("bsummary",c(list(x = var), object$Summary_para)))))

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
  options(digits=3)
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


#' Personal Preference Based Bayesian Summary. Function in brmsmargins but change default CIType to "ETI".
#'
#' Returns a summary of a posterior distribution for a single
#' parameter / value. It is based on personal preference. Notably, it does not
#' only use \code{bayestestR::describe_posterior}, an excellent function,
#' because of the desire to also describe the percentage of the full posterior
#' distribution that is at or exceeding the value of a
#' Minimally Important Difference (MID). MIDs are used in clinical studies with outcome
#' measures where there are pre-defined differences that are considered clinically
#' important, which is distinct from the ROPE or general credible intervals capturing
#' uncertainty.
#'
#' @param x The posterior distribution of a parameter
#' @param CI A numeric value indicating the desired width of the credible interval.
#'   Defaults to \code{0.99} currently, but this is subject to change.
#'   a 99% interval was chosen as the default as there have been recent arguments
#'   made in the realm of meta science that there are, essentially, too many
#'   false positives and that many of the \dQuote{findings} in science are not able
#'   to be replicated.
#'   In any case, users should ideally specify a desired CI width, and not rely on
#'   defaults.
#' @param CIType A character string indicating the type of credible interval, passed on
#'   to the \code{\link[bayestestR]{ci}} function as the method for CIs.
#' @param ROPE Either left as \code{NULL}, the default, or a numeric vector of
#'   length 2, specifying the lower and upper thresholds for the
#'   Region of Practical Equivalence (ROPE).
#' @param MID Either left as \code{NULL}, the default, or a numeric vector of
#'   length 2, specifying the lower and upper thresholds for a
#'   Minimally Important Difference (MID). Unlike the ROPE, percentages for
#'   the MID are calculated as at or exceeding the bounds specified by this
#'   argument, whereas the ROPE is the percentage of the posterior at or inside
#'   the bounds specified.
#' @return A \code{data.table} with the mean, \code{M}
#' \describe{
#'   \item{M}{the mean of the posterior samples}
#'   \item{Mdn}{the median of the posterior samples}
#'   \item{LL}{the lower limit of the credible interval}
#'   \item{UL}{the upper limit of the credible interval}
#'   \item{PercentROPE}{the percentage of posterior samples falling into the ROPE}
#'   \item{PercentMID}{the percentage of posterior samples falling at or beyond the MID}
#'   \item{CI}{the width of the credible interval used}
#'   \item{CIType}{the type of credible interval used (e.g., highest density interval)}
#'   \item{ROPE}{a label describing the values included in the ROPE}
#'   \item{MID}{a label describing the values included in the MID}
#' }
#' @export
#' @importFrom bayestestR ci
#' @importFrom data.table data.table
#' @importFrom stats median
#' @references
#' Kruschke, J. K. (2018).
#' \doi{10.1177/2515245918771304}
#' \dQuote{Rejecting or accepting parameter values in Bayesian estimation}
#' @examples
#'
#' bsummary(rnorm(1000))
#'
#' bsummary(rnorm(1000), ROPE = c(-.5, .5), MID = c(-1, 1))
bsummary <- function(x, CI = 0.99, CIType = "ETI", ROPE = NULL, MID = NULL) {
  if (isFALSE(is.numeric(x))) {
    stop(sprintf("to be summarized x must be numeric, but %s class was found",
                 paste(class(x), collapse = "; ")))
  }
  ropes <- brmsmargins:::.percent(x, window = ROPE, within = TRUE)
  mids <- brmsmargins:::.percent(x, window = MID, within = FALSE)

  m <- mean(x, na.rm = TRUE)
  mdn <- median(x, na.rm = TRUE)
  cis <- bayestestR::ci(x, ci = CI, method = CIType)
  out <- data.table(
    M = as.numeric(m),
    Mdn = as.numeric(mdn),
    LL = as.numeric(cis$CI_low),
    UL = as.numeric(cis$CI_high),
    PercentROPE = as.numeric(ropes$Percent),
    PercentMID = as.numeric(mids$Percent),
    CI = as.numeric(CI),
    CIType = CIType,
    ROPE = ropes$Label,
    MID = mids$Label)

  return(out)
}



