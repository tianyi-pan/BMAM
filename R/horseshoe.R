#' Title Define a model formula with horseshoe prior
#'
#' @param shrinkage.term term with horseshoe prior.
#' @param nonshrinkage.term terms with horseshoe prior.
#' @param ... arguments in horseshoe() see ?horseshoe
#'
#' @return a list containing formula used for brms and prior.
#' @import rlang
#' @import brms
#' @export
#'
brms.horseshoe <- function(shrinkage.term, nonshrinkage.term, ...){
  shrinkage.term.len <- length(shrinkage.term)
  nonshrinkage.term.len <- length(nonshrinkage.term)

  stopifnot(shrinkage.term.len == 3 || nonshrinkage.term.len == 3)

  ## bf for brms
  if(length(shrinkage.term) == 3){
    formula <- bf(eval(expr(!!shrinkage.term[[2]] ~ a + b)) , nl = TRUE) +
      lf(eval(expr(a ~ 0+!!shrinkage.term[[shrinkage.term.len]])), center = TRUE) +
      lf(eval(expr(b ~ !!nonshrinkage.term[[nonshrinkage.term.len]])), cmc = FALSE)
  }else{
    formula <- bf(eval(expr(!!nonshrinkage.term[[2]] ~ a + b)) , nl = TRUE) +
      lf(eval(expr(a ~ !!shrinkage.term[[shrinkage.term.len]])), center = TRUE) +
      lf(eval(expr(b ~ 0+!!nonshrinkage.term[[nonshrinkage.term.len]])), cmc = FALSE)
  }
  # if("prior" %in% names(list(...))) stop("horseshoe prior has been set")
  horseshoe_prior <- horseshoe(...)
  prior <- set_prior(horseshoe_prior, class = "b", nlpar = "a")

  return(list(formula = formula,
         prior = prior))
}
