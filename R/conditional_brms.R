conditional_brms <- function(object, data, CI, ...){
  yhat <- fitted(
    object = object, newdata = data,
    re_formula = NA, scale = "linear",
    probs =  abs(c(0,1) - (1-CI)/2),
    ...)
  as.data.frame(yhat)
}
