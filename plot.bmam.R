plot.bmam <- function(object, x, n.variable){
  length(object$Predicted_Summary$M) / n.variable
  dfplot <- data.frame(x = beaverspred$time[101:200],
                       fitted = mc$Predicted_Summary$M[101:200],
                       uci = mc$Predicted_Summary$UL[101:200],
                       lci = mc$Predicted_Summary$LL[101:200],
                       Method = rep(c("BMAM"), each = 100))
  
}