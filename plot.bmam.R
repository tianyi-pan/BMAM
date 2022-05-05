#' @title plot BMAM model
#' 
#' @param object a bmam model object
#' @param showplot TRUE: display the plot.  
#'
#' @return a list containing ggplot objects. 
#' @import ggplot2
#'
#' TO DO: change object as a S3 class!
plot.bmam <- function(object, showplot = TRUE){
  preddat <- object$Preddat # pred data in object 
  
  plot_var <- unique(preddat$varname)
  gg <- vector("list", length = length(plot_var))
  
  for(i in seq_along(plot_var)){
    var <- plot_var[i]
    index <- which(preddat$varname == var)
    dfplot <- data.frame(x = preddat[[var]][index],
                         fitted = mc$Predicted_Summary$M[index],
                         uci = mc$Predicted_Summary$UL[index],
                         lci = mc$Predicted_Summary$LL[index],
                         Method = rep(c("BMAM"), each = length(index)))
    gg[[i]] <- ggplot(data=dfplot,aes(x=x,y=fitted))+
      geom_ribbon(aes(ymin=lci,ymax=uci),alpha=0.2, size = 0.3)+
      geom_line(aes(colour = Method), size = 1)+
      # ylim(myrange)+
      xlab(var)+
      ylab(expression(f~(X)))
      # ggtitle("Female") +
      # theme(text = element_text(size = GGPLOTTEXTSIZE))
  }
  if(showplot){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
    for (gg_i in gg) {
      dev.hold()
      plot(gg_i)
      dev.flush()      
    }
    devAskNewPage(oask)
  }
  # return(gg)
  # else{
  #   return(gg)
  # }
  invisible(gg)
  # gg
}
