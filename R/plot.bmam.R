#' @title plot() Method for Objects of Class 'bmam'
#' 
#' @param object Objects of Class 'bmam'
#' @param compared.model Other model compared with BMAM. 
#'                       supported models: 1. mam
#'                                         2. gam
#'                                         3. brms gam
#' @param display TRUE: display the plot.  
#'
#' @return a list containing ggplot objects. 
#' @import ggplot2
#'
plot.bmam <- function(object, compared.model, display = TRUE){

  preddat <- object$Preddat # pred data in object 
  
  plot_var <- unique(preddat$varname) # smooth term
  


  gg <- list( # list to store ggplot object
    BMAM = vector("list", length = length(plot_var)),
    Compared.Model = vector("list", length = length(plot_var)),
    Both = vector("list", length = length(plot_var))
  )

  

  for(i in seq_along(plot_var)){
    var <- plot_var[i]
    index <- which(preddat$varname == var)
    

    dfplotM <- data.frame(x = preddat[[var]][index],
                          fitted = mc$Predicted_Summary$M[index],
                          uci = mc$Predicted_Summary$UL[index],
                          lci = mc$Predicted_Summary$LL[index],
                          Method = rep(c("BMAM"), each = length(index)))
    gg$BMAM[[i]] <- ggplot(data=dfplotM,aes(x=x,y=fitted))+
      geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method, colour= Method),alpha=0.2, size = 0.3)+
      geom_line(aes(colour = Method), size = 1)+
      xlab(var)+
      ylab(expression(f~(X)))+
      scale_colour_manual(values = c("coral"), breaks = c("BMAM"))+
      scale_fill_manual(values = c("coral"), breaks = c("BMAM"))+
      ggtitle("BMAM")
    
    
    
    if(!missingArg(compared.model)){
      ## check whether the model is supported. 
      
      if(!is.null(compared.model$mam)){ 
        ## MAM model 
        mam.uci <- compared.model$mam$fitted+1.96*compared.model$mam$fitted_se
        mam.lci <- compared.model$mam$fitted-1.96*compared.model$mam$fitted_se
        
        dfplotC <- data.frame(x = preddat[[var]][index],
                              fitted = compared.model$mam$fitted[index],
                              uci = mam.uci[index],
                              lci = mam.lci[index],
                              Method = rep(c("MAM"), each = length(index)))
      
      }else if(class(compared.model)[1] == "gam"){
        ## GAM model
        gamfit <- predict(compared.model,newdata=preddat,se.fit=TRUE)
        gamfitted <- gamfit$fit
        gamfitted_se <- gamfit$se.fit
        
        gam.uci <- gamfitted+1.96*gamfitted_se
        gam.lci <- gamfitted-1.96*gamfitted_se
        
        dfplotC <- data.frame(x = preddat[[var]][index],
                              fitted = gamfitted[index],
                              uci = gam.uci[index],
                              lci = gam.lci[index],
                              Method = rep(c("GAM"), each = length(index)))
      
      }else if(class(compared.model)[1] == "brmsfit"){
        ## brms model 
        brms_posterior <- fitted(compared.model, newdata = preddat, 
                                  summary = FALSE, scale = "linear")
        brms_summary <- as.data.table(do.call(rbind,
                                              apply(brms_posterior, 2, function(post.) do.call("bsummary",c(list(x = post.), object$Summary_para)))))
        dfplotC <- data.frame(x = preddat[[var]][index],
                              fitted = brms_summary$M[index],
                              uci = brms_summary$UL[index],
                              lci = brms_summary$LL[index],
                              Method = rep(c("BGAM"), each = length(index)))
      
      }else{
        stop("the model is not supported yet.")
      }
      
      gg$Compared.Model[[i]] <- ggplot(data=dfplotC,aes(x=x,y=fitted))+
        geom_ribbon(aes(ymin=lci,ymax=uci, fill=Method, colour= Method),alpha=0.2, size = 0.3)+
        geom_line(aes(colour = Method), size = 1)+
        xlab(var)+
        ylab(expression(f~(X)))+
        scale_colour_manual(values = c("deepskyblue"), breaks = c(dfplotC$Method[1]))+
        scale_fill_manual(values = c("deepskyblue"), breaks = c(dfplotC$Method[1]))+
        ggtitle(as.character(dfplotC$Method[1]))
      
      dfplotCompare <- rbind(dfplotM,dfplotC)
      gg$Both[[i]] <- ggplot(data=dfplotCompare,aes(x=x,y=fitted, group=Method))+
        geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method, colour= Method),alpha=0.2, size = 0.3)+
        geom_line(aes(colour = Method), size = 1)+
        xlab(var)+
        ylab(expression(f~(X))) + 
        scale_colour_manual(values = c("coral", "deepskyblue"), breaks = c("BMAM",dfplotC$Method[1]))+
        scale_fill_manual(values = c("coral", "deepskyblue"), breaks = c("BMAM",dfplotC$Method[1]))+
        ggtitle("Comparision")
    }
  }
  
  if(display){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
    for (gg_model in gg){
      for (gg_i in gg_model) {
        if(!is.null(gg_i)){
          dev.hold()
          plot(gg_i)
          dev.flush()               
        }
      }    
    }

    devAskNewPage(oask)
  }
  invisible(gg)
}
