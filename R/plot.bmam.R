#' @title plot() Method for Objects of Class 'bmamfit'
#'
#' @param object Objects of Class 'bmamfit'
#' @param compared.model Other model compared with BMAM.
#'                       supported models: 1. mam
#'                                         2. gam
#'                                         3. brms gam
#' @param display Whether or not to display the plot. Default: TRUE
#' @param conditional Whether or not to plot the conditional model. Default: TRUE
#' @param smooth.function A list. True values of the smooth functions.
#' @return a list containing ggplot objects.
#' @import ggplot2
#' @export
#'
plot.bmamfit <- function(object, compared.model, conditional = TRUE, display = TRUE, smooth.function){

  ### set up  #################
  theme_set(theme_bw())
  theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  preddat <- object$Preddat # pred data in object
  if(is.null(preddat$varname)) stop("Provide name of each smooth function in column varname in preddat")

  plot_var <- unique(preddat$varname) # smooth term


  ## check arguments
  if(object$Centered){
    if((!missingArg(compared.model)) || (!missingArg(smooth.function)))
      message("BMAM is centered.")
  }


  ## list to store ggplot object
  gg <- list(
    BMAM = vector("list", length = length(plot_var)),
    Conditional = list(Conditional.Model = vector("list", length = length(plot_var)),
                       Comparison = vector("list", length = length(plot_var))),
    Compared.Model = list(Compared.Model = vector("list", length = length(plot_var)),
                          Comparison = vector("list", length = length(plot_var)))
  )



  ### call ggplot to draw the plots #####################
  for(i in seq_along(plot_var)){

    ## i th smooth function

    var <- plot_var[i]
    index <- which(preddat$varname == var)

    ## marginal model
    dfplotM <- data.frame(x = preddat[[var]][index],
                          fitted = object$Predicted_Summary$M[index],
                          uci = object$Predicted_Summary$UL[index],
                          lci = object$Predicted_Summary$LL[index],
                          Method = rep(c("BMAM"), each = length(index)))
    gg$BMAM[[i]] <- ggplot(data=dfplotM,aes(x=x,y=fitted))+
      geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method),alpha=0.2, size = 0.3, colour = NA)+
      geom_line(aes(colour = Method), size = 1)+
      xlab(var)+
      ylab(expression(f~(X)))+
      scale_colour_manual(values = c("coral"), breaks = c("BMAM"))+
      scale_fill_manual(values = c("coral"), breaks = c("BMAM"))+
      ggtitle("BMAM")




    ## add true value to comparison
    if(!missingArg(smooth.function)){
      stopifnot(length(plot_var) == length(smooth.function))
      fun <- smooth.function[[i]]

      truevalue <- fun(preddat[[var]])

      dfplotT <- data.frame(x = preddat[[var]][index],
                            fitted = truevalue[index],
                            uci = truevalue[index],
                            lci = truevalue[index],
                            Method = rep(c("True Value"), each = length(index)))
      dfplotBoth <- rbind(dfplotT, dfplotM)
    }



    ## compared model
    if(!missingArg(compared.model)){
      if (!exists("dfplotBoth")) dfplotBoth <- dfplotM
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




      gg$Compared.Model$Compared.Model[[i]] <- ggplot(data=dfplotC,aes(x=x,y=fitted))+
        geom_ribbon(aes(ymin=lci,ymax=uci, fill=Method),alpha=0.2, size = 0.3, colour = NA)+
        geom_line(aes(colour = Method), size = 1)+
        xlab(var)+
        ylab(expression(f~(X)))+
        scale_colour_manual(values = c("deepskyblue"), breaks = c(dfplotC$Method[1]))+
        scale_fill_manual(values = c("deepskyblue"), breaks = c(dfplotC$Method[1]))+
        ggtitle(as.character(dfplotC$Method[1]))



      dfplotBoth <- rbind(dfplotBoth,dfplotC)
    }




    ## settings for ggplot
    if(!missingArg(smooth.function)){
      values <- c("black", "coral")
      breaks <- c("True Value", "BMAM")
    }else{
      values <- c("coral")
      breaks <- c("BMAM")
    }
    if(!missingArg(compared.model)){
      values <- c(values, "deepskyblue")
      breaks <- c(breaks, dfplotC$Method[1])
    }


    if(exists("dfplotBoth")){
      gg$Compared.Model$Comparison[[i]] <- ggplot(data=dfplotBoth,aes(x=x,y=fitted, group=Method))+
        geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method),alpha=0.2, size = 0.3, colour = NA)+
        geom_line(aes(colour = Method), size = 1)+
        xlab(var)+
        ylab(expression(f~(X))) +
        scale_colour_manual(values = values, breaks = breaks)+
        scale_fill_manual(values = values, breaks = breaks)+
        ggtitle("BMAM v.s. Compared.Model")
      remove(dfplotBoth)
    }

    ## conditional
    if(conditional){

      conditional_predicted <- object$Conditional$Predicted ## results of conditional model
      dfplotConditional <- data.frame(x = preddat[[var]][index],
                                      fitted = conditional_predicted$M[index],
                                      uci = conditional_predicted$UL[index],
                                      lci = conditional_predicted$LL[index],
                                      Method = rep(c("Conditional Model"), each = length(index)))
      dfplotBoth <- rbind(dfplotM, dfplotConditional)
      if(!missingArg(smooth.function)) dfplotBoth <- rbind(dfplotT, dfplotBoth)

      ## settings for ggplot
      if(!missingArg(smooth.function)){
        values <- c("black", "coral", "deepskyblue")
        breaks <- c("True Value", "BMAM", dfplotConditional$Method[1])
        labels <- c("True Value", "Marginal", "Conditional")
      }else{
        values <- c("coral", "deepskyblue")
        breaks <- c("BMAM", dfplotConditional$Method[1])
        labels <- c("Marginal", "Conditional")
      }

      gg$Conditional$Conditional.Model[[i]] <- ggplot(data=dfplotConditional,aes(x=x,y=fitted))+
        geom_ribbon(aes(ymin=lci,ymax=uci, fill=Method),alpha=0.2, size = 0.3, colour = NA)+
        geom_line(aes(colour = Method), size = 1)+
        xlab(var)+
        ylab(expression(f~(X)))+
        scale_colour_manual(values = c("deepskyblue"), breaks = c(dfplotConditional$Method[1]))+
        scale_fill_manual(values = c("deepskyblue"), breaks = c(dfplotConditional$Method[1]))+
        ggtitle(as.character(dfplotConditional$Method[1]))


      gg$Conditional$Comparison[[i]] <- ggplot(data=dfplotBoth,aes(x=x,y=fitted, group=Method))+
        geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method),alpha=0.2, size = 0.3, colour = NA)+
        geom_line(aes(colour = Method), size = 1)+
        xlab(var)+
        ylab(expression(f~(X))) +
        scale_colour_manual(values = values, breaks = breaks, labels = labels)+
        scale_fill_manual(values = values, breaks = breaks, labels = labels)+
        ggtitle("BMAM: Marginal v.s. Conditional Model")

      remove(dfplotBoth)
    }
  }



  ### display ##################################
  if(display){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))

    ## choose which plot will be shown
    if(!missingArg(compared.model)){
      gg_display <- list(gg$Compared.Model$Comparison)
      # if(conditional){
      #   gg_display <- list(gg$Conditional$Comparison, gg$Compared.Model$Comparison)
      # }else{
      #   gg_display <- list(gg$Compared.Model$Comparison)
      # }
    } else{

      if(conditional){
        gg_display <- list(gg$Conditional$Comparison)
      }else{
        gg_display <- list(gg$BMAM)
      }

    }
    for (gg_model in gg_display){
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
