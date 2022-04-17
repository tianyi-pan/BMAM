#' Marginal Coefficients from a 'brms' Model
#' modified from brmsmargins::marginalcoef
#' Calculate marginal coefficients from a \code{brms}
#' generalized linear mixed model using the method proposed by Hedeker (2018).
#'
#' @param object A fitted brms model object that includes random effects. Required.
#' @param summarize A logical value, whether or not to
#'   calculate summaries of the posterior predictions.
#'   Defaults to \code{TRUE}.
#' @param posterior A logical value whether or not to
#'   save and return the posterior samples. Defaults
#'   to \code{FALSE} as the assumption is a typical
#'   use case is to return the summaries only.
#' @param index An optional integer vector, giving the posterior draws
#'   to be used in the calculations. If omitted, defaults to all
#'   posterior draws.
#' @param backtrans A character string indicating the type of
#'   back transformation to be applied. Can be one of
#'   \dQuote{response} meaning to use the response scale,
#'   \dQuote{linear} or \dQuote{identity} meaning to use the linear predictor scale,
#'   or a specific back transformation desired, from a possible list of
#'   \dQuote{invlogit}, \dQuote{exp}, \dQuote{square}, or \dQuote{inverse}.
#'   Custom back transformations should only be needed if, for example,
#'   the outcome variable was transformed prior to fitting the model.
#'   default: \dQuote{response}
#' @param k An integer providing the number of random draws to use for
#'   integrating out the random effects. Only relevant when \code{effects}
#'   is \dQuote{integrateoutRE}.
#' @param ... Additional arguments passed to \code{fitted()}
#' @return A list with \code{Summary} and \code{Posterior}.
#'   Some of these may be \code{NULL} depending on the arguments used.
#' @references
#' Hedeker, D., du Toit, S. H., Demirtas, H. & Gibbons, R. D. (2018)
#' \doi{10.1111/biom.12707}
#' \dQuote{A note on marginalization of regression parameters from mixed models of binary outcomes}
#' @importFrom data.table as.data.table
#' @importFrom stats formula
#' @importFrom posterior as_draws_df ndraws
#' @importFrom brms make_standata
#' @importFrom methods missingArg
#' @export
marginalcoef <- function(object, preddat, summarize = TRUE, posterior = FALSE, index,
                         backtrans = c("response", "linear", "identity",
                                       "invlogit", "exp", "square", "inverse"),
                         centered = FALSE,
                         k = 100L, ...) {
  ## check smooth term
  smooth <- !is.null(brmsterms(object$formula)$dpars$mu$sm)
  
  ## checks and assertions
  brmsmargins:::.assertbrmsfit(object)

  if (isFALSE(brmsmargins:::is.random(object))) {
    stop("object must have random effects to use marginalcoef()")
  }

  ## assert the assumed family / distribution is a supported one
  brmsmargins:::.assertfamily(object)
  ## assert the link function used is a supported one
  brmsmargins:::.assertlink(object)
  ## assert that all random effects in the model are Gaussian
  brmsmargins:::.assertgaussian(object)

  # number of draws. generate the corresponding index: 1 to the number.
  if (isTRUE(missingArg(index))) {
    index <- seq_len(ndraws(object))
  }
  
  # Convert a Link Function Name to a List
  backtrans <- match.arg(backtrans)
  links <- brmsmargins:::.links(
    link = brmsmargins:::.extractlink(object, NULL),
    effects = "integrateoutRE", backtrans = backtrans)
  
  # get the dataset in the model.
  mf <- model.frame(object)
  
  # see prediction function. 
  # get the predictive posterior distribution of \mu 
  lambda <- brmsmargins::prediction(
    object, data = mf,
    summarize = FALSE, posterior = TRUE, index = index,
    effects = "integrateoutRE", backtrans = backtrans,
    k = k, raw = TRUE)
  
  # convert \mu to \eta = g(\mu). links$fun: link function
  # get lambda_pa in equation 6 in Hedeker's paper
  y <- links$fun(t(lambda$Posterior))
  ## dim(y) <- number of observations * number of draws(in MCMC)
  
  
  # design matrix in dataset
  standata <- make_standata(formula(object), data = mf)
  
  if(smooth){
    data_names <- names(standata) # get the names of data
    
    # ## Z: design matrix for random effect
    # Z_name <- data_names[grep(pattern = "Z_\\d_\\d", data_names)]
    # Z <- do.call(cbind,standata[Z_name])
    
    ## Zs: basis function for smooth term (ncol = k-2)
    Zs_name <- data_names[grep(pattern = "Zs_\\d_\\d", data_names)]
    Zs <- do.call(cbind, standata[Zs_name])
    ## set names for Zs
    Zs_name_list <- mapply(function(i,j)paste(j,seq_len(i),sep="_alpha_"), 
                           lapply(standata[Zs_name],ncol), # number of basis function
                           Zs_name,
                           SIMPLIFY = FALSE)
    colnames(Zs) <- as.character(unlist(Zs_name_list))
    
    
    ## Xs: basis function for smooth term, without penalty (ncol = 1)
    Xs_name <- data_names[grep(pattern = "Xs", data_names)]
    Xs <- do.call(cbind, standata[Xs_name])
    ## set names for Xs
    Xs_name_list <- mapply(function(i,j)paste(j,seq_len(i),"alpha",sep="_"), 
                           lapply(standata[Xs_name],ncol), # number of basis function
                           Xs_name)
    colnames(Xs) <- as.character(unlist(Xs_name_list))
    
    
    ## X: linear term, for example intercept + x1 + x2 + x1:x2
    X <- standata$X
    
    ## design matrix for 
    B <- cbind(X, Xs, Zs)
  }else{
    B <- standata$X # in GLMM, B=X, linear term, for example 1 + x1 + x2 + x1:x2
  }
  

  # calculate beta_pa using equation 6 in Hedeker's paper
  # heagerty's method, solve 
  beta <- lmcpp(B, y)
  rownames(beta) <- colnames(B)
  out <- list(
    Summary = NULL,
    Posterior = NULL,
    DesignMatrix = NULL,
    Bname = NULL,
    Smooth = NULL,
    Predicted_Summary = NULL,
    Predicted = NULL)

  if (isTRUE(summarize)) { # does not work for 
    out$Summary <- as.data.table(do.call(rbind, apply(beta, 1, bsummary, ...)))
    out$Summary[, Label := colnames(B)]
  }

  if(smooth){
    Bname <- mapply(c, as.list(Xs_name_list), Zs_name_list) # variable names for basis function. 
    ## estimate the smooth function by B*\alpha
    if(!is.list(Bname)){ # convert Bname to a list
      Bname <- lapply(seq_len(ncol(Bname)), function(i) Bname[,i])
    }
    smooth_estimates <- lapply(Bname, function(name)B[,name] %*% beta[name, ])
    names(smooth_estimates) <- paste0("f",seq_along(Bname))
    
    out$DesignMatrix <- B
    out$Bname <- Bname
    out$Smooth <- smooth_estimates
    
    if(!missingArg(preddat)){

      # fe_formula <- object$formula$formula[[3]][[2]][[2]]
      # fe_formula <- as.formula(paste0("y~", paste(as.character(fe_formula)[as.character(fe_formula)!="+"], collapse='+')))
      # 

      object <- restructure(object)
      prep <- prepare_predictions(object, newdata = preddat,check_response = FALSE, re_formula = NA)
      
      pred_Xs <- prep$dpars$mu$sm$fe$Xs
      
      pred_Zs <- sapply(prep$dpars$mu$sm$re, function(re.)re.$Zs)
      pred_Zs <- do.call(cbind, pred_Zs)
      
      pred_X <- prep$dpars$mu$fe$X
      
      
      # pred_standata <- make_standata(fe_formula, data = preddat)
      
      # pred_data_names <- names(pred_standata) # get the names of data

      # pred_Zs_name <- pred_data_names[grep(pattern = "Zs_\\d_\\d", pred_data_names)]
      # pred_Zs <- do.call(cbind, pred_standata[pred_Zs_name])
      
      ## set names for Zs
      # Zs_name_list <- mapply(function(i,j) paste(j,seq_len(i),sep="_alpha_"), 
      #                        lapply(pred_standata[pred_Zs_name],ncol), # number of basis function
      #                        pred_Zs_name,
      #                        SIMPLIFY = FALSE)
      # colnames(pred_Zs) <- as.character(unlist(Zs_name_list))
      
      
      ## Xs: basis function for smooth term, without penalty (ncol = 1)
      # pred_Xs_name <- pred_data_names[grep(pattern = "Xs", pred_data_names)]
      # pred_Xs <- do.call(cbind, pred_standata[pred_Xs_name])
      # ## set names for Xs
      # Xs_name_list <- mapply(function(i,j) paste(j,seq_len(i),"alpha",sep="_"), 
      #                        lapply(pred_standata[pred_Xs_name],ncol), # number of basis function
      #                        pred_Xs_name)
      # colnames(pred_Xs) <- as.character(unlist(Xs_name_list))
      
      
      
      ## X: linear term, for example intercept + x1 + x2 + x1:x2
      # pred_X <- pred_standata$X
      ## design matrix
      
      
      
      pred_B <- cbind(pred_X, pred_Xs, pred_Zs)
      
      if(centered) pred_B <- sweep(pred_B,2,colMeans(pred_B),'-') # Return centered smooths.
      
      
      
      # Bname <- mapply(c, as.list(Xs_name_list), Zs_name_list) # variable names for basis function.
      ## estimate the smooth function by B*\alpha
      # if(!is.list(Bname)){ # convert Bname to a list
      #   Bname <- lapply(seq_len(ncol(Bname)), function(i) Bname[,i])
      # }
      # smooth_pred <- lapply(Bname, function(name)pred_B[,name] %*% beta[name, ])
      # names(smooth_pred) <- paste0("f",seq_along(Bname))
      # out$smooth_pred <- smooth_pred
      
      ## only consider column without 0. 
      
      Predicted <- pred_B %*% beta
      out$Predicted <- Predicted
      out$Predicted_Summary <- as.data.table(do.call(rbind, apply(Predicted, 1, bsummary, ...)))
    }
    if (isTRUE(posterior)) {
      out$Posterior <- beta
    }
  }

  return(out)
}
