#' Marginal Posterior Predictions from a 'brms' Model
#'
#' Calculate marginal predictions from a \code{brms} model.
#' Marginal predictions average over the input data for each posterior draw.
#' Marginal predictions for models with random effects will integrate
#' over random effects.
#'
#' @param object A fitted brms model object. Required.
#' @param data A data frame or data table passed to \code{fitted()}
#'   as the new data to be used for predictions. Required.
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
#' @param dpar Parameter passed on the \code{dpar}
#'   argument of \code{fitted()} in brms. Defaults to \code{NULL}
#'   indicating the mean or location parameter typically.
#' @param resample An integer indicating the number of
#'   bootstrap resamples of the posterior predictions to
#'   use when calculating summaries. Defaults to \code{0L}.
#'   See documentation from [.averagePosterior()] for more details.
#' @param resampleseed A seed for random number generation. Defaults to \code{FALSE},
#'   which means no seed is set.
#'   Only used if \code{resample} is a positive, non-zero integer.
#'   See documentation from [.averagePosterior()] for more details.
#' @param effects A character string indicating the type of
#'   prediction to be made. Can be one of
#'   \dQuote{fixedonly} meaning only use fixed effects,
#'   \dQuote{includeRE} meaning that random effects should be
#'   included in the predictions, or
#'   \dQuote{integrateoutRE} meaning that random effects should be
#'    integrated out / over in the predictions.
#' @param backtrans A character string indicating the type of
#'   back transformation to be applied. Can be one of
#'   \dQuote{response} meaning to use the response scale,
#'   \dQuote{linear} or \dQuote{identity} meaning to use the linear predictor scale,
#'   or a specific back transformation desired, from a possible list of
#'   \dQuote{invlogit}, \dQuote{exp}, \dQuote{square}, or \dQuote{inverse}.
#'   Custom back transformations should only be needed if, for example,
#'   the outcome variable was transformed prior to fitting the model.
#' @param k An integer providing the number of random draws to use for
#'   integrating out the random effects. Only relevant when \code{effects}
#'   is \dQuote{integrateoutRE}.
#' @param raw A logical value indicating whether to return the raw output or
#'   to average over the Monte Carlo samples. Defaults to \code{FALSE}.
#'   Setting it to \code{TRUE} can be useful if you want not only the
#'   full posterior distribution but also the \code{k} Monte Carlo samples
#'   used for the numerical integration. This cannot be used with
#'   \code{summarize = TRUE}.
#' @param ... Additional arguments passed to \code{fitted()}
#' @return A list with \code{Summary} and \code{Posterior}.
#'   Some of these may be \code{NULL} depending on the arguments used.
#' @references
#' Pavlou, M., Ambler, G., Seaman, S., & Omar, R. Z. (2015)
#' \doi{10.1186/s12874-015-0046-6}
#' \dQuote{A note on obtaining correct marginal predictions from a random intercepts model for binary outcomes}
#' and
#' Skrondal, A., & Rabe-Hesketh, S. (2009)
#' \doi{10.1111/j.1467-985X.2009.00587.x}
#' \dQuote{Prediction in multilevel generalized linear models}
#' @importFrom data.table as.data.table
#' @importFrom stats fitted formula
#' @importFrom posterior as_draws_df ndraws
#' @importFrom brms standata
#' @importFrom methods missingArg
#' @export
prediction <- function(object, data, summarize = TRUE, posterior = FALSE,
                       index, dpar = NULL, resample = 0L, resampleseed = FALSE,
                       effects = c("fixedonly", "includeRE", "integrateoutRE"),
                       fullbayesian = TRUE, 
                       backtrans = c("response", "linear", "identity",
                                     "invlogit", "exp", "square", "inverse"),
                       k = 100L, raw = FALSE, ...) {
  ## checks and assertions
  brmsmargins:::.assertbrmsfit(object)
  brmsmargins:::.assertdpar(object, dpar = dpar)

  if (isFALSE(brmsmargins:::is.random(object))) {
    if (isFALSE(effects == "fixedonly")) {
      stop("object does not have random effects: must use \"effects = 'fixedonly'\"")
    }
  }

  effects <- match.arg(effects, several.ok = FALSE)
  backtrans <- match.arg(backtrans, several.ok = FALSE)

  if (isTRUE(effects == "integrateoutRE")) {
    ## assert the assumed family / distribution is a supported one
    brmsmargins:::.assertfamily(object)
    ## assert the link function used is a supported one
    brmsmargins:::.assertlink(object, dpar = dpar)
    ## assert that all random effects in the model are Gaussian
    brmsmargins:::.assertgaussian(object)
  }

  if (isTRUE(missingArg(index))) {
    index <- seq_len(ndraws(object))
  }

  links <- brmsmargins:::.links(
    link = brmsmargins:::.extractlink(object, dpar),
    effects = effects, backtrans = backtrans)

  ## set whether fitted() should include RE (NULL) or not (NA)
  ## see help for ?fitted.brmsfit for more details ?fitted, which is alias of posterior_epred()
  # posterior_epred. Expectation of posterior predictive distribution.
  # posterior_predict
  # rstanarm::posterior_predict() brms::posterior_predict().
  
  
  # 1. draw parameters from posterior distributional (just the posterior sample in fit1, 4000 draw)
  # 2. given a set of parameters and new X, draw a y from likelihood(add_predicted_draws), 
  #   or calculate the mean value of likelihood(add_epred_draws). add_epred_draws has a lower variance. 
  #   But they are similar when calculate the average over different draws. (4000 draw for per value of X)
  
  if (isTRUE(effects %in% c("fixedonly", "integrateoutRE"))) {
    useRE <- NA
  } else if (isTRUE(effects == "includeRE")) {
    useRE <- NULL
  }

  ## generate all predictions (if fixedonly or includeRE)
  ## or generate just the fixed effects predictions (if integrateRE)
  
  
  # possible: Might be wrong!
  # In marginalcoefs, useRE <- NA. Do not consider random effect when doing posterior predictive distribution
  # It can be termed as the posterior predictive distribution and random effects are zero.
  # when integrate the random effects, we should draw random effects from posterior distribution,
  # and then add them into the predicted value of y.  For One y, several drawed random effects, take average over these draws.
  # This procedure can integrate out the random effects.
  # see margins.R 
  
  
  ## or generate just the fixed effects predictions (if integrateoutRE)
  yhat <- fitted(
    object = object, newdata = data,
    re_formula = useRE, # useRE = NA if integrateoutRE (in marginaleffects and AMEs)
    scale = links$scale, dpar = dpar,
    draw_ids = index, summary = FALSE)
  yhat <- links$useifun(yhat) #  a matrix

  if (isTRUE(effects == "integrateoutRE")) {
    if (isTRUE(links$ilink != "identity")) {
      # posterior sample from model. 
      post <- as.data.table(as_draws_df(object))[index, ]
      
      # data
      dtmp <- standata(object, newdata = data, check_response = FALSE)

      # random effects
      re <- as.data.table(object$ranef)
      
      if (is.null(dpar)) {
        usedpar <- ""
      } else {
        usedpar <- dpar
      }
      
      # usually all random effects. 
      re <- re[dpar == usedpar]

      
      blocks <- unique(re$id)
      nblocks <- length(blocks)

      d2 <- sd <- L <- vector("list", nblocks)
      # r <- vector("list",nblocks)
      J <- vector("list",nblocks)
      z <- vector("list",nblocks)
      
      
      # if (fullbayesian){
      # for (i in seq_len(nblocks)) {
      #   useblock <- blocks[i]
      #   usere <- re[id == useblock]
      #   num <- max(usere$cn)
      #   usecoef <- usere$coef
      #   usegroup <- unique(usere$group)
      #   d2[[i]] <- brmsmargins:::.buildZ(data = dtmp, block = useblock, number = num, dpar = dpar)
      #   J[[i]] <- .buildJ(dtmp, useblock)[[1]]
      #   # r[[i]] <- .buildr(data = post, id = J[[i]], usecoef,usegroup)  # TO DO: change dim
      #   sd[[i]] <- brmsmargins:::.buildSD(data = post, ranef = usere, block = useblock, dpar = dpar)
      #   L[[i]] <- brmsmargins:::.buildL(data = post, block = useblock, number = num)
      #   z[[i]] <- .buildz(data = post, id = J[[i]], usecoef,usegroup)
      #   # Z[[i]] <- sapply(J[[i]], function(row.){
      #   #   # r_index <- which(unique(J[[i]])==row.) 
      #   #   # r <- .buildr(data = post, id = row., usecoef,usegroup)
      #   #   z <- .buildz(data = post, id = row., usecoef,usegroup)
      #   #   # d2[[i]][row.,] * r[[i]][[r_index]] 
      #   #   do.call(cbind,r) %*% d2[[i]][row.,] 
      #   # })
      # }
      
      # ## TO Do: change dim
      # Z.sum <- do.call("+",Z)
      # ## integrate Z1+Z2?
      # inte <- function(yhat_row) ## integrate out random effect (full bayeisan).
      #   colMeans(links$ifun(yhat_row + Z.sum))
      # yhat.df <- as.data.frame(t(yhat))
      # yhat <- t(sapply(yhat.df,inte))
      

      # # integrae for Z1 and Z2
      ## sampling might not be independent. Randomly select k samples from Z1 and the other k samples from Z2, 
      ## then sum them up as Z.sum
      ## have a look at the dist. of r
      
      ## how to deal with two random effects. 1. See the brmsmargins, 2.try model with one random effect. 
      # inte <- function(yhat_row) ## integrate out random effect (full bayeisan).
      #   colMeans(links$ifun(yhat_row + Z[[1]]))
      # yhat.df <- as.data.frame(t(yhat))
      # yhat <- t(sapply(yhat.df,inte))
      # 
      # # yhat <- links$ifun(yhat)
      
      ## how to fasten the function? better function in R or Cpp
      
      # %>% apply(1, sum)
      
      ## !! to do: change cpp code!! 
      ## use z instead of N(0,1), sample from z not generate data from N(0,1) 
       # yhat <- integratere(d = d2, sd = sd, L = L, k = k,
       #                    yhat = yhat, backtrans = links$useilinknum)
       # 
    # }else{
      for (i in seq_len(nblocks)) {
        useblock <- blocks[i]
        usere <- re[id == useblock]
        num <- max(usere$cn)
        d2[[i]] <- brmsmargins:::.buildZ(data = dtmp, block = useblock, number = num, dpar = dpar)
        sd[[i]] <- brmsmargins:::.buildSD(data = post, ranef = usere, block = useblock, dpar = dpar)
        L[[i]] <- brmsmargins:::.buildL(data = post, block = useblock, number = num)
        
        usecoef <- usere$coef
        usegroup <- unique(usere$group)
        J[[i]] <- .buildJ(dtmp, useblock)[[1]]
        z[[i]] <- .buildz(data = post, id = J[[i]], usecoef,usegroup)
        
        names(d2)[i] <- names(sd)[i] <- names(L)[i] <- names(z)[i] <- names(J)[i] <- 
          sprintf("Block%d", useblock)
      }
      # cpp file
      
      if(fullbayesian){
        yhat <- integratere_fullbayesian(d2, sd, L, J, post, usecoef,usegroup, k, yhat, links$ifun)
        ## TO DO: write it as cpp code
      }else{
        yhat <- integratere(d = d2, sd = sd, L = L, k = k,
                            yhat = yhat, backtrans = links$useilinknum)  
      }
      
      # }
      
    }
  }
  # raw = true in marginalcoef.R
  if (isTRUE(raw)) {
    if (isTRUE(summarize)) {
      message("summarize cannot be TRUE when raw = TRUE, setting to FALSE")
      summarize <- FALSE
    }
    if (isFALSE(posterior)) {
      message("posterior cannot be FALSE when raw = TRUE, setting to TRUE")
      posterior <- TRUE
    }
  } else {
    ## average across rows
    ## either using row wise means, or row wise bootstrapped means
    yhat <- brmsmargins:::.averagePosterior(yhat, resample = resample, seed = resampleseed)
  }
  # in marginalceofs. nrow(yhat) = num of samples, ncol(yhat) = num of observations.
  out <- list(
    Summary = NULL,
    Posterior = NULL)

  if (isTRUE(summarize)) {
    out$Summary <- bsummary(yhat, ...)
  }
  if (isTRUE(posterior)) {
    out$Posterior <- yhat
  }

  return(out)
}
