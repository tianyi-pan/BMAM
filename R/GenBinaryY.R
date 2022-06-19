### Auxiliary functions ####################
##' @title Generate binary response variable of GLMM model based on a known Marginalized model
##' @useDynLib bmam DeconvolveGH_CALL1
Etai2Deltai1 <-function(etai, gamma, sigma, q.points, Z, W){
  expit <- function(aa){exp(aa)/(1+exp(aa))}
  n <- length(etai)

  gam   <- gamma
  if (length(gam)==1){ gam<-rep(gamma, n) }
  sig   <- sigma
  if (length(sig)==1){ sig<-rep(sigma, n) }

  if (length(sig) != n | length(gam) !=n) { stop("Error in etai.2.deltai") }

  deltai <- rep(NA, n)
  deltai <-.Call("DeconvolveGH_CALL1", etai, gam, sig, Z, W)
  deltai
}

expit <- function(aa) {
  exp_aa <- exp(aa)
  exp_aa/(1+exp_aa)
}

get.GH <- function(q, scale_abscissa = sqrt(2), scale_weight=1/sqrt(pi)) {
  rule = gaussHermiteData(q)
  if(scale_abscissa!=1) rule$x = rule$x*scale_abscissa
  if(scale_weight!=1)   rule$w = rule$w*scale_weight
  list(z=rule$x,w=rule$w)
}


### load c files ###########
## in shell: R CMD SHLIB src.c
# dyn.load("src/DeconvolveGH.so")



##' @title Generate binary Y
##' @param lv.formula lv formula for GLMM model
##' @param t.formula transform (lag term) for GLMM model
##' @param beta coefficients for marginalized model
##' @param Sigma covariance matrix for random effects
##' @param gamma coefficients for lag term
##' @param id id for clusters
##' @param data a dataframe containing explanatory variables and id
##' @param q number of nodes
##' @param Yname name of generated reponse variable in the returned dataframe
##' @return a dataframe contains explanatory variables, id and generated binary
##'   response variable
##' @author Tianyi Pan changing the GenBinaryY function in MMLB
##'   package(https://github.com/mercaldo/MMLB).
##' @note the function reduce the \f$X\gamma\f$ into one-dim normal distributed
##'   variable, and then using Gauss-Hermite Quadrature to solve delta. After
##'   getting the solution, generate data based on \f$X\gamma\f$ but not the
##'   one-dim term. This function need include the R file Etai2Deltai1.R,
##'   get.GH.R, expit.R and C code source from src1.so. See more details in
##'   GenBinY_Summary.md and 1 dim in LME.md
##' @import fastGHQuad
##' @export


GenBinaryY <- function(mean.formula, lv.formula = NULL, t.formula = NULL,
                        beta = NULL, Sigma = NULL, gamma = NULL, id, data,
                        q = 10, Yname='Y')
{  # gamma is for transition Y_{t} and Y_{t-1}. Sigma is for random effects
  # Sigma: Variance matrix of random effects.
  if (is.null(lv.formula) & is.null(t.formula)) {
    stop("Specify association model (both lv.formula and t.formula arguments cannot be NULL).")
  }
  if(is.null(beta)) {
    stop('Specify beta values for marginal mean model.')
  }
  if (is.null(Sigma) & is.null(gamma)) {
    stop("Specify Sigma and/or gamma values (both Sigma and gamma arguments cannot be NULL).")
  }

  terms = unique(c(all.vars(mean.formula), all.vars(lv.formula),  # extract variables from formula
                   all.vars(t.formula), as.character(substitute(id))))
  data0 = data # original data
  data  = data[, terms] # the dataset only containing the variables in formula
  # if (any(is.na(data))) data = na.omit(data)

  id = data$id = data[, as.character(substitute(id))] # substitute(id) get name

  x = model.matrix(mean.formula, model.frame(mean.formula, data)) # design matrix for fixed effects
  if(ncol(x) != length(beta)) { # check the dim of beta
    stop('Issue with beta and design matrix associated with mean model.')
  }
  etai <- x %*% cbind(beta, NULL)

  x.t = x.lv = matrix(0, ncol = 1, nrow = nrow(data))
  Sigma.tmp = gamma.tmp = 0
  gam = sig = x.t

  if (!is.null(t.formula)) {
    if(is.null(gamma)) stop('Specify both t.formula and gamma.')
    x.t = model.matrix(t.formula, model.frame(t.formula, data))
    gamma.tmp = gamma
    if(ncol(x.t) != length(gamma.tmp)) {
      stop('Issue with gamma and design matrix associated with association model.')
    }
    gam = x.t %*% cbind(gamma.tmp, NULL)
  }

  if (!is.null(lv.formula)) {
    # *************** Modifications start **********
    if(is.null(Sigma)) stop('Specify both lv.formula and Sigma')
    x.lv = model.matrix(lv.formula, model.frame(lv.formula, data))  # design matrix for random effects

    if(length(Sigma)==1){
      # Only one random effect. (In most cases, only a random slope)
      Sigma.tmp = as.matrix(Sigma)
    }else if(!is.matrix(Sigma)){
      # Only variance for each random effects, with covariance
      Sigma.tmp = diag(Sigma)
    }else{
      stopifnot(nrow(Sigma) == ncol(Sigma)) # check sigma is a square matrix
      Sigma.tmp = Sigma
    }
    # print(Sigma.tmp)
    if(ncol(x.lv) != nrow(Sigma.tmp)) {
      stop('Issue with Sigma and design matrix associated with association model.')
    }
    chol.Sigma <- t(chol(Sigma.tmp)) # cholesky decomposition for Sigma.
    # Sigma.tmp == chol.Sigma %*% t(chol.Sigma)
    sig = x.lv %*% chol.Sigma

    # convert sig to one column
    # Assume we have 2 columns (random intercept and one random slope)
    # denote the first line of sig as a, second line of sig as b
    # in hedeker's method, z T \theta = a \theta_1  + b \theta_2. because theta ~ N(0,I)
    # the i the element: a_i\theta_1 ~ N(0,a_i^2), b_i\theta_2 ~ N(0,b_i^2)
    # the i the element: a_i\theta_1  + b_i\theta_2 ~ N(0,a_i^2+b+i^2)
    # N(0,a_i^2+b_i^2) = sqrt(a_i^2+b_i^2) N(0,1)

    # sig = sqrt(apply(sig^2,1,sum))
    # ************ Modifications stop ***********
  }
  lps <- data.frame(id, etai,gam,sig) # design matrix (random effect and fixed effect)
  colnames(lps)[1:3] <- c('id','etai','gam')



  lps.s = split(lps, lps$id) # split the design matrix by id (cluster)
  n.id  = length(lps.s) # number of id
  YY    = vector('list', n.id) # number of vector = number of id
  # delta_vector    = vector('list', n.id) # add delta
  # mu_vector    = vector('list', n.id) # add mu
  # Zg_vector    = vector('list', n.id) # add Z%*%gamma

  tmp = get.GH(q) # Computes Gauss-Hermite quadrature rule of requested order using Golub-Welsch algorithm.
  W   = tmp$w
  Z   = tmp$z

  for(ix in seq(n.id)) { # id = each of id
    lp.tmp = lps.s[[ix]] # get the matrix of id=ix
    nr.lp  = nrow(lp.tmp) # number units in id=ix

    etai.ix = lp.tmp$etai
    gam.ix  = lp.tmp$gam
    sig.ix  = lp.tmp[,-c(1,2,3)]

    sig.ix_1 = sqrt(apply(sig.ix^2,1,sum))
    deltai = Etai2Deltai1(etai=etai.ix, gamma=gam.ix, sigma=sig.ix_1, Z=Z, W=W)
    # print(etai.ix)
    # print(",")
    # print(deltai)
    # return(0)

    d <- ncol(sig.ix)
    z <- mvtnorm::rmvnorm(1, mean=rep(0,d), sigma=diag(d))

    MuC    = numeric(nr.lp)
    Y      = numeric(nr.lp)
    # Zg = numeric(nr.lp) # add Z%*%gamma

    # MuC[1] = expit(deltai[1]+sig[1]*z) # Assume lagged value = 0

    sig.ix <- as.matrix(sig.ix)

    MuC[1] = expit(deltai[1] + sig.ix[1,] %*% t(z))


    Y[1]   = as.integer(runif(1)<MuC[1])
    # Zg[1] = sig.ix[1, ] %*% t(z) # add Z%*%gamma

    if(nr.lp > 1) {
      for ( it in seq(2, nr.lp) ) {
        MuC[it] = expit( deltai[it] + gam.ix[it]*Y[it-1] + sig.ix[it, ] %*% t(z) )
        Y[it]   = as.integer( runif(1)<MuC[it] )
        # Zg[it] = sig.ix[it, ] %*% t(z) # add Z%*%gamma
      }
    } else {
      Y <- Y #Y <- NA
    }
    YY[[ix]] <- Y
    # mu_vector[[ix]] <- MuC # add mu
    # delta_vector[[ix]] <- deltai # add delta
    # Zg_vector[[ix]] <- Zg # add Z%*%gamma
  }
  Y <- unlist(YY)
  data0[,Yname] <- Y

  # Delta <- unlist(delta_vector) #add delta
  # data0[,"delta"] <- Delta # add delta
  # data0[,"mu"] <- unlist(mu_vector) # add mu
  # data0[,"ZGamma"] <- unlist(Zg_vector) # add Z%*%gamma

  data0
}




