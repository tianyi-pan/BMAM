## TO DO: make it flexible to support more models.
## functions in mam paper (https://github.com/awstringer1/mam-paper-code)

get_data_mam <- function(SSmat,x1,x2,x3,K,Nk,f1,f2,beta3){
  if(SSmat[2,2]==0){
    V <- as.data.frame(rnorm(K,0,sqrt(SSmat[1,1])))
    colnames(V) <- "intercepts"
    V$id <- 1:K
  }else{
    V <- mvtnorm::rmvnorm(K,sigma = SSmat) # Intercept and slope
    V <- as.data.frame(V)
    colnames(V) <- c("intercepts","slopes")
    V$id <- 1:K
  }
  dat <- data.frame(id = rep(1:K,each=Nk))
  dat <- merge(dat,V,by="id",all.x=TRUE)
  dat$x1 <- x1
  dat$x2 <- x2
  dat$x3 <- x3
  dat$fx1 <- f1(x1)
  dat$fx2 <- f2(x2)

  ## approximation from Hedeker
  if(SSmat[2,2]==0){
    Ui <- dat[,c("intercepts")]
    vars <- SSmat[1,1]+ ##
      +(pi^2)/3 # (pi^2)/3 is approximate due to logistic distn
  }else{
    Z <- cbind(1,x3) # ranefs with design matrix
    Ui <- apply(Z*dat[,c("intercepts","slopes")],1,sum)
    vars <- colSums(((chol(SSmat)%*%t(Z)))^2)+ ## faster than #diag(Z%*%SS%*%t(Z))
      +(pi^2)/3 # (pi^2)/3 is approximate due to logistic distn
  }
  estar <- rlogis(K*Nk) # approx: rnorm(K*Nk,0,sqrt((pi^2)/3))
  dat <- within(dat,{y <- as.numeric(pnorm(estar+Ui,0,sqrt(vars)) < mam::ilogit(fx1 + fx2 +               ## non linear fixed effects
                                                                                  beta3*x3 +                ## linear fixed effects
                                                                                  0))  })  ## random effects on other side
  return(dat)
}


#' Generate simulated data
#'
#' @param K Number of clusters
#' @param Nk Number of units within clusters
#'
#' @return a list containing generated data and true functions.
#' @export SimData
SimData <- function(K, Nk){
  f1temp <- function(x){(1.5*sin(pi*x)-2*x)}

  f1 <- function(x) f1temp(x)-f1temp(0)
  # f1 <- function(x) x

  f2temp <- function(x){5*(dnorm(2*(x)))}

  f2 <- function(x) f2temp(x)-f2temp(0)
  # f2 <- function(x) x

  beta3 <- 0

  sigma0 <- 2 #2 # random intercept sd
  sigma3 <- 1 # random slope sd
  rho <- 0.5 # correlation between random intercept and slope

  trueSigma <- matrix(c(sigma0^2,
                        rho*sigma0*sigma3,rho*sigma0*sigma3,
                        sigma3^2),nrow=2,ncol=2)

  ## generate pred data
  x1 <- round(runif(K*Nk,-1,1),3)
  x3 <- round(runif(K*Nk,-1,1),3)
  x2 <- round(runif(K*Nk,-1,1),3)
  dat <- get_data_mam(trueSigma,x1,x2,x3,K,Nk,f1,f2,beta3)
  return(list(data = dat, f = list(f1,f2)))
}

