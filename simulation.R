### compare with frequentist MAM ##############################
library(mam)
library(mgcv)
library(gamm4)
library(tidyverse)
library(brms)
library(brmsmargins)
library(data.table)
## generate data
get_data_mam <- function(SSmat,x1,x2,x3,K,Nk){
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


f1temp <- function(x){(1.5*sin(pi*x)-2*x)}
f1 <- function(x) f1temp(x)-f1temp(0)
f2temp <- function(x){5*(dnorm(2*(x)))}
f2 <- function(x) f2temp(x)-f2temp(0)
beta3 <- 0

sigma0 <- 2 #2 # random intercept sd
sigma3 <- 1 # random slope sd
rho <- 0.5 # correlation between random intercept and slope

K = 200
Nk = 10
sigma0 = 2
sigma3 = 1
rho = 0.5
# itr = 1:B


trueSigma <- matrix(c(sigma0^2,
                      rho*sigma0*sigma3,rho*sigma0*sigma3,
                      sigma3^2),nrow=2,ncol=2)

## generate data
x1 <- round(runif(K*Nk,-1,1),3)
x3 <- round(runif(K*Nk,-1,1),3)
x2 <- round(runif(K*Nk,-1,1),3)
dat <- get_data_mam(trueSigma,x1,x2,x3,K,Nk)

gridlen=100
dat2pred <- data.frame(x1 = c(seq(-1,1,length=gridlen),rep(0,2*gridlen)),
                       x2 = c(rep(0,gridlen),seq(-1,1,length=gridlen),rep(0,gridlen)),
                       x3 = c(rep(0,2*gridlen),seq(-1,1,length=gridlen)))

dat2pred$fx1 <- f1(dat2pred$x1)
dat2pred$fx2 <- f2(dat2pred$x2)

View(dat2pred)

model <- brm(bf(y ~  x3 + s(x1) + s(x2) + (1+x3|id)),
             data = dat, family = "bernoulli", cores = 2, seed = 17,
             warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")

source("marginalcoef.R")
mc <- marginalcoef(object = model, preddat = dat2pred, CI = 0.95, posterior = T)

# mc <- marginalcoef(object = model, preddat = dat, CI = 0.95, posterior = T)



themam <- mam::mam(smooth = list(s(x1),s(x2)),
                   re = y ~ (1+x3|id),
                   fe = ~ x3,
                   dat = dat,
                   margdat = dat,
                   preddat = dat2pred ,
                   control = mam_control(
                     method = 'trust',
                     varmethod = 1,
                     verbose = FALSE,
                     retcond = TRUE))

1.96*as.numeric(themam$mam$fitted_se)



## predicted 
dat2pred_results <- cbind(dat2pred, mc$Predicted_Summary$M, themam$mam$fitted)
names(dat2pred_results)[6:7] <- c("BMAM", "MAM")
dat2pred_results_plot <- filter(dat2pred_results,x2==0,x3==0)
colors <- c("BMAM" = "red", "MAM" = "blue", "Fx1" = "black")
ggplot(dat2pred_results_plot, aes(x=x1)) + 
  geom_line(aes(y=BMAM, colour = "BMAM")) + 
  geom_line(aes(y=MAM, colour = "MAM")) + 
  geom_line(aes(y=fx1, colour = "Fx1")) + 
  scale_color_manual(values = colors) + labs(color = "", x="X1", y = "f(X1)")

dat2pred_results_plot <- filter(dat2pred_results,x1==0,x3==0)
colors <- c("BMAM" = "red", "MAM" = "blue", "Fx1" = "black")
ggplot(dat2pred_results_plot, aes(x=x2)) + 
  geom_line(aes(y=BMAM, colour = "BMAM")) + 
  geom_line(aes(y=MAM, colour = "MAM")) + 
  geom_line(aes(y=fx2, colour = "Fx1")) + 
  scale_color_manual(values = colors) + labs(color = "", x="X2", y = "f(X2)")
