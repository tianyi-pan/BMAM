### compare with frequentist MAM 
### load packages #####################
library(mam)
library(mgcv)
library(gamm4)
library(tidyverse)
library(brms)
library(brmsmargins)
library(data.table)

### generate data #######################
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

K = 50
Nk = 10
sigma0 = 2
sigma3 = 1
rho = 0.5
# itr = 1:B


trueSigma <- matrix(c(sigma0^2,
                      rho*sigma0*sigma3,rho*sigma0*sigma3,
                      sigma3^2),nrow=2,ncol=2)

## generate pred data
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




### Model #################### 

## BMAM
if(TRUE){
  # don't run. Very slow. use load("data/simu_brms.rds")
  model_brms <- brm(bf(y ~  x3 + s(x1) + s(x2) + (1+x3|id)),
                    data = dat, family = "bernoulli", cores = 2, seed = 17,
                    warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")
  # save(model_brms, file = "data/simu_brms.rds")
}
# load("data/simu_brms.rds")

source("marginalcoef.R")
source("prediction.R")
source("builder.R")
source("integretere_fullbayesian.R")
mc.T <- marginalcoef(object = model_brms, fullbayesian = T, preddat = dat2pred,CIType="ETI", CI = 0.95, posterior = T)
save(mc.T, file = "mc_simulation.rds")

mc.F <- marginalcoef(object = model_brms, fullbayesian = F, preddat = dat2pred,CIType="ETI", CI = 0.95, posterior = T)

# mc <- marginalcoef(object = model, preddat = dat, CI = 0.95, posterior = T)



## MAM
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




### plot ##################
## plot set up
mc <- mc.T
library(gridExtra)
library(patchwork)
GGPLOTTEXTSIZE <- 15
PLOTWIDTH <- PLOTHEIGHT <- 7
myrange=c(-3,2.5)
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## data prepare
themam.uci <- themam$mam$fitted+1.96*themam$mam$fitted_se
themam.lci <- themam$mam$fitted-1.96*themam$mam$fitted_se


dfplotfx1 <- data.frame(x = rep(dat2pred$x1[1:100],2),
                        fitted = c(mc$Predicted_Summary$M[1:100],themam$mam$fitted[1:100]),
                        uci = c(mc$Predicted_Summary$UL[1:100], themam.uci[1:100]),
                        lci = c(mc$Predicted_Summary$LL[1:100], themam.lci[1:100]),
                        Method = rep(c("BMAM", "MAM"), each = 100))

dfplotfx2 <- data.frame(x = rep(dat2pred$x2[101:200],2),
                      fitted = c(mc$Predicted_Summary$M[101:200],themam$mam$fitted[101:200]),
                      uci = c(mc$Predicted_Summary$UL[101:200], themam.uci[101:200]),
                      lci = c(mc$Predicted_Summary$LL[101:200], themam.lci[101:200]),
                      Method = rep(c("BMAM", "MAM"), each = 100))



## plot
ggx1 <- ggplot(data=dfplotfx1,aes(x=x,y=fitted, group=Method))+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method, colour= Method),alpha=0.2, size = 0.3)+
  geom_line(aes(colour = Method), size = 1)+
  geom_line(aes(y=rep(dat2pred$fx1[1:100],2)), size = 1, lty = "dashed")+
  ylim(myrange)+
  xlab(expression(X[1]))+
  ylab(expression(f~(X[1])))+
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggx1
ggx2 <- ggplot(data=dfplotfx2,aes(x=x,y=fitted, group=Method))+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method, colour= Method),alpha=0.2, size = 0.3)+
  geom_line(aes(colour = Method), size = 1)+
  geom_line(aes(y=rep(dat2pred$fx2[101:200],2)), size = 1, lty = "dashed")+
  ylim(myrange)+
  xlab(expression(X[2]))+
  ylab(expression(f~(X[2])))+
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggx2



combined <- ggx1 + ggx2 & theme(legend.position = "right")
gg_combined <- combined + plot_layout(guides = "collect")  # combine two plots
gg_combined


## save plots
ggsave(filename = file.path(paste0('figures/fullbayes(wrong)-simulation-combined', 'K = ', as.character(K),
                                   'Nk = ', as.character(Nk), '.pdf')),
       plot = gg_combined,
       width=2*PLOTWIDTH,height=PLOTHEIGHT)



