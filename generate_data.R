### include functions 
source("R/GenBinaryY.R")
### load packages 
library(mam)
library(mgcv)
library(gamm4)
library(tidyverse)
library(brms)
library(brmsmargins)
library(data.table)
# source("R/bmam.R")
# source("R/generate_pred.R")
# source("R/SimData.R")
# source("R/summary.R")
# source("R/plot.bmam.R")
# source("R/conditional_brms.R")


### Case 1 #####################

## setting
K <- 30
Nk <- 5

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


set.seed(54321)
x1 <- round(runif(K*Nk,-1,1),3)
x3 <- round(runif(K*Nk,-1,1),3)
x2 <- round(runif(K*Nk,-1,1),3)


fx1 <- f1(x1)
fx2 <- f2(x2)
beta <- c(0,1,1,beta3)



mean.formula <- ~fx1+fx2+x3
lv.formula <- ~1+x3
Sigma <- matrix(c(sigma0^2,
                  rho*sigma0*sigma3,rho*sigma0*sigma3,
                  sigma3^2),nrow=2,ncol=2)


nclust <- rep(Nk,K) # number of units in each cluster
id <- rep(seq(K), nclust) # id of each unit
data  <- data.frame(id, x1, fx1, x2, fx2, x3) # generate data (without y)
data  <- data[order(data$id),] # order data by id and time


data_gen <- GenBinaryY(mean.formula = mean.formula, lv.formula = lv.formula, 
                         beta=beta, Sigma=Sigma, id=id, data=data, q=120, 
                         Yname = "y")



## check: compare with SimData()
set.seed(54321)
data_simdata <-SimData(K,Nk)


library(brms)
model_gen <- brm(bf(y ~  x3 + s(x1) + s(x2) + (1+x3|id)),
                  data = data_gen, family = "bernoulli", cores = 4, seed = 4321,
                  warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")

bmam_gen <- bmam(object = model_gen, centered = F,
                 k=100, CIType="ETI", CI = 0.95)


model_simdata <- brm(bf(y ~  x3 + s(x1) + s(x2) + (1+x3|id)),
                 data = data_simdata$data, family = "bernoulli", cores = 4, seed = 4321,
                 warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")

bmam_simdata <- bmam(object = model_simdata, centered = F,
                     k=100, CIType="ETI", CI = 0.95)


p_gen <- plot(bmam_gen, display = FALSE)
p_simdata <- plot(bmam_simdata, display = FALSE)

g_gen <- plot(bmam_gen, smooth.function = data_simdata$f, display = F)
g_simdata <- plot(bmam_simdata,smooth.function = data_simdata$f, display = F)

library(gridExtra)
g <- grid.arrange(g_gen$Conditional$Comparison[[2]]+ylim(-8,4), 
             g_simdata$Conditional$Comparison[[2]]+ylim(-8,4),
             nrow = 1)

ggsave(g, filename = paste0("figures/compare2gen2_",as.character(K),"_",as.character(Nk),"true.png"),
       width = 10, height = 6)

View(bmam_gen$Preddat)
View(bmam_simdata$Preddat)

## mam 
mam_gen <- mam::mam(smooth = list(s(x1),s(x2)),
                   re = y ~ (1+x3|id),
                   fe = ~ x3,
                   dat = data_gen,
                   margdat = data_gen,
                   preddat = bmam_gen$Preddat,
                   control = mam_control(
                     method = 'trust',
                     varmethod = 1,
                     verbose = FALSE,
                     retcond = TRUE))

dat <- data_simdata$data
mam_simdata <- mam::mam(smooth = list(s(x1),s(x2)),
                    re = y ~ (1+x3|id),
                    fe = ~ x3,
                    dat = dat,
                    margdat = dat,
                    preddat = bmam_simdata$Preddat,
                    control = mam_control(
                      method = 'trust',
                      varmethod = 1,
                      verbose = FALSE,
                      retcond = TRUE))






plot(bmam_simdata, compared.model = mam_simdata)

plot(bmam_gen, compared.model = mam_gen)




### a mixture of two Normal distributions (Bie's paper) ###########


## setting
K <- 10
Nk <- 6

data_simdata <-SimData(K,Nk)



f1temp <- function(x){(1.5*sin(pi*x)-2*x)}

f1 <- function(x) f1temp(x)-f1temp(0)
# f1 <- function(x) x

f2temp <- function(x){5*(dnorm(2*(x)))}

f2 <- function(x) f2temp(x)-f2temp(0)
# f2 <- function(x) x

beta3 <- 1
beta4 <- 2 # coefficient for Xe



# set.seed(32)
x1 <- round(runif(K*Nk,-1,1),3)
x3 <- round(runif(K*Nk,-1,1),3)
x2 <- round(runif(K*Nk,-1,1),3)

nclust <- rep(Nk,K) # number of units in each cluster
trt <- rep(1e-4,K)
trt[1:(K*0.5)] <- 1
xe <- rep(trt, nclust) # binary exposure


fx1 <- f1(x1)
fx2 <- f2(x2)
beta <- c(0,1,1,beta3,beta4)
# beta <- c(0,1,1,beta3)



mean.formula <- ~fx1+fx2+x3+xe
# mean.formula <- ~fx1+fx2+x3
lv.formula <- ~1+x3

sigma0 <- 2 #2 # random intercept sd
sigma3 <- 0.75 # random slope sd
rho <- 0 # correlation between random intercept and slope

Sigma1 <- matrix(c(sigma0^2,
                  rho*sigma0*sigma3,rho*sigma0*sigma3,
                  sigma3^2),nrow=2,ncol=2)
Sigma2 <- matrix(c(sigma3^2,
                    0*sigma0*sigma3,0*sigma0*sigma3,
                    sigma0^2),nrow=2,ncol=2)


id <- rep(seq(K), nclust) # id of each unit
data  <- data.frame(id, x1, fx1, x2, fx2, x3, xe) # generate data (without y)
# data  <- data.frame(id, x1, fx1, x2, fx2, x3) # generate data (without y)
data  <- data[order(data$id),] # order data by id and time


data_gen1 <- GenBinaryY(mean.formula = mean.formula, lv.formula = lv.formula, 
                       beta=beta, Sigma=Sigma1, id=id, data=data[which(xe==1),], q=120, 
                       Yname = "y")
data_gen2 <- GenBinaryY(mean.formula = mean.formula, lv.formula = lv.formula, 
                        beta=beta, Sigma=Sigma2, id=id, data=data[which(xe<1),], q=120, 
                        Yname = "y")

data_gen <- rbind(data_gen1, data_gen2)




## bmam
model_gen <- brm(bf(y ~  x3 + xe + s(x1) + s(x2) + (1+x3|id)),
                 data = data_gen, family = "bernoulli", cores = 4, seed = 4321,
                 warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")

bmam_gen <- bmam(object = model_gen, centered = F,
                 k=100, CIType="ETI", CI = 0.95)

plot(bmam_gen, smooth.function = data_simdata$f)


## mam 
mam_gen <- mam::mam(smooth = list(s(x1),s(x2)),
                    re = y ~ (1+x3|id),
                    fe = ~ x3 + xe,
                    dat = data_gen,
                    margdat = data_gen,
                    preddat = bmam_gen$Preddat,
                    control = mam_control(
                      method = 'trust',
                      varmethod = 1,
                      verbose = FALSE,
                      retcond = TRUE))





g <- plot(bmam_gen, smooth.function = data_simdata$f,
          compared.model = mam_gen, display = F)

library(gridExtra)
g_combined <- grid.arrange(g$Compared.Model$Comparison[[1]]+ylim(-8,4), 
                           g$Compared.Model$Comparison[[2]]+ylim(-8,4),
                           nrow = 1)


ggsave(g_combined, filename = paste0("figures/mis_",as.character(K),"_",as.character(Nk),".png"),
       width = 10, height = 5)


## compare the results with mam, especially in small sample size



