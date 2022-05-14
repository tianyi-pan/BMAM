### compare with frequentist MAM 
### load packages #####################
library(mam)
library(mgcv)
library(gamm4)
library(tidyverse)
library(brms)
library(brmsmargins)
library(data.table)
source("R/bmam.R")
source("R/prediction.R")
source("R/generate_pred.R")
source("R/SimData.R")
source("R/summary.R")
source("R/plot.bmam.R")
source("R/conditional_brms.R")
### generate data #######################
set.seed(4321)
simdata <- SimData(100,10)
dat <- simdata$data
fun <- simdata$f



# gridlen <- 100
# dat2pred <- data.frame(x1 = c(seq(-1,1,length=gridlen),rep(0,2*gridlen)),
#                        x2 = c(rep(0,gridlen),seq(-1,1,length=gridlen),rep(0,gridlen)),
#                        x3 = c(rep(0,2*gridlen),seq(-1,1,length=gridlen)))
# 
# dat2pred$fx1 <- f1(dat2pred$x1)
# dat2pred$fx2 <- f2(dat2pred$x2)



### Model #################### 

## BMAM
if(FALSE){
  # don't run. Very slow. use load("data/simu_brms.rds")
  model_brms <- brm(bf(y ~  x3 + s(x1) + s(x2) + (1+x3|id)),
                    data = dat, family = "bernoulli", cores = 4, seed = 4321,
                    warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")
  # save(model_brms, file = "data/simu_brms.rds")
}
load("data/simu_brms.rds")


bmam <- bmam(object = model_brms,
                     k=100, CIType="ETI", CI = 0.95)

## preddat 
dat2pred <- bmam$Preddat

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

## bgam
if(FALSE){
  bgam <- brm(bf(y ~  x3 + s(x1) + s(x2)),
                    data = dat, family = "bernoulli", cores = 4, seed = 4321,
                    warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")
  # save(bgam, file = "data/simu_bgam.rds")
}
load("data/simu_bgam.rds")

### plot ##################
## plot set up

library(gridExtra)
library(patchwork)
library(ggplot2)
GGPLOTTEXTSIZE <- 15
PLOTWIDTH <- PLOTHEIGHT <- 7
myrange=c(-3,2.5)
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## plot function
plot(bmam, compared.model = themam, smooth.function = fun)

## data prepare
themam.uci <- themam$mam$fitted+1.96*themam$mam$fitted_se
themam.lci <- themam$mam$fitted-1.96*themam$mam$fitted_se


dfplotfx1 <- data.frame(x = rep(dat2pred$x1[1:100],2),
                        fitted = c(mc$Predicted_Summary$M[1:100], themam$mam$fitted[1:100]),
                        uci = c(mc$Predicted_Summary$UL[1:100], themam.uci[1:100]),
                        lci = c(mc$Predicted_Summary$LL[1:100], themam.lci[1:100]),
                        Method = rep(c("BMAM", "MAM"), each = 100))



dfplotfx2 <- data.frame(x = rep(dat2pred$x2[101:200],2),
                        fitted = c(mc$Predicted_Summary$M[101:200], themam$mam$fitted[101:200]),
                        uci = c(mc$Predicted_Summary$UL[101:200], themam.uci[101:200]),
                        lci = c(mc$Predicted_Summary$LL[101:200], themam.lci[101:200]),
                        Method = rep(c("BMAM", "MAM"), each = 100))


## plot
ggx1 <- ggplot(data=dfplotfx1,aes(x=x,y=fitted, group=Method))+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method, colour= Method),alpha=0.1, size = 0.3)+
  geom_line(aes(colour = Method), size = 1)+
  geom_line(aes(y=rep(dat2pred$fx1[1:100],3)), size = 1, lty = "dashed")+
  ylim(myrange)+
  xlab(expression(X[1]))+
  ylab(expression(f~(X[1])))+
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggx1

ggx2 <- ggplot(data=dfplotfx2,aes(x=x,y=fitted, group=Method))+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method, colour= Method),alpha=0.2, size = 0.3)+
  geom_line(aes(colour = Method), size = 1)+
  geom_line(aes(y=rep(dat2pred$fx2[101:200],3)), size = 1, lty = "dashed")+
  ylim(myrange)+
  xlab(expression(X[2]))+
  ylab(expression(f~(X[2])))+
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggx2



combined <- ggx1 + ggx2 & theme(legend.position = "right")
gg_combined <- combined + plot_layout(guides = "collect")  # combine two plots
gg_combined


## save plots
ggsave(filename = file.path(paste0('figures/simulation-combined', 'K = ', as.character(K),
                                   'Nk = ', as.character(Nk), '.pdf')),
       plot = gg_combined,
       width=2*PLOTWIDTH,height=PLOTHEIGHT)



