## import libraries
library(mam)
library(mgcv)
library(ggplot2)
library(mgcv)
library(brms)
library(data.table)
library(brmsmargins)
library(dplyr)
source("marginalcoef.R")
source("prediction.R")
source("builder.R")
source("generate_pred.R")
source("integretere_fullbayesian.R")

### data ################

data(beavers)
# Create Sex indicator variables. Change 0's to small nonzero value
# due to internal behaviour of lme4 and Matrix packages which convert data zeroes
# to structural zeroes. See https://github.com/lme4/lme4/issues/671
beavers$sex0 <- as.numeric(beavers$sex=="M");beavers$sex0[beavers$sex0==0] <- 0.0000001
beavers$sex1 <- as.numeric(beavers$sex=="F");beavers$sex1[beavers$sex1==0] <- 0.0000001

# Create the prediction data
beaverspred <- data.frame(time=rep(seq(min(beavers$time),max(beavers$time),length=100),2),
                          timeM=c(seq(min(beavers$time),max(beavers$time),length=100),rep(0,100)),
                          timeF=c(rep(0,100),seq(min(beavers$time),max(beavers$time),length=100)),
                          year=factor(levels(beavers$year)[1],levels=levels(beavers$year)),
                          sex=factor(rep(c("M","F"),each=100),levels=levels(beavers$sex)))


### Model ###############

## BMAM
if(FALSE){ 
  ## Don't run. Very slow. use load("data/beavers_brms.rds")
  brms_model <- brm(bf(y ~ year + s(timeF) + s(timeM) + (sex0-1|id)+(sex1-1|id)),
                    data = beavers, family = "bernoulli", cores = 2, seed = 17,
                    warmup = 1000, iter = 2000, chains = 4,backend = "cmdstanr")
  # save(brms_model,file="data/beavers_brms.rds")
}

load("data/beavers_brms.rds")


beaverspred <- generate_pred(brms_model, length = 100)
mc <- marginalcoef(object = brms_model, fullbayesian = F, centered = TRUE, 
                   preddat = beaverspred, CI = 0.95, CIType="ETI", posterior = T)

mc <- marginalcoef(object = brms_model, fullbayesian = F, centered = TRUE, 
                   CI = 0.95, CIType="ETI", posterior = T)


# save(mc, file="mc.rds")
load("mc.rds")
## type of CI: Can be 'ETI' (default), 'HDI', 'BCI', 'SPI' or 'SI'.
## https://easystats.github.io/bayestestR/reference/ci.html


## MAM
bv.mam <- mam(
  smooth = list(s(timeF),s(timeM)),
  re = y ~ (sex0-1|id)+(sex1-1|id),
  fe=~year,#+sex
  dat = beavers,
  margdat = beavers,
  preddat = beaverspred,
  control = mam_control(
    centered=TRUE,
    method = 'trust',
    varmethod = 1,
    verbose = TRUE,
    retcond = TRUE
  )
)




### plot ##################
## plot set up
library(gridExtra)
library(patchwork)
GGPLOTTEXTSIZE <- 15
PLOTWIDTH <- PLOTHEIGHT <- 7
myrange=c(-2,3.2)
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## data prepare
bv.mam.uci <- bv.mam$mam$fitted+1.96*bv.mam$mam$fitted_se
bv.mam.lci <- bv.mam$mam$fitted-1.96*bv.mam$mam$fitted_se

dfplotF <- data.frame(time = rep(beaverspred$timeM[101:200],2),
                      fitted = c(mc$Predicted_Summary$M[101:200],bv.mam$mam$fitted[101:200]),
                      uci = c(mc$Predicted_Summary$UL[101:200], bv.mam.uci[101:200]),
                      lci = c(mc$Predicted_Summary$LL[101:200], bv.mam.lci[101:200]),
                      Method = rep(c("BMAM", "MAM"), each = 100))

dfplotM <- data.frame(time = rep(beaverspred$timeF[1:100],2),
                      fitted = c(mc$Predicted_Summary$M[1:100],bv.mam$mam$fitted[1:100]),
                      uci = c(mc$Predicted_Summary$UL[1:100], bv.mam.uci[1:100]),
                      lci = c(mc$Predicted_Summary$LL[1:100], bv.mam.lci[1:100]),
                      Method = rep(c("BMAM", "MAM"), each = 100))


## plot
ggF <- ggplot(data=dfplotF,aes(x=time,y=fitted, group=Method))+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method, colour= Method),alpha=0.2, size = 0.3)+
  geom_line(aes(colour = Method), size = 1)+
  ylim(myrange)+
  xlab("Day of Year")+
  ylab(expression(f^M~(time)))+
  ggtitle("Female") +
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggF

ggM <- ggplot(data=dfplotM,aes(x=time,y=fitted, group=Method))+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Method, colour= Method),alpha=0.2, size = 0.3)+
  geom_line(aes(colour = Method), size = 1)+
  ylim(myrange)+
  xlab("Day of Year")+
  ylab(expression(f^M~(time)))+
  ggtitle("Male") +
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggM


ggM_combined <- ggM + ylab("")
combined <- ggF + ggM_combined & theme(legend.position = "right")
gg_combined <- combined + plot_layout(guides = "collect")  # combine two plots
gg_combined


## save plots
ggsave(filename = file.path('figures/beaver-time-female.pdf'),plot = ggF,width=PLOTWIDTH,height=PLOTHEIGHT)
ggsave(filename = file.path('figures/beaver-time-male.pdf'),plot = ggM,width=PLOTWIDTH,height=PLOTHEIGHT)
ggsave(filename = file.path('figures/beaver-time-combined.pdf'),plot = gg_combined,width=2*PLOTWIDTH,height=PLOTHEIGHT)




