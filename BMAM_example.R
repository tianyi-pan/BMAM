library(mgcv)
library(brms)
library(data.table)
dat <- gamSim(6,n=200,scale=.2,dist="poisson")
names(dat)[11] <- "id"
## brms
# https://fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
model_brms <- brm(bf(y ~  s(x0,k=8) + s(x1,k=6) + (1+x0|id)),
                  data = dat, family = poisson(), cores = 4, seed = 17,
                  warmup = 1000, iter = 2000, chains = 4, refresh=0, backend = "cmdstanr")
# model_brms <- brm(bf(y ~ s(x0,k=8) + s(x1,k=6)),
#                   data = dat, family = poisson(), cores = 4, seed = 17,
#                   warmup = 1000, iter = 2000, chains = 4, backend = "cmdstanr")

# model_brms
# tmp <- as_draws_df(model_brms)
# View(tmp)
tmp <- model.frame(model_brms)
formula(model_brms)
model_brms
save(model_brms, file="model_brms.rds")
library(brmsmargins)
source("marginalcoef.R")

mc <- marginalcoef(model_brms, CI = 0.95, posterior = T)
View(mc)



