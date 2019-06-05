## compare code to mediation package with survreg
## == challenge == 
## this comparison is difficult because the mediation
## package will only return effects on mean difference scale
## for parametric survreg functions in survival. the weibull
## distribution in survreg is a special case of coxph, but
## but the mean survival is sensitive to tails of the distribution
## where there is little data and it is difficult to estimate baseline
## hazard in the coxph
## but this measure is very sensitive to 
rm(list=ls())
setwd("/Users/jplong/Desktop/set/funcs/tests")
library(mediation)
library(survival)
library(pryr)
source("../funcs_sim.R")
source("../funcs.R")
set.seed(1234)

n <- 10000
family <- "cox"
sim_params <- OneDSim(n,n_snp_noise=0,n_gs_noise=0,n_co=0,family=family)
dat <- SimulateData(sim_params)


ds <- data.frame(as.matrix(dat$y),dat$snp,dat$gs)
summary(ds$time)
hist(ds$time)
head(ds)

model.m <- lm(gs1~SNP1,data=ds)
model.y <- survreg(Surv(time, status) ~ SNP1 + gs1, data=ds,dist="weibull")
summary(model.y)
out <- predict(model.y)
summary(out)


fit <- mediate(model.m, model.y,treat="SNP1",mediator="gs1",sims=1000)




rmean <- max(ds$time[ds$status==1])
rmean
fit_my <- ComputeDirecteQTL(dat)
fit_my$directfit
direct <- ComputeEffectSNP(dat,fit_my,"direct",risk_scale="diff",rmean=rmean)
indirect <- ComputeEffectSNP(dat,fit_my,"indirect",risk_scale="diff",rmean=rmean)
total <- ComputeEffectSNP(dat,fit_my,"total",risk_scale="diff",rmean=rmean)
direct
indirect
total

# direct effects
direct
fit$z0
##fit$z.avg

# indirect effects
indirect
##indirect
##[1] 0.005257949
fit$d1
##fit$d0
##fit$d.avg

# total effects
total
fit$tau.coef
fit$d1 + fit$z0

table(ds$status)

