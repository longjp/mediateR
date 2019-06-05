## comparison of mediation package in R
## to some of my code
rm(list=ls())
setwd("~/Desktop/set/funcs/tests")
library(mediation)
source("../funcs.R")

## returns simulation parameters
## SNP1 and gs: both have direct causes
# n = sample size
#
# family = "gaussian" or "binomial"
OneDSim <- function(n,n_snp_noise=0,n_gs_noise=0,n_co=1,family="gaussian"){
  eqtl <- matrix(0,nrow=1+n_snp_noise,ncol=1+n_gs_noise)
  snpnames <- paste0("SNP",1:nrow(eqtl))
  gsnames <- paste0("gs",1:ncol(eqtl))
  rownames(eqtl) <- snpnames
  colnames(eqtl) <- gsnames
  eqtl[1,1] <- 1
  ## coefficients linking snp with gs
  eqtl_model <- eqtl
  eqtl_model[,] <- 0
  eqtl_model[1,1] <- -1.5
  ## covariates
  conames <- paste0("c",1:n_co)
  co_gs <- matrix(-1,nrow=n_co,ncol=ncol(eqtl))
  if(n_co > 0){
    rownames(co_gs) <- conames
  }
  colnames(co_gs) <- colnames(eqtl)
  ## direct effects
  co_direct <- rep(-1,n_co)
  if(n_co > 0){
    names(co_direct) <- conames
  }
  snp_direct <- c(.5,rep(0,n_snp_noise))
  names(snp_direct) <- snpnames
  gs_direct <- c(.8,rep(0,n_gs_noise))
  names(gs_direct) <- gsnames
  ## constant parameters
  const_gs <- rep(0,ncol(eqtl))
  const_direct <- 2
  sim_params <- list(n=n,eqtl=eqtl,eqtl_model=eqtl_model,co_gs=co_gs,const_gs=const_gs,
                     snp_direct=snp_direct,gs_direct=gs_direct,co_direct=co_direct,const_direct=const_direct,
                     snp_prob=0.3,
                     linkage="indep",family=family)
  return(sim_params)
}


## simulate data
SimulateData <- function(params){
  ## for now either simulate old linkage structure or do independent bernoullis 
  if(params$linkage=="empty" & nrow(params$eqtl)==3){
    s1 <- rnorm(params$n,mean=-.3)
    s2 <- sqrt(.5)*s1 + sqrt(.5)*rnorm(params$n,mean=-.3)
    s3 <- sqrt(1/3)*s1 + sqrt(1/3)*s2 + sqrt(1/3)*rnorm(params$n,mean=-.3)
    snp <- cbind(1*(s1>0),1*(s2>0),1*(s3>0))
  } else {
    snp <- matrix(rbinom(params$n*nrow(params$eqtl),size=1,prob=params$snp_prob),
                  ncol=nrow(params$eqtl))
  }
  colnames(snp) <- rownames(params$eqtl)
  ## simulate covariates
  co <- matrix(rnorm(params$n*nrow(params$co_gs)),nrow=params$n)
  colnames(co) <- rownames(params$co_gs)
  ## create gene expressions from snp, eqtl_model
  gs <- snp%*%params$eqtl_model + co%*%params$co_gs + matrix(rnorm(params$n*ncol(params$eqtl),sd=.8),nrow=params$n)
  colnames(gs) <- colnames(params$eqtl)
  ## simulate response y:
  xbeta <- colSums(t(snp)*params$snp_direct) + colSums(t(co)*params$co_direct) + colSums(t(gs)*params$gs_direct)
  if(params$family=="gaussian"){
    y <-  xbeta + rnorm(params$n)
  }
  if(params$family=="binomial"){
    y <- rbinom(n=length(xbeta),size=1,prob=1/(1 + exp(-xbeta)))
  }
  return(list(y=y,gs=gs,snp=snp,co=co,family=params$family,eqtl=params$eqtl))
}




######## simulate gaussian model

sim_params <- OneDSim(500000,n_snp_noise=0,n_gs_noise=0,n_co=0)
dat <- SimulateData(sim_params)


## snp direct =1, gs direct =2, snp on gs is 1
d <- data.frame(y=dat$y,gs=dat$gs[,1],snp=dat$snp[,1])
model.y <- lm(y~.,data=d)
model.m <- lm(gs~snp,data=d)
model.y
model.m
fit <- mediate(model.m,model.y,sims=10,treat="snp",mediator="gs")

## use my code
fit_my <- ComputeDirecteQTL(dat)
direct <- ComputeEffectSNP(dat,fit_my,"direct")
indirect <- ComputeEffectSNP(dat,fit_my,"indirect")
total <- ComputeEffectSNP(dat,fit_my,"total")

## compare output
# direct effects
direct
fit$z0
# indirect effects
indirect
fit$d0
# total effects
total
fit$tau.coef





######## simulate logistic model
sim_params <- OneDSim(500000,n_snp_noise=0,n_gs_noise=0,n_co=0,
                      family="binomial")
dat <- SimulateData(sim_params)

## snp direct =1, gs direct =2, snp on gs is 1
d <- data.frame(y=dat$y,gs=dat$gs[,1],snp=dat$snp[,1])
model.y <- glm(y~.,data=d,family=binomial)
model.m <- lm(gs~snp,data=d)
model.y
model.m
fit <- mediate(model.m, model.y, sims=10, treat="snp", mediator="gs")

## use my code
fit_my <- ComputeDirecteQTL(dat)
str(fit_my)
fit_my$snp_direct
fit_my$gs_direct
fit_my$eqtl_model


direct <- ComputeEffectSNP(dat,fit_my,"direct",risk_scale="diff")
indirect <- ComputeEffectSNP(dat,fit_my,"indirect",risk_scale="diff")
total <- ComputeEffectSNP(dat,fit_my,"total",risk_scale="diff")

## compare output
# direct effects
direct
fit$z0
##fit$z1
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
