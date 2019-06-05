## the following script roughly reproduces Figure 2 D
## in Gaynor 2018 "Mediation analysis for common binary outcomes"
## the results of this script suggest our code is working
rm(list=ls())
setwd("/Users/jplong/Desktop/set/funcs/tests")
library(boot)
source("../funcs.R")
source("../funcs_sim.R")
source("funcs_gaynor.R")

### issues with gaynor code
##  - article says simulating N(0.3A,0.75) but simulating N(0.3*0.4,0.75) for C
##  - sourcing the wrong script
##  - what does Args[[1]] store?
##  - should the line   const <-  seq(-5,1.5,0.005)[forInd] have a -3 instead of a -5
##  - why is the Re function needed?
  

## simulation parameters
snp_direct <- 0.4
names(snp_direct) <- "SNP1"
gs_direct <- 0.5
names(gs_direct) <- "gs1"
co_direct <- 0.25
names(co_direct) <- "co1"
eqtl <- matrix(1,ncol=1,nrow=1)
rownames(eqtl) <- names(snp_direct)
colnames(eqtl) <- names(gs_direct)
eqtl_model <- eqtl
eqtl_model[1,1] <- 0.5
co_gs <- matrix(0.4,ncol=1,nrow=1)
rownames(co_gs) <- names(co_direct)
colnames(co_gs) <- names(gs_direct)
const_gs <- 0.1
var_gs <- 0.75^2
family <- "binomial"



## loop across k
## loop across number of samples
ks <- seq(from=-3,to=1.5,length.out=100)
N <- 1000
n <- 500

## for storing results
probit_approx <- matrix(0,nrow=length(ks),ncol=N)
num_approx <- matrix(0,nrow=length(ks),ncol=N)
prevs_true <- matrix(0,nrow=length(ks),ncol=2)

for(ii in 1:length(ks)){
  print(ii)
  k <- ks[ii]
  dat <- SimulatePaper(1e6,k)
  prevs_true[ii,1] <- dat$prev
  dat2 <- ReformatData(dat)
  sim_params <- list(n=n,eqtl=eqtl,eqtl_model=eqtl_model,co_gs=co_gs,const_gs=const_gs,var_gs=var_gs,
                     snp_direct=snp_direct,gs_direct=gs_direct,co_direct=co_direct,const_direct=k,family=family)
  prevs_true[ii,2] <- ComputeEffectSNP(dat2,sim_params,"direct")
  for(jj in 1:N){
    dat <- SimulatePaper(n,k)
    dat2 <- ReformatData(dat)
    probit_approx[ii,jj] <- EstimateDirect(dat)
    fit <- ComputeDirecteQTL(dat2)
    num_approx[ii,jj] <- ComputeEffectSNP(dat2,fit,"direct")
  }
}

probit_approx_mean <- rowMeans(probit_approx)
num_approx_mean <- rowMeans(num_approx)

## see Figure 4b
ylim <- range(c(abs(prevs_true[,2]-probit_approx_mean),abs(prevs_true[,2]-num_approx_mean)))
plot(prevs_true[,1],abs(prevs_true[,2]-probit_approx_mean),ylim=ylim)
points(prevs_true[,1],abs(prevs_true[,2]-num_approx_mean),col='red')
abline(v=0.2)
abline(h=0.045)
## we have systematically lower bias, due to not making the probit approximation

## Notes
## 1.  The rare disease assumption would estimate direct effect at exp(theta1) = exp(0.4) = 1.49.
##     This is not a bad estimate, as the true values seem around 1.47
## 2. Two sources of bias in estimates: probit approximation and MLE parameter bias. which is causing
##    the error