#### functions for simulating data and choosing
#### parameters for simulation

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
  if(is.null(params$dag)){
    err <- matrix(rnorm(params$n*ncol(params$eqtl),sd=sqrt(params$var_gs)),
                  nrow=params$n,byrow=TRUE)
    gs <- snp%*%params$eqtl_model + co%*%params$co_gs + err
  } else {
    dag <- params$dag
    Q <- solve(diag(1,nrow(dag$path)) - dag$path)
    cov_gs <- Q%*%diag(dag$var)%*%t(Q)
    err <- MASS::mvrnorm(params$n,rep(0,nrow(cov_gs)),cov_gs)
    mu <- snp%*%params$eqtl_model + co%*%params$co_gs
    mu <- t(Q%*%t(mu))
    gs <- mu + err
  }
  colnames(gs) <- colnames(params$eqtl)
  ## simulate response y:
  xbeta <- colSums(t(snp)*params$snp_direct) + colSums(t(co)*params$co_direct) + colSums(t(gs)*params$gs_direct)
  if(params$family=="gaussian"){
    y <-  xbeta + rnorm(params$n)
  }
  if(params$family=="binomial"){
    y <- rbinom(n=length(xbeta),size=1,prob=1/(1 + exp(-xbeta)))
  }
  ## simulate exponential survival times
  if(params$family=="cox"){
    tt <- rexp(length(xbeta),rate=(exp(xbeta)/1e4))
    death <- 1*(tt<=median(tt))
    tt[tt>median(tt)] <- median(tt)
    # wshape <- 10
    # tt <- rweibull(length(xbeta),shape=wshape,scale=25*exp(-xbeta/wshape))
    # cc <- rweibull(length(xbeta),shape=wshape,scale=25*exp(-xbeta/wshape))
    # cc[sample(c(TRUE,FALSE),replace=TRUE,size=length(cc))] <- Inf
    # death <- 1*(tt<=cc)
    # tt <- pmin(tt,cc)
    y <- survival::Surv(tt,death)
  }
  return(list(y=y,gs=gs,snp=snp,co=co,eqtl=params$eqtl,family=params$family))
}



## returns simulation parameters assuming a simple eqtl structure
# n = sample size
# nsnp = number of snps
# ngs = number of gene sets
# family = "gaussian" or "binomial" or "cox"
QuickSim <- function(n,nsnp,ngs,family,snp_prob=0.5){
  if(ngs < 2){
    stop("ngs (number of gene sets) must be at least 2")
  }
  if(nsnp < 3){
    stop("nsnp (number of snps) must be at least 3")
  }
  eqtl <- matrix(rbinom(nsnp*ngs,size=1,prob=.1),nrow=nsnp,ncol=ngs)
  snpnames <- paste0("SNP",1:nrow(eqtl))
  gsnames <- paste0("gs",1:ncol(eqtl))
  rownames(eqtl) <- snpnames
  colnames(eqtl) <- gsnames
  eqtl[1,1] <- 1
  eqtl[1,2] <- 0
  eqtl[2,1] <- 1
  eqtl[2,2] <- 1
  eqtl[3,1] <- 0
  eqtl[3,2] <- 1
  if(nsnp >=4){
    eqtl[4:nsnp,1:2] <- 0
  }
  ## coefficients linking snp with gs
  eqtl_model <- matrix(1/2*sample(c(-1,1),size=nrow(eqtl)*ncol(eqtl),replace=TRUE),
                       nrow=nrow(eqtl),ncol=ncol(eqtl))
  eqtl_model <- eqtl_model*eqtl
  ## fix eqtl model for relevant subgraph
  eqtl_model[1,1] <- 1/2
  eqtl_model[1,2] <- 0
  eqtl_model[2,1] <- -1/2
  eqtl_model[2,2] <- 1/2
  eqtl_model[3,1] <- 0
  eqtl_model[3,2] <- 1/2
  ## simulate empy matrices for covariates
  co_gs <- matrix(-1,nrow=0,ncol=ncol(eqtl))
  colnames(co_gs) <- colnames(eqtl)
  co_direct <- rep(-1,0)
  ## direct effects
  snp_direct <- c(rep(1,3),rep(0,nrow(eqtl)-3))
  names(snp_direct) <- snpnames
  gs_direct <- c(rep(3,2),rep(0,ncol(eqtl)-2))
  names(gs_direct) <- gsnames
  ## constant parameters
  const_gs <- rep(0,ncol(eqtl))
  const_direct <- 0
  ## gs variance
  var_gs <- rep(1,ncol(eqtl))
  sim_params <- list(n=n,eqtl=eqtl,eqtl_model=eqtl_model,co_gs=co_gs,const_gs=const_gs,var_gs=var_gs,
                     snp_direct=snp_direct,gs_direct=gs_direct,co_direct=co_direct,const_direct=const_direct,
                     snp_prob=snp_prob,
                     linkage="indep",family=family)
  return(sim_params)
}



## returns simulation parameters assuming a simple eqtl structure
# n = sample size
# nsnp = number of snps
# ngs = number of gene sets
# family = "gaussian" or "binomial"
QuickSim2 <- function(n,nsnp,ngs,family,snp_prob=0.5){
  if(ngs < 2){
    stop("ngs (number of gene sets) must be at least 2")
  }
  if(nsnp < 30){
    stop("nsnp (number of snps) must be at least 30")
  }
  ##eqtl <- matrix(rbinom(nsnp*ngs,size=1,prob=.2),nrow=nsnp,ncol=ngs)
  eqtl <- matrix(1,nrow=nsnp,ncol=ngs)
  snpnames <- paste0("SNP",1:nrow(eqtl))
  gsnames <- paste0("gs",1:ncol(eqtl))
  rownames(eqtl) <- snpnames
  colnames(eqtl) <- gsnames
  # eqtl[1:10,1] <- 1
  # eqtl[1,2] <- 0
  # eqtl[2,1] <- 1
  # eqtl[2,2] <- 1
  # eqtl[3,1] <- 0
  # eqtl[3,2] <- 1
  # if(nsnp >=4){
  #   eqtl[4:nsnp,1:2] <- 0
  # }
  ## coefficients linking snp with gs
  eqtl_model <- matrix(1/2*sample(c(-1,1,0,0,0,0),size=nrow(eqtl)*ncol(eqtl),replace=TRUE),
                       nrow=nrow(eqtl),ncol=ncol(eqtl))
  eqtl_model[1:10,1] <- 1/2*sample(c(-1,1),size=10,replace=TRUE)
  eqtl_model[11:20,1] <- 0
  eqtl_model <- eqtl_model*eqtl
  ## fix eqtl model for relevant subgraph
  eqtl_model[11:nsnp,1] <- 0
  ## simulate empy matrices for covariates
  co_gs <- matrix(-1,nrow=0,ncol=ncol(eqtl))
  colnames(co_gs) <- colnames(eqtl)
  co_direct <- rep(-1,0)
  ## direct effects
  snp_direct <- c(rep(0,10),rep(1,10),rep(0,nsnp-20))
  names(snp_direct) <- snpnames
  gs_direct <- c(1,rep(0,ncol(eqtl)-1))
  names(gs_direct) <- gsnames
  ## constant parameters
  const_gs <- rep(0,ncol(eqtl))
  const_direct <- 0
  ## gs variance
  var_gs <- rep(1,ncol(eqtl))
  sim_params <- list(n=n,eqtl=eqtl,eqtl_model=eqtl_model,co_gs=co_gs,const_gs=const_gs,var_gs=var_gs,
                     snp_direct=snp_direct,gs_direct=gs_direct,co_direct=co_direct,const_direct=const_direct,
                     snp_prob=snp_prob,
                     linkage="indep",family=family)
  return(sim_params)
}


## returns simulation parameters assuming a simple eqtl structure
# n = sample size
# ngs = number of gene sets
QuickSimMultipleMediator <- function(n,ngs,family,
                                     var_gs=10,snp_direct=0.5,gs_direct=.5,snp_prob=0.2){
  ngs_min <- 5
  if(ngs < ngs_min){
    stop(paste0("ngs (number of gene sets) must be at least ",ngs_min))
  }
  eqtl <- matrix(1,ncol=ngs,nrow=1)
  snpnames <- paste0("SNP",1:nrow(eqtl))
  gsnames <- paste0("gs",1:ncol(eqtl))
  rownames(eqtl) <- snpnames
  colnames(eqtl) <- gsnames
  ## coefficients linking snp with gs
  eqtl_model <- matrix(c(rep(1.0,ngs_min),rep(0,ncol(eqtl)-ngs_min)),nrow=nrow(eqtl),ncol=ncol(eqtl))
  colnames(eqtl_model) <- gsnames
  rownames(eqtl_model) <- snpnames
  ## simulate empy matrices for covariates
  co_gs <- matrix(-1,nrow=0,ncol=ncol(eqtl))
  colnames(co_gs) <- colnames(eqtl)
  co_direct <- rep(-1,0)
  ## direct effects
  names(snp_direct) <- snpnames
  gs_direct <- c(rep(gs_direct,ngs_min),rep(0,ncol(eqtl)-ngs_min))
  names(gs_direct) <- gsnames
  ## constant parameters
  const_gs <- rep(0,ncol(eqtl))
  const_direct <- 0
  ## gs variance
  var_gs <- rep(var_gs,ncol(eqtl))
  sim_params <- list(n=n,eqtl=eqtl,eqtl_model=eqtl_model,co_gs=co_gs,const_gs=const_gs,var_gs=var_gs,
                     snp_direct=snp_direct,gs_direct=gs_direct,co_direct=co_direct,const_direct=const_direct,
                     snp_prob=snp_prob,
                     linkage="indep",family=family)
  return(sim_params)
}

## returns simulation parameters
## SNP1 and gs: both have direct causes
# n = sample size
#
# family = "gaussian" or "binomial" or "cox"
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
  eqtl_model[1,1] <- 1
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
  snp_direct <- c(1,rep(0,n_snp_noise))
  names(snp_direct) <- snpnames
  gs_direct <- c(2,rep(0,n_gs_noise))
  names(gs_direct) <- gsnames
  ## constant parameters
  const_gs <- rep(0,ncol(eqtl))
  const_direct <- 0
  ## gs variance
  var_gs <- rep(1,ncol(eqtl))
  sim_params <- list(n=n,eqtl=eqtl,eqtl_model=eqtl_model,co_gs=co_gs,const_gs=const_gs,var_gs=var_gs,
                     snp_direct=snp_direct,gs_direct=gs_direct,co_direct=co_direct,const_direct=const_direct,
                     snp_prob=0.5,
                     linkage="indep",family=family)
  return(sim_params)
}
