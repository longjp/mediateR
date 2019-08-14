#' Generate simulation data.
#'
#' \code{SimulateData} creates simulation data set using specified
#' input parameters. Output list is \code{dat} argument for fitting
#' path coefficients and direct and indirect effects.
#'
#' Output list used as \code{dat} argument for \code{ComputePath}
#' and \code{ComputEffectxx}. Has named elements:
#'
#' \itemize{
#'
#'  \item    \code{y} :   response vector of length n, either double or Surv object
#'
#'  \item     \code{mm} :   n x q mediator matrix
#'
#'  \item    \code{xx} :   n x p treamtment matrix
#'
#'  \item     \code{co} :   n x r matrix of confounders
#'
#'  \item   \code{path} :   p x q matrix of 0s and 1s denoting paths to be fit.
#'                        path[i,j]=0 indicates a priori
#'                        knowledge of no path xx[i] and mm[j]
#'
#'  \item \code{family} : 'gaussian', 'binomial', or 'cox' indicating
#'                        whether a linear, logistic, or cox model
#'                        should be fit with the y on xx,mm,co regression
#'  }
#'
#' @param params Simulation parameters (see SimpleSim for format)
#' @return List containing simulation data.
#' @examples
#' set.seed(1234)
#' params <- SimpleSim()
#' dat <- SimulateData(params)
#' @export
SimulateData <- function(params){
  ## for now either simulate old linkage structure or do independent bernoullis
  if(params$linkage=="empty" & nrow(params$path)==3){
    s1 <- rnorm(params$n,mean=-.3)
    s2 <- sqrt(.5)*s1 + sqrt(.5)*rnorm(params$n,mean=-.3)
    s3 <- sqrt(1/3)*s1 + sqrt(1/3)*s2 + sqrt(1/3)*rnorm(params$n,mean=-.3)
    xx <- cbind(1*(s1>0),1*(s2>0),1*(s3>0))
  } else {
    xx <- matrix(rbinom(params$n*nrow(params$path),size=1,prob=params$xx_prob),
                  ncol=nrow(params$path))
  }
  colnames(xx) <- rownames(params$path)
  ## simulate covariates
  co <- matrix(rnorm(params$n*nrow(params$co_mm)),nrow=params$n)
  colnames(co) <- rownames(params$co_mm)
  ## create gene expressions from xx, path_model
  if(is.null(params$dag)){
    err <- matrix(rnorm(params$n*ncol(params$path),sd=sqrt(params$var_mm)),
                  nrow=params$n,byrow=TRUE)
    mm <- xx%*%params$path_model + co%*%params$co_mm + err
  } else {
    dag <- params$dag
    Q <- solve(diag(1,nrow(dag$path)) - dag$path)
    cov_mm <- Q%*%diag(dag$var)%*%t(Q)
    err <- MASS::mvrnorm(params$n,rep(0,nrow(cov_mm)),cov_mm)
    mu <- xx%*%params$path_model + co%*%params$co_mm
    mu <- t(Q%*%t(mu))
    mm <- mu + err
  }
  colnames(mm) <- colnames(params$path)
  ## simulate response y:
  xbeta <- (colSums(t(xx)*params$xx_direct) + colSums(t(co)*params$co_direct) +
              colSums(t(mm)*params$mm_direct))
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
  return(list(y=y,mm=mm,xx=xx,co=co,path=params$path,family=params$family))
}

#' Generate simulation parameters.
#'
#' \code{SimpleSim} creates simulation parameters which are then used
#'   to generate data with \code{SimulateData}. Mostly used for testing.
#'
#' @param n Sample Size
#' @param family Either 'gaussian', 'binomial', or 'cox' depending on
#'   desired response model.
#' @return List that the parameterizes simulation. Typically used as
#'   argument in subsequent call to \code{SimulateData}
#' @examples
#' params <- SimpleSim()
#' names(params)
#' @export
SimpleSim <- function(n=500,family="gaussian"){
  xxnames <- "x"
  mmnames <- "m"
  xx_direct <- 2
  names(xx_direct) <- xxnames
  mm_direct <- 3
  names(mm_direct) <- mmnames
  ## create path variables
  path <- matrix(1,nrow=1,ncol=1)
  rownames(path) <- xxnames
  colnames(path) <- mmnames
  path_model <- matrix(2.5,nrow=1,ncol=1)
  rownames(path_model) <- xxnames
  colnames(path_model) <- mmnames
  co_mm <- matrix(-0.5,nrow=1,ncol=1)
  colnames(co_mm) <- colnames(path)
  rownames(co_mm) <- rownames(path)
  co_direct <- 1.5
  names(co_direct) <- "co"
  ## constant parameters
  const_mm <- pi/2
  const_direct <- pi
  ## mm variance
  var_mm <- 1
  xx_prob <- 0.5
  sim_params <- list(n=n,path=path,path_model=path_model,
                     co_mm=co_mm,const_mm=const_mm,var_mm=var_mm,
                     xx_direct=xx_direct,mm_direct=mm_direct,
                     co_direct=co_direct,const_direct=const_direct,
                     xx_prob=xx_prob,
                     linkage="indep",family=family)
  return(sim_params)
}





## returns simulation parameters assuming a simple path structure
# n = sample size
# nxx = number of xxs
# nmm = number of gene sets
# family = "gaussian" or "binomial" or "cox"
QuickSim <- function(n,nxx,nmm,family,xx_prob=0.5,const_direct=0){
  if(nmm < 2){
    stop("nmm (number of gene sets) must be at least 2")
  }
  if(nxx < 3){
    stop("nxx (number of xxs) must be at least 3")
  }
  path <- matrix(rbinom(nxx*nmm,size=1,prob=.1),nrow=nxx,ncol=nmm)
  xxnames <- paste0("x",1:nrow(path))
  mmnames <- paste0("m",1:ncol(path))
  rownames(path) <- xxnames
  colnames(path) <- mmnames
  path[1,1] <- 1
  path[1,2] <- 0
  path[2,1] <- 1
  path[2,2] <- 1
  path[3,1] <- 0
  path[3,2] <- 1
  if(nxx >=4){
    path[4:nxx,1:2] <- 0
  }
  ## coefficients linking xx with mm
  path_model <- matrix(1/2*sample(c(-1,1),size=nrow(path)*ncol(path),
                                  replace=TRUE),
                       nrow=nrow(path),ncol=ncol(path))
  path_model <- path_model*path
  ## fix path model for relevant subgraph
  path_model[1,1] <- 1/2
  path_model[1,2] <- 0
  path_model[2,1] <- -1/2
  path_model[2,2] <- 1/2
  path_model[3,1] <- 0
  path_model[3,2] <- 1/2
  ## simulate empy matrices for covariates
  co_mm <- matrix(-1,nrow=0,ncol=ncol(path))
  colnames(co_mm) <- colnames(path)
  co_direct <- rep(-1,0)
  ## direct effects
  xx_direct <- c(rep(1,3),rep(0,nrow(path)-3))
  names(xx_direct) <- xxnames
  mm_direct <- c(rep(3,2),rep(0,ncol(path)-2))
  names(mm_direct) <- mmnames
  ## constant parameters
  const_mm <- rep(0,ncol(path))
  ## mm variance
  var_mm <- rep(1,ncol(path))
  sim_params <- list(n=n,path=path,path_model=path_model,
                     co_mm=co_mm,const_mm=const_mm,var_mm=var_mm,
                     xx_direct=xx_direct,mm_direct=mm_direct,
                     co_direct=co_direct,const_direct=const_direct,
                     xx_prob=xx_prob,
                     linkage="indep",family=family)
  return(sim_params)
}



## returns simulation parameters assuming a simple path structure
# n = sample size
# nxx = number of xxs
# nmm = number of gene sets
# family = "gaussian" or "binomial"
QuickSim2 <- function(n,nxx,nmm,family,xx_prob=0.5){
  if(nmm < 2){
    stop("nmm (number of gene sets) must be at least 2")
  }
  if(nxx < 30){
    stop("nxx (number of xxs) must be at least 30")
  }
  ##path <- matrix(rbinom(nxx*nmm,size=1,prob=.2),nrow=nxx,ncol=nmm)
  path <- matrix(1,nrow=nxx,ncol=nmm)
  xxnames <- paste0("x",1:nrow(path))
  mmnames <- paste0("m",1:ncol(path))
  rownames(path) <- xxnames
  colnames(path) <- mmnames
  # path[1:10,1] <- 1
  # path[1,2] <- 0
  # path[2,1] <- 1
  # path[2,2] <- 1
  # path[3,1] <- 0
  # path[3,2] <- 1
  # if(nxx >=4){
  #   path[4:nxx,1:2] <- 0
  # }
  ## coefficients linking xx with mm
  path_model <- matrix(1/2*sample(c(-1,1,0,0,0,0),size=nrow(path)*ncol(path),
                                  replace=TRUE),
                       nrow=nrow(path),ncol=ncol(path))
  path_model[1:10,1] <- 1/2*sample(c(-1,1),size=10,replace=TRUE)
  path_model[11:20,1] <- 0
  path_model <- path_model*path
  ## fix path model for relevant subgraph
  path_model[11:nxx,1] <- 0
  ## simulate empy matrices for covariates
  co_mm <- matrix(-1,nrow=0,ncol=ncol(path))
  colnames(co_mm) <- colnames(path)
  co_direct <- rep(-1,0)
  ## direct effects
  xx_direct <- c(rep(0,10),rep(1,10),rep(0,nxx-20))
  names(xx_direct) <- xxnames
  mm_direct <- c(1,rep(0,ncol(path)-1))
  names(mm_direct) <- mmnames
  ## constant parameters
  const_mm <- rep(0,ncol(path))
  const_direct <- 0
  ## mm variance
  var_mm <- rep(1,ncol(path))
  sim_params <- list(n=n,path=path,path_model=path_model,
                     co_mm=co_mm,const_mm=const_mm,var_mm=var_mm,
                     xx_direct=xx_direct,mm_direct=mm_direct,
                     co_direct=co_direct,const_direct=const_direct,
                     xx_prob=xx_prob,
                     linkage="indep",family=family)
  return(sim_params)
}


#' Generate simulation parameters assuming a simple path structure
#'
#' \code{SimpleSim} creates simulation parameters which are then used
#'   to generate data with \code{SimulateData}. Mostly used for testing.
#'
#' @param n Sample Size
#' @param nmm Number of mediators
#' @param family Response type
#' @param var_mm Variance of the mediators
#' @param xx_direct Path coefficient form x to response
#' @param mm_direct Path coefficient from mediators to response
#' @param xx_prob Probability x is 1.
#' @return List that the parameterizes simulation. Typically used as
#'   argument in subsequent call to \code{SimulateData}
#' @examples
#' params <- QuickSimMultipleMediator(500,5,"cox")
#' names(params)
#' @export
QuickSimMultipleMediator <- function(n,nmm,family,
                                     var_mm=10,xx_direct=0.5,
                                     mm_direct=.5,xx_prob=0.2){
  nmm_min <- 5
  if(nmm < nmm_min){
    stop(paste0("nmm (number of mediators) must be at least ",nmm_min))
  }
  path <- matrix(1,ncol=nmm,nrow=1)
  xxnames <- paste0("x",1:nrow(path))
  mmnames <- paste0("m",1:ncol(path))
  rownames(path) <- xxnames
  colnames(path) <- mmnames
  ## coefficients linking xx with mm
  path_model <- matrix(c(rep(1.0,nmm_min),rep(0,ncol(path)-nmm_min)),
                       nrow=nrow(path),ncol=ncol(path))
  colnames(path_model) <- mmnames
  rownames(path_model) <- xxnames
  ## simulate empy matrices for covariates
  co_mm <- matrix(-1,nrow=0,ncol=ncol(path))
  colnames(co_mm) <- colnames(path)
  co_direct <- rep(-1,0)
  ## direct effects
  names(xx_direct) <- xxnames
  mm_direct <- c(rep(mm_direct,nmm_min),rep(0,ncol(path)-nmm_min))
  names(mm_direct) <- mmnames
  ## constant parameters
  const_mm <- rep(0,ncol(path))
  const_direct <- 0
  ## mm variance
  var_mm <- rep(var_mm,ncol(path))
  sim_params <- list(n=n,path=path,path_model=path_model,
                     co_mm=co_mm,const_mm=const_mm,var_mm=var_mm,
                     xx_direct=xx_direct,mm_direct=mm_direct,
                     co_direct=co_direct,const_direct=const_direct,
                     xx_prob=xx_prob,
                     linkage="indep",family=family)
  return(sim_params)
}

## returns simulation parameters
## xx1 and mm: both have direct causes
# n = sample size
#
# family = "gaussian" or "binomial" or "cox"
OneDSim <- function(n,n_xx_noise=0,n_mm_noise=0,n_co=1,family="gaussian"){
  path <- matrix(0,nrow=1+n_xx_noise,ncol=1+n_mm_noise)
  xxnames <- paste0("x",1:nrow(path))
  mmnames <- paste0("m",1:ncol(path))
  rownames(path) <- xxnames
  colnames(path) <- mmnames
  path[1,1] <- 1
  ## coefficients linking xx with mm
  path_model <- path
  path_model[,] <- 0
  path_model[1,1] <- 1
  ## covariates
  conames <- paste0("c",1:n_co)
  co_mm <- matrix(-1,nrow=n_co,ncol=ncol(path))
  if(n_co > 0){
    rownames(co_mm) <- conames
  }
  colnames(co_mm) <- colnames(path)
  ## direct effects
  co_direct <- rep(-1,n_co)
  if(n_co > 0){
    names(co_direct) <- conames
  }
  xx_direct <- c(1,rep(0,n_xx_noise))
  names(xx_direct) <- xxnames
  mm_direct <- c(2,rep(0,n_mm_noise))
  names(mm_direct) <- mmnames
  ## constant parameters
  const_mm <- rep(0,ncol(path))
  const_direct <- 0
  ## mm variance
  var_mm <- rep(1,ncol(path))
  sim_params <- list(n=n,path=path,path_model=path_model,
                     co_mm=co_mm,const_mm=const_mm,var_mm=var_mm,
                     xx_direct=xx_direct,mm_direct=mm_direct,
                     co_direct=co_direct,const_direct=const_direct,
                     xx_prob=0.5,
                     linkage="indep",family=family)
  return(sim_params)
}
