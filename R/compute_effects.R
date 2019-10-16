#########
######### COMPUTATION FUNCTIONS


#' Computes path coefficients.
#'
#' \code{ComputePath} computes the coefficients for for all paths linking
#' response, mediators, treatment and potential confounders.
#'
#' @param dat List containing data (see SimulateData for format)
#' @param reg Should path coefficient estimates be regularized. Experimental.
#'   Defaults to FALSE.
#' @param mmn Should residuals from mm | xx regression be returned.
#'   If FALSE return NULL. Must be true if estimating direct/indirect effects
#'   with graph structure.
#' @return List containing path coefficients, variance estimates, residuals.
#' @examples
#' params <- SimpleSim()
#' dat <- SimulateData(params)
#' fit <- ComputePath(dat)
#' print(fit$xx_direct)
#' print(params$xx_direct)
#' @export
ComputePath <- function(dat,reg=FALSE,mmn=FALSE){
  ## unpack list
  xx <- dat$xx
  mm <- dat$mm
  y <- dat$y
  path <- dat$path
  co <- dat$co
  family <- dat$family
  if(is.null(co)){
    co <- matrix(0,nrow=nrow(xx),ncol=0)
  }
  ## argument checks
  if(nrow(xx)!=nrow(mm)){
    stop("mm and xx must have same number of rows")
  }
  if(length(y)!=nrow(mm)){
    stop("mm and y must have same number of rows")
  }
  if(!(family %in% c("gaussian","binomial","cox"))){
    stop("family must be either gaussian (for OLS), binomial (for logistic), or cox (with right censored survival y)")
  }
  if((nrow(xx) < ncol(xx) + ncol(mm) + ncol(co)) & !reg){
    stop("reg must be TRUE if p > n")
  }
  if(sum(is.na(co)) + sum(is.na(xx)) + sum(is.na(mm)) + sum(is.na(y)) > 0.5){
    stop("cannot handle NAs in dat")
  }
  ## create matrix for storing mm residuals, if requested
  if(mmn){
    mm_resid <- matrix(0,nrow=nrow(mm),ncol=ncol(mm))
  } else {
    mm_resid <- NULL
  }
  n_co <- ncol(co)
  ## store xx-mm coefficients in list before transforming to path_model
  coeffs <- vector("list",length=ncol(mm))
  ## store xx-mm constant coefficients and variances of fits
  const_mm <- rep(0,ncol(mm))
  var_mm <- rep(0,ncol(mm))
  ## regress y on everything to get direct effects
  x <- cbind(xx,mm,co)
  alpha <- 0.99 ## balance of L1 versus L2 used in glmnet
  if(reg){
    out <- glmnet::cv.glmnet(x,y,family=family,alpha=alpha)
    if(out$lambda[1]==out$lambda.1se){
      warning("need to choose larger lambdas for computing direct effects")
    }
    temp <- coef(out)[,1]
  } else{
    if(family=="cox"){
      ## fitter occasionally fails in simulations, we catch error and infer 0 coefficients
      out <- NULL
      try(out <- survival::coxph(y~x),silent=TRUE)
      if(is.null(out)){
        out <- "coxph fitting FAIL, likely due to high dimensionality. proceeding as if all coefficients 0"
        temp <- rep(0,ncol(x))
      }
      if(class(out)=="coxph"){
        temp <- coef(out)
      }
      names(temp) <- colnames(x)
      ## fitter warns when wald test statistic not useful, but we don't care
      if(sum(is.na(temp)>0)){
        temp[is.na(temp)] <- 0
        warning("coxph returned NA coefficient estimates. NAs imputed to 0. likely cause is singularity of design. consider using regularized regression.")
      }
    } else {
      out <- glm(y~x,family=family)
      temp <- coef(out)
      names(temp) <- c("intercept",colnames(x))
    }
  }
  directfit <- out
  ## cox survival model does not have intercept
  if(family!="cox"){
    direct <- temp[2:length(temp)]
    const_direct <- temp[1]
  }
  if(family=="cox"){
    direct <- temp
    const_direct <- 0
  }
  names(direct) <- c(colnames(xx),colnames(mm),colnames(co))
  for(ii in 1:ncol(mm)){
    x <- cbind(xx[,path[,ii]==1,drop=FALSE],co)
    y <- mm[,ii]
    ## compute regularized path relations when requested reg=TRUE
    ## AND there are at least 2 path xxs
    if(reg & ncol(x)>=2){
      out <- glmnet::cv.glmnet(x,y,family='gaussian',alpha=alpha)
      if(out$lambda[1]==out$lambda.1se){
        warning(paste0("need to choose larger lambdas for mm",ii))
      }
      temp <- coef(out)[,1]
      rs <- predict(out,newx=x) - y
    } else{
      if(ncol(x)>0){
        out <- lm(y~x)
        temp <- coef(out)
        rs <- predict(out) - y
      } else {
        temp <- mean(y)
        rs <- temp - y
      }
      names(temp) <- c("intercept",colnames(x))
    }
    const_mm[ii] <- temp[1]
    if(mmn){
      mm_resid[,ii] <- rs
    }
    var_mm[ii] <- mean(rs^2)
    if(ncol(x) > 0){
      coeffs[[ii]] <- temp[2:length(temp)]
    }
  }
  ## store coeffs in path_model and co_mm
  path_model <- path
  path_model[,] <- 0
  co_mm <- matrix(0,nrow=n_co,ncol=ncol(path))
  colnames(co_mm) <- colnames(mm)
  rownames(co_mm) <- colnames(co)
  for(ii in 1:ncol(mm)){
    if(!is.null(coeffs[[ii]])){
      te <- coeffs[[ii]]
      if(n_co>0){
        co_mm[,ii] <- te[(length(te)-n_co+1):length(te)]
      }
      if(length(te) > n_co){
        te <- te[1:(length(te)-n_co)]
        path_model[names(te),ii] <- te
      }
    }
  }
  ## break direct into xx, mm, co
  names(direct) <- c(colnames(xx),colnames(mm),colnames(co))
  xx_direct <- direct[1:nrow(path_model)]
  mm_direct <- direct[(nrow(path_model)+1):(length(direct)-n_co)]
  if(n_co > 0){
    co_direct <- direct[(length(direct)-n_co+1):length(direct)]
  } else {
    co_direct <- NULL
  }
  if(mmn){
    cov_mm <- cov(mm_resid)
  } else {
    cov_mm <- NULL
  }
  return(list(xx_direct=xx_direct,mm_direct=mm_direct,co_direct=co_direct,
              const_direct=const_direct,path_model=path_model,co_mm=co_mm,
              const_mm=const_mm,var_mm=var_mm,mm_resid=mm_resid,cov_mm=cov_mm,
              directfit=directfit))
}

#' Computes direct and indirect effects for linear models.
#'
#' \code{ComputeEffectsLinear} computes the direct and indirect effects
#' for linear models. These are simple functions of the path coefficients.
#' Argument can be path coefficients fitted to data, in format output by
#' \code{ComputePath} or simulation parameters output by function such as
#' \code{SimpleSim}.
#'
#' @param fit List with path coefficients.
#' @return Matrix of direct, indirect, and total effects for each xx.
#' @examples
#' set.seed(1234)
#' params <- SimpleSim()
#' dat <- SimulateData(params)
#' fit <- ComputePath(dat)
#' ComputeEffectsLinear(params) ## exact
#' ComputeEffectsLinear(fit) ## based on sample
#' @export
ComputeEffectsLinear <- function(fit){
  ## indirect effects are product: sum_{mediator} (xx -> mediator) x (mediator -> y)
  xx_indirect <- colSums(t(fit$path_model)*fit$mm_direct)
  ## add indirect effects of mediator sets (by definition 0)
  indirect <- c(xx_indirect,rep(0,ncol(fit$path_model)))
  direct <- c(fit$xx_direct,fit$mm_direct)
  eff <- cbind(direct,indirect,direct+indirect)
  rownames(eff) <- c(rownames(fit$path_model),colnames(fit$path_model))
  colnames(eff) <- c("direct","indirect","total")
  return(eff)
}

Computexxp <- function(dat,fit,ii,dox,doxpa,mmn,rmean){
  ## simulate mm data
  xx_s <- dat$xx
  xx_s[,ii] <- doxpa
  if(mmn){
    err <- MASS::mvrnorm(nrow(dat$mm),rep(0,ncol(dat$mm)),fit$cov_mm)
  } else {
    err <- matrix(rnorm(prod(dim(dat$mm)),mean=0,sd=sqrt(fit$var_mm)),
                  nrow=nrow(dat$mm),byrow=TRUE)
  }
  mm_s <- xx_s %*% fit$path_model + dat$co %*% fit$co_mm + err
  mm_s <- t(t(mm_s) + fit$const_mm)
  ## set xx to dox
  xx_s[,ii] <- dox
  x <- cbind(1,xx_s,mm_s,dat$co)
  direct <- c(fit$const_direct,fit$xx_direct,fit$mm_direct,fit$co_direct)
  preds <- colSums(t(x)*direct)
  ## transform to logistic scale
  if(dat$family=="binomial"){
    preds <- vapply(preds,function(x){1 / (1 + exp(-x))},c(0))
  }
  ## compute mean restricted survival
  if(dat$family=="cox"){
    ## get baseline predictions
    x <- cbind(1,dat$xx,dat$mm,dat$co)
    preds0 <- colSums(t(x)*direct)
    ## get baseline hazard using linear predictors and response
    H0 <- ComputeBaselineSurvival(dat$y,preds0)
    s <- H0$s
    ti <- H0$ti
    ## compute restricted mean using baseline survival s,
    ## times ti, and predictors pred
    ## 8.8.4 in klein
    s <- s[ti<rmean]
    ti <- ti[ti<rmean]
    ti <- c(ti,rmean)
    dti <- diff(ti)
    preds <- vapply(preds,function(pred){sum(dti*s^(exp(pred)))},c(0))
  }
  return(preds)
}

## compute baseline survival function, see 8.8 in
## survival analysis by klein, equations 8.8.1-8.8.3
ComputeBaselineSurvival <- function(y,preds0){
  ymat <- as.matrix(y)
  ti <- ymat[,1]
  del <- ymat[,2]
  ## compute 8.8.1
  ## order by times with deaths as tie breakers
  ords <- order(ti,1-del)
  ti <- ti[ords]
  del <- del[ords]
  preds0 <- preds0[ords]
  r <- rev(cumsum(rev(exp(preds0))))
  r <- r[del==1]
  ti <- ti[del==1]
  di <- table(ti)
  nd <- !duplicated(ti)
  r <- r[nd]
  ti <- ti[nd]
  di <- di[as.character(ti)]
  ## compute 8.8.2
  r <- cumsum(di / r)
  ## compute 8.8.3
  s <- exp(-r)
  ti <- c(0,ti)
  s <- c(1,s)
  names(s) <- NULL
  if(mean(s==1)==1){
    warning("estimated survival function identically 1. likely cause is numerical instability due to poor parameter estimates")
  }
  return(list(ti=ti,s=s))
}

#' Computes direct, indirect, or total effects of xx
#'
#' \code{ComputeEffectxx} computes either direct, indirect, or total effect
#' of each xx variable on response y when transition x from xp to xpp.
#'
#' @param dat List containing data (see SimulateData for format)
#' @param fit List with path coefficients.
#' @param effect Either 'direct', 'indirect', or 'total'
#' @param xp The initial value of x.
#' @param xpp The new value of x.
#' @param risk_scale Either 'diff' (for risk difference) or 'ratio' (for odds).
#'  see details for more information.
#' @param mmn Should mediator network be simulated. If FALSE assume
#'   conditionally independent mediators.
#' @return Vector of effects.
#' @examples
#' params <- SimpleSim()
#' dat <- SimulateData(params)
#' fit <- ComputePath(dat)
#' d <- ComputeEffectxx(dat,fit,"direct")
#' i <- ComputeEffectxx(dat,fit,"indirect")
#' to <- ComputeEffectxx(dat,fit,"total")
#' ## effects computed by estimating integrals
#' print(c(d,i,to))
#' ## effects using model linearity assumption
#' ComputeEffectsLinear(fit) ## based on sample
#' ## exact effects using simulation coefficients
#' ComputeEffectsLinear(params) ## exact
#' @export
ComputeEffectxx <- function(dat,fit,effect,xp=0,xpp=1,risk_scale=NULL,mmn=FALSE,rmean=NULL){
  if(is.null(fit$cov_mm) & mmn){
    stop("fit$cov_mm must be non-null when mmn=TRUE. either set mmn=FALSE or specify mmn=TRUE when calling ComputePaths to produce fit object.")
  }
  if(dat$family=="cox" & is.null(rmean)){
    stop("In ComputeEffectxx, rmean must be non NULL for dat$family=cox")
  }
  ## occasionally coxph fit fails. in this case class(fit) will
  ## not be coxph. we return 0 effects
  ## most useful in simulations
  if(dat$family=="cox" & class(fit$directfit)!="coxph"){
    return(rep(0,nrow(dat$path)))
  }

  ## by default gaussian family is computed on risk difference scale
  ## and binomial family is computed on odds ratio, but computing binomial
  ## on risk difference could also make sense and can be forced with risk_scale="diff"
  if(is.null(risk_scale)){
    if(dat$family=="gaussian" | dat$family=="cox"){
      risk_scale <- "diff"
    }
    if(dat$family=="binomial"){
      risk_scale <- "ratio"
    }
  }
  if(risk_scale=="ratio" & dat$family=="gaussian"){
    stop("Effect risk scale is ratio but model is gaussian linear. Incompatible.")
  }
  if(is.null(dat$co)){
    dat$co <- matrix(0,nrow=nrow(dat$xx),ncol=0)
    fit$co_direct <- NULL
    fit$co_mm <- matrix(0,nrow=0,ncol=ncol(dat$path))
  }
  ## effects are differences of e(i,j) where i,j are xp or xpp, see paper for details
  ## here we determine which to compute
  if(effect=="direct"){
    dox1 <- xpp
    doxpa1 <- xp
    dox0 <- xp
    doxpa0 <- xp
  }
  if(effect=="total"){
    dox1 <- xpp
    doxpa1 <- xpp
    dox0 <- xp
    doxpa0 <- xp
  }
  if(effect=="indirect"){
    dox1 <- xpp
    doxpa1 <- xpp
    dox0 <- xpp
    doxpa0 <- xp
  }
  total_effects <- rep(0,nrow(dat$path))
  for(ii in 1:nrow(dat$path)){
    preds0 <- Computexxp(dat,fit,ii,dox0,doxpa0,mmn,rmean)
    preds1 <- Computexxp(dat,fit,ii,dox1,doxpa1,mmn,rmean)
    ## report on the risk difference scale
    if(risk_scale=="diff"){
      total_effects[ii] <- mean(preds1) - mean(preds0)
    }
    ## report on the log odds scale
    if(risk_scale=="ratio"){
      total_effects[ii] <- (mean(preds1) / (1-mean(preds1))) / (mean(preds0) / (1-mean(preds0)))
    }
  }
  return(total_effects)
}


## subsets data, selecting only ix rows
## useful when bootstrap sampling
SubsetDat <- function(dat,ix){
  dat$y <- dat$y[ix]
  dat$mm <- dat$mm[ix,,drop=FALSE]
  dat$xx <- dat$xx[ix,,drop=FALSE]
  dat$co <- dat$co[ix,,drop=FALSE]
  return(dat)
}
