#########
######### COMPUTATION FUNCTIONS

## computes path coefficients
##    arguments
##         dat$xx   :   n x p matrix with xx values
##          dat$mm   :   n x q matrix with gene set values
##        dat$path   :   p x q 0,1 matrix with path[i,j] = 1 if path between xx i and mm j
##          dat$co   :   n x r matrix with covariates, NULL if no covariates
##           dat$y   :   length n response vector
##      dat$family   :   gaussian for continuous y, binomial for binary y
##             reg   :   should direct effect estimates be regularized
##             mmr   :   should residuals from mm | xx regression be returned, if FALSE return NULL
##
##      value
##        list with path_model, direct, and covariates-gene set coefficients
ComputePath <- function(dat,reg=FALSE,mmr=FALSE){
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
  ## create matrix for storing mm residuals, if requested
  if(mmr){
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
      ## fitter warns when wald test statistic not useful, but we don't care
      out <- survival::coxph(y~x)
      temp <- coef(out)
      names(temp) <- colnames(x)
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
    if(mmr){
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
  return(list(xx_direct=xx_direct,mm_direct=mm_direct,co_direct=co_direct,const_direct=const_direct,
              path_model=path_model,co_mm=co_mm,const_mm=const_mm,var_mm=var_mm,mm_resid=mm_resid,directfit=directfit))
}

## input direct and path relations, outputs direct and indirect effects table
## WARNING: only appropriate for linear model where effects are sums of products
ComputeDirectIndirect <- function(path,xx_direct,mm_direct){
  ## indirect effects are product: sum_{gene set} (xx -> gene set) x (gene_set -> y)
  xx_indirect <- colSums(t(path)*mm_direct)
  ## add indirect effects of gene sets (by definition 0)
  indirect <- c(xx_indirect,rep(0,ncol(path)))
  ## put output in data frame
  direct <- c(xx_direct,mm_direct)
  eff <- cbind(direct,indirect,direct+indirect)
  rownames(eff) <- c(rownames(path),colnames(path))
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
  ## order by times with deaths as tie breakers
  ords <- order(ti,1-del)
  ti <- ti[ords]
  del <- del[ords]
  preds0 <- preds0[ords]
  di <- table(ti[del==1])
  r <- rev(cumsum(rev(exp(preds0))))
  r <- r[del==1]
  ti <- ti[del==1]
  nd <- !duplicated(ti)
  r <- r[nd]
  ti <- ti[nd]
  di <- di[as.character(ti)] ## TODO: dangerous?
  ##print(di)
  ##print(r)
  r <- cumsum(di / r)
  ##print(r)
  s <- exp(-r)
  ##print(s)
  ti <- c(0,ti)
  s <- c(1,s)
  names(s) <- NULL
  ##print(s)
  if(mean(s==1)==1){
    warning("estimated survival function identically 1. likely cause is numerical instability due to poor parameter estimates")
  }
  return(list(ti=ti,s=s))
}


## computes direct, indirect, or total effects of xx
##    arguments
##         dat   :   list of xx, mm, path, co, y, family
##         fit   :   output from ComputePath (path coefficients) or
##                   the true coefficients if we are approximating the true value
##  risk_scale   :   "diff" for risk difference, "ratio" for odds ratio
##         mmn   :   should gene set network be simulated
##
##      value
##        vector with all xx effects
ComputeEffectxx <- function(dat,fit,effect,risk_scale=NULL,mmn=FALSE,rmean=NULL){
  if(dat$family=="cox" & is.null(rmean)){
    stop("In ComputeEffectxx, rmean must be non NULL for dat$family=cox")
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
  if(is.null(dat$co)){
    dat$co <- matrix(0,nrow=nrow(dat$xx),ncol=0)
    fit$co_direct <- NULL
    fit$co_mm <- matrix(0,nrow=0,ncol=ncol(dat$path))
  }
  ## effects are functions of p_xy - p_ij where x,y,i,j are in 0,1
  ## here we determine which to compute
  if(effect=="direct"){
    dox1 <- 1
    doxpa1 <- 0
    dox0 <- 0
    doxpa0 <- 0
  }
  if(effect=="total"){
    dox1 <- 1
    doxpa1 <- 1
    dox0 <- 0
    doxpa0 <- 0
  }
  if(effect=="indirect"){
    dox1 <- 1
    doxpa1 <- 1
    dox0 <- 1
    doxpa0 <- 0
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
