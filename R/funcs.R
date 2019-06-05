#########
######### COMPUTATION FUNCTIONS

## computes path coefficients
##    arguments
##         dat$snp   :   n x p matrix with snp values
##          dat$gs   :   n x q matrix with gene set values
##        dat$eqtl   :   p x q 0,1 matrix with eqtl[i,j] = 1 if eqtl between snp i and gs j
##          dat$co   :   n x r matrix with covariates, NULL if no covariates
##           dat$y   :   length n response vector
##      dat$family   :   gaussian for continuous y, binomial for binary y
##             reg   :   should direct effect estimates be regularized
##             gsr   :   should residuals from gs | snp regression be returned, if FALSE return NULL
##
##      value
##        list with eqtl_model, direct, and covariates-gene set coefficients
ComputeDirecteQTL <- function(dat,reg=FALSE,gsr=FALSE){
  ## unpack list
  snp <- dat$snp
  gs <- dat$gs
  y <- dat$y
  eqtl <- dat$eqtl
  co <- dat$co
  family <- dat$family
  if(is.null(co)){
    co <- matrix(0,nrow=nrow(snp),ncol=0)
  }
  ## argument checks
  if(nrow(snp)!=nrow(gs)){
    stop("gs and snp must have same number of rows")
  }
  if(length(y)!=nrow(gs)){
    stop("gs and y must have same number of rows")
  }
  if(!(family %in% c("gaussian","binomial","cox"))){
    stop("family must be either gaussian (for OLS), binomial (for logistic), or cox (with right censored survival y)")
  }
  if((nrow(snp) < ncol(snp) + ncol(gs) + ncol(co)) & !reg){
    stop("reg must be TRUE if p > n")
  }
  ## create matrix for storing gs residuals, if requested
  if(gsr){
    gs_resid <- matrix(0,nrow=nrow(gs),ncol=ncol(gs))
  } else {
    gs_resid <- NULL
  }
  n_co <- ncol(co)
  ## store snp-gs coefficients in list before transforming to eqtl_model
  coeffs <- vector("list",length=ncol(gs))
  ## store snp-gs constant coefficients and variances of fits
  const_gs <- rep(0,ncol(gs))
  var_gs <- rep(0,ncol(gs))
  ## regress y on everything to get direct effects
  x <- cbind(snp,gs,co)
  alpha <- 0.99 ## balance of L1 versus L2 used in glmnet
  if(reg){
    out <- cv.glmnet(x,y,family=family,alpha=alpha)
    if(out$lambda[1]==out$lambda.1se){
      warning("need to choose larger lambdas for computing direct effects")
    }
    temp <- coef(out)[,1]
  } else{
    if(family=="cox"){
      ## fitter warns when wald test statistic not useful, but we don't care
      out <- coxph(y~x)
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
  names(direct) <- c(colnames(snp),colnames(gs),colnames(co))
  for(ii in 1:ncol(gs)){
    x <- cbind(snp[,eqtl[,ii]==1,drop=FALSE],co)
    y <- gs[,ii]
    ## compute regularized eqtl relations when requested reg=TRUE
    ## AND there are at least 2 eqtl SNPs
    if(reg & ncol(x)>=2){
      out <- cv.glmnet(x,y,family='gaussian',alpha=alpha)
      if(out$lambda[1]==out$lambda.1se){
        warning(paste0("need to choose larger lambdas for gs",ii))
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
    const_gs[ii] <- temp[1]
    if(gsr){
      gs_resid[,ii] <- rs
    }
    var_gs[ii] <- mean(rs^2)
    if(ncol(x) > 0){
      coeffs[[ii]] <- temp[2:length(temp)]
    }
  }
  ## store coeffs in eqtl_model and co_gs
  eqtl_model <- eqtl
  eqtl_model[,] <- 0
  co_gs <- matrix(0,nrow=n_co,ncol=ncol(eqtl))
  colnames(co_gs) <- colnames(gs)
  rownames(co_gs) <- colnames(co)
  for(ii in 1:ncol(gs)){
    if(!is.null(coeffs[[ii]])){
      te <- coeffs[[ii]]
      if(n_co>0){
        co_gs[,ii] <- te[(length(te)-n_co+1):length(te)]
      }
      if(length(te) > n_co){
        te <- te[1:(length(te)-n_co)]
        eqtl_model[names(te),ii] <- te
      }
    }
  }
  ## break direct into snps, gs, co
  names(direct) <- c(colnames(snp),colnames(gs),colnames(co))
  snp_direct <- direct[1:nrow(eqtl_model)]
  gs_direct <- direct[(nrow(eqtl_model)+1):(length(direct)-n_co)]
  if(n_co > 0){
    co_direct <- direct[(length(direct)-n_co+1):length(direct)]
  } else {
    co_direct <- NULL
  }
  return(list(snp_direct=snp_direct,gs_direct=gs_direct,co_direct=co_direct,const_direct=const_direct,
              eqtl_model=eqtl_model,co_gs=co_gs,const_gs=const_gs,var_gs=var_gs,gs_resid=gs_resid,directfit=directfit))
}

## input direct and eqtl relations, outputs direct and indirect effects table
## WARNING: only appropriate for linear model where effects are sums of products
ComputeDirectIndirect <- function(eqtl,snp_direct,gs_direct){
  ## indirect effects are product: sum_{gene set} (SNP -> gene set) x (gene_set -> y)
  snp_indirect <- colSums(t(eqtl)*gs_direct)
  ## add indirect effects of gene sets (by definition 0)
  indirect <- c(snp_indirect,rep(0,ncol(eqtl)))
  ## put output in data frame
  direct <- c(snp_direct,gs_direct)
  eff <- cbind(direct,indirect,direct+indirect)
  rownames(eff) <- c(rownames(eqtl),colnames(eqtl))
  colnames(eff) <- c("direct","indirect","total")
  return(eff)
}


ComputeSNPp <- function(dat,fit,ii,dox,doxpa,gsn,rmean){
  ## simulate gs data
  snp_s <- dat$snp
  snp_s[,ii] <- doxpa
  if(gsn){
    err <- mvrnorm(nrow(dat$gs),rep(0,ncol(dat$gs)),fit$cov_gs)
  } else {
    err <- matrix(rnorm(prod(dim(dat$gs)),mean=0,sd=sqrt(fit$var_gs)),
                  nrow=nrow(dat$gs),byrow=TRUE)
  }
  gs_s <- snp_s %*% fit$eqtl_model + dat$co %*% fit$co_gs + err
  gs_s <- t(t(gs_s) + fit$const_gs)
  ## set snp to dox
  snp_s[,ii] <- dox
  x <- cbind(1,snp_s,gs_s,dat$co)
  direct <- c(fit$const_direct,fit$snp_direct,fit$gs_direct,fit$co_direct)
  preds <- colSums(t(x)*direct)
  ## transform to logistic scale
  if(dat$family=="binomial"){
    preds <- vapply(preds,function(x){1 / (1 + exp(-x))},c(0))
  }
  ## compute mean restricted survival
  if(dat$family=="cox"){
    ## get baseline predictions
    x <- cbind(1,dat$snp,dat$gs,dat$co)
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


## computes direct, indirect, or total effects of snp
##    arguments
##         dat   :   list of snp, gs, eqtl, co, y, family
##         fit   :   output from ComputeDirecteQTL (path coefficients) or
##                   the true coefficients if we are approximating the true value
##  risk_scale   :   "diff" for risk difference, "ratio" for odds ratio
##         gsn   :   should gene set network be simulated
##
##      value
##        vector with all SNP effects
ComputeEffectSNP <- function(dat,fit,effect,risk_scale=NULL,gsn=FALSE,rmean=NULL){
  if(dat$family=="cox" & is.null(rmean)){
    stop("In ComputeEffectSNP, rmean must be non NULL for dat$family=cox")
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
    dat$co <- matrix(0,nrow=nrow(dat$snp),ncol=0)
    fit$co_direct <- NULL
    fit$co_gs <- matrix(0,nrow=0,ncol=ncol(dat$eqtl))
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
  total_effects <- rep(0,nrow(dat$eqtl))
  for(ii in 1:nrow(dat$eqtl)){
    preds0 <- ComputeSNPp(dat,fit,ii,dox0,doxpa0,gsn,rmean)
    preds1 <- ComputeSNPp(dat,fit,ii,dox1,doxpa1,gsn,rmean)
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


## determines covariance of dag based on path coefficients
## and variances
ComputeDAGCov <- function(dag){
  Q <- solve(diag(1,nrow(dag$path)) - dag$path)
  return(Q%*%diag(dag$var)%*%t(Q))
}



## subsets data, selecting only ix rows
## useful when bootstrap sampling
SubsetDat <- function(dat,ix){
  dat$y <- dat$y[ix]
  dat$gs <- dat$gs[ix,,drop=FALSE]
  dat$snp <- dat$snp[ix,,drop=FALSE]
  dat$co <- dat$co[ix,,drop=FALSE]
  return(dat)
}


#########
######### PLOTTING AND DISPLAY FUNCTIONS

## returns subgraph (as eqtl matrix) which has causal influence on y
FindSubgraph <- function(params){
  temp <- rowSums(params$eqtl_model[,params$gs_direct!=0,drop=FALSE]!=0)!=0
  snp_to_keep <- temp | params$snp_direct!=0
  gs_to_keep <- params$gs_direct!=0
  return(params$eqtl_model[snp_to_keep,gs_to_keep,drop=FALSE])
}

## creates GraphNEL object which can then
## be plotted
MakeGraphNELObject <- function(eqtl,snp_direct,gs_direct,gs_path=NULL){
  if(is.null(gs_path)){
    gs_path <- matrix(0,nrow=length(gs_direct),ncol=length(gs_direct))
  }
  l <- c("h",rownames(eqtl),colnames(eqtl),"y")
  edL <- vector("list",length=length(l))
  names(edL) <- l
  ## draw arrows from h to snps
  edL[[1]] <- list(edges=2:(nrow(eqtl)+1),weights=rep(1,nrow(eqtl)))
  ## snp arrows
  for(ii in 2:(nrow(eqtl)+1)){
    ix <- which(eqtl[ii-1,]!=0) + nrow(eqtl)+1
    names(ix) <- NULL
    if(snp_direct[ii-1]!=0){
      ix <- c(ix,sum(dim(eqtl))+2)
    }
    edL[[ii]] <- list(edges=ix,weights=rep(1,length(ix)))
  }
  ## gs arrows
  for(ii in (nrow(eqtl)+2):(nrow(eqtl)+ncol(eqtl)+1)){
    jj <- ii - nrow(eqtl) - 1
    if(gs_direct[jj]!=0){
      edL[[ii]] <- list(edges=nrow(eqtl)+ncol(eqtl)+2,weights=1)
    } else {
      edL[[ii]] <- list(edges=c(),weights=c())
    }
    if(length(which(gs_path[,jj]!=0))!=0){
      edL[[ii]]$edges <- c(edL[[ii]]$edges,which(path[,jj]!=0) + nrow(eqtl) + 1)
      edL[[ii]]$weights <- rep(1,length(edL[[ii]]$edges))
    }
  }
  edL[[length(edL)]] <- list(edges=c(),weights=c())
  gR <- graphNEL(nodes=l, edgeL=edL, edgemode="directed")
  return(gR)
}


## makes edge attributes for graphNEL plotting using eqtl_model (gs, SNP coeffs) and direct effects
MakeedgeAttrs <- function(eqtl_model,snp_direct,gs_direct,path=NULL){
  ew <- as.vector(eqtl_model)
  names(ew) <- unlist(lapply(colnames(eqtl_model),function(x){paste0(rownames(eqtl_model),"~",x)}))
  ew <- round(ew,2)
  ef <- c(snp_direct,gs_direct)
  names(ef) <- paste0(c(names(snp_direct),names(gs_direct)),"~y")
  ef <- round(ef,2)
  ew <- c(ew,ef)
  if(!is.null(path)){
    enet <- as.vector(t(path))
    names(enet) <- unlist(lapply(colnames(eqtl_model),function(x){paste0(colnames(eqtl_model),"~",x)}))
    enet <- round(enet,2)
    enet <- enet[enet!=0]
    ew <- c(ew,enet)
  }
  eAttrs <- list()
  eAttrs$label <- ew
  return(eAttrs)
}

## find effects which are non-zero in both true and estimates
FindInterestingEffs <- function(eff_comb){
  direct_both0 <- (eff_comb[,1]==0) & (eff_comb[,1]==eff_comb[,2])
  indirect_both0 <- (eff_comb[,3] == 0) & (eff_comb[,3]==eff_comb[,4])
  indirect_both0[is.na(indirect_both0)] <- TRUE
  to_remove <- direct_both0 & indirect_both0
  return(eff_comb[!to_remove,])
}
