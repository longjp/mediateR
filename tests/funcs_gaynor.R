####### Modification of code from github repo Gaynor "Mediation analysis for common binary outcomes"

## simulate following section 3.1
SimulatePaper <- function(n,k){
  mu_a <- 0.4
  mu_c <- mu_a*0.3
  sig <- 0.75
  beta <- matrix(c(0.1,0.5,0.4),ncol=1)
  theta <- matrix(c(k,0.4,0.5,0.25),ncol=1)
  ##prev <- inv.logit(trueT0 + trueT1*trueA + trueT2*(trueB0 + trueB1*trueA + trueB2*(trueC)) + trueT4*(trueC)) 
  prev <- inv.logit(theta[1] + theta[2]*mu_a + theta[3]*(beta[1] + beta[2]*mu_a + beta[3]*(mu_c)) + theta[4]*(mu_c))
  a <- matrix(rnorm(n,mean=mu_a,sd=sig),ncol=1)
  con <- matrix(rnorm(n,mean=mu_c,sd=sig),ncol=1)
  m <- cbind(1,a,con)%*%beta + rnorm(n,mean=0,sd=sig)
  p <- 1 / (1 + exp(-cbind(1,a,m,con)%*%theta))
  y <- rbinom(n,size=1,prob=p)
  return(list(y=y,con=as.vector(con),a=as.vector(a),m=as.vector(m),prev=prev))
}

## from a,m,con,y to snp,gs,co,y data format
ReformatData <- function(dat){
  snp <- matrix(dat$a,ncol=1)
  colnames(snp) <- "SNP1"
  gs <- matrix(dat$m,ncol=1)
  colnames(gs) <- "gs1"
  co <- matrix(dat$con,ncol=1)
  colnames(co) <- "co1"
  eqtl <- matrix(1,nrow=1,ncol=1)
  colnames(eqtl) <- colnames(gs)
  rownames(eqtl) <- colnames(snp)
  return(list(y=dat$y,snp=snp,gs=gs,co=co,eqtl=eqtl,family="binomial"))
}


directEffectProposed   <-function(s,t0,t1,t2,t4,c,sigma,b0,b1,b2){ 
  return( oddCalc(pnorm(    (s*t0 + s*t1 + s*t4*c + s*t2*(b0 + b2*c)) / (sqrt(1 + (s*t2*sigma)^2))  )) /
            oddCalc(pnorm(    (s*t0 +        s*t4*c + s*t2*(b0 + b2*c)) / (sqrt(1 + (s*t2*sigma)^2))  ))) }
indirectEffectProposed <-function(s,t0,t1,t2,t4,c,sigma,b0,b1,b2){ 
  return( oddCalc(pnorm(    (s*t0 + s*t1 + s*t4*c + s*t2*(b0 + b1 + b2*c)) / (sqrt(1 + (s*t2*sigma)^2))  )) / 
            oddCalc(pnorm(    (s*t0 + s*t1 + s*t4*c + s*t2*(b0 +      b2*c)) / (sqrt(1 + (s*t2*sigma)^2))  ))) }
oddCalc<-function(p){ return((p/(1-p))) }
mse <- function(sm){ mean(sm$residuals^2) }
erfinv <-  function(x) { qnorm((1 + x) /2) / sqrt(2) }


EstimateDirect <- function(dat){
  y <- dat$y
  a <- dat$a
  m <- dat$m
  con <- dat$con
  cCon<- median(con)
  logitModel <- glm(y ~ a + m  + con, na.action=na.omit, family=binomial(link="logit"))
  probitModel <- glm(y ~ a + m  + con, na.action=na.omit, family=binomial(link="probit"))
  sEst <- median(coef(probitModel)/coef(logitModel), na.rm = TRUE)
  linearModel <- lm(m ~ a + con, na.action=na.omit)
  logitModelCoef <- coef(logitModel); linearModelCoef <- coef(linearModel)
  de <- directEffectProposed(s= sEst,t0=logitModelCoef[1],
                             t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                             c=cCon,sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],
                             b2=linearModelCoef[3])
  return(de)
}


