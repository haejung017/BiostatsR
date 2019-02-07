
## Hazard, cumhazard, survival and density function under PE time-dep covariate
marg_fun_ORtvc<- function(mpars, TDfunc, Y, data){
  
  #Y=Y=data$time
  G=data$mgene
  
  # screening indicator
  Sc.Ind=data$gender
  
  # time at OR
  Yor=data$ov.censortime
  
  lamb<-mpars[1]
  rho<-mpars[2]
  beta.g<-mpars[3]
  beta.or<-mpars[4]
  
  beta.st<-mpars[5]
  #beta.nofst<-mpars[6]
  
  inds0=which(Y<=Yor)
  inds1=which(Y>Yor )
  
     haz=vector(length=length(Y))
     H=vector(length=length(Y))
  
     xbeta0 =beta.g*G  + beta.st*Sc.Ind 
  
     # hazard and cumhaz for non Ovarian remove
    haz[inds0]=lamb*rho*exp(xbeta0[inds0])*Y[inds0]^(rho-1)
    H[inds0]=lamb*exp(xbeta0[inds0])*Y[inds0]^rho
    
    # hazard and cumhaz for non Ovarian remove
    haz[inds1]=lamb*rho*exp(xbeta0[inds1])*Y[inds1]^(rho-1)*exp(beta.or)
    H[inds1]=lamb*exp(xbeta0[inds1])*( Yor[inds1]^rho + exp(beta.or)*(Y[inds1]^rho - Yor[inds1]^rho) )
    
    mf<-list(S=exp(-H), f=haz*exp(-H) , haz=haz, H=H)
 
  return( mf )
}



#### tv log-likelihood function
loglik.Acop.ORtvc=function(pars,data, agemin, TDfunc, cop.fam="Clayton"){
  
  
  data = data[data$currentage>=agemin,]
  data$ov.censortime=data$ov.censortime-agemin
  
  status=data$status; 
  
  alpha<-exp(pars[1])
  lamb<-exp(pars[2])
  rho<-exp(pars[3])
  beta.g<-pars[4]
  beta.or<-pars[5]
  beta.st<-pars[6]

  mpars=c(lamb,rho, beta.g=beta.g, beta.or=beta.or, beta.st=beta.st)

  if(alpha<0 || rho<0 || lamb<0) 
    return(Inf)
  
  
  print(c(alpha, mpars))
  
  mdist=marg_fun_ORtvc(mpars=mpars,TDfunc=TDfunc, Y=data$time, data=data)
  
  haz=mdist$haz
  H=mdist$H
  Sy=mdist$S
  
  
  iGenSy<-invGenerator(x=Sy, alpha=alpha, cop.fam=cop.fam)
  fiGenSy <- aggregate(iGenSy, by=list(data$famID), sum)[,2]
  nObs <- aggregate(status, by=list(data$famID), sum)[,2]   # vector of number of observed events
  
  # the two terms of the log-likelihood
  term1=log(haz) - H  - log( -vDGenerator(x=iGenSy, alpha=alpha, cop.fam=cop.fam, order=1) ) 
  term2=log( (-1)^(nObs)*vDGenerator(x=fiGenSy, alpha=alpha, cop.fam=cop.fam, order=nObs) )
  
  ## corrected term for proband
  
  # for effected proband
  datap<- data[data$proband==1,]
  datap$currentage<- datap$currentage-agemin
  statusp <- datap$affect  # proband disease status at the study entry
  
  Hp<-marg_fun_ORtvc(mpars=mpars, TDfunc=TDfunc, Y=datap$currentage, data=datap)$H
  
  logasc0 <- - Hp                   # log survival until study entry for probands
  logasc1 <- log(1-exp(-Hp) )       # log penetrace at study entry for probands
  
  slogasc = sum(logasc0[statusp==0],na.rm=TRUE) + sum(logasc1[statusp==1],na.rm=TRUE)
  cloglik<- sum(status*term1, na.rm=T) + sum(term2, na.rm=T)   -  slogasc 
  
  return(-cloglik)
  
}

## Fit copula models
fit.cop.ORtvc<- function(pars, data, TDfunc, cop.fam, agemin){

  
  out <- optim(pars, loglik.Acop.ORtvc, data=data, agemin=agemin, TDfunc=TDfunc, cop.fam=cop.fam, #method = "Nelder-Mead")
               hessian=TRUE, control = list(trace = 1, maxit = 2500) )
  
  # Results
  pars.est <- out$par
  pars.est[1:3]<-exp(pars.est[1:3])
  print(out$convergence)
  Var <- try(solve(out$hessian), TRUE)
  #Var <- try(ginv(out$hessian), TRUE)
  se <- sqrt(diag(Var))
  se.exp <- exp(out$par)*se
  SE <- c(se.exp[1:3], se[4:6])

  # AIC
  AIC <- 2*(length(out$par)-(-out$value))
  #AIC<-round(AIC,digits=1)
  # WALD & Pvalue 
  wald <-  pars.est/SE
  #wald<-round(wald,digits=2)
  pval<- 2*(1-pnorm(abs(wald)))
  
  ## estimates of beta20 and beta21
  
  alpha <- pars.est[1]
  tau<-copPar2KT(alpha,cop.fam)
  
  mle<-matrix(NA,ncol=4,nrow=7+2 )
  rownames(mle)<-c("alpha","lambda","rho","betag","betasor","betascr","phi","tau","AIC")
  colnames(mle)<-c("Est","SE","Stat","Pvalue")
  mle[,1]=c(pars.est,tau,AIC)
  mle[,2]=c(SE,NA,NA)
  mle[,3]=c(wald,NA,NA)
  mle[,4]=c(pval,NA,NA)
  
  return( mle )
}
