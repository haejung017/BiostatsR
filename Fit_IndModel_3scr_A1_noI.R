
# ficher's z-transformation
alpha_2_eta<-function(alpha) log( (1+alpha)/(1-alpha))/2
eta_2_alpha<-function(eta) (exp(2*eta)-1)/(exp(2*eta)+1)

# alpha_2_eta(0)
# eta_2_alpha(0)

#### fit Copula model
fit.Ind<-function(pars, dataset, TDf_or, TDf_sc,  agemin, phi0=1e-5,  phi=c(3,2,1e-5), wgt=FALSE){
  
  
  ee0=phi0
  ee1=phi[1]
  ee2=phi[2]
  ee3=phi[3]
  
  phi_init0=log(ee0)
  phi_init1=log(ee1)
  phi_init2=log(ee2)
  phi_init3=log(ee3)
  phi_sc=c(phi_init1, phi_init2, phi_init3)
  
  a00=0.001
  a01=0.001
  a02=0.001
  a03=0.001
  a0=c(0.001, 0.001, 0.001)
  
  pars1=pars
  if(TDf_or=="PE" && TDf_sc=="ED")
    pars1=c(pars, phi_sc)
  if(TDf_or=="PE" && TDf_sc=="CO")
    pars1=c(pars, phi_init1, a01, phi_init2, a02, phi_init3, a03 )
  if(TDf_or=="ED" && TDf_sc=="PE")
    pars1=c(pars, phi_init0)
  if(TDf_or=="ED" && TDf_sc=="ED")
    pars1=c(pars,phi_init0, phi_sc)
  if(TDf_or=="ED" && TDf_sc=="CO")
    pars1=c(pars,phi_init0, phi_init1, a01, phi_init2, a02, phi_init3, a03 )
  
  if(TDf_or=="CO" && TDf_sc=="PE")
    pars1=c(pars, phi_init0, a00)
  if(TDf_or=="CO" && TDf_sc=="ED")
    pars1=c(pars,phi_init0, a00, phi_sc)
  if(TDf_or=="CO" && TDf_sc=="CO")
    pars1=c(pars,phi_init0, a00, phi_init1, a01, phi_init2, a02, phi_init3, a03 )
  
  
  
  
  out <- try( optim(pars1, loglik.Ind, data=dataset, agemin=agemin, TDf_or=TDf_or, TDf_sc=TDf_sc, wgt=wgt,
                    hessian=TRUE, control = list(trace = 1, maxit = 10000) ) , TRUE)
  
  if(inherits(out ,'try-error')){
    warning(as.matrix(out)) 
    out <- optim(pars1, loglik.Ind, data=dataset, agemin=agemin, TDf_or=TDf_or, TDf_sc=TDf_sc, wgt=wgt, method="Nelder-Mead",
                 control = list(trace = 1, maxit = 10000) )
    out$hessian <- hessian(loglik.Ind, out$par, data=dataset, agemin=agemin, TDf_or=TDf_or, TDf_sc=TDf_sc, wgt=wgt)
  }
  
  print(out$convergence)
  # print(out$hessian)
  #out$hessian[out$hessian==0]=1e-7
  # Var1 <- try(solve(out$hessian), TRUE)
  Var1 <- try(ginv(out$hessian), TRUE)
  if(inherits(Var1 ,'try-error')){
    warning(as.matrix(Var1))
    Var1 <- try(solve(out$hessian), TRUE)
    if(inherits(Var1 ,'try-error')){
      warning(as.matrix(Var1))
      Var1 <-NA_real_
    }
  }
  
  #print(Var1)
  print(diag(Var1))
  
  Var<-Var1
  # Results
  pars.est <- out$par
  pars.est[1:2]<-exp(pars.est[1:2])
  dg<-abs(diag(Var))
  se <- sqrt(dg)
  se.exp <- exp(out$par)*se
  # se.cg <- eta_2_alpha(out$par)*se
  SE <- c(se.exp[1:2], se[3:7])
  
  
  pars.est1<-c(pars.est[1:7], phi0=NA, a00=NA, phi1=NA, a01=NA, phi2=NA, a02=NA, phi3=NA, a03=NA, LRTphi=NA, alpha=NA, tau=NA, df=NA, LRT=NA,   AIC=NA, BIC=NA, LogL=NA)
  SE1<-c(SE[1:7], phi0=NA, a00=NA, phi1=NA, a01=NA, phi2=NA, a02=NA, phi3=NA, a03=NA,  LRTphi=NA, alpha=NA, tau=NA, df=NA, LRT=NA,  AIC=NA,  BIC=NA, LogL=NA)
  
  LogL<- -out$value
  AIC <- 2*(length(out$par)-LogL)
  dd=sum(dataset$status, na.rm=T)
  
  BIC=log(dd)*length(out$par) -2*LogL
  
  # WALD & Pvalue 
  wald <-  pars.est1/SE1
  pval<-ifelse(wald>0, 2*(1-pnorm(wald)), 2*(pnorm(wald)) )
  
  if(TDf_or=="PE" && TDf_sc=="ED"){             # two terms of the likelihood under frailty model
    k=7
    
    pars.est1[k+3]<-exp(pars.est[k+1])
    SE1[k+3]<-se.exp[k+1]
    wald[k+3]=pars.est1[k+3]/SE1[k+3]
    pval[k+3]<- 1-pnorm(wald[k+3])
    
    pars.est1[k+5]<-exp(pars.est[k+2])
    SE1[k+5]<-se.exp[k+2]
    wald[k+5]=pars.est1[k+5]/SE1[k+5]
    pval[k+5]<- 1-pnorm(wald[k+5])
    
    pars.est1[k+7]<-exp(pars.est[k+3])
    SE1[k+7]<-se.exp[k+3]
    wald[k+7]=pars.est1[k+7]/SE1[k+7]
    pval[k+7]<- 1-pnorm(wald[k+7])
  }
  
  if(TDf_or=="PE" && TDf_sc=="CO"){             # two terms of the likelihood under frailty model
    k=7
    pars.est1[k+3]<-exp(pars.est[k+1])
    pars.est1[k+4]<-pars.est[k+2]
    pars.est1[k+5]<-exp(pars.est[k+3])
    pars.est1[k+6]<-pars.est[k+4]
    pars.est1[k+7]<-exp(pars.est[k+5])
    pars.est1[k+8]<-pars.est[k+6]
    
    SE1[(k+3)]<-se.exp[(k+1)]
    SE1[(k+4)]<-se[(k+2)]
    
    SE1[(k+5)]<-se.exp[(k+3)]
    SE1[(k+6)]<-se[(k+4)]
    
    SE1[(k+7)]<-se.exp[(k+5)]
    SE1[(k+8)]<-se[(k+6)]
    
    
    wald[c(10,12,14)] <-   pars.est1[c(10,12,14)]/ SE1[c(10,12,14)]
    pval[c(10,12,14)]<- 1-pnorm(wald[c(10,12,14)])
    
    wald[c(11,13,15)] <-   pars.est1[c(11,13,15)]/ SE1[c(11,13,15)]
    pval[c(11,13,15)]<- ifelse(wald[c(11,13,15)]>0, 2*(1-pnorm(wald[c(11,13,15)])), 2*(pnorm(wald[c(11,13,15)])) )
  }
  
  if(TDf_or=="ED" && TDf_sc=="PE"){             # two terms of the likelihood under frailty model
    pars.est1[8]<-exp(pars.est[length(pars.est)])
    SE1[8]<-se.exp[length(pars.est)]
    wald[8]=pars.est1[8]/SE1[8]
    pval[8]<- 1-pnorm(wald[8])
  }
  
  if(TDf_or=="ED" && TDf_sc=="ED"){             # two terms of the likelihood under frailty model
    k=7
    pars.est1[k+1]<-exp(pars.est[k+1])
    SE1[k+1]<-se.exp[k+1]
    wald[k+1]=pars.est1[k+1]/SE1[k+1]
    pval[k+1]<- 1-pnorm(wald[k+1])
    
    pars.est1[k+3]<-exp(pars.est[k+2])
    SE1[k+3]<-se.exp[k+2]
    wald[k+3]=pars.est1[k+3]/SE1[k+3]
    pval[k+3]<- 1-pnorm(wald[k+3])
    
    pars.est1[k+5]<-exp(pars.est[k+3])
    SE1[k+5]<-se.exp[k+3]
    wald[k+5]=pars.est1[k+5]/SE1[k+5]
    pval[k+5]<- 1-pnorm(wald[k+5])
    
    pars.est1[k+7]<-exp(pars.est[k+4])
    SE1[k+7]<-se.exp[k+4]
    wald[k+7]=pars.est1[k+7]/SE1[k+7]
    pval[k+7]<- 1-pnorm(wald[k+7])
    
  }
  
  if(TDf_or=="ED" && TDf_sc=="CO"){             # two terms of the likelihood under frailty model
    k=7
    pars.est1[k+1]<-exp(pars.est[k+1])
    SE1[k+1]<-se.exp[k+1]
    wald[k+1]=pars.est1[k+1]/SE1[k+1]
    pval[k+1]<- 1-pnorm(wald[k+1])
    
    pars.est1[k+3]<-exp(pars.est[k+2])
    SE1[k+3]<-se.exp[k+2]
    wald[k+3]=pars.est1[k+3]/SE1[k+3]
    pval[k+3]<- 1-pnorm(wald[k+3])
    
    pars.est1[k+4]<-pars.est[k+3]
    SE1[k+4]<-se[k+3]
    wald[k+4]=pars.est1[k+4]/SE1[k+4]
    pval[k+4]<- ifelse(wald[k+4]>0, 2*(1-pnorm(wald[k+4])), 2*(pnorm(wald[k+4])) )
    
    pars.est1[k+5]<-exp(pars.est[k+4])
    SE1[k+5]<-se.exp[k+4]
    wald[k+5]=pars.est1[k+5]/SE1[k+5]
    pval[k+5]<- 1-pnorm(wald[k+5])
    
    pars.est1[k+6]<-pars.est[k+5]
    SE1[k+6]<-se[k+5]
    wald[k+6]=pars.est1[k+6]/SE1[k+6]
    pval[k+6]<-ifelse(wald[k+6]>0, 2*(1-pnorm(wald[k+6])), 2*(pnorm(wald[k+6])) )
    
    pars.est1[k+7]<-exp(pars.est[k+6])
    SE1[k+7]<-se.exp[k+6]
    wald[k+7]=pars.est1[k+7]/SE1[k+7]
    pval[k+7]<- 1-pnorm(wald[k+7]) 
    
    pars.est1[k+8]<-pars.est[k+7]
    SE1[k+8]<-se[k+7]
    wald[k+8]=pars.est1[k+8]/SE1[k+8]
    pval[k+8]<-ifelse(wald[k+8]>0, 2*(1-pnorm(wald[k+8])), 2*(pnorm(wald[k+8])) )
    
  }
  if(TDf_or=="CO" && TDf_sc=="PE"){             # two terms of the likelihood under frailty model
    pars.est1[8]<-exp(pars.est[length(pars.est)-1])
    SE1[8]<-se.exp[length(pars.est)-1]
    wald[8]=pars.est1[8]/SE1[8]
    pval[8]<- 1-pnorm(wald[8])
    
    pars.est1[9]<-pars.est[length(pars.est)]
    SE1[9]<-se[length(pars.est)]
    wald[9]=pars.est1[9]/SE1[9]
    pval[9]<- ifelse(wald[9]>0, 2*(1-pnorm(wald[9])), 2*(pnorm(wald[9])) )
  } 
  if(TDf_or=="CO" && TDf_sc=="ED"){             # two terms of the likelihood under frailty model
    k=7
    pars.est1[8]<-exp(pars.est[k+1])
    SE1[8]<-se.exp[k+1]
    wald[8]=pars.est1[8]/SE1[8]
    pval[8]<- 1-pnorm(wald[8])
    
    pars.est1[9]<-pars.est[k+2]
    SE1[9]<-se[k+2]
    wald[9]=pars.est1[9]/SE1[9]
    pval[9]<- ifelse(wald[9]>0, 2*(1-pnorm(wald[9])), 2*(pnorm(wald[9])) )
    
    pars.est1[k+3]<-exp(pars.est[k+3])
    SE1[k+3]<-se.exp[k+3]
    wald[k+3]=pars.est1[k+3]/SE1[k+3]
    pval[k+3]<- 1-pnorm(wald[k+3])
    
    pars.est1[k+5]<-exp(pars.est[k+4])
    SE1[k+5]<-se.exp[k+4]
    wald[k+5]=pars.est1[k+5]/SE1[k+5]
    pval[k+5]<- 1-pnorm(wald[k+5])
    
    pars.est1[k+7]<-exp(pars.est[k+5])
    SE1[k+7]<-se.exp[k+5]
    wald[k+7]=pars.est1[k+7]/SE1[k+7]
    pval[k+7]<- 1-pnorm(wald[k+7])
  }
  
  if(TDf_or=="CO" && TDf_sc=="CO"){             # two terms of the likelihood under frailty model
    k=7
    pars.est1[8]<-exp(pars.est[k+1])
    SE1[8]<-se.exp[k+1]
    wald[8]=pars.est1[8]/SE1[8]
    pval[8]<- 1-pnorm(wald[8])
    
    pars.est1[9]<-pars.est[k+2]
    SE1[9]<-se[k+2]
    wald[9]=pars.est1[9]/SE1[9]
    pval[9]<- ifelse(wald[9]>0, 2*(1-pnorm(wald[9])), 2*(pnorm(wald[9])) )
    
    pars.est1[k+3]<-exp(pars.est[k+3])
    SE1[k+3]<-se.exp[k+3]
    wald[k+3]=pars.est1[k+3]/SE1[k+3]
    pval[k+3]<- 1-pnorm(wald[k+3])
    
    pars.est1[k+4]<-exp(pars.est[k+4])
    SE1[k+4]<-se[k+4]
    wald[k+4]=pars.est1[k+4]/SE1[k+4]
    pval[k+4]<- ifelse(wald[k+4]>0, 2*(1-pnorm(wald[k+4])), 2*(pnorm(wald[k+4])) )
    
    pars.est1[k+5]<-exp(pars.est[k+5])
    SE1[k+5]<-se.exp[k+5]
    wald[k+5]=pars.est1[k+5]/SE1[k+5]
    pval[k+5]<- 1-pnorm(wald[k+5])
    
    pars.est1[k+6]<-pars.est[k+6]
    SE1[k+6]<-se[k+6]
    wald[k+6]=pars.est1[k+6]/SE1[k+6]
    pval[k+6]<-ifelse(wald[k+6]>0, 2*(1-pnorm(wald[k+6])), 2*(pnorm(wald[k+6])) )
    
    pars.est1[k+7]<-exp(pars.est[k+7])
    SE1[k+7]<-se[k+7]
    wald[k+7]=pars.est1[k+7]/SE1[k+7]
    pval[k+7]<-2*(1-pnorm(abs(wald[k+7])))
    
    pars.est1[k+8]<-pars.est[k+8]
    SE1[k+8]<-se[k+8]
    wald[k+8]=pars.est1[k+8]/SE1[k+8]
    pval[k+8]<-ifelse(wald[k+8]>0, 2*(1-pnorm(wald[k+8])), 2*(pnorm(wald[k+8])) )
  }
  
  pars.est1[length(pars.est1)-2]<-AIC
  pars.est1[length(pars.est1)-1]<-BIC
  pars.est1[length(pars.est1)]<-LogL
  
  mle<-matrix(NA,ncol=3,nrow=length(pars.est1) )
  colnames(mle)<-c("Est","SE","Pvalue")
  mle[,1]=pars.est1
  mle[,2]=SE1
  mle[,3]=pval
  
  rownames(mle)=c("Lambda","Rho","gene","Ooph", "scr1","scr2","scr3",
                  "phi0", "a00", "phi1", "a01", "phi2", "a02", "phi3", "a03","LRTphi", "alpha","tau","df","LRT", "AIC","BIC","LogL")
  
  return( mle )
}


#### TD Log-likelihood function under copula model  
loglik.Ind=function(pars, data, agemin, TDf_or="PE", TDf_sc="PE", wgt=FALSE){
  
  # the data
  data = data[data$currentage>=agemin,]
  Y=data$time - agemin
  #G=data$mgene
  Ybr=data$br.censortime-agemin
  Yor=data$ov.censortime-agemin
  Yst1=data$st1-agemin
  Yst2=data$st2-agemin
  Yst3=data$st3-agemin
  
  cp = data$carrp.geno                                         #carrier probabilty
  
  status=data$status 
  wt<-rep(1,length(Y))
  if(wgt==TRUE)
    wt<-1/data$weight
  
  
  # BRI <- ifelse(Y-Ybr> 0,1,0)
  # Y[BRI==1]<-Ybr[BRI==1]
  # status[BRI==1]=0
  
  #parameters
  lamb<-exp(pars[1])
  rho<-exp(pars[2])
  beta.g<-pars[3]
  beta.or<-pars[4]
  #beta.gor<-pars[5]
  beta.sc<-pars[5:7]
  #beta.gsc<-c(0,0,0)
  # beta.orsc<-pars[12:14]
  
  
  if(TDf_or=="PE" && TDf_sc=="PE" ){
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc)
  } else  if(TDf_or=="PE" && TDf_sc=="ED" ){
    phi1=exp(pars[length(pars)-2])
    phi2=exp(pars[length(pars)-1])
    phi3=exp(pars[length(pars)])
    
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc,  phi1, phi2, phi3)
    
    
  } else  if(TDf_or=="PE" & TDf_sc=="CO" ){
    
    phi1=exp(pars[length(pars)-5])
    a01=pars[length(pars)-4]
    phi2=exp(pars[length(pars)-3])
    a02=pars[length(pars)-2]
    phi3=exp(pars[length(pars)-1])
    a03=pars[length(pars)]
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc,  phi1, a01, phi2, a02, phi3, a03)
    
  } else  if(TDf_or=="ED" && TDf_sc=="PE" ){
    phi0=exp(pars[length(pars)])
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc,  phi0)
   # if(phi0<e1 || phi0>e2) return( Inf )
  } else  if(TDf_or=="ED" && TDf_sc=="ED" ){
    phi0=exp(pars[length(pars)-3])
    phi1=exp(pars[length(pars)-2])
    phi2=exp(pars[length(pars)-1])
    phi3=exp(pars[length(pars)])
    
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc,  phi0, phi1, phi2, phi3)
    
  } else  if(TDf_or=="ED" && TDf_sc=="CO" ){
    phi0=exp(pars[length(pars)-6])
    phi1=exp(pars[length(pars)-5])
    a01=pars[length(pars)-4]
    phi2=exp(pars[length(pars)-3])
    a02=pars[length(pars)-2]
    phi3=exp(pars[length(pars)-1])
    a03=pars[length(pars)]
    
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc,  phi0, phi1, a01, phi2, a02, phi3, a03)
    
  } else  if(TDf_or=="CO" && TDf_sc=="PE" ){
    phi0=exp(pars[length(pars)-1])
    a00=pars[length(pars)]
    
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc,  phi0, a00)
    
  } else  if(TDf_or=="CO" && TDf_sc=="ED" ){
    phi0=exp(pars[length(pars)-4])
    a00=pars[length(pars)-3]
    phi1=exp(pars[length(pars)-2])
    phi2=exp(pars[length(pars)-1])
    phi3=exp(pars[length(pars)])
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc,  phi0, a00, phi1, phi2, phi3)

  } else  if(TDf_or=="CO" && TDf_sc=="CO" ){
    
    phi0=exp(pars[length(pars)-7])
    a00=pars[length(pars)-6]
    phi1=exp(pars[length(pars)-5])
    a01=pars[length(pars)-4]
    phi2=exp(pars[length(pars)-3])
    a02=pars[length(pars)-2]
    phi3=exp(pars[length(pars)-1])
    a03=pars[length(pars)]
    mpars=c(lamb, rho,beta.g, beta.or, beta.sc,  phi0, a00, phi1, a01, phi2, a02, phi3, a03)
  }
  
  
  #marginal distributions 
  mdist=marg_fun(mpars=mpars, Y=Y, Yor=Yor, Yst1=Yst1, Yst2=Yst2, Yst3=Yst3,  Ybr=Ybr, TDf_or=TDf_or, TDf_sc=TDf_sc, cp=cp)
  loghaz=mdist$logh
  H=mdist$H
  
  print( mpars )
  term1= wt*loghaz 
  term2=-wt*H
  # the naive log-likelihood function
  loglik<- sum(status*term1, na.rm=T) + sum(term2, na.rm=T) 
  
  
  ## corrected term bases on affected and unaffected  proband
  ip<- which(data$proband==1)
  cage<- data$currentage[ip]-agemin
  timep<- data$time[ip]-agemin
  Gp=data$mgene[ip]
  Yst1p=data$st1p[ip]-agemin
  Yst2p=data$st2p[ip]-agemin
  Yst3p=data$st3p[ip]-agemin
  Ybrp=data$br.censortimep[ip]-agemin
  Yorp=data$ov.censortimep[ip]-agemin
  nbSCp=data$screen_number_p
  statusp <- data$affect[ip]          # proband disease status at the study entry
  wtp=wt[ip]
  cpp=cp[ip]
  # BRIp=ifelse(cage<=Ybrp,0,1)
  # cage[BRIp==1]<-Ybrp[BRIp==1]
  # statusp[BRIp==1]=0
  
  ### For unaffected proband
  ip0<-which(statusp==0)
  cage0<- cage[ip0]
  Gp0=Gp[ip0]
  Yst1p0=Yst1p[ip0]
  Yst2p0=Yst2p[ip0]
  Yst3p0=Yst3p[ip0]
  Ybrp0=Ybrp[ip0]
  Yorp0=Yorp[ip0]
  wtp0=wtp[ip0]
  cpp0=cpp[ip0]
  Hp0<-marg_fun(mpars=mpars, Y=cage0, Yor=Yorp0, Yst1=Yst1p0, Yst2=Yst2p0, Yst3=Yst3p0, Ybr=Ybrp0, TDf_or=TDf_or, TDf_sc=TDf_sc, cp=cpp0)$H
  #Hp0[Hp0<=0]=1e-10
  logasc0 <- - Hp0                   # log survival until study entry for unaffacted probands
  logasc0 <- sum(wtp0*logasc0,na.rm=TRUE)
  
  ### For affected proband
  ip1<-which(statusp==1 | statusp==2)
  cage1<- cage[ip1]
  Gp1=Gp[ip1]
  #Yscp1=Yscp[statusp==1 |statusp==2]
  Ybrp1=Ybrp[ip1]
  Yorp1=Yorp[ip1]
  Yst1p1=Yst1p[ip1]
  Yst2p1=Yst2p[ip1]
  Yst3p1=Yst3p[ip1]
  
  
  wtp1=wtp[ip1]
  cpp1=cpp[ip1]
  #Hp1<-marg_funp1(mpars=mpars, Y=cage1, Yor=Yorp1, TDfunc=TDfunc, cp=cpp1)$H
  Hp1<-marg_fun(mpars=mpars, Y=cage1, Yor=Yorp1, Yst1=Yst1p1, Yst2=Yst2p1, Yst3=Yst3p1, Ybr=Ybrp1, TDf_or=TDf_or, TDf_sc=TDf_sc, cp=cpp1)$H
  #print(Hp1)
  #Hp1[Hp1<=0]=1e-10
  logasc1 <- log(1-exp(-Hp1) )                   # log penetrace at study entry for affected probands
  logasc1 <- sum(wtp1*logasc1,na.rm=TRUE)
  
  slogasc = logasc0 + logasc1
  
  # corrected log-likelihood
  cloglik<- loglik  -  slogasc
  
  # cloglik<- loglik
  return(-cloglik)
}

# # test between ED and PE
#  mpars=c(log(I.PE.PE[1:2,1]),I.PE.PE[3:10,1])
# loglik.Ind(pars=mpars, data, agemin, TDf_or="PE", TDf_sc="PE", wgt=FALSE)
# loglik.Ind(pars=c(mpars, log(1e-7)), data, agemin, TDf_or="ED", TDf_sc="PE", wgt=FALSE)
# loglik.Ind(pars=c(mpars,log(1e-7), log(1e-7), log(1e-7)), data, agemin, TDf_or="PE", TDf_sc="ED", wgt=FALSE)
# loglik.Ind(pars=c(mpars,log(1e-7), 0, log(1e-7), 0, log(1e-7), 0), data, agemin, TDf_or="PE", TDf_sc="CO", wgt=FALSE)
# loglik.Ind(pars=c(mpars, log(1e-7), log(1e-7), log(1e-7), log(1e-7) ), data, agemin, TDf_or="ED", TDf_sc="ED", wgt=FALSE)

# loglik.Ind(pars=c(mpars, log(1e-7)), data, agemin, TDf_or="ED", TDf_sc="PE", wgt=FALSE)
# loglik.Ind(pars=c(mpars,log(1e-7), log(1e-7), log(1e-7), log(1e-7)), data, agemin, TDf_or="ED", TDf_sc="ED", wgt=FALSE)
# loglik.Ind(pars=c(mpars,log(1e-7),log(1e-7), 0, log(1e-7), 0, log(1e-7), 0), data, agemin, TDf_or="ED", TDf_sc="CO", wgt=FALSE)
# loglik.Ind(pars=c(mpars, log(1e-7),0), data, agemin, TDf_or="CO", TDf_sc="PE", wgt=FALSE)




