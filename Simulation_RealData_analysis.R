b1 <- read.csv("C:/Users/jay/Dropbox/Jay/New Era/CompetingRisk Logliklihood/Data/dfcomp_ScreenUpto5_March9.csv")
b1=carrierprobgeno(method="data",b1)
table(b1$mgene,useNA  = "always")
impute=function(b1){
  set.seed(99)
for(i in 1:nrow(b1)){
  if(is.na(b1$mgene[i]))b1[i,"mgene"]<-sample(c(0,1),1,prob=c(1-b1$carrp.geno[i],b1$carrp.geno[i]),replace=TRUE)
}
  b1
}  
b1=impute(b1)
b1[b1$status==4,"status"]<-0
b1[b1$status==3,"status"]<-2
table(b1$status)
table(b1$mgene,useNA  = "always")
b1$TVCstatus=b1$gender
b1$tvc=b1$st1

Comp_simulation3_RD_one <- function(data = data, initparm=c(0.009,2.4424,1,1.5),missing.method = "data", base.dist = "Weibull", tvctype="PE"){
  
  #initial parameters
  base.parms <- initparm[1:4]
  vbeta <- initparm[5:7]
  theta = theta0 = c(log(base.parms),vbeta) 
  
  
  datacopy = data
  fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
  df1 <- aggregate(datacopy$status==1, by=list(datacopy$famID), FUN=sum)[,2]
  df2 <- aggregate(datacopy$status==2, by=list(datacopy$famID), FUN=sum)[,2]
  data$df1 <- rep(df1, fsize)
  data$df2 <- rep(df2, fsize)
  fp <- c(initparm[8:9])
  
  theta = theta0 = c(theta, log(fp))
  
  
  if(tvctype=="ED")theta <- c(theta,log(initparm[10]))
  if(tvctype=="CO")theta <- c(theta,log(initparm[10]),initparm[11])
  # data preparation
  # EM
  data.cooked = carrierprobgeno(method = missing.method, data=data)  
  
  est0 <- est <- theta
  
  est1 =        optim(est, loglik_Comp_simulation3_RD_one, data=data.cooked,  
                      agemin=16, frailty=TRUE,tvctype=tvctype,
                      hessian=FALSE, control=list(maxit=100000))
  # est2 =        optim(log(fp), loglik_Comp_simulation3_RD_second, data=data.cooked,  
  #                     agemin=16, frailty=TRUE,tvctype=tvctype,stage1=est1$par,
  #                     hessian=FALSE, control=list(maxit=100000))
  
  
  if(tvctype=="TI")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:9]))
  if(tvctype=="PE")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:9]))
  if(tvctype=="ED")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:10]))
  if(tvctype=="CO")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:10]),est1$par[11])
  
  
  return(mle)
}
Comp_simulation3_RD <- function(data = data, initparm=c(0.009,2.4424,1,1.5),missing.method = "data", base.dist = "Weibull", tvctype="PE"){
  
  #initial parameters
  base.parms <- initparm[1:4]
  vbeta <- initparm[5:7]
  theta = theta0 = c(log(base.parms),vbeta) 
  
  
    datacopy = data
    fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
    df1 <- aggregate(datacopy$status==1, by=list(datacopy$famID), FUN=sum)[,2]
    df2 <- aggregate(datacopy$status==2, by=list(datacopy$famID), FUN=sum)[,2]
    data$df1 <- rep(df1, fsize)
    data$df2 <- rep(df2, fsize)
    fp <- c(initparm[8:9])
    
  
  
  if(tvctype=="ED")theta <- c(theta,log(initparm[10]))
  if(tvctype=="CO")theta <- c(theta,log(initparm[10]),initparm[11])
  # data preparation
  # EM
  data.cooked = carrierprobgeno(method = missing.method, data=data)  
  
  #  theta = theta0 = c(theta, log(fp))
  est0 <- est <- theta
  
  est1 =        optim(est, loglik_Comp_simulation3_RD, data=data.cooked,  
                      agemin=16, frailty=FALSE,tvctype=tvctype,
                      hessian=FALSE, control=list(maxit=100000))
  est2 =        optim(log(fp), loglik_Comp_simulation3_RD_second, data=data.cooked,  
                      agemin=16, frailty=TRUE,tvctype=tvctype,stage1=est1$par,
                      hessian=FALSE, control=list(maxit=100000))
  
 
    if(tvctype=="TI")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:9]))
    if(tvctype=="PE")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est2$par[1:2]))
    if(tvctype=="ED")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est2$par[1:2]),exp(est1$par[8]))
    if(tvctype=="CO")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est2$par[1:2]),exp(est1$par[8]),est1$par[9])
  
  
  return(mle)
}
loglik_Comp_simulation3_RD_one=function(theta, data, agemin, frailty,tvctype){
  
  #data = data[data$currentage>=agemin,]
  # base parameters
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  
  
  # vbeta for breast  
  beta.sc1_1 = theta[5]
  beta.gen1 = theta[6]
  
  # vbeta for ovarian
  #beta.sc1_2 =  theta[7]
  beta.gen2 = theta[7]
  
  if(tvctype=="ED"){
    eta <- exp(theta[10])
  }else if(tvctype=="CO"){
    eta <- c(exp(theta[10]),theta[11])
  }
  
  # frailty parameter
  if(frailty==TRUE){
    fp <- exp(theta[8:9])
  }else{
    fp=NULL
  }
  
  
  # Y, delta, carrier prob, proband indicator, sc number
  time0 = data$time-agemin
  status = data$status
  ip = which(data$proband==1)
  #ex1 = data$carrp.geno 
  #gender = data$gender
  mgene=data$mgene
  st1=data$tvc-agemin
  gender=data$TVCstatus
  st1[gender==0]<-time0[gender==0]
  
  # baseline hazard(weibull), baseline cumulative hazard(weibull)
  bhaz1 = (lambda1^rho1)*rho1*time0^(rho1-1)
  bhaz2 = (lambda2^rho2)*rho2*time0^(rho2-1)
  
  logh1 = log(bhaz1) + mgene*beta.gen1
  logh2 = log(bhaz2) + mgene*beta.gen2
  
  
  # creating screen,gene,oophorectomy interaction vectors for hazard
  g=function(t,tsc,eta){exp(-eta*(t-tsc))}
  
  if(tvctype=="TI"){scvec1=gender*beta.sc1_1}
  if(tvctype=="PE"){scvec1=gender*beta.sc1_1}
  if(tvctype=="ED"){scvec1=gender*(g(time0,st1,eta[1])*beta.sc1_1)}
  if(tvctype=="CO"){scvec1=gender*(g(time0,st1,eta[1])*beta.sc1_1+eta[2])}
  
  # sum of log-hazard
  sum1 = sum(( logh1 + scvec1)[status==1], na.rm=TRUE) 
  sum2 = sum(  logh2[status==2], na.rm=TRUE)              
  
  if(tvctype=="PE"){cH1=cH(st1[gender==1],time0[gender==1],c(lambda1,rho1,beta.sc1_1),"PE")}
  if(tvctype=="ED"){cH1=cH(st1[gender==1],time0[gender==1],c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}
  if(tvctype=="CO"){cH1=cH(st1[gender==1],time0[gender==1],c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}
  
  # sum of survival until event time for all the individuals in the data
  CH_bc=vector()
  if(tvctype!="TI"){
    CH_bc[gender==0] <- exp(beta.gen1*mgene[gender==0])*((lambda1*time0[gender==0])^rho1)
    CH_bc[gender==1] <- exp(beta.gen1*mgene[gender==1])*((lambda1*st1[gender==1])^rho1+cH1)
  }else{
    CH_bc <- exp(beta.gen1*mgene)*((lambda1*time0)^rho1)*exp(scvec1)  
  }
  
  CH_ov <- exp(beta.gen2*mgene)*((lambda2*time0)^rho2)
  
  
  k <- CH_bc+CH_ov
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam1 <- (CH_bc)
    Hfam2 <- (CH_ov)
    
    df1 <- data$df1[data$proband==1]
    df2 <- data$df2[data$proband==1]
    
    Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
    Hfam2 <- aggregate(Hfam2,by=list(data$famID),FUN=sum)[,2]
    
    
    sum4.1 <- sum((lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1])), na.rm=T)
    sum4.2 <- sum((lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2])), na.rm=T)
    
    loglik = sum1+sum2 + sum4.1+ sum4.2
  }else{
    sum4 <- sum(-k,na.rm = TRUE)
    loglik = sum1+sum2+sum4 # numerator in loglikelihood
  }
  
  # Ascertainment correction by design="pop+"
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.gen2)
  cageptotal <- data$currentage[ip]-agemin
  statusp <- data$affect[ip]  # proband disease status at the study entry
  genderptotal = data$screenp[ip]
  st1ptotal=data$st1p[ip]-agemin
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  
  if(tvctype=="PE"|tvctype=="TI")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp)
  if(tvctype=="ED"|tvctype=="CO")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp,eta)
  print(round(dest,4))
  
  cagep=cageptotal[statusp==1|statusp==2]
  genderp=genderptotal[statusp==1|statusp==2]
  st1p=st1ptotal[statusp==1|statusp==2]
  # sum of survival until event time for all the individuals in the data
  screenp0<-which(genderp==0)
  screenp1<-which(genderp==1)
  
  if(tvctype=="PE"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1),"PE")}#;cH2=cH_PE(st1,time0,c(lambda2,rho2,beta.sc1_2))}
  if(tvctype=="ED"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2]),"ED")}
  if(tvctype=="CO"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2],eta[4]),"CO")
  
  Hp1=vector()
  if(tvctype!="TI"){
    Hp1[screenp0]<- exp(beta.gen1)*((lambda1*cagep[screenp0])^rho1)
    Hp1[screenp1]<- exp(beta.gen1)*((lambda1*st1p[screenp1])^rho1+cH1p)
  }else{
    Hp1[screenp0]<- exp(beta.gen1)*((lambda1*cagep[screenp0])^rho1)
    Hp1[screenp1]<- exp(beta.gen1)*((lambda1*cagep[screenp1])^rho1)*exp(beta.sc1_1)  
  }
  
  Hp2<- exp(beta.gen2)*((lambda2*cagep)^rho2)
  
  if(frailty==TRUE){
    slogasc12=sum(log(1- ((1+Hp1/fp[1])^(-fp[1]))*((1+Hp2/fp[2])^(-fp[2]))),na.rm=TRUE)
  }else{
    slogasc12=sum(log(1- exp(-Hp1-Hp2)),na.rm=TRUE)  
  }
  cagep0=cageptotal[statusp==0]
  genderp0=genderptotal[statusp==0]
  st1p0=st1ptotal[statusp==0]
  # sum of survival until event time for all the individuals in the data
  screenp00<-which(genderp0==0)
  screenp10<-which(genderp0==1)
  
  if(tvctype=="PE"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1),"PE")}#;cH2=cH_PE(st1,time0,c(lambda2,rho2,beta.sc1_2))}
  if(tvctype=="ED"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2]),"ED")}
  if(tvctype=="CO"){cH1p0=cH(st1p0[screenp10],cagep0[screenp10],c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2],eta[4]),"CO")
  
  Hp10=vector()
  if(tvctype!="TI"){
    Hp10[screenp00]<- exp(beta.gen1)*((lambda1*cagep0[screenp00])^rho1)
    Hp10[screenp10]<- exp(beta.gen1)*((lambda1*st1p0[screenp10])^rho1+cH1p0)
  }else{
    Hp10[screenp00]<- exp(beta.gen1)*((lambda1*cagep0[screenp00])^rho1)
    Hp10[screenp10]<- exp(beta.gen1)*((lambda1*cagep0[screenp10])^rho1)*exp(beta.sc1_1)
  }
  
  Hp20<- exp(beta.gen2)*((lambda2*cagep0)^rho2)
  
  if(frailty==TRUE){
    slogasc0=sum(   log(  ((1+Hp10/fp[1])^(-fp[1]))*((1+Hp20/fp[2])^(-fp[2]))  )    ,na.rm=TRUE)
  }else{
    slogasc0=sum(   -Hp10-Hp20,na.rm=TRUE)  
  }
  
  #print(round(c(cest,fp),3))
  # cagep1.0 <- cagep[statusp==1&genderp==0]
  # cagep1.1 <- cagep[statusp==1&genderp==1]
  # st1p1<-st1p[statusp==1&genderp==1]
  # 
  # if(length(cagep1.0)!=0){
  # logasc1.0 = log(sapply(1:length(cagep1.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep1.0[index])))
  # }else{logasc1.0=0}
  # if(length(cagep1.1)!=0){
  # logasc1.1 = log(sapply(1:length(cagep1.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p1[index]))+
  #                 sapply(1:length(cagep1.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=1,tvc=st1p1[index],tt=tvctype,eta=eta,lower=st1p1[index],upper=cagep1.1[index])))
  # }else{logasc1.1=0}  
  # logasc1 <- sum(logasc1.0,na.rm=TRUE)+sum(logasc1.1,na.rm=TRUE)
  # 
  # # subsetting the probands according to his affection status at study entry
  # cagep2.0 <- cagep[statusp==2&genderp==0]
  # cagep2.1 <- cagep[statusp==2&genderp==1]
  # st1p2<-st1p[statusp==2&genderp==1]
  # 
  # if(length(cagep2.0)!=0){
  # logasc2.0 = log(sapply(1:length(cagep2.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep2.0[index])))
  # }else{logasc2.0=0}
  # if(length(cagep2.1)!=0){
  # logasc2.1 = log(sapply(1:length(cagep2.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p2[index]))+
  #             sapply(1:length(cagep2.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=2,tvc=st1p2[index],tt=tvctype,eta=eta,lower=st1p2[index],upper=cagep2.1[index])))
  # }else{logasc2.1=0}
  # logasc2 <- sum(logasc2.0,na.rm=TRUE)+sum(logasc2.1,na.rm=TRUE)
  
  slogasc = slogasc12 + slogasc0
  likelihood  <- loglik - slogasc
  
  return(-likelihood)
} 
loglik_Comp_simulation3_RD=function(theta, data, agemin, frailty,tvctype){
  
  data = data[data$currentage>=agemin,]
  # base parameters
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  
  
  # vbeta for breast  
  beta.sc1_1 = theta[5]
  beta.gen1 = theta[6]
  
  # vbeta for ovarian
  #beta.sc1_2 =  theta[7]
  beta.gen2 = theta[7]
  
  if(tvctype=="ED"){
    eta <- exp(theta[8])
  }else if(tvctype=="CO"){
    #eta <- c(exp(theta[10]),theta[11])
    eta <- c(exp(theta[8]),theta[9])
  }
  
  # frailty parameter
  if(frailty==TRUE){
    fp <- exp(theta[8:9])
  }else{
    fp=NULL
  }
  
  
  # Y, delta, carrier prob, proband indicator, sc number
  time0 = data$time-agemin
  status = data$status
  ip = which(data$proband==1)
  ex1 = data$carrp.geno 
  #gender = data$gender
  st1=data$tvc-agemin
  gender=data$TVCstatus
  st1[gender==0]<-time0[gender==0]
  
  # baseline hazard(weibull), baseline cumulative hazard(weibull)
  bhaz1 = (lambda1^rho1)*rho1*time0^(rho1-1)
  bhaz2 = (lambda2^rho2)*rho2*time0^(rho2-1)
  
  logh1.0 = log(bhaz1) 
  logh1.1 = log(bhaz1) + beta.gen1
  logh2.0 = log(bhaz2) 
  logh2.1 = log(bhaz2) + beta.gen2
  
  
  # creating screen,gene,oophorectomy interaction vectors for hazard
  g=function(t,tsc,eta){exp(-eta*(t-tsc))}
  if(tvctype=="TI")scvec1=gender*beta.sc1_1
  if(tvctype=="PE")scvec1=gender*beta.sc1_1#;scvec2=gender*beta.sc1_2
  if(tvctype=="ED")scvec1=gender*(g(time0,st1,eta[1])*beta.sc1_1)#;scvec2=gender*g(time0,st1,eta[2])*beta.sc1_2
  if(tvctype=="CO")scvec1=gender*(g(time0,st1,eta[1])*beta.sc1_1+eta[2])#;scvec2=gender*(g(time0,st1,eta[2])*beta.sc1_2+eta[4])
  
  # sum of log-hazard
  sum1 = sum((      ex1*(logh1.1 + scvec1 ) +
                      (1-ex1)*(logh1.0 + scvec1))[status==1], na.rm=TRUE) 
  sum2 = sum((      ex1*(logh2.1) +
                      (1-ex1)*(logh2.0))[status==2], na.rm=TRUE)              
  
  if(tvctype=="PE"){cH1=cH(st1,time0,c(lambda1,rho1,beta.sc1_1),"PE")}#;cH2=cH_PE(st1,time0,c(lambda2,rho2,beta.sc1_2))}
  if(tvctype=="ED"){cH1=cH(st1,time0,c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2]),"ED")}
  if(tvctype=="CO"){cH1=cH(st1,time0,c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2],eta[4]),"CO")
  
  # sum of survival until event time for all the individuals in the data
  if(tvctype!="TI"){
    CH0_bc <-              (1-ex1)*((lambda1*st1)^rho1+cH1)
    CH1_bc <- exp(beta.gen1)*(ex1)*((lambda1*st1)^rho1+cH1)
  }else{
    CH0_bc <-              (1-ex1)*((lambda1*time0)^rho1)*exp(scvec1)
    CH1_bc <- exp(beta.gen1)*(ex1)*((lambda1*time0)^rho1)*exp(scvec1)  
  }
  
  CH0_ov <-              (1-ex1)*((lambda2*time0)^rho2)
  CH1_ov <- exp(beta.gen2)*(ex1)*((lambda2*time0)^rho2)
  
  
  k <- CH0_bc + CH1_bc + CH0_ov + CH1_ov
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam1 <- (CH0_bc+CH1_bc)
    Hfam2 <- (CH0_ov+CH1_ov)
    
    df1 <- data$df1[data$proband==1]
    df2 <- data$df2[data$proband==1]
    
    Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
    Hfam2 <- aggregate(Hfam2,by=list(data$famID),FUN=sum)[,2]
    
    
    sum4.1 <- sum((lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1])), na.rm=T)
    sum4.2 <- sum((lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2])), na.rm=T)
    
    loglik = sum1+sum2 + sum4.1+ sum4.2
  }else{
    sum4 <- sum(-k,na.rm = TRUE)
    loglik = sum1+sum2+sum4 # numerator in loglikelihood
  }
  
  # Ascertainment correction by design="pop+"
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.gen2)
  cageptotal <- data$currentage[ip]-agemin
  statusp <- data$affect[ip]  # proband disease status at the study entry
  genderptotal = data$screenp[ip]
  st1ptotal=data$st1p[ip]-agemin
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  
  if(tvctype=="PE"|tvctype=="TI")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp)
  if(tvctype=="ED"|tvctype=="CO")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp,eta)
  print(round(dest,4))
  
  cagep=cageptotal[statusp==1|statusp==2]
  genderp=genderptotal[statusp==1|statusp==2]
  st1p=st1ptotal[statusp==1|statusp==2]
  # sum of survival until event time for all the individuals in the data
  screenp0<-which(genderp==0)
  screenp1<-which(genderp==1)
  
  if(tvctype=="PE"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1),"PE")}#;cH2=cH_PE(st1,time0,c(lambda2,rho2,beta.sc1_2))}
  if(tvctype=="ED"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2]),"ED")}
  if(tvctype=="CO"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2],eta[4]),"CO")
  
  Hp1=vector()
  if(tvctype!="TI"){
    Hp1[screenp0]<- exp(beta.gen1)*((lambda1*cagep[screenp0])^rho1)
    Hp1[screenp1]<- exp(beta.gen1)*((lambda1*st1p[screenp1])^rho1+cH1p)
  }else{
    Hp1[screenp0]<- exp(beta.gen1)*((lambda1*cagep[screenp0])^rho1)
    Hp1[screenp1]<- exp(beta.gen1)*((lambda1*cagep[screenp1])^rho1)*exp(beta.sc1_1)  
  }
  
  Hp2<- exp(beta.gen2)*((lambda2*cagep)^rho2)
  
  if(frailty==TRUE){
  slogasc12=sum(log(1- ((1+Hp1/fp[1])^(-fp[1]))*((1+Hp2/fp[2])^(-fp[2]))),na.rm=TRUE)
  }else{
  slogasc12=sum(log(1- exp(-Hp1-Hp2)),na.rm=TRUE)  
  }
  cagep0=cageptotal[statusp==0]
  genderp0=genderptotal[statusp==0]
  st1p0=st1ptotal[statusp==0]
  # sum of survival until event time for all the individuals in the data
  screenp00<-which(genderp0==0)
  screenp10<-which(genderp0==1)

  if(tvctype=="PE"){cH1p0=cH(st1p0[screenp10],cagep0[screenp10],c(lambda1,rho1,beta.sc1_1),"PE")}#;cH2=cH_PE(st1,time0,c(lambda2,rho2,beta.sc1_2))}
  if(tvctype=="ED"){cH1p0=cH(st1p0[screenp10],cagep0[screenp10],c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2]),"ED")}
  if(tvctype=="CO"){cH1p0=cH(st1p0[screenp10],cagep0[screenp10],c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2],eta[4]),"CO")

  Hp10=vector()
  if(tvctype!="TI"){
    Hp10[screenp00]<- exp(beta.gen1)*((lambda1*cagep0[screenp00])^rho1)
    Hp10[screenp10]<- exp(beta.gen1)*((lambda1*st1p0[screenp10])^rho1+cH1p0)
  }else{
    Hp10[screenp00]<- exp(beta.gen1)*((lambda1*cagep0[screenp00])^rho1)
    Hp10[screenp10]<- exp(beta.gen1)*((lambda1*cagep0[screenp10])^rho1)*exp(beta.sc1_1)
  }

  Hp20<- exp(beta.gen2)*((lambda2*cagep0)^rho2)
  
  if(frailty==TRUE){
  slogasc0=sum(   log(  ((1+Hp10/fp[1])^(-fp[1]))*((1+Hp20/fp[2])^(-fp[2]))  )    ,na.rm=TRUE)
  }else{
  slogasc0=sum(   -Hp10-Hp20,na.rm=TRUE)  
  }
  
  #print(round(c(cest,fp),3))
  # cagep1.0 <- cagep[statusp==1&genderp==0]
  # cagep1.1 <- cagep[statusp==1&genderp==1]
  # st1p1<-st1p[statusp==1&genderp==1]
  # 
  # if(length(cagep1.0)!=0){
  # logasc1.0 = log(sapply(1:length(cagep1.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep1.0[index])))
  # }else{logasc1.0=0}
  # if(length(cagep1.1)!=0){
  # logasc1.1 = log(sapply(1:length(cagep1.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p1[index]))+
  #                 sapply(1:length(cagep1.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=1,tvc=st1p1[index],tt=tvctype,eta=eta,lower=st1p1[index],upper=cagep1.1[index])))
  # }else{logasc1.1=0}  
  # logasc1 <- sum(logasc1.0,na.rm=TRUE)+sum(logasc1.1,na.rm=TRUE)
  # 
  # # subsetting the probands according to his affection status at study entry
  # cagep2.0 <- cagep[statusp==2&genderp==0]
  # cagep2.1 <- cagep[statusp==2&genderp==1]
  # st1p2<-st1p[statusp==2&genderp==1]
  # 
  # if(length(cagep2.0)!=0){
  # logasc2.0 = log(sapply(1:length(cagep2.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep2.0[index])))
  # }else{logasc2.0=0}
  # if(length(cagep2.1)!=0){
  # logasc2.1 = log(sapply(1:length(cagep2.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p2[index]))+
  #             sapply(1:length(cagep2.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=2,tvc=st1p2[index],tt=tvctype,eta=eta,lower=st1p2[index],upper=cagep2.1[index])))
  # }else{logasc2.1=0}
  # logasc2 <- sum(logasc2.0,na.rm=TRUE)+sum(logasc2.1,na.rm=TRUE)
  
  slogasc = slogasc12 + slogasc0
  likelihood  <- loglik - slogasc
  
  return(-likelihood)
} 
loglik_Comp_simulation3_RD_second=function(theta, data, agemin, frailty,tvctype,stage1){
  
  data = data[data$currentage>=agemin,]
  # base parameters
  lambda1 = exp(stage1[1])
  rho1 = exp(stage1[2])
  lambda2 = exp(stage1[3])
  rho2 = exp(stage1[4])
  
  
  # vbeta for breast  
  beta.sc1_1 = stage1[5]
  beta.gen1 = stage1[6]
  
  # vbeta for ovarian
  #beta.sc1_2 =  theta[7]
  beta.gen2 = stage1[7]
  
  if(tvctype=="ED"){
    eta <- exp(stage1[8])
  }else if(tvctype=="CO"){
    #eta <- c(exp(theta[10]),theta[11])
    eta <- c(exp(stage1[8]),stage1[9])
  }
  
  # frailty parameter
  if(frailty==TRUE){
    fp <- exp(theta[1:2])
  }
  
  
  # Y, delta, carrier prob, proband indicator, sc number
  time0 = data$time-agemin
  status = data$status
  ip = which(data$proband==1)
  ex1 = data$carrp.geno 
  #gender = data$gender
  st1=data$tvc-agemin
  gender=data$TVCstatus
  st1[gender==0]<-time0[gender==0]
  
  # baseline hazard(weibull), baseline cumulative hazard(weibull)
  bhaz1 = (lambda1^rho1)*rho1*time0^(rho1-1)
  bhaz2 = (lambda2^rho2)*rho2*time0^(rho2-1)
  
  logh1.0 = log(bhaz1) 
  logh1.1 = log(bhaz1) + beta.gen1
  logh2.0 = log(bhaz2) 
  logh2.1 = log(bhaz2) + beta.gen2
  
  
  # creating screen,gene,oophorectomy interaction vectors for hazard
  g=function(t,tsc,eta){exp(-eta*(t-tsc))}
  if(tvctype=="TI")scvec1=gender*beta.sc1_1
  if(tvctype=="PE")scvec1=gender*beta.sc1_1#;scvec2=gender*beta.sc1_2
  if(tvctype=="ED")scvec1=gender*(g(time0,st1,eta[1])*beta.sc1_1)#;scvec2=gender*g(time0,st1,eta[2])*beta.sc1_2
  if(tvctype=="CO")scvec1=gender*(g(time0,st1,eta[1])*beta.sc1_1+eta[2])#;scvec2=gender*(g(time0,st1,eta[2])*beta.sc1_2+eta[4])
  
  # sum of log-hazard
  sum1 = sum((      ex1*(logh1.1 + scvec1 ) +
                      (1-ex1)*(logh1.0 + scvec1))[status==1], na.rm=TRUE) 
  sum2 = sum((      ex1*(logh2.1) +
                      (1-ex1)*(logh2.0))[status==2], na.rm=TRUE)              
  
  if(tvctype=="PE"){cH1=cH(st1,time0,c(lambda1,rho1,beta.sc1_1),"PE")}#;cH2=cH_PE(st1,time0,c(lambda2,rho2,beta.sc1_2))}
  if(tvctype=="ED"){cH1=cH(st1,time0,c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2]),"ED")}
  if(tvctype=="CO"){cH1=cH(st1,time0,c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2],eta[4]),"CO")
  
  # sum of survival until event time for all the individuals in the data
  if(tvctype!="TI"){
    CH0_bc <-              (1-ex1)*((lambda1*st1)^rho1+cH1)
    CH1_bc <- exp(beta.gen1)*(ex1)*((lambda1*st1)^rho1+cH1)
  }else{
    CH0_bc <-              (1-ex1)*((lambda1*time0)^rho1)*exp(scvec1)
    CH1_bc <- exp(beta.gen1)*(ex1)*((lambda1*time0)^rho1)*exp(scvec1)  
  }
  
  CH0_ov <-              (1-ex1)*((lambda2*time0)^rho2)
  CH1_ov <- exp(beta.gen2)*(ex1)*((lambda2*time0)^rho2)
  
  
  k <- CH0_bc + CH1_bc + CH0_ov + CH1_ov
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam1 <- (CH0_bc+CH1_bc)
    Hfam2 <- (CH0_ov+CH1_ov)
    
    df1 <- data$df1[data$proband==1]
    df2 <- data$df2[data$proband==1]
    
    Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
    Hfam2 <- aggregate(Hfam2,by=list(data$famID),FUN=sum)[,2]
    
    
    sum4.1 <- sum((lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1])), na.rm=T)
    sum4.2 <- sum((lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2])), na.rm=T)
    
    loglik = sum1+sum2 + sum4.1+ sum4.2
  }else{
    sum4 <- sum(-k,na.rm = TRUE)
    loglik = sum1+sum2+sum4 # numerator in loglikelihood
  }
  
  # Ascertainment correction by design="pop+"
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.gen2)
  cageptotal <- data$currentage[ip]-agemin
  statusp <- data$affect[ip]  # proband disease status at the study entry
  genderptotal = data$screenp[ip]
  st1ptotal=data$st1p[ip]-agemin
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  
  if(tvctype=="PE"|tvctype=="TI")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp)
  if(tvctype=="ED"|tvctype=="CO")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp,eta)
  print(round(dest,4))
  
  cagep=cageptotal[statusp==1|statusp==2]
  genderp=genderptotal[statusp==1|statusp==2]
  st1p=st1ptotal[statusp==1|statusp==2]
  # sum of survival until event time for all the individuals in the data
  screenp0<-which(genderp==0)
  screenp1<-which(genderp==1)
  
  if(tvctype=="PE"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1),"PE")}#;cH2=cH_PE(st1,time0,c(lambda2,rho2,beta.sc1_2))}
  if(tvctype=="ED"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2]),"ED")}
  if(tvctype=="CO"){cH1p=cH(st1p[screenp1],cagep[screenp1],c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2],eta[4]),"CO")
  
  Hp1=vector()
  if(tvctype!="TI"){
    Hp1[screenp0]<- exp(beta.gen1)*((lambda1*cagep[screenp0])^rho1)
    Hp1[screenp1]<- exp(beta.gen1)*((lambda1*st1p[screenp1])^rho1+cH1p)
  }else{
    Hp1[screenp0]<- exp(beta.gen1)*((lambda1*cagep[screenp0])^rho1)
    Hp1[screenp1]<- exp(beta.gen1)*((lambda1*cagep[screenp1])^rho1)*exp(beta.sc1_1)  
  }
  
  Hp2<- exp(beta.gen2)*((lambda2*cagep)^rho2)
  
  if(frailty==TRUE){
    slogasc12=sum(log(1- ((1+Hp1/fp[1])^(-fp[1]))*((1+Hp2/fp[2])^(-fp[2]))),na.rm=TRUE)
  }else{
    slogasc12=sum(log(1- exp(-Hp1-Hp2)),na.rm=TRUE)  
  }
  cagep0=cageptotal[statusp==0]
  genderp0=genderptotal[statusp==0]
  st1p0=st1ptotal[statusp==0]
  # sum of survival until event time for all the individuals in the data
  screenp00<-which(genderp0==0)
  screenp10<-which(genderp0==1)
  
  if(tvctype=="PE"){cH1p0=cH(st1p0[screenp10],cagep0[screenp10],c(lambda1,rho1,beta.sc1_1),"PE")}#;cH2=cH_PE(st1,time0,c(lambda2,rho2,beta.sc1_2))}
  if(tvctype=="ED"){cH1p0=cH(st1p0[screenp10],cagep0[screenp10],c(lambda1,rho1,beta.sc1_1,eta[1]),"ED")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2]),"ED")}
  if(tvctype=="CO"){cH1p0=cH(st1p0[screenp10],cagep0[screenp10],c(lambda1,rho1,beta.sc1_1,eta[1],eta[2]),"CO")}#;cH2=cH(st1,time0,c(lambda2,rho2,beta.sc1_2,eta[2],eta[4]),"CO")
  
  Hp10=vector()
  if(tvctype!="TI"){
    Hp10[screenp00]<- exp(beta.gen1)*((lambda1*cagep0[screenp00])^rho1)
    Hp10[screenp10]<- exp(beta.gen1)*((lambda1*st1p0[screenp10])^rho1+cH1p0)
  }else{
    Hp10[screenp00]<- exp(beta.gen1)*((lambda1*cagep0[screenp00])^rho1)
    Hp10[screenp10]<- exp(beta.gen1)*((lambda1*cagep0[screenp10])^rho1)*exp(beta.sc1_1)
  }
  
  Hp20<- exp(beta.gen2)*((lambda2*cagep0)^rho2)
  
  if(frailty==TRUE){
    slogasc0=sum(   log(  ((1+Hp10/fp[1])^(-fp[1]))*((1+Hp20/fp[2])^(-fp[2]))  )    ,na.rm=TRUE)
  }else{
    slogasc0=sum(   -Hp10-Hp20,na.rm=TRUE)  
  }
  
  #print(round(c(cest,fp),3))
  # cagep1.0 <- cagep[statusp==1&genderp==0]
  # cagep1.1 <- cagep[statusp==1&genderp==1]
  # st1p1<-st1p[statusp==1&genderp==1]
  # 
  # if(length(cagep1.0)!=0){
  # logasc1.0 = log(sapply(1:length(cagep1.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep1.0[index])))
  # }else{logasc1.0=0}
  # if(length(cagep1.1)!=0){
  # logasc1.1 = log(sapply(1:length(cagep1.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p1[index]))+
  #                 sapply(1:length(cagep1.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=1,tvc=st1p1[index],tt=tvctype,eta=eta,lower=st1p1[index],upper=cagep1.1[index])))
  # }else{logasc1.1=0}  
  # logasc1 <- sum(logasc1.0,na.rm=TRUE)+sum(logasc1.1,na.rm=TRUE)
  # 
  # # subsetting the probands according to his affection status at study entry
  # cagep2.0 <- cagep[statusp==2&genderp==0]
  # cagep2.1 <- cagep[statusp==2&genderp==1]
  # st1p2<-st1p[statusp==2&genderp==1]
  # 
  # if(length(cagep2.0)!=0){
  # logasc2.0 = log(sapply(1:length(cagep2.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep2.0[index])))
  # }else{logasc2.0=0}
  # if(length(cagep2.1)!=0){
  # logasc2.1 = log(sapply(1:length(cagep2.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p2[index]))+
  #             sapply(1:length(cagep2.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=2,tvc=st1p2[index],tt=tvctype,eta=eta,lower=st1p2[index],upper=cagep2.1[index])))
  # }else{logasc2.1=0}
  # logasc2 <- sum(logasc2.0,na.rm=TRUE)+sum(logasc2.1,na.rm=TRUE)
  
  slogasc = slogasc12 + slogasc0
  likelihood  <- loglik - slogasc
  
  return(-likelihood)
} 

trueparm= c(0.01,2.5, 0.01,2.5, 2,2,2, 3.5,5,  0.1,0.2) ;mrate=0#CO
                                              
realdataestimate1=Comp_simulation3_RD_one(b1,trueparm,tvctype="CO")#survivedproband  corrected
realdataestimate1=Comp_simulation3_RD_one(b1,trueparm,tvctype="CO")

realdataestimate2=Comp_simulation3_RD(b1,trueparm,tvctype="ED")#survivedproband  corrected
realdataestimatePE=Comp_simulation3_RD(b1,trueparm,tvctype="PE")#
realdataestimateCO=Comp_simulation3_RD(b1,trueparm,tvctype="CO")#

1                                                                        0.008
2                                                                        2.300
3                                                                        0.007
4                                                                        2.932
5                                                                        1.872
6                                                                        1.858
7                                                                        1.224
8                                                                       10.000
9                                                                        3.240
10                                                                       0.278

