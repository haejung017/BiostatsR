# First Run everything after line 15 in the file

# Run 
# RobustSE, delta method
rse=matrix(ncol=9,nrow=100)
pse1=matrix(ncol=16,nrow=100)
pse2=matrix(ncol=16,nrow=100)
for (i in 1:100){
obj=variance(i)
print(i)
rse[i,]=obj[[1]]
pse1[i,]=obj[[2]]
pse2[i,]=obj[[3]]
}
# Run 
# RobustSE, delta method
qrse=matrix(ncol=9,nrow=100)
qpse1=matrix(ncol=16,nrow=100)
qpse2=matrix(ncol=16,nrow=100)
for (i in 101:200){
  obj=variance(i)
  print(i)
  i=i-100
  qrse[i,]=obj[[1]]
  qpse1[i,]=obj[[2]]
  qpse2[i,]=obj[[3]]
}
# Run 
# RobustSE, delta method
qqrse=matrix(ncol=9,nrow=100)
qqpse1=matrix(ncol=16,nrow=100)
qqpse2=matrix(ncol=16,nrow=100)
for (i in 201:300){
  obj=variance(i)
  print(i)
  i=i-200
  qqrse[i,]=obj[[1]]
  qqpse1[i,]=obj[[2]]
  qqpse2[i,]=obj[[3]]
}
# RobustSE, delta method
qqqrse=matrix(ncol=9,nrow=100)
qqqpse1=matrix(ncol=16,nrow=100)
qqqpse2=matrix(ncol=16,nrow=100)
for (i in 301:400){
  obj=variance(i)
  print(i)
  i=i-300
  qqqrse[i,]=obj[[1]]
  qqqpse1[i,]=obj[[2]]
  qqqpse2[i,]=obj[[3]]
}
# RobustSE, delta method
qqqqrse=matrix(ncol=9,nrow=100)
qqqqpse1=matrix(ncol=16,nrow=100)
qqqqpse2=matrix(ncol=16,nrow=100)
for (i in 401:500){
  obj=variance(i)
  print(i)
  i=i-400
  qqqqrse[i,]=obj[[1]]
  qqqqpse1[i,]=obj[[2]]
  qqqqpse2[i,]=obj[[3]]
}













# HessianSE, delta method
drse=matrix(ncol=9,nrow=100)
dpse1=matrix(ncol=16,nrow=100)
dpse2=matrix(ncol=16,nrow=100)
for (i in 1:100){
  obj=variance2(i)
  print(i)
  drse[i,]=obj[[1]]
  dpse1[i,]=obj[[2]]
  dpse2[i,]=obj[[3]]
}





x=rbind(rse,qrse,qqrse,qqqrse,qqqqrse)
y=rbind(pse1,qpse1,qqpse1,qqqpse1,qqqqpse1)
z=rbind(pse2,qpse2,qqpse2,qqqpse2,qqqqpse2)

a=rbind(drse,qdrse)
b=rbind(qdpse1,qdpse1)
c=rbind(qdpse2,qdpse2)
# HessianSE, delta method
qdrse=matrix(ncol=9,nrow=100)
qdpse1=matrix(ncol=16,nrow=100)
qdpse2=matrix(ncol=16,nrow=100)
for (i in 101:200){
  obj=variance2(i)
  print(i)
  i=i-100
  qdrse[i,]=obj[[1]]
  qdpse1[i,]=obj[[2]]
  qdpse2[i,]=obj[[3]]
}

write.excel(
rbind(
cbind(colMeans(x)),
cbind(colMeans(y)),
cbind(colMeans(z))
))
write.excel(
rbind(
cbind(colMeans(a)),
cbind(colMeans(b)),
cbind(colMeans(c))
))
# [[1]] parameters SE
# [[2]] penetrance by 70 SE for event 1, ESE is 

petest=excelresult(PE_i,list(RSE=x,PSE1=y,PSE2=z),"PE",10)
write.excel(petest$cov)
write.excel(petest$pen1)
write.excel(petest$pen2)
#0.0086 for NS,NC
#0.015 for S,NC
#0.0231 for NS,C
#0.0233 for S,C


library(pracma)
library(MASS)
loglik_variance_frailty=function(theta,data, agemin, frailty=TRUE,tvctype){
  
  data = data[data$currentage>=agemin,]
  theta = theta
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
    eta <- exp(stage1[8])
  }else if(tvctype=="CO"){
    #eta <- c(exp(theta[10]),theta[11])
    eta <- c(exp(stage1[8]),stage1[9])
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
  proband=data$proband
  #ex1 = data$carrp.geno 
  gender=data$TVCstatus
  mgene=data$mgene
  st1=data$tvc-agemin
  genderpp=data$TVCstatusp
  st1[gender==0&proband==0]<-time0[gender==0&proband==0]
  
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
  
  #### UNI
  loglik <- ifelse(status==1, logh1 + scvec1, ifelse(status==2, logh2, 0))
  
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
    Hfam1 <- CH_bc
    Hfam2 <- CH_ov
    
    df1 <- data$df1[data$proband==1]
    df2 <- data$df2[data$proband==1]
    
    Hfam1 <- aggregate(Hfam1, by=list(data$famID), FUN=sum)[,2]
    Hfam2 <- aggregate(Hfam2, by=list(data$famID), FUN=sum)[,2]
    
    
    sum4.1 <- (lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1]))
    sum4.2 <- (lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2]))
    
    ####  UNI
    loglikfam = aggregate(loglik, by=list(data$famID), sum)[,2]
    loglikfam = loglikfam + sum4.1 + sum4.2 
    #loglik = sum1+sum2 + sum4.1+ sum4.2
  }else{
    sum4 <- sum(-k,na.rm = TRUE)
    loglik <- loglik - k
    #loglik = sum1+sum2+sum4 # numerator in loglikelihood
    loglikfam = aggregate(loglik, by=list(data$famID), sum)[,2]
  }
  
  # Ascertainment correction by design="pop+"
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.gen2)
  cagep <- data$currentage[ip]-agemin
  statusp <- data$status[ip]  # proband disease status at the study entry
  genderp = genderpp[ip]
  st1p=st1[ip]
  st1p[genderp==0]<-cagep[genderp==0]
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  
  if(tvctype=="PE"|tvctype=="TI")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp)
  if(tvctype=="ED"|tvctype=="CO")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp,eta)
  #print(round(dest,4))
  
  #-----------------------------------------------------------------------------------
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
    ### UNI
    logasc= log(1- ((1+Hp1/fp[1])^(-fp[1]))*((1+Hp2/fp[2])^(-fp[2])))
    loglikfam = loglikfam - logasc
    
  }else{
    logasc = log(1- exp(-Hp1-Hp2))
    loglikfam = loglikfam - logasc
  }
  #-----------------------------------------------------------------------------------
  
  # #-----------------------------------------------------------------------------------
  # cagep1.0 <- cagep[statusp==1&genderp==0]
  # cagep1.1 <- cagep[statusp==1&genderp==1]
  # st1p1<-st1p[statusp==1&genderp==1]
  # 
  # if(length(cagep1.0)!=0){
  #   logasc1.0 = log(sapply(1:length(cagep1.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep1.0[index])))
  # }else{logasc1.0=0}
  # if(length(cagep1.1)!=0){
  #   logasc1.1 = log(sapply(1:length(cagep1.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p1[index]))+
  #                     sapply(1:length(cagep1.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=1,tvc=st1p1[index],tt=tvctype,eta=eta,lower=st1p1[index],upper=cagep1.1[index])))
  # }else{logasc1.1=0}
  # logasc1 <- sum(logasc1.0,na.rm=TRUE)+sum(logasc1.1,na.rm=TRUE)
  # 
  # # subsetting the probands according to his affection status at study entry
  # cagep2.0 <- cagep[statusp==2&genderp==0]
  # cagep2.1 <- cagep[statusp==2&genderp==1]
  # st1p2<-st1p[statusp==2&genderp==1]
  # 
  # if(length(cagep2.0)!=0){
  #   logasc2.0 = log(sapply(1:length(cagep2.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep2.0[index])))
  # }else{logasc2.0=0}
  # if(length(cagep2.1)!=0){
  #   logasc2.1 = log(sapply(1:length(cagep2.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p2[index]))+
  #                     sapply(1:length(cagep2.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=2,tvc=st1p2[index],tt=tvctype,eta=eta,lower=st1p2[index],upper=cagep2.1[index])))
  # }else{logasc2.1=0}
  # logasc2 <- sum(logasc2.0,na.rm=TRUE)+sum(logasc2.1,na.rm=TRUE)
  # slogasc = logasc1 + logasc2
  # #-----------------------------------------------------------------------------------
  
  #likelihood  <- loglik - slogasc
  sloglikfam <- sum(loglikfam)
  return(-sloglikfam)
}
loglik_variance_vec_frailty=function(theta,data, agemin, frailty=TRUE,tvctype){
  
  data = data[data$currentage>=agemin,]
  theta = theta
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
    #eta <- c(exp(theta[10]),theta[11])
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
  proband=data$proband
  #ex1 = data$carrp.geno 
  gender=data$TVCstatus
  mgene=data$mgene
  st1=data$tvc-agemin
  genderpp=data$TVCstatusp
  st1[gender==0&proband==0]<-time0[gender==0&proband==0]
  
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
  
  #### UNI
  loglik <- ifelse(status==1, logh1 + scvec1, ifelse(status==2, logh2, 0))
  
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
    Hfam1 <- CH_bc
    Hfam2 <- CH_ov
    
    df1 <- data$df1[data$proband==1]
    df2 <- data$df2[data$proband==1]
    
    Hfam1 <- aggregate(Hfam1, by=list(data$famID), FUN=sum)[,2]
    Hfam2 <- aggregate(Hfam2, by=list(data$famID), FUN=sum)[,2]
    
    
    sum4.1 <- (lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1]))
    sum4.2 <- (lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2]))
    
    ####  UNI
    loglikfam = aggregate(loglik, by=list(data$famID), sum)[,2]
    loglikfam = loglikfam + sum4.1 + sum4.2 
    #loglik = sum1+sum2 + sum4.1+ sum4.2
  }else{
    sum4 <- sum(-k,na.rm = TRUE)
    loglik <- loglik - k
    #loglik = sum1+sum2+sum4 # numerator in loglikelihood
    loglikfam = aggregate(loglik, by=list(data$famID), sum)[,2]
  }
  
  # Ascertainment correction by design="pop+"
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.gen2)
  cagep <- data$currentage[ip]-agemin
  statusp <- data$status[ip]  # proband disease status at the study entry
  genderp = genderpp[ip]
  st1p=st1[ip]
  st1p[genderp==0]<-cagep[genderp==0]
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  
  if(tvctype=="PE"|tvctype=="TI")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp)
  if(tvctype=="ED"|tvctype=="CO")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp,eta)
  #print(round(dest,4))
  
  #-----------------------------------------------------------------------------------
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
    ### UNI
    logasc= log(1- ((1+Hp1/fp[1])^(-fp[1]))*((1+Hp2/fp[2])^(-fp[2])))
    loglikfam = loglikfam - logasc
    
  }else{
    logasc = log(1- exp(-Hp1-Hp2))
    loglikfam = loglikfam - logasc
  }
  #-----------------------------------------------------------------------------------
  
  # #-----------------------------------------------------------------------------------
  # cagep1.0 <- cagep[statusp==1&genderp==0]
  # cagep1.1 <- cagep[statusp==1&genderp==1]
  # st1p1<-st1p[statusp==1&genderp==1]
  # 
  # if(length(cagep1.0)!=0){
  #   logasc1.0 = log(sapply(1:length(cagep1.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep1.0[index])))
  # }else{logasc1.0=0}
  # if(length(cagep1.1)!=0){
  #   logasc1.1 = log(sapply(1:length(cagep1.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p1[index]))+
  #                     sapply(1:length(cagep1.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=1,tvc=st1p1[index],tt=tvctype,eta=eta,lower=st1p1[index],upper=cagep1.1[index])))
  # }else{logasc1.1=0}
  # logasc1 <- sum(logasc1.0,na.rm=TRUE)+sum(logasc1.1,na.rm=TRUE)
  # 
  # # subsetting the probands according to his affection status at study entry
  # cagep2.0 <- cagep[statusp==2&genderp==0]
  # cagep2.1 <- cagep[statusp==2&genderp==1]
  # st1p2<-st1p[statusp==2&genderp==1]
  # 
  # if(length(cagep2.0)!=0){
  #   logasc2.0 = log(sapply(1:length(cagep2.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep2.0[index])))
  # }else{logasc2.0=0}
  # if(length(cagep2.1)!=0){
  #   logasc2.1 = log(sapply(1:length(cagep2.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p2[index]))+
  #                     sapply(1:length(cagep2.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=2,tvc=st1p2[index],tt=tvctype,eta=eta,lower=st1p2[index],upper=cagep2.1[index])))
  # }else{logasc2.1=0}
  # logasc2 <- sum(logasc2.0,na.rm=TRUE)+sum(logasc2.1,na.rm=TRUE)
  # slogasc = logasc1 + logasc2
  # #-----------------------------------------------------------------------------------
  
  #likelihood  <- loglik - slogasc
  #sloglikfam <- sum(loglikfam)
  return(-loglikfam)
}



variance=function(i){
tvctype="PE"
data<- read.csv( paste0("simdata",i,".csv"))
estbas<-read.csv( paste0("est",i,".csv"))
# data <- read.csv("simdata1.csv")
# estbas=read.csv("est1.csv")
#estcov=read.csv("esto1.csv")
estbas=as.vector(t(estbas))
#estcov=as.vector(t(estcov))
output = estbas#c(estbas[1:4],estcov[-c(1:4)])
#st1est=output[-c(8,9)]
#st2est=output[c(8,9)]

#### UNI
size<- aggregate(data$status==1, by=list(data$famID), length)[,2]
data$df1 <- rep(aggregate(data$status==1, by=list(data$famID), sum)[,2],size)
data$df2 <- rep(aggregate(data$status==2, by=list(data$famID), sum)[,2],size)
## 

          H=hessian(loglik_variance_frailty,output,data=data, agemin=16, frailty=TRUE,tvctype=tvctype) 
     grad <- jacobian(loglik_variance_vec_frailty,output,data=data,agemin=16,frailty=TRUE,tvctype=tvctype)
   Jscore <- t(grad)%*%grad
RobustVar <- solve(H)%*%(Jscore)%*%solve(H)
 RobustSE <- sqrt(diag(RobustVar))
 
 # RobustVar <- solve(H)
 # RobustSE <- sqrt(diag(RobustVar))

# Delta method
dpenSE1=NULL
for(age in c(25,35,45,55)){
  FD1=jacobian(pen4comb,output,type=1,t=age,C=c(0,0))
  FD2=jacobian(pen4comb,output,type=1,t=age,C=c(1,0))
  FD3=jacobian(pen4comb,output,type=1,t=age,C=c(0,1))
  FD4=jacobian(pen4comb,output,type=1,t=age,C=c(1,1))
  f1=FD1%*%RobustVar%*%t(FD1) #penetrance se by age NS,NC
  f2=FD2%*%RobustVar%*%t(FD2) #penetrance se by age S,NC
  f3=FD3%*%RobustVar%*%t(FD3) #penetrance se by age NS,C
  f4=FD4%*%RobustVar%*%t(FD4) #penetrance se by age S,C
  se=sqrt(c(f1,f2,f3,f4))
  dpenSE1=c(dpenSE1,se)
}
dpenSE2=NULL
for(age in c(25,35,45,55)){
  FD1=jacobian(pen4comb,output,type=2,t=age,C=c(0,0))
  FD2=jacobian(pen4comb,output,type=2,t=age,C=c(1,0))
  FD3=jacobian(pen4comb,output,type=2,t=age,C=c(0,1))
  FD4=jacobian(pen4comb,output,type=2,t=age,C=c(1,1))
  f1=FD1%*%RobustVar%*%t(FD1) #penetrance se by age NS,NC
  f2=FD2%*%RobustVar%*%t(FD2) #penetrance se by age S,NC
  f3=FD3%*%RobustVar%*%t(FD3) #penetrance se by age NS,C
  f4=FD4%*%RobustVar%*%t(FD4) #penetrance se by age S,C
  se=sqrt(c(f1,f2,f3,f4))
  dpenSE2=c(dpenSE2,se)
}
 
 
# #Other method
# sest <- try(mvrnorm(n=1000, mu=output, Sigma=RobustVar),silent=TRUE)
# if(inherits(sest ,'try-error')){
#   warning(as.vector(sest))
#   RobustSE=rep(NA_real_, length(output))
#   penSE1=penSE2=rep(NA_real_,16)
# }else{
# type=1
# penSE1=c(
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=25),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=25),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=25),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=25),na.rm=TRUE)),
# 
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=35),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=35),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=35),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=35),na.rm=TRUE)),
# 
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=45),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=45),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=45),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=45),na.rm=TRUE)),
# 
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=55),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=55),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=55),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=55),na.rm=TRUE))
# )
# type=2
# penSE2=c(
# sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=25),na.rm=TRUE)),
# sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=25),na.rm=TRUE)),
# sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=25),na.rm=TRUE)),
# sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=25),na.rm=TRUE)),
# 
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=35),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=35),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=35),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=35),na.rm=TRUE)),
# 
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=45),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=45),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=45),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=45),na.rm=TRUE)),
# 
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=55),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=55),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=55),na.rm=TRUE)),
#   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=55),na.rm=TRUE))
# )
# }
return(list(RobustSE,dpenSE1,dpenSE2))
#return(list(RobustSE,penSE1,penSE2))
}
variance2=function(i){
  tvctype="PE"
  data<- read.csv( paste0("simdata",i,".csv"))
  estbas<-read.csv( paste0("est",i,".csv"))
  #estcov=read.csv("esto1.csv")
  estbas=as.vector(t(estbas))
  #estcov=as.vector(t(estcov))
  output = estbas#c(estbas[1:4],estcov[-c(1:4)])
  #st1est=output[-c(8,9)]
  #st2est=output[c(8,9)]
  
  #### UNI
  size<- aggregate(data$status==1, by=list(data$famID), length)[,2]
  data$df1 <- rep(aggregate(data$status==1, by=list(data$famID), sum)[,2],size)
  data$df2 <- rep(aggregate(data$status==2, by=list(data$famID), sum)[,2],size)
  ## 
  
  H=hessian(loglik_variance_frailty,output,data=data, agemin=16, frailty=TRUE,tvctype=tvctype) 
  # grad <- jacobian(loglik_variance_vec_frailty,output,data=data,agemin=16,frailty=TRUE,tvctype=tvctype)
  # Jscore <- t(grad)%*%grad
  # RobustVar <- solve(H)%*%(Jscore)%*%solve(H)
  # RobustSE <- sqrt(diag(RobustVar))
  
  RobustVar <- solve(H)
  RobustSE <- sqrt(diag(RobustVar))
  
  # Delta method
  dpenSE1=NULL
  for(age in c(25,35,45,55)){
    FD1=jacobian(pen4comb,output,type=1,t=age,C=c(0,0))
    FD2=jacobian(pen4comb,output,type=1,t=age,C=c(1,0))
    FD3=jacobian(pen4comb,output,type=1,t=age,C=c(0,1))
    FD4=jacobian(pen4comb,output,type=1,t=age,C=c(1,1))
    f1=FD1%*%RobustVar%*%t(FD1) #penetrance se by age NS,NC
    f2=FD2%*%RobustVar%*%t(FD2) #penetrance se by age S,NC
    f3=FD3%*%RobustVar%*%t(FD3) #penetrance se by age NS,C
    f4=FD4%*%RobustVar%*%t(FD4) #penetrance se by age S,C
    se=sqrt(c(f1,f2,f3,f4))
    dpenSE1=c(dpenSE1,se)
  }
  dpenSE2=NULL
  for(age in c(25,35,45,55)){
    FD1=jacobian(pen4comb,output,type=2,t=age,C=c(0,0))
    FD2=jacobian(pen4comb,output,type=2,t=age,C=c(1,0))
    FD3=jacobian(pen4comb,output,type=2,t=age,C=c(0,1))
    FD4=jacobian(pen4comb,output,type=2,t=age,C=c(1,1))
    f1=FD1%*%RobustVar%*%t(FD1) #penetrance se by age NS,NC
    f2=FD2%*%RobustVar%*%t(FD2) #penetrance se by age S,NC
    f3=FD3%*%RobustVar%*%t(FD3) #penetrance se by age NS,C
    f4=FD4%*%RobustVar%*%t(FD4) #penetrance se by age S,C
    se=sqrt(c(f1,f2,f3,f4))
    dpenSE2=c(dpenSE2,se)
  }
  
  
  # #Other method
  # sest <- try(mvrnorm(n=1000, mu=output, Sigma=RobustVar),silent=TRUE)
  # if(inherits(sest ,'try-error')){
  #   warning(as.vector(sest))
  #   RobustSE=rep(NA_real_, length(output))
  #   penSE1=penSE2=rep(NA_real_,16)
  # }else{
  # type=1
  # penSE1=c(
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=25),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=25),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=25),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=25),na.rm=TRUE)),
  # 
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=35),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=35),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=35),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=35),na.rm=TRUE)),
  # 
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=45),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=45),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=45),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=45),na.rm=TRUE)),
  # 
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=55),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=55),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=55),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=55),na.rm=TRUE))
  # )
  # type=2
  # penSE2=c(
  # sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=25),na.rm=TRUE)),
  # sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=25),na.rm=TRUE)),
  # sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=25),na.rm=TRUE)),
  # sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=25),na.rm=TRUE)),
  # 
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=35),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=35),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=35),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=35),na.rm=TRUE)),
  # 
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=45),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=45),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=45),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=45),na.rm=TRUE)),
  # 
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,0),t=55),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,0),t=55),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(0,1),t=55),na.rm=TRUE)),
  #   sqrt(var(apply(sest, 1 , pen4comb,type=type, C=c(1,1),t=55),na.rm=TRUE))
  # )
  # }
  return(list(RobustSE,dpenSE1,dpenSE2))
  #return(list(RobustSE,penSE1,penSE2))
}


# Penetrance computation functions
pen4comb=function(est,C,t,type){
  
  #print(type)
  est[c(1:4,8:9)] <- sapply(est[c(1:4,8:9)],  exp)
  
  pencomp=function(){
    if(C[1]==0){
      pen=penf(integrand,X=C,est=est,frailty=TRUE,type=type,t=t)
    }else if(C[1]==1){
      pen=penf(integrand,X=C,est=est,frailty=TRUE,type=type,t=20)+
        penf.tsc(integrand.tsc,X=C,est=est,frailty=TRUE,type=type,t=t)
    }
    pen
  }
  
  pen=try(pencomp(), silent=TRUE)
  
  if(inherits(pen ,'try-error')){
    warning(as.vector(pen))
    pen <-NA_real_
  } else {
    pen <- pen
  }
  pen
}
penf <- function(f, X, est, frailty,type,  t){  
  integrate(integrand, X=X, est=est, frailty=frailty,type=type,lower=0,upper=t)$value
}
integrand <- function(u, X, est, type, frailty){
  u <- u
  fp=c(est[8],est[9])
  index=X
  xbeta1 =  index[2]*est[6]
  xbeta2 =  index[2]*est[7]
  if(type==1){
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
  }else{
  haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)  
  }
  H1 = - (est[1]*u)^est[2]*exp(xbeta1)
  H2 = - (est[3]*u)^est[4]*exp(xbeta2)
  
  if(frailty==TRUE){
    if(type==1){
      f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
    }else{
      f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))  
    }
  } else{
    f = haz1*exp(H1+H2) # h * S
  }
  return(f)
}
penf.tsc <- function(f, X, est, frailty,type,t){  
  integrate(integrand.tsc, X=X, est=est, frailty=frailty,type=type,lower=20,upper=t)$value
}
integrand.tsc <- function(u, X, est, type, frailty){
  u <- u
  fp=c(est[8],est[9])
  index=X
  xbeta1 =  index[2]*est[6]
  xbeta2 =  index[2]*est[7]
  if(type==1){
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)*exp(est[5])
  }else{
  haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)  
  }
  H1 = - (sapply(u,cH_PE,a=20,parm=c(est[1:2],est[5]))+(est[1]*20)^est[2])*exp(xbeta1)
  H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
  
  if(frailty==TRUE){
    if(type==1){
    f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
    }else{
    f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))  
    }
  } else{
    f = haz1*exp(H1+H2) # h * S
  }
  return(f)
}
##  if we integrate over one TDC
## under Permanent Exposure (PE) TVC
cH_PE<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  
  bH0 <- function(u){  # Cumulative hazard function
    (lamb*u)^rho
  }
  
  
  return( exp(bz)*(bH0(b) - bH0 (a) ) )
  
}

# combine all these function in one for subject j
cHj<-function(x, parm,TDfunc){
  a=x[1]; b=x[2]
  if(TDfunc=="PE") res<-cH_PE(a, b, parm)
  else if(TDfunc=="ED") res<-cH_ED(a, b, parm)
  else if(TDfunc=="CO") res<-cH_CO(a, b, parm)
  else stop("TDfunc should be PE or ED or CO")
  res
}
#for several subjects :  use lapply function
cH<-function(A, B, parm, TDfunc){
  unlist( lapply(1:length(A), function(i) cHj(x=c(A[i],B[i]), parm=parm, TDfunc=TDfunc) ) )
}

numericGradient <- function(f, t0, eps=1e-6, ...) {
  ### numeric gradient of a vector-valued function
  ### f    : function, return Nval x 1 vector of values
  ### t0   : NPar x 1 vector of parameters
  ### return:
  ### NvalxNPar matrix, gradient
  NPar <- length(t0)
  NVal <- length(f0 <- f(t0, ...))
  grad <- matrix(NA, NVal, NPar)
  row.names(grad) <- names(f0)
  colnames(grad) <- names(t0)
  for(i in 1:NPar) {
    t2 <- t1 <- t0
    t1[i] <- t0[i] - eps/2
    t2[i] <- t0[i] + eps/2
    grad[,i] <- (f(t2, ...) - f(t1, ...))/eps
  }
  return(grad)
}

