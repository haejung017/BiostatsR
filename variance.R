# library(mnormt)
# #library(survival)
# library(pbivnorm)
library(pracma)

#example
# data=readRDS("simdata") 
# st1est=readRDS("st1est")#stage 1 est, baseline+covariate
# st2est=readRDS("st2est")#stage 2 est, two frailty parameters

#source("loglik_variance.R")#stage 2 est, two frailty parameters
#source("TVCtools.R")#stage 2 est, two frailty parameters
variance=function(){
# dataR=readRDS("simdata_nof.rds")
# output_final4=readRDS("output_simdata_nof.rds")
# output_final4[i,-c(8,9)]
# i=1
#setotal=NULL
tvctype="CO"

#data=dataR[[i]]
data <- read.csv("simdata.csv")
output<-read.csv("est.csv")
output=as.vector(t(output))
st1est=output[-c(8,9)]
st2est=output[c(8,9)]

#### UNI
size<- aggregate(data$status==1, by=list(data$famID), length)[,2]
data$df1 <- rep(aggregate(data$status==1, by=list(data$famID), sum)[,2],size)
data$df2 <- rep(aggregate(data$status==2, by=list(data$famID), sum)[,2],size)


#If frailty is FALSE, it is independent model likelihood
#If frailty is TRUE, it is frailty model likelihood

#loglik_variance(stage1=st1est,h.Y=st2est, data=data, agemin=16, frailty=FALSE,tvctype="PE")
#loglik_variance(stage1=st1est,h.Y=st2est, data=data, agemin=16, frailty=TRUE,tvctype="PE")
## 

# Compute Matrix A -----------------------------------------
n=length(st1est)
n2=length(st2est)
nt=n+n2

Amat = Bmat = matrix(0,ncol=n+n2,nrow=n+n2) 
Amat[1:n,1:n]=hessian(loglik_variance,st1est,h.Y=st2est, #upper left
                      data=data, agemin=16, frailty=FALSE,tvctype=tvctype) 
Amat[(n+1):nt,(n+1):nt]=hessian(loglik_variance,st2est,stage1=st1est, #lower right
                                data=data, agemin=16, frailty=TRUE,tvctype=tvctype)
#Amat[(n+1):nt,1:n]=cbind(dd)%*%d #wrong

s2.Y=function(stage1,h.Y,data){jacobian(loglik_variance,h.Y,stage1=stage1,data=data, frailty=TRUE,agemin=16,tvctype=tvctype)} 
#s2.Y(stage1, h.Y, data)

##
### UNI: grad -> jacobian
  aa<-jacobian(s2.Y,st1est,h.Y=st2est, data=data) #lower left
  Amat[(n+1):nt,1:n] = aa
  #Amat[1:n, (n+1):nt] = t(aa)
  

# Compute Matrix B -----------------------------------------
L1T1=grad(loglik_variance,st1est,h.Y=st2est,data=data, agemin=16, frailty=FALSE,tvctype=tvctype)
L2T2=grad(loglik_variance,st2est,stage1=st1est,data=data, agemin=16, frailty=TRUE,tvctype=tvctype)
B1=cbind(L1T1)%*%L1T1;B2=cbind(L1T1)%*%L2T2
B3=cbind(L2T2)%*%L1T1;B4=cbind(L2T2)%*%L2T2
Bmat[1:n,1:n]=B1
Bmat[(n+1):nt,(n+1):nt]=B4
Bmat[(n+1):nt,1:n]=B3
Bmat[1:n,(n+1):nt]=B2
#-------------------------------------------------------------
A=solve(Amat)
varcov=A%*%Bmat%*%t(A)
se=sqrt(diag(varcov))
#setotal=rbind(setotal,se)
#print(i)
#print(colMeans(setotal))
return(se)
}

#[1] 0.036167440 0.021828297 0.005614409 0.039668520 1.355541895 0.183674173
#[7] 0.240956818 0.358769970 0.557147405 1.046353303 0.372551383

#saveRDS(sqrt(diag(varcov)),"v1")

# d=grad(loglik_variance,st1est,h.Y=st2est,data=data, agemin=16, frailty=TRUE,tvctype="CO")
# dd=grad(loglik_variance,st2est,stage1=st1est,data=data, agemin=16, frailty=TRUE,tvctype="CO")
# 
# ttm = output_final8[,9]
# boxplot(ttm)
# sd(ttm)
# mean(ttm)
# median(ttm)
# abline(h=1.0647)

loglik_variance=function(h.Y,stage1,data, agemin, frailty=TRUE,tvctype){
  
  data = data[data$currentage>=agemin,]
  theta = h.Y
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
  print(round(dest,4))
  
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

######## integrale term in the cum hazard function 

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

## under Exponential Decay (ED) TVC
cH_ED<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  phi<-parm[4]
  
  bh0<-function(u){
    (lamb*rho)*(lamb*u)^(rho-1)
  }
  
  integrand <- function(u) { bh0(u)*exp(bz*exp(-(u-a)*phi)) }
  int <- try( integrate( integrand, lower=a, upper=b), silent = TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <-NA_real_
  } else {
    integrated <-int$value
  }
  
  return(integrated)
  
}


## under Cox and Oaks (CO) TVC
cH_CO<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  phi<-parm[4]
  a0<-parm[5]
  
  bh0<-function(u){
    (lamb*rho)*(lamb*u)^(rho-1)
  }
  
  integrand <- function(u) { bh0( u)*exp(a0+bz*exp(-(u-a)*phi)) }
  int <- try( integrate( integrand, lower=a, upper=b), silent = TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <-NA_real_
  } else {
    integrated <-int$value
  }
  
  return(integrated)
  
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

