# Example of computing penetrance SE by age 70 for noncarrier, nonscreened
# Empirical SE 0.0096
# run everything after line 6 first
variance()

# compute the gradient inside the integral
variance=function(){
tvctype="PE"
data <- read.csv("simdata1.csv")
estbas=read.csv("est1.csv")
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

# Var(theta)
RobustVar <- solve(H)

# gradient(theta)
gra=vector()
for(j in 1:length(output)){
  #compute the partial derivate of the theta[j] separately
  gra[j]=pen4comb(output,type=1,t=55,C=c(0,0),j)
}

#delta method 
penetranceSE=t(gra)%*%RobustVar%*%cbind(gra)
return(penetranceSE)
}

# Penetrance computation functions
library(pracma)
library(MASS)
pen4comb=function(est,C,t,type,j){
  #print(type)
  est[c(1:4,8:9)] <- sapply(est[c(1:4,8:9)],  exp)
  pen=penf(integrand,X=C,est=est,frailty=TRUE,type=type,t=t,j=j)
  pen
}
#compute outer integral with the function "penf"
penf <- function(f, X, est, frailty,type,  t,j){  
  integrate(integrand, X=X, est=est, frailty=frailty,type=type,j=j,lower=0,upper=t)$value
}
# integrate function needs to be vectorized function of u
# integrand involves computing partial derivative of theta[j]
# integrand returns   dF/dtheta[j]|u=u
integrand <- function(u, X, est, type, frailty,j){
  u <- u
  #u= c(10, 0.260934714828283, 19.7390652851717, 1.34936633311015, 
  #     18.6506336668898, 3.20590431700976, 16.7940956829902, 5.66604605870753, 
  #     14.3339539412925, 8.51125661018369, 11.4887433898163)
  index=X
  func=function(est,index,type,frailty,u){
    fp=c(est[8],est[9])
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
        density = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
      }else{
        density = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))  
      }
    } else{
      density = haz1*exp(H1+H2) # h * S
    }
    return(density)
  }
  # below returns dfunc/dtheta[j]|u=u
  out=jacobian(f=func,x0=est,index=index,type=type,frailty=frailty,u=u)[,j]
  return(out)
}
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

