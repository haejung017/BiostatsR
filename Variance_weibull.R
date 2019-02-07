
s1.Y=function(h.Y,theta,h.eps,data2) {grad(NegLogLik.corr.wt,theta, data=data2, h.Y=h.Y, h.eps=h.eps)}
s2.Y=function(theta,h.Y,data1) {grad(NegLogLik.h.Y.wt,h.Y,theta=theta,data=data1)}
s1.eps=function(h.eps,theta,h.Y,data2) {grad(NegLogLik.corr.wt,theta, data=data2, h.Y=h.Y, h.eps=h.eps)}
s2.eps=function(theta,h.eps,data3) {grad(NegLkogLik.h.eps.wt,h.eps,theta=theta,data=data3)}



compute.second.derivatives = function(theta,h.Y,h.eps,data1,data2,data3)
{
  n=length(theta)
  result = matrix(0,ncol=(n+2),nrow=(n+2))
  result[1:n,1:n]=hessian(NegLogLik.corr.wt,theta, data=data2, h.Y=h.Y, h.eps=h.eps)
  result[1:n,(n+1)]=jacobian(s1.Y,h.Y,theta=theta,h.eps=h.eps,data2=data2)
  result[(n+1),1:n]=grad(s2.Y,theta,h.Y=h.Y,data1=data1)
  result[1:n,(n+2)]=jacobian(s1.eps,h.eps,theta=theta,h.Y=h.Y,data2=data2)
  result[(n+2),1:n]=t(result[1:n,(n+2)])
  result[(n+1),(n+1)]=hessian(NegLogLik.h.Y.wt,h.Y,theta=theta,data=data1)
  result[(n+2),(n+2)]=hessian(NegLogLik.h.eps.wt,h.eps,theta=theta,data=data3)
  result
}


NegLogLik.corr.wt.vec = function(theta,data,h.Y,h.eps,famid)
{
  
  Y1 = data$Y1
  Y2 = data$Y2
  
  eps1 = data$eps1
  eps2 = data$eps2
  
  G1 = data$G1
  G2 = data$G2
  
  kinship = data$kinship
  npair = data$npair
  
  al1=exp(theta[1])
  lam1=exp(theta[2])
  al2=exp(theta[3])
  lam2=exp(theta[4])
  beta1=theta[5]
  beta2=theta[6]
  
  S1 = exp(-lam1*Y1^al1*exp(beta1*G1)-lam2*Y1^al2*exp(beta2*G1))
  S2 = exp(-lam1*Y2^al1*exp(beta1*G2)-lam2*Y2^al2*exp(beta2*G2))
  
  
  f11 = (lam1*al1*exp(beta1*G1)*Y1^(al1-1))*S1 
  f21 = (lam2*al2*exp(beta2*G1)*Y1^(al2-1))*S1 
  f12 = (lam1*al1*exp(beta1*G2)*Y2^(al1-1))*S2 
  f22 = (lam2*al2*exp(beta2*G2)*Y2^(al2-1))*S2 
  
  f1 = f11+f21
  f2 = f12+f22
  
  P1 = f11/f1
  P2 = f12/f2
  
  P11 = BiGauCDF(u1=P1,u2=P2,par=h.eps*kinship)
  P12 = P1 - P11 
  P21 = P2 - P11 
  P22 = 1 - (P11+P12+P21)
  
  corr.term = log(1-BiGauCDF(S1, S2, par=h.Y*kinship))
  
  cond11 = (eps1==1) & (eps2==1)
  cond12 = (eps1==1) & (eps2==2)
  cond21 = (eps1==2) & (eps2==1)
  cond22 = (eps1==2) & (eps2==2)
  
  cond10 = (eps1==1) & (eps2==0)
  cond20 = (eps1==2) & (eps2==0)
  cond01 = (eps1==0) & (eps2==1)
  cond02 = (eps1==0) & (eps2==2)
  
  loglik = 0
  
  cond1 = (eps1>0) & (eps2>0)
  cond2 = (eps1>0) & (eps2==0)
  cond3 = (eps1==0) & (eps2>0)
  
  loglik[cond1] = log(BiCopPDF(u1=S1[cond1],u2=S2[cond1],family=1,par=h.Y*kinship[cond1])) + log(f1[cond1]) + log(f2[cond1]) 
  loglik[cond2] = log(BiCopHfunc1(u1=S1[cond2],u2=S2[cond2],family=1,par=h.Y*kinship[cond2])) 
  loglik[cond3] = log(BiCopHfunc2(u1=S1[cond3],u2=S2[cond3],family=1,par=h.Y*kinship[cond3]))
  
  loglik[cond11] = loglik[cond11] + log(P11[cond11])
  loglik[cond12] = loglik[cond12] + log(P12[cond12])
  loglik[cond21] = loglik[cond21] + log(P21[cond21])
  loglik[cond22] = loglik[cond22] + log(P22[cond22])
  loglik[cond10] = loglik[cond10] + log(f11[cond10])
  loglik[cond20] = loglik[cond20] + log(f21[cond20])
  loglik[cond01] = loglik[cond01] + log(f12[cond01])
  loglik[cond02] = loglik[cond02] + log(f22[cond02])
  ret<- -(loglik-corr.term)/npair
  
  out <- rep(0,length(famid))
  out[famid%in%data$fam.id] <-aggregate(ret, list(data$fam.id),sum)[,2]  
  out
  
}

NegLogLik.h.Y.wt.vec = function(h,theta,data,famid)
{
  Y1 = data$Y1
  Y2 = data$Y2
  
  eps1 = data$eps1
  eps2 = data$eps2
  
  G1 = data$G1
  G2 = data$G2
  
  kinship = data$kinship
  npair = data$npair
  
  al1=exp(theta[1])
  lam1=exp(theta[2])
  al2=exp(theta[3])
  lam2=exp(theta[4])
  beta1=theta[5]
  beta2=theta[6]
  
  S1 = exp(-lam1*Y1^al1*exp(beta1*G1)-lam2*Y1^al2*exp(beta2*G1))
  S2 = exp(-lam1*Y2^al1*exp(beta1*G2)-lam2*Y2^al2*exp(beta2*G2))
  
  cond00 = (eps1==0) & (eps2==0)
  cond11 = (eps1!=0) & (eps2!=0)
  
  cond10 = (eps1!=0) & (eps2==0)
  cond01 = (eps1==0) & (eps2!=0)
  
  loglik = 0 
  loglik[cond00] = log(BiGauCDF(u1=S1[cond00],u2=S2[cond00],par=h*kinship[cond00]))
  
  loglik[cond11] = log(BiCopPDF(u1=S1[cond11],u2=S2[cond11],family=1,par=h*kinship[cond11])) 
  
  loglik[cond10] = log(BiCopHfunc1(u1=S1[cond10],u2=S2[cond10],family=1,par=h*kinship[cond10])) 
  
  loglik[cond01] = log(BiCopHfunc2(u1=S1[cond01],u2=S2[cond01],family=1,par=h*kinship[cond01])) 
  
  ret<- -loglik/npair
  out <- rep(0,length(famid))
  out[famid%in%data$fam.id] <-aggregate(ret, list(data$fam.id),sum)[,2]  
  out
  
}

NegLogLik.h.eps.wt.vec = function(h,theta,data, famid)
{
  Y1 = data$Y1
  Y2 = data$Y2
  
  eps1 = data$eps1
  eps2 = data$eps2
  
  G1 = data$G1
  G2 = data$G2
  
  kinship = data$kinship
  npair = data$npair
  
  al1=exp(theta[1])
  lam1=exp(theta[2])
  al2=exp(theta[3])
  lam2=exp(theta[4])
  beta1=theta[5]
  beta2=theta[6]
  
  haz11 = lam1*al1*Y1^(al1-1)*exp(beta1*G1)
  haz21 = lam2*al2*Y1^(al2-1)*exp(beta2*G1)
  haz12 = lam1*al1*Y2^(al1-1)*exp(beta1*G2)
  haz22 = lam2*al2*Y2^(al2-1)*exp(beta2*G2)
  
  P1 = haz11 / (haz11+haz21)  # P(eps1=1)
  P2 = haz12 / (haz12+haz22)  # P(eps2=1)
  
  P11 = BiGauCDF(u1=P1,u2=P2,par=h*kinship) # P(eps1=1,eps2=1)
  P12 = P1 - P11 # P(eps1=1,eps2=2)
  P21 = P2 - P11 # P(eps1=2,eps2=1)
  P22 = 1 - (P11+P12+P21)
  
  cond11 = (eps1==1) & (eps2==1)
  cond12 = (eps1==1) & (eps2==2)
  
  cond21 = (eps1==2) & (eps2==1)
  cond22 = (eps1==2) & (eps2==2)
  
  loglik = 0
  loglik[cond11] = log(P11[cond11])
  loglik[cond12] = log(P12[cond12])
  loglik[cond21] = log(P21[cond21])
  loglik[cond22] = log(P22[cond22])
  
  ret <- -loglik/npair
  out <- rep(0,length(famid))
  out[famid%in%data$fam.id] <-aggregate(ret, list(data$fam.id),sum)[,2]  
  out
  
}
variance <- function(theta, h.Y, h.eps, data1, data2, data3){
  famid<-unique(data1$fam.id)
  ee <- compute.second.derivatives(theta, h.Y, h.eps, data1, data2, data3)
  D=jacobian(NegLogLik.corr.wt.vec,theta, data=data2, h.Y=h.Y, h.eps=h.eps, famid=famid)
  D=cbind(D,jacobian(NegLogLik.h.Y.wt.vec, h.Y, theta=theta, data=data1,famid=famid))
  D=cbind(D,jacobian(NegLogLik.h.eps.wt.vec, h.eps, theta=theta, data=data3,famid=famid))
  B=t(D)%*%D
  A=solve(ee)
  varcov=t(A)%*%B%*%A
  return(list(cov=varcov, se=sqrt(diag(varcov))))
}