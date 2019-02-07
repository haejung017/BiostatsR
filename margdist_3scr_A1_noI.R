

###########################
########################## log-hazards under PH model with two time-varying covariates
##########################

loghaz<- function(mpars, Y, Yor, Yst1, Yst2, Yst3,  Ybr, TDf_or, TDf_sc, cp){
  
  lamb<-mpars[1]
  rho<-mpars[2]
  beta.g<-mpars[3]
  beta.or<-mpars[4]
  #beta.gor<-mpars[4]
  beta.sc<-mpars[5:7]
  beta.gsc<-c(0,0,0)
  
  a00=0
  a01=0
  a02=0
  a03=0
  phi0=1e-9
  phi1=1e-9
  phi2=1e-9
  phi3=1e-9
  
  if(TDf_or=="PE" && TDf_sc=="ED" ){
    phi1=mpars[length(mpars)-2]
    phi2=mpars[length(mpars)-1]
    phi3=mpars[length(mpars)]
  } else  if(TDf_or=="PE" && TDf_sc=="CO" ){
    
    phi1=mpars[length(mpars)-5]
    a01=mpars[length(mpars)-4]
    phi2=mpars[length(mpars)-3]
    a02=mpars[length(mpars)-2]
    phi3=mpars[length(mpars)-1]
    a03=mpars[length(mpars)]
    
  } else  if(TDf_or=="ED" && TDf_sc=="PE" ){
    phi0=mpars[length(mpars)]
  } else  if(TDf_or=="ED" && TDf_sc=="ED" ){
    phi0=mpars[length(mpars)-3]
    phi1=mpars[length(mpars)-2]
    phi2=mpars[length(mpars)-1]
    phi3=mpars[length(mpars)]
    
  } else  if(TDf_or=="ED" && TDf_sc=="CO" ){
    phi0=mpars[length(mpars)-6]
    phi1=mpars[length(mpars)-5]
    a01=mpars[length(mpars)-4]
    phi2=mpars[length(mpars)-3]
    a02=mpars[length(mpars)-2]
    phi3=mpars[length(mpars)-1]
    a03=mpars[length(mpars)]
    
  } else  if(TDf_or=="CO" && TDf_sc=="PE" ){
    phi0=mpars[length(mpars)-1]
    a00=mpars[length(mpars)]
  } else  if(TDf_or=="CO" && TDf_sc=="ED" ){
    phi0=mpars[length(mpars)-4]
    a00=mpars[length(mpars)-3]
    phi1=mpars[length(mpars)-2]
    phi2=mpars[length(mpars)-1]
    phi3=mpars[length(mpars)]
  } else  if(TDf_or=="CO" && TDf_sc=="CO" ){
    
    phi0=mpars[length(mpars)-7]
    a00=mpars[length(mpars)-6]
    phi1=mpars[length(mpars)-5]
    a01=mpars[length(mpars)-4]
    phi2=mpars[length(mpars)-3]
    a02=mpars[length(mpars)-2]
    phi3=mpars[length(mpars)-1]
    a03=mpars[length(mpars)]
    
  }
  
  beta.gor=0
  beta.orsc<-c(0,0,0)
  
  BRI <- ifelse(Y-Ybr> 0, 1, 0)
  Y[BRI==1]=Ybr[BRI==1]
  
  # Oopharectomy PE
  
  ORI <- ifelse(Y-Yor> 0, 1, 0)
  Zor<-ifelse(ORI==0,0,1)
  if(TDf_or=="ED" || TDf_or=="CO")
    Zor<-ifelse(ORI==0,0,exp(-phi0*(Y-Yor)))
  
  # screening ED or PE
  SCI1 <- ifelse(!is.na(Yst1), 1, 0)
  
  SCI2 <- ifelse(!is.na(Yst2), 1, 0)
  
  SCI3 <- ifelse(!is.na(Yst3), 1, 0)
  
  
  
  # screening ED or PE
  SCI1 <- ifelse(!is.na(Yst1), 1, 0)
  Zst1<-ifelse(SCI1==0,0,1)
  Zst2<-ifelse(SCI2==0,0,1)
  Zst3<-ifelse(SCI3==0,0,1)
  if(TDf_sc=="ED" || TDf_sc=="CO"){
    Zst1<-ifelse(SCI1==0,0,exp(-phi1*(Y-Yst1)))
    Zst2<-ifelse(SCI2==0,0,exp(-phi2*(Y-Yst2)))
    Zst3<-ifelse(SCI3==0,0,exp(-phi3*(Y-Yst3)))
  }
  nbSC<-apply(cbind(SCI1, SCI2, SCI3),1, sum)
  
  h0<-function(u){
    (lamb*rho)*(lamb*u)^(rho-1)
  }
  
  
  # The log of the hazard 
  loghg0=vector(length=length(ORI))
  loghg1=vector(length=length(ORI))
  
  
  
  ###################################  No Oopharectomy ORI==0
  
  
  ########### nbSC=0 & ORI=0 & BRI=0 
  id000=which(nbSC==0 & ORI==0)
  loghg1[id000]=log(h0(Y[id000])) + beta.g
  loghg0[id000]=log(h0(Y[id000]))
  
  ########### nbSC=1 & ORI=0 & BRI=0 
  id100=which(nbSC==1 & ORI==0)
  loghg1[id100]=log(h0(Y[id100])) + beta.g + (beta.sc[1] + beta.gsc[1])*Zst1[id100] + a01
  loghg0[id100]=log(h0(Y[id100])) +  (beta.sc[1])*Zst1[id100] + a01 
  
  ########### nbSC=2 & ORI=0 & BRI=0 
  id200=which(nbSC==2 & ORI==0)
  loghg1[id200]=log(h0(Y[id200])) + beta.g + (beta.sc[2] + beta.gsc[2])*Zst2[id200] + a02
  loghg0[id200]=log(h0(Y[id200])) +  (beta.sc[2])*Zst2[id200] + a02
  
  ########### nbSC=3 & ORI=0 & BRI=0 
  id300=which(nbSC==3 & ORI==0) 
  loghg1[id300]=log(h0(Y[id300])) + beta.g + (beta.sc[3] + beta.gsc[3])*Zst3[id300] + a03
  loghg0[id300]=log(h0(Y[id300])) + (beta.sc[3])*Zst3[id300] + a03
  
  
  ###################################  Oopharectomy ORI==1
  
  ########### nbSC=0 & ORI=1 & BRI=0  
  id010=which(nbSC==0 & ORI==1)
  loghg1[id010]=log(h0(Y[id010])) + beta.g + (beta.or+beta.gor)*Zor[id010] + a00
  loghg0[id010]=log(h0(Y[id010])) + (beta.or)*Zor[id010] + a00
  
  ########### nbSC=1 & ORI=1 & BRI=0 
  id110=which(nbSC==1 & ORI==1)
  loghg1[id110]=log(h0(Y[id110])) + beta.g + (beta.or+beta.gor)*Zor[id110] + (beta.sc[1] + beta.gsc[1])*Zst1[id110] + a00 + a01
  loghg0[id110]=log(h0(Y[id110])) + (beta.or)*Zor[id110] + (beta.sc[1])*Zst1[id110]+ a00 + a01
  
  
  ########### nbSC=2 & ORI=1 & BRI=0 
  id210=which(nbSC==2 & ORI==1)
  loghg1[id210]=log(h0(Y[id210])) + beta.g + (beta.or+beta.gor)*Zor[id210] + (beta.sc[2] + beta.gsc[2])*Zst2[id210] + a00 + a02
  loghg0[id210]=log(h0(Y[id210])) + (beta.or)*Zor[id210] + (beta.sc[2])*Zst2[id210] + a00 + a02
  
  
  
  ########### nbSC=3 & ORI=1 & BRI=0 
  id310=which(nbSC==3 & ORI==1)
  loghg1[id310]=log(h0(Y[id310])) + beta.g + (beta.or+beta.gor)*Zor[id310] + (beta.sc[3] + beta.gsc[3])*Zst3[id310] + a00 + a03
  loghg0[id310]=log(h0(Y[id310])) + (beta.or)*Zor[id310] + (beta.sc[3])*Zst3[id310] + a00 + a03
  
  logh=(1-cp)*loghg0 + cp*loghg1
  
  
  return( list(logh=logh) )
}



###########################
########################## cumulative hazard under PH model with two time-varying covariates
##########################


cumHaz<- function(mpars, Y, Yor, Yst1, Yst2, Yst3,  Ybr, TDf_or, TDf_sc, cp){
  
  lamb<-mpars[1]
  rho<-mpars[2]
  beta.g<-mpars[3]
  beta.or<-mpars[4]
  #beta.gor<-mpars[4]
  beta.sc<-mpars[5:7]
  beta.gsc<-c(0,0,0)
  
  a00=0
  a01=0
  a02=0
  a03=0
  phi0=1e-9
  phi1=1e-9
  phi2=1e-9
  phi3=1e-9
  
  if(TDf_or=="PE" && TDf_sc=="ED" ){
    phi1=mpars[length(mpars)-2]
    phi2=mpars[length(mpars)-1]
    phi3=mpars[length(mpars)]
  } else  if(TDf_or=="PE" && TDf_sc=="CO" ){
    
    phi1=mpars[length(mpars)-5]
    a01=mpars[length(mpars)-4]
    phi2=mpars[length(mpars)-3]
    a02=mpars[length(mpars)-2]
    phi3=mpars[length(mpars)-1]
    a03=mpars[length(mpars)]
    
  } else  if(TDf_or=="ED" && TDf_sc=="PE" ){
    phi0=mpars[length(mpars)]
  } else  if(TDf_or=="ED" && TDf_sc=="ED" ){
    phi0=mpars[length(mpars)-3]
    phi1=mpars[length(mpars)-2]
    phi2=mpars[length(mpars)-1]
    phi3=mpars[length(mpars)]
    
  } else  if(TDf_or=="ED" && TDf_sc=="CO" ){
    phi0=mpars[length(mpars)-6]
    phi1=mpars[length(mpars)-5]
    a01=mpars[length(mpars)-4]
    phi2=mpars[length(mpars)-3]
    a02=mpars[length(mpars)-2]
    phi3=mpars[length(mpars)-1]
    a03=mpars[length(mpars)]
    
  } else  if(TDf_or=="CO" && TDf_sc=="PE" ){
    phi0=mpars[length(mpars)-1]
    a00=mpars[length(mpars)]
  } else  if(TDf_or=="CO" && TDf_sc=="ED" ){
    phi0=mpars[length(mpars)-4]
    a00=mpars[length(mpars)-3]
    phi1=mpars[length(mpars)-2]
    phi2=mpars[length(mpars)-1]
    phi3=mpars[length(mpars)]
  } else  if(TDf_or=="CO" && TDf_sc=="CO" ){
    
    phi0=mpars[length(mpars)-7]
    a00=mpars[length(mpars)-6]
    phi1=mpars[length(mpars)-5]
    a01=mpars[length(mpars)-4]
    phi2=mpars[length(mpars)-3]
    a02=mpars[length(mpars)-2]
    phi3=mpars[length(mpars)-1]
    a03=mpars[length(mpars)]
    
  }
  
  beta.gor=0
  beta.orsc<-c(0,0,0)
  
  BRI <- ifelse(Y-Ybr> 0, 1, 0)
  
  # Oopharectomy PE
  
  ORI <- ifelse(Y-Yor> 0, 1, 0)
  
  
  # screening ED or PE
  SCI1 <- ifelse(!is.na(Yst1), 1, 0)
  SCI2 <- ifelse(!is.na(Yst2), 1, 0)
  SCI3 <- ifelse(!is.na(Yst3), 1, 0)
  
  nbSC<-apply(cbind(SCI1, SCI2, SCI3),1, sum)
  
  
  
  
  h0<-function(u){
    (lamb*rho)*(lamb*u)^(rho-1)
  }
  
  
  H0 <- function(u){  # Cumulative hazard function
    (lamb*u)^rho
  }
  
  # The cum hazard
  Int0=function(A,B,bz)  Intfunc(A=A, B=B, parm=c(lamb, rho, bz, phi0, a00), TDfunc=TDf_or)
  
  Int1=function(A,B,bz) Intfunc(A=A, B=B, parm=c(lamb, rho, bz, phi1, a01), TDfunc=TDf_sc)
  Int2=function(A,B,bz)  Intfunc(A=A, B=B,  parm=c(lamb, rho, bz, phi2, a02), TDfunc=TDf_sc)
  Int3=function(A,B,bz)  Intfunc(A=A, B=B,  parm=c(lamb, rho, bz, phi3, a03), TDfunc=TDf_sc)
  
  # integrate from Yor to B ( Ysc or Y)
  
  # Int01=function(A,B, Yor, Yst1, bz0, bz1, bz01) Vectorize( apply(as.matrix(1:length(A)),1, function(i) 
  #   cH2_a_b(A[i], B[i], tor=Yor[i], tsc=Yst1[i],  lamb=lamb, rho=rho, bz0=bz0, bz1=bz1, phi0=phi0,  phi1=phi1, TDf_or=TDf_or, TDf_sc=TDf_sc) ) )
  # Int02=function(A,B,Yor, Yst2, bz0, bz1, bz01) Vectorize( apply(as.matrix(1:length(A)),1, function(i) 
  #   cH2_a_b(A[i], B[i], tor=Yor[i], tsc=Yst2[i],  lamb=lamb, rho=rho, bz0=bz0, bz1=bz1, phi0=phi0,  phi1=phi2, TDf_or=TDf_or, TDf_sc=TDf_sc) ) )
  # Int03=function(A,B,Yor, Yst3, bz0, bz1, bz01) Vectorize( apply(as.matrix(1:length(A)),1, function(i) 
  #   cH2_a_b(A[i], B[i], tor=Yor[i], tsc=Yst3[i],  lamb=lamb, rho=rho, bz0=bz0, bz1=bz1, phi0=phi0,  phi1=phi3, TDf_or=TDf_or, TDf_sc=TDf_sc) ) )
  
  # integral for subjects with oopharectomy and screening
  Int01=function(A,B, Yor, Yst1, bz0, bz1, bz01=0) 
    Int_or1sc(A=A, B=B, Tor=Yor,  Tsc=Yst1, parm=c(lamb, rho), bz0=bz0, phi0=phi0, a00=a00,   bz=bz1,  phi=phi1, a0=a01, TDf_or=TDf_or, TDf_sc=TDf_sc)
  
  Int02=function(A,B, Yor, Yst2, bz0, bz2, bz02=0) 
    Int_or1sc(A=A, B=B, Tor=Yor,  Tsc=Yst2, parm=c(lamb, rho), bz0=bz0, phi0=phi0, a00=a00,   bz=bz2,  phi=phi2, a0=a02, TDf_or=TDf_or, TDf_sc=TDf_sc)
  
  Int03=function(A,B, Yor, Yst3, bz0, bz3, bz03=0) 
    Int_or1sc(A=A, B=B, Tor=Yor,  Tsc=Yst3, parm=c(lamb, rho), bz0=bz0, phi0=phi0, a00=a00,   bz=bz3,  phi=phi3, a0=a03, TDf_or=TDf_or, TDf_sc=TDf_sc)
  
  
  
  ## Now we calculate the cum hazard function for all subjects
  
  H=vector(length=length(ORI))
  Hg0=vector(length=length(ORI))
  Hg1=vector(length=length(ORI))
  
  xbg0=0
  xbg1=beta.g
  ######################################  
  ######################################    Subjects without Mastectomy BRI==0
  ###################################### 
  
  ############ 1)   No Oopharectomy ORI==0
  
  
  ########### nbSC=0 & ORI=0 & BRI=0 
  id000=which(nbSC==0 & ORI==0 & BRI==0)
  Hg1[id000]=exp( xbg1 )*H0(Y[id000])
  Hg0[id000]=exp( xbg0 )*H0(Y[id000])
  
  ########### nbSC=1 & ORI=0 & BRI=0 
  id100=which(nbSC==1 & ORI==0 & BRI==0)
  if(length(id100)>0){
  Hg1[id100]=exp(xbg1)*( H0(Yst1[id100]) + Int1(Yst1[id100], Y[id100], beta.sc[1] + beta.gsc[1]) )  
  Hg0[id100]=exp(xbg0)*( H0(Yst1[id100]) + Int1(Yst1[id100], Y[id100], beta.sc[1]) )
  }
  ########### nbSC=2 & ORI=0 & BRI=0 
  id200=which(nbSC==2 & ORI==0 & BRI==0) 
  if(length(id200)>0){
  Hg1[id200]=exp(xbg1)*( H0(Yst1[id200]) + Int1(Yst1[id200], Yst2[id200], beta.sc[1] + beta.gsc[1]) + Int2(Yst2[id200], Y[id200], beta.sc[2] + beta.gsc[2]) ) 
  Hg0[id200]=exp(xbg0)*( H0(Yst1[id200]) + Int1(Yst1[id200], Yst2[id200], beta.sc[1] ) + Int2(Yst2[id200], Y[id200], beta.sc[2] ) ) 
  }
  ########### nbSC=3 & ORI=0 & BRI=0 
  id300=which(nbSC==3 & ORI==0 & BRI==0) 
  if(length(id300)>0){
  Hg1[id300]=exp(xbg1)*( H0(Yst1[id300]) + Int1(Yst1[id300], Yst2[id300], beta.sc[1] + beta.gsc[1]) + Int2(Yst2[id300], Yst3[id300], beta.sc[2] + beta.gsc[2]) +
                           Int3(Yst3[id300], Y[id300], beta.sc[3] + beta.gsc[3]) )
  
  Hg0[id300]=exp(xbg0)*( H0(Yst1[id300]) + Int1(Yst1[id300], Yst2[id300], beta.sc[1]) + Int2(Yst2[id300], Yst3[id300], beta.sc[2]) + Int3(Yst3[id300], Y[id300], beta.sc[3]) )
  }
  
  
  ############ 2)  Oopharectomy ORI==1
  
  ########### nbSC=0 & ORI=1 & BRI=0  
  id010=which(nbSC==0 & ORI==1 & BRI==0)
  if(length(id010)>0){
  Hg1[id010]=exp(xbg1)*( H0(Yor[id010]) + Int0(Yor[id010], Y[id010], beta.or+beta.gor) )
  Hg0[id010]=exp(xbg0)*( H0(Yor[id010]) + Int0(Yor[id010], Y[id010], beta.or) )
  }
  ########### nbSC=1 & ORI=1 & BRI=0 
  id110a=which(nbSC==1 & ORI==1 & BRI==0 & Yor<Yst1)
  if(length(id110a)>0){
  Hg1[id110a]=exp( xbg1 )*( H0(Yor[id110a]) +  Int0(Yor[id110a], Yst1[id110a], beta.or+beta.gor) + 
                              Int01(Yst1[id110a], Y[id110a], Yor[id110a], Yst1[id110a],beta.or+beta.gor, beta.sc[1] + beta.gsc[1], beta.orsc[1]  ) )
  Hg0[id110a]=exp( xbg0 )*( H0(Yor[id110a]) + Int0(Yor[id110a], Yst1[id110a], beta.or) + 
                              Int01(Yst1[id110a], Y[id110a],Yor[id110a], Yst1[id110a], beta.or, beta.sc[1] , beta.orsc[1]  ) )
  }
  id110b=which(nbSC==1 & ORI==1 & BRI==0 & Yor>=Yst1)
  if(length(id110b)>0){
  Hg1[id110b]=exp(xbg1)*( H0(Yst1[id110b]) + Int1(Yst1[id110b], Yor[id110b], beta.sc[1] + beta.gsc[1]) + 
                            Int01(Yor[id110b], Y[id110b], Yor[id110b], Yst1[id110b],  beta.or+beta.gor, beta.sc[1] + beta.gsc[1], beta.orsc[1]) ) 
  Hg0[id110b]= exp(xbg0)*( H0(Yst1[id110b]) + Int1(Yst1[id110b], Yor[id110b], beta.sc[1]) + 
                             Int01(Yor[id110b], Y[id110b], Yor[id110b], Yst1[id110b],  beta.or, beta.sc[1], beta.orsc[1]) ) 
  }
  
  ########### nbSC=2 & ORI=1 & BRI=0 
  # a)
  id210a=which(nbSC==2 & ORI==1 & BRI==0 &  Yor<Yst1)
  if(length(id210a)>0){
  Hg1[id210a]=exp( xbg1 )*( H0(Yor[id210a]) +  Int0(Yor[id210a], Yst1[id210a], beta.or+beta.gor) +
                              Int01(Yst1[id210a], Yst2[id210a], Yor[id210a], Yst1[id210a], beta.or+beta.gor,  beta.sc[1] + beta.gsc[1], beta.orsc[1]) +
                              Int02(Yst2[id210a], Y[id210a], Yor[id210a], Yst2[id210a], beta.or+beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
  
  Hg0[id210a]=exp( xbg0 )*( H0(Yor[id210a]) + Int0(Yor[id210a], Yst1[id210a], beta.or) +
                              Int01(Yst1[id210a], Yst2[id210a], Yor[id210a], Yst1[id210a], beta.or,  beta.sc[1], beta.orsc[1]) +
                              Int02(Yst2[id210a], Y[id210a], Yor[id210a], Yst2[id210a],    beta.or,  beta.sc[2], beta.orsc[2]) )
  }
  # b)
  id21b=which(nbSC==2 & ORI==1 & BRI==0 & Yor>=Yst1 & Yor<Yst2)
  if(length(id21b)>0){
  Hg1[id21b]=exp( xbg1 )*( H0(Yst1[id21b]) + Int1(Yst1[id21b], Yor[id21b], beta.sc[1] + beta.gsc[1]) +
                             Int01(Yor[id21b], Yst2[id21b], Yor[id21b], Yst1[id21b],  beta.or+beta.gor, beta.sc[1] + beta.gsc[1], beta.orsc[1]) +
                             Int02(Yst2[id21b], Y[id21b], Yor[id21b], Yst2[id21b], beta.or+beta.gor, beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
  
  Hg0[id21b]=exp( xbg0 )*( H0(Yst1[id21b]) + Int1(Yst1[id21b], Yor[id21b], beta.sc[1]) +
                             Int01(Yor[id21b], Yst2[id21b], Yor[id21b], Yst1[id21b],  beta.or, beta.sc[1], beta.orsc[1]) +
                             Int02(Yst2[id21b], Y[id21b],  Yor[id21b], Yst2[id21b], beta.or, beta.sc[2], beta.orsc[2]) )
  }
  # c)
  id21c=which(nbSC==2 & ORI==1 & BRI==0 & Yor>=Yst2)
  if(length(id21c)>0){
  Hg1[id21c]=exp( xbg1 )*( H0(Yst1[id21c]) + Int1(Yst1[id21c], Yst2[id21c], beta.sc[1] + beta.gsc[1]) +
                             Int2(Yst2[id21c], Yor[id21c], beta.sc[2] + beta.gsc[2]) +
                             Int02(Yor[id21c], Y[id21c], Yor[id21c], Yst2[id21c], beta.or +beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
  
  
  Hg0[id21c]=exp( xbg0 )*( H0(Yst1[id21c]) + Int1(Yst1[id21c], Yst2[id21c], beta.sc[1]) +
                             Int2(Yst2[id21c], Yor[id21c], beta.sc[2]) +
                             Int02(Yor[id21c], Y[id21c], Yor[id21c], Yst2[id21c], beta.or,  beta.sc[2], beta.orsc[2]) )
  
  }
  ########### nbSC=3 & ORI=1 & BRI=0 
  # a)
  id310a=which(nbSC==3 & ORI==1 & BRI==0 &  Yor<Yst1)
  if(length(id310a)>0){
  Hg1[id310a]=exp( xbg1 )*( H0(Yor[id310a]) +  Int0(Yor[id310a], Yst1[id310a], beta.or+beta.gor) +
                              Int01(Yst1[id310a], Yst2[id310a], Yor[id310a], Yst1[id310a], beta.or+beta.gor,  beta.sc[1] + beta.gsc[1], beta.orsc[1]) +
                              Int02(Yst2[id310a], Yst3[id310a], Yor[id310a], Yst2[id310a],  beta.or+beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) +
                              Int03(Yst3[id310a], Y[id310a],  Yor[id310a], Yst3[id310a], beta.or +beta.gor,  beta.sc[3] + beta.gsc[3], beta.orsc[3]) )
  
  Hg0[id310a]=exp( xbg0 )*( H0(Yor[id310a]) +  Int0(Yor[id310a], Yst1[id310a], beta.or) +
                              Int01(Yst1[id310a], Yst2[id310a], Yor[id310a], Yst1[id310a], beta.or,  beta.sc[1], beta.orsc[1]) +
                              Int02(Yst2[id310a], Yst3[id310a], Yor[id310a], Yst2[id310a], beta.or,  beta.sc[2], beta.orsc[2]) +
                              Int03(Yst3[id310a], Y[id310a], Yor[id310a], Yst3[id310a], beta.or,  beta.sc[3], beta.orsc[3]) )
  }
  
  # b)
  id31b=which(nbSC==3 & ORI==1 & BRI==0 & Yor>=Yst1 & Yor<Yst2)
  if(length(id31b)>0){
  Hg1[id31b]=exp( xbg1 )*( H0(Yst1[id31b]) + Int1(Yst1[id31b], Yor[id31b], beta.sc[1] + beta.gsc[1]) +
                             Int01(Yor[id31b], Yst2[id31b], Yor[id31b], Yst1[id31b], beta.or+beta.gor, beta.sc[1] + beta.gsc[1], beta.orsc[1]) +
                             Int02(Yst2[id31b], Yst3[id31b], Yor[id31b], Yst2[id31b], beta.or+beta.gor, beta.sc[2] + beta.gsc[2], beta.orsc[2]) +
                             Int03(Yst3[id31b], Y[id31b], Yor[id31b], Yst3[id31b],    beta.or+beta.gor, beta.sc[3] + beta.gsc[3], beta.orsc[3]) )
  
  Hg0[id31b]=exp( xbg0 )*( H0(Yst1[id31b]) + Int1(Yst1[id31b], Yor[id31b], beta.sc[1]) +
                             Int01(Yor[id31b], Yst2[id31b], Yor[id31b], Yst1[id31b],  beta.or, beta.sc[1], beta.orsc[1]) +
                             Int02(Yst2[id31b], Yst3[id31b], Yor[id31b], Yst2[id31b], beta.or, beta.sc[2], beta.orsc[2]) +
                             Int03(Yst3[id31b], Y[id31b], Yor[id31b], Yst3[id31b],   beta.or, beta.sc[3], beta.orsc[3]) )
  }
  # c)
  id31c=which(nbSC==3 & ORI==1 & BRI==0 & Yor>=Yst2 & Yor<Yst3 )
  if(length(id31c)>0){
  Hg1[id31c]=exp( xbg1 )*( H0(Yst1[id31c]) + Int1(Yst1[id31c], Yst2[id31c], beta.sc[1] + beta.gsc[1]) +
                             Int2(Yst2[id31c], Yor[id31c], beta.sc[2] + beta.gsc[2]) +
                             Int02(Yor[id31c], Yst3[id31c], Yor[id31c], Yst2[id31c], beta.or +beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) + 
                             Int03(Yst3[id31c], Y[id31c],  Yor[id31c], Yst3[id31c],  beta.or  +beta.gor,  beta.sc[3] + beta.gsc[3], beta.orsc[3] ) )
  
  Hg0[id31c]=exp( xbg0 )*( H0(Yst1[id31c]) + Int1(Yst1[id31c], Yst2[id31c], beta.sc[1]) +
                             Int2(Yst2[id31c], Yor[id31c], beta.sc[2]) +
                             Int02(Yor[id31c], Yst3[id31c], Yor[id31c], Yst2[id31c], beta.or,  beta.sc[2], beta.orsc[2]) + 
                             Int03(Yst3[id31c], Y[id31c],  Yor[id31c], Yst3[id31c],  beta.or,  beta.sc[3], beta.orsc[3] ) )
  }
  
  # c)
  id31d=which(nbSC==3 & ORI==1 & BRI==0 & Yor>=Yst3)
  if(length(id31d)>0){
  Hg1[id31d]=exp( xbg1 )*( H0(Yst1[id31d]) + Int1(Yst1[id31d], Yst2[id31d], beta.sc[1] + beta.gsc[1]) +
                             Int2(Yst2[id31d], Yst3[id31d], beta.sc[2] + beta.gsc[2]) +
                             Int3(Yst3[id31d], Yor[id31d], beta.sc[3] + beta.gsc[3]) +
                             Int03(Yor[id31d], Y[id31d], Yor[id31d], Yst3[id31d],  beta.or+beta.gor,  beta.sc[3] + beta.gsc[3], beta.orsc[3] ) )
  
  Hg0[id31d]=exp( xbg0 )*( H0(Yst1[id31d]) + Int1(Yst1[id31d], Yst2[id31d], beta.sc[1]) +
                             Int2(Yst2[id31d], Yst3[id31d], beta.sc[2]) +
                             Int3(Yst3[id31d], Yor[id31d], beta.sc[3]) +
                             Int03(Yor[id31d], Y[id31d], Yor[id31d], Yst3[id31d],   beta.or,  beta.sc[3], beta.orsc[3] ) )
  }
  
  
  ######################################  
  ######################################    Subjects with Mastectomy BRI==1
  ###################################### 
  
  ############ 1)   No Oopharectomy ORI==0
  
  
  
  ########### nbSC=0 & ORI=0 & BRI=1 
  id001=which(nbSC==0 & ORI==0 & BRI==1)
  Hg1[id001]=exp( xbg1 )*H0(Ybr[id001])
  Hg0[id001]=exp( xbg0 )*H0(Ybr[id001])
  
  ########### nbSC=1 & ORI=0 & BRI=1 
  
  id101a=which(SCI1==1 & ORI==0 & BRI==1 & Ybr<Yst1)
  Hg1[id101a]=exp( xbg1 )*H0(Ybr[id101a])
  Hg0[id101a]=exp( xbg0 )*H0(Ybr[id101a])
  
  id101b=which(SCI1==1 & ORI==0 & BRI==1 & Ybr>=Yst1)
  if(length(id101b)>0){
  Hg1[id101b]=exp(xbg1)*( H0(Yst1[id101b]) + Int1(Yst1[id101b], Ybr[id101b], beta.sc[1] + beta.gsc[1]) )
  Hg0[id101b]=exp(xbg0)*( H0(Yst1[id101b]) + Int1(Yst1[id101b], Ybr[id101b], beta.sc[1] ) )
  }
  ########### nbSC=2 & ORI=0 & BRI=1 
  
  id201a=which(nbSC==2 & ORI==0 & BRI==1 & Ybr<Yst1)
  Hg1[id201a]=exp(xbg1)*( H0(Ybr[id201a]) )
  Hg0[id201a]=exp(xbg0)*( H0(Ybr[id201a]) )
  
  id201b=which(nbSC==2 & ORI==0 & BRI==1 & Ybr>=Yst1 & Ybr<Yst2)
  if(length(id201b)>0){
  Hg1[id201b]=exp(xbg1)*( H0(Yst1[id201b]) + Int1(Yst1[id201b], Ybr[id201b], beta.sc[1] + beta.gsc[1]) )
  Hg0[id201b]=exp(xbg0)*( H0(Yst1[id201b]) + Int1(Yst1[id201b], Ybr[id201b], beta.sc[1] ) )
  }
  id201c=which(nbSC==2 & ORI==0 & BRI==1 & Ybr>=Yst2)
  if(length(id201c)>0){
  Hg1[id201c]=exp(xbg1)*( H0(Yst1[id201c]) + Int1(Yst1[id201c], Yst2[id201c], beta.sc[1] + beta.gsc[1]) + Int2(Yst2[id201c], Ybr[id201c], beta.sc[2] + beta.gsc[2]) )
  Hg0[id201c]=exp(xbg0)*( H0(Yst1[id201c]) + Int1(Yst1[id201c], Yst2[id201c], beta.sc[1] ) + Int2(Yst2[id201c], Ybr[id201c], beta.sc[2] ) )
  }
  ########### nbSC=3 & ORI=0 & BRI=1 
  
  id301a=which(nbSC==3 & ORI==0 & BRI==1 & Ybr<Yst1) 
  Hg1[id301a]=exp(xbg1)*( H0(Ybr[id301a]) ) 
  Hg0[id301a]=exp(xbg0)*( H0(Ybr[id301a]) ) 
  
  id301b=which(nbSC==3 & ORI==0 & BRI==1 & Ybr>=Yst1 & Ybr<Yst2)
  if(length(id301b)>0){
  Hg1[id301b]=exp(xbg1)*( H0(Yst1[id301b]) + Int1(Yst1[id301b], Ybr[id301b], beta.sc[1] + beta.gsc[1]) ) 
  Hg0[id301b]=exp(xbg0)*( H0(Yst1[id301b]) + Int1(Yst1[id301b], Ybr[id301b], beta.sc[1] ) ) 
  }
  id301c=which(nbSC==3 & ORI==0 & BRI==1 & Ybr>=Yst2 & Ybr<Yst3)
  if(length(id301c)>0){
  Hg1[id301c]=exp(xbg1)*( H0(Yst1[id301c]) + Int1(Yst1[id301c], Yst2[id301c], beta.sc[1] + beta.gsc[1]) + Int2(Yst2[id301c], Ybr[id301c], beta.sc[2] + beta.gsc[2]) ) 
  Hg0[id301c]=exp(xbg0)*( H0(Yst1[id301c]) + Int1(Yst1[id301c], Yst2[id301c], beta.sc[1] ) + Int2(Yst2[id301c], Ybr[id301c], beta.sc[2] ) ) 
  }
  id301d=which(nbSC==3 & ORI==0 & BRI==1 & Ybr>=Yst3) 
  if(length(id301d)>0){
  Hg1[id301d]=exp(xbg1)*( H0(Yst1[id301d]) + Int1(Yst1[id301d], Yst2[id301d], beta.sc[1] + beta.gsc[1]) + 
                            Int2(Yst2[id301d], Yst3[id301d], beta.sc[2] + beta.gsc[2]) +
                            Int3(Yst3[id301d], Ybr[id301d], beta.sc[3] + beta.gsc[3]) ) 
  
  Hg0[id301d]=exp(xbg0)*( H0(Yst1[id301d]) + Int1(Yst1[id301d], Yst2[id301d], beta.sc[1]) + 
                            Int2(Yst2[id301d], Yst3[id301d], beta.sc[2]) +
                            Int3(Yst3[id301d], Ybr[id301d], beta.sc[3]) ) 
  }
  
  ############ 2)   Oopharectomy ORI==1
  
  ########### nbSC=0
  
  id011a=which(nbSC==0 & ORI==1 & BRI==1 & Ybr<Yor)
  Hg1[id011a]=exp( xbg1 )*H0(Ybr[id011a])
  Hg0[id011a]=exp( xbg0 )*H0(Ybr[id011a])
  
  id011b=which(nbSC==0 & ORI==1 & BRI==1 & Ybr>=Yor)
  if(length(id011b)>0){
  Hg1[id011b]=exp(xbg1)*( H0(Yor[id011b]) + Int0(Yor[id011b], Ybr[id011b], beta.or+beta.gor) ) 
  Hg0[id011b]=exp(xbg0)*( H0(Yor[id011b]) + Int0(Yor[id011b], Ybr[id011b], beta.or) )
  }
  ########### nbSC=1
  
  # nbSC=1 & ORI=1 BRI=1 Yor<Yst1 & Ybr<Yor
  id111a1=which(SCI1==1 & ORI==1 & BRI==1 & Yor<Yst1 & Ybr<Yor)
  Hg1[id111a1]=exp( xbg1 )*H0(Ybr[id111a1])
  Hg0[id111a1]=exp( xbg0 )*H0(Ybr[id111a1])
  
  # nbSC=1 & ORI=1 BRI=1 Yor<Yst1 & Ybr>=Yor & Ybr<Yst1
  id111a2=which(SCI1==1 & ORI==1 & BRI==1 & Yor<Yst1 & Ybr>=Yor & Ybr<Yst1 )
  if(length(id111a2)>0){
  Hg1[id111a2]=exp(xbg1)*( H0(Yor[id111a2]) + Int0(Yor[id111a2], Ybr[id111a2], beta.or+beta.gor) )
  Hg0[id111a2]=exp(xbg0)*( H0(Yor[id111a2]) + Int0(Yor[id111a2], Ybr[id111a2], beta.or) )
  }
  # nbSC=1 & ORI=1 BRI=1 Yor<Yst1 & Ybr>=Yst1 
  id111a3=which(SCI1==1 & ORI==1 & BRI==1 & Yor<Yst1 & Ybr>=Yst1)
  if(length(id111a3)>0){
    Hg1[id111a3]=exp( xbg1 )*( H0(Yor[id111a3]) + Int0(Yor[id111a3], Ybr[id111a3], beta.or+beta.gor) +
                                 Int01(Yst1[id111a3], Ybr[id111a3], Yor[id111a3], Yst1[id111a3], beta.or+beta.gor, beta.sc[1] + beta.gsc[1],  beta.orsc[1] ) )
    
    Hg0[id111a3]=exp( xbg0 )*( H0(Yor[id111a3]) +  Int0(Yor[id111a3], Ybr[id111a3], beta.or) +
                                 Int01(Yst1[id111a3], Ybr[id111a3],Yor[id111a3], Yst1[id111a3], beta.or, beta.sc[1],  beta.orsc[1] ) )
  }
  
  ## nbSC=1 & ORI=1 & BRI=1 & Yor>=Yst1
  
  # nbSC=1 & ORI=1 & BRI=1 & Yor>=Yst1 & Ybr<Yst1
  id111b1=which(SCI1==1 & ORI==1 & BRI==1 & Yor>=Yst1 & Ybr<Yst1)
  Hg1[id111b1]=exp( xbg1 )*H0(Ybr[id111b1])
  Hg0[id111b1]=exp( xbg0 )*H0(Ybr[id111b1])
  
  # nbSC=1 & ORI=1 & BRI=1 & Yor>=Yst1 & Ybr>=Yst1 & Ybr<Yor
  id111b2=which(SCI1==1 & ORI==1 & BRI==1 & Yor>=Yst1 & Ybr>=Yst1 & Ybr<Yor)
  if(length(id111b2)>0){
  Hg1[id111b2]=exp(xbg1)*( H0(Yst1[id111b2]) + Int1(Yst1[id111b2], Ybr[id111b2], beta.sc[1] + beta.gsc[1]) )
  Hg0[id111b2]=exp(xbg0)*( H0(Yst1[id111b2]) + Int1(Yst1[id111b2], Ybr[id111b2], beta.sc[1] ) )
  }
  # nbSC=1 & ORI=1 & BRI=1 & Yor>=Yst1 & Ybr>=Yor
  id111b3=which(SCI1==1 & ORI==1 & BRI==1 & Yor>=Yst1 & Ybr>=Yor)
  if(length(id111b3)>0){
    Hg1[id111b3]=exp(xbg1)*( H0(Yst1[id111b3]) + Int1(Yst1[id111b3], Yor[id111b3], beta.sc[1] + beta.gsc[1]) + 
                               Int01(Yor[id111b3], Ybr[id111b3],Yor[id111b3], Yst1[id111b3], beta.or+beta.gor,  beta.sc[1] + beta.gsc[1],  beta.orsc[1] ) )
    
    Hg0[id111b3]= exp(xbg0)*( H0(Yst1[id111b3]) + Int1(Yst1[id111b3], Yor[id111b3], beta.sc[1]) + 
                                Int01(Yor[id111b3], Ybr[id111b3],Yor[id111b3], Yst1[id111b3], beta.or,  beta.sc[1],  beta.orsc[1] ) )
  }
  
  
  ########### nbSC=2
  
  # # a) 1
  id211a1=which(nbSC==2 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr<Yor)
  Hg1[id211a1]=exp( xbg1 )*( H0(Ybr[id211a1]) )
  Hg0[id211a1]=exp( xbg0 )*( H0(Ybr[id211a1]) )
  
  # # a) 2
  id211a2=which(nbSC==2 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr>=Yor & Ybr <Yst1)
  if(length(id211a2)>0){
  Hg1[id211a2]=exp( xbg1 )*( H0(Yor[id211a2]) + Int0(Yor[id211a2], Ybr[id211a2], beta.or+beta.gor) )
  Hg0[id211a2]=exp( xbg0 )*( H0(Yor[id211a2]) + Int0(Yor[id211a2], Ybr[id211a2], beta.or) )
  }
  # 
  # # a) 3
  id211a3=which(nbSC==2 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr>=Yst1 & Ybr<Yst2)
  if(length(id211a3)>0){
    Hg1[id211a3]=exp( xbg1 )*( H0(Yor[id211a3]) + Int0(Yor[id211a3], Yst1[id211a3], beta.or+beta.gor) +
                                 Int01(Yst1[id211a3], Ybr[id211a3],Yor[id211a3], Yst1[id211a3], beta.or+beta.gor,  beta.sc[1] + beta.gsc[1],  beta.orsc[1]) )
    
    Hg0[id211a3]=exp( xbg0 )*( H0(Yor[id211a3]) + Int0(Yor[id211a3], Yst1[id211a3], beta.or) +
                                 Int01(Yst1[id211a3], Ybr[id211a3],Yor[id211a3], Yst1[id211a3], beta.or,  beta.sc[1],  beta.orsc[1]) )
  }
  # # a) 4
  id211a4=which(nbSC==2 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr>=Yst2)
  if(length(id211a4)>0){
    Hg1[id211a4]=exp( xbg1 )*( H0(Yor[id211a4]) +  Int0(Yor[id211a4], Yst1[id211a4], beta.or+beta.gor) +
                                 Int01(Yst1[id211a4], Yst2[id211a4],Yor[id211a4], Yst1[id211a4], beta.or+beta.gor, beta.sc[1] + beta.gsc[1],  beta.orsc[1]) +
                                 Int02(Yst2[id211a4], Ybr[id211a4], Yor[id211a4], Yst2[id211a4], beta.or+beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
    
    Hg0[id211a4]=exp( xbg0 )*( H0(Yor[id211a4]) +  Int0(Yor[id211a4], Yst1[id211a4], beta.or) +
                                 Int01(Yst1[id211a4], Yst2[id211a4],Yor[id211a4], Yst1[id211a4], beta.or, beta.sc[1],  beta.orsc[1]) +
                                 Int02(Yst2[id211a4], Ybr[id211a4], Yor[id211a4], Yst2[id211a4],beta.or,  beta.sc[2], beta.orsc[2]) )
  }
  # b) 1
  id21b1=which(nbSC==2 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr<Yst1)
  
  Hg1[id21b1]=exp( xbg1 )*( H0(Ybr[id21b1]) )
  Hg0[id21b1]=exp( xbg0 )*( H0(Ybr[id21b1]) )
  
  # b) 2                         
  id21b2=which(nbSC==2 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr>=Yst1 & Ybr<Yor)
  if(length(id21b2)>0){
  Hg1[id21b2]=exp( xbg1 )*( H0(Yst1[id21b2]) + Int1(Yst1[id21b2], Ybr[id21b2], beta.sc[1] + beta.gsc[1]) )
  Hg0[id21b2]=exp( xbg0 )*( H0(Yst1[id21b2]) + Int1(Yst1[id21b2], Ybr[id21b2], beta.sc[1]) )
  }
  # b) 3                          
  id21b3=which(nbSC==2 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr>=Yor & Ybr<Yst2)
  if(length(id21b3)>0){
    Hg1[id21b3]=exp( xbg1 )*( H0(Yst1[id21b3]) + Int1(Yst1[id21b3], Yor[id21b3], beta.sc[1] + beta.gsc[1]) +
                                Int01(Yor[id21b3], Ybr[id21b3], Yor[id21b3], Yst1[id21b3], beta.or + beta.gor, beta.sc[1] + beta.gsc[1],  beta.orsc[1])  )
    
    Hg0[id21b3]=exp( xbg0 )*( H0(Yst1[id21b3]) + Int1(Yst1[id21b3], Yor[id21b3], beta.sc[1]) +
                                Int01(Yor[id21b3], Ybr[id21b3],Yor[id21b3], Yst1[id21b3], beta.or,  beta.sc[1],  beta.orsc[1])  )
  }
  # b) 4                          
  id21b4=which(nbSC==2 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr>=Yst2)
  if(length(id21b4)>0){
    Hg1[id21b4]=exp( xbg1 )*( H0(Yst1[id21b4]) + Int1(Yst1[id21b4], Yor[id21b4], beta.sc[1] + beta.gsc[1]) +
                                Int01(Yor[id21b4], Yst2[id21b4], Yor[id21b4], Yst1[id21b4], beta.or+beta.gor, beta.sc[1] + beta.gsc[1], beta.orsc[1]) +
                                Int02(Yst2[id21b4], Ybr[id21b4], Yor[id21b4], Yst2[id21b4], beta.or+beta.gor, beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
    
    Hg0[id21b4]=exp( xbg0 )*( H0(Yst1[id21b4]) + Int1(Yst1[id21b4], Yor[id21b4], beta.sc[1]) +
                                Int01(Yor[id21b4], Yst2[id21b4],Yor[id21b4], Yst1[id21b4], beta.or,  beta.sc[1], beta.orsc[1]) +
                                Int02(Yst2[id21b4], Ybr[id21b4], Yor[id21b4], Yst2[id21b4], beta.or,  beta.sc[2], beta.orsc[2]) )
  }
  
  # c)
  id21c1=which(nbSC==2 & ORI==1 & BRI==1 & Yor>=Yst2 & Ybr<Yst1)
  Hg1[id21c1]=exp( xbg1 )*( H0(Yst1[id21c1]) )
  Hg0[id21c1]=exp( xbg0 )*( H0(Yst1[id21c1]) )
  
  id21c2=which(nbSC==2 & ORI==1 & BRI==1 & Yor>=Yst2 & Ybr>=Yst1 & Ybr<Yst2)
  if(length(id21c2)>0){
  Hg1[id21c2]=exp( xbg1 )*( H0(Yst1[id21c2]) + Int1(Yst1[id21c2], Ybr[id21c2], beta.sc[1] + beta.gsc[1])  )
  Hg0[id21c2]=exp( xbg0 )*( H0(Yst1[id21c2]) + Int1(Yst1[id21c2], Ybr[id21c2], beta.sc[1])  )
  }
  
  id21c3=which(nbSC==2 & ORI==1 & BRI==1 & Yor>=Yst2 & Ybr>=Yst2 & Ybr<Yor)
  if(length(id21c3)>0){
  Hg1[id21c3]=exp( xbg1 )*( H0(Yst1[id21c3]) + Int1(Yst1[id21c3], Yst2[id21c3], beta.sc[1] + beta.gsc[1]) +
                              Int2(Yst2[id21c3], Ybr[id21c3], beta.sc[2] + beta.gsc[2]) )
  Hg0[id21c3]=exp( xbg0 )*( H0(Yst1[id21c3]) + Int1(Yst1[id21c3], Yst2[id21c3], beta.sc[1]) +
                              Int2(Yst2[id21c3], Ybr[id21c3], beta.sc[2]) )
  }
  id21c4=which(nbSC==2 & ORI==1 & BRI==1 & Yor>=Yst2 & Ybr>=Yor)
  if(length(id21c4)>0){
    Hg1[id21c4]=exp( xbg1 )*( H0(Yst1[id21c4]) + Int1(Yst1[id21c4], Yst2[id21c4], beta.sc[1] + beta.gsc[1]) +
                                Int2(Yst2[id21c4], Yor[id21c4], beta.sc[2] + beta.gsc[2]) +
                                Int02(Yor[id21c4], Y[id21c4], Yor[id21c4], Yst2[id21c4], beta.or+beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
    
    
    Hg0[id21c4]=exp( xbg0 )*( H0(Yst1[id21c4]) + Int1(Yst1[id21c4], Yst2[id21c4], beta.sc[1]) +
                                Int2(Yst2[id21c4], Yor[id21c4], beta.sc[2]) +
                                Int02(Yor[id21c4], Y[id21c4],Yor[id21c4], Yst2[id21c4], beta.or,  beta.sc[2],  beta.orsc[2]) )
    
  }
  
  
  ########### nbSC=3
  
  # # a)
  # ## 1)
  id311a1=which(nbSC==3 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr<Yor)
  if(length(id311a1)>0){
    Hg1[id311a1]=exp( xbg1 )*( H0(Ybr[id311a1]) )
    Hg0[id311a1]=exp( xbg0 )*( H0(Ybr[id311a1]) )
  }
  # ## 2)
  id311a2=which(nbSC==3 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr>=Yor & Ybr <Yst1)
  if(length(id311a2)>0){
  Hg1[id311a2]=exp( xbg1 )*( H0(Yor[id311a2]) + Int0(Yor[id311a2], Ybr[id311a2], beta.or+beta.gor) )
  Hg0[id311a2]=exp( xbg0 )*( H0(Yor[id311a2]) + Int0(Yor[id311a2], Ybr[id311a2], beta.or) )
  }
  #
  # ## 3)
  id311a3=which(nbSC==3 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr>=Yst1 & Ybr <Yst2)
  if(length(id311a3)>0){
    Hg1[id311a3]=exp( xbg1 )*( H0(Yor[id311a3]) + Int0(Yor[id311a3], Yst1[id311a3], beta.or+beta.gor) +
                                 Int01(Yst1[id311a3], Ybr[id311a3],Yor[id311a3], Yst1[id311a3], beta.or+beta.gor,  beta.sc[1] + beta.gsc[1],  beta.orsc[1]) )
    
    Hg0[id311a3]=exp( xbg0 )*( H0(Yor[id311a3]) + Int0(Yor[id311a3], Yst1[id311a3], beta.or) +
                                 Int01(Yst1[id311a3], Ybr[id311a3], Yor[id311a3], Yst1[id311a3], beta.or,  beta.sc[1],  beta.orsc[1]) )
  }
  # ## 4)
  id311a4=which(nbSC==3 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr>=Yst2  & Ybr <Yst3)
  if(length(id311a4)>0){
    Hg1[id311a4]=exp( xbg1 )*( H0(Yor[id311a4]) +  Int0(Yor[id311a4], Yst1[id311a4], beta.or+beta.gor) +
                                 Int01(Yst1[id311a4], Yst2[id311a4], Yor[id311a4], Yst1[id311a4], beta.or+beta.gor, beta.sc[1] + beta.gsc[1],  beta.orsc[1]) +
                                 Int02(Yst2[id311a4], Ybr[id311a4], Yor[id311a4], Yst2[id311a4], beta.or+beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
    
    Hg0[id311a4]=exp( xbg0 )*( H0(Yor[id311a4]) +  Int0(Yor[id311a4], Yst1[id311a4], beta.or) +
                                 Int01(Yst1[id311a4], Yst2[id311a4],Yor[id311a4], Yst1[id311a4], beta.or, beta.sc[1],  beta.orsc[1]) +
                                 Int02(Yst2[id311a4], Ybr[id311a4],Yor[id311a4], Yst2[id311a4], beta.or,  beta.sc[2], beta.orsc[2]) )
  }
  # ## 5)
  id311a5=which(nbSC==3 & ORI==1 & BRI==1 &  Yor<Yst1  & Ybr>=Yst3)
  if(length(id311a5)>0){
    Hg1[id311a5]=exp( xbg1 )*( H0(Yor[id311a5]) +  Int0(Yor[id311a5], Yst1[id311a5], beta.or+beta.gor) +
                                 Int01(Yst1[id311a5], Yst2[id311a5],Yor[id311a5], Yst1[id311a5], beta.or+beta.gor, beta.sc[1] + beta.gsc[1],  beta.orsc[1]) +
                                 Int02(Yst2[id311a5], Yst3[id311a5],Yor[id311a5], Yst2[id311a5], beta.or+beta.gor, beta.sc[2] + beta.gsc[2],  beta.orsc[2]) +
                                 Int03(Yst3[id311a5], Ybr[id311a5], Yor[id311a5], Yst3[id311a5], beta.or+beta.gor,  beta.sc[3] + beta.gsc[3], beta.orsc[3]) )
    
    Hg0[id311a5]=exp( xbg0 )*( H0(Yor[id311a5]) +  Int0(Yor[id311a5], Yst1[id311a5], beta.or) +
                                 Int01(Yst1[id311a5], Yst2[id311a5],Yor[id311a5], Yst1[id311a5], beta.or, beta.sc[1],  beta.orsc[1]) +
                                 Int02(Yst2[id311a5], Yst3[id311a5],Yor[id311a5], Yst2[id311a5], beta.or, beta.sc[2],  beta.orsc[2]) +
                                 Int03(Yst3[id311a5], Ybr[id311a5],Yor[id311a5], Yst3[id311a5], beta.or,  beta.sc[3],  beta.orsc[3]) )
  }
  # b)
  # ## 1)
  id31b1=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr<Yst1)
  Hg1[id31b1]=exp( xbg1 )*( H0(Ybr[id31b1]) )
  Hg0[id31b1]=exp( xbg0 )*( H0(Ybr[id31b1]) )
  #
  # ## 2)
  id31b2=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr>=Yst1 & Ybr<Yor)
  if(length(id31b2)>0){
  Hg1[id31b2]=exp( xbg1 )*( H0(Yst1[id31b2]) + Int1(Yst1[id31b2], Ybr[id31b2], beta.sc[1] + beta.gsc[1]) )
  Hg0[id31b2]=exp( xbg0 )*( H0(Yst1[id31b2]) + Int1(Yst1[id31b2], Ybr[id31b2], beta.sc[1]) )
  }
  # ## 3)
  id31b3=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr>=Yor & Ybr<Yst2)
  if(length(id31b3)>0){
    Hg1[id31b3]=exp( xbg1 )*( H0(Yst1[id31b3]) + Int1(Yst1[id31b3], Yor[id31b3], beta.sc[1] + beta.gsc[1]) +
                                Int01(Yor[id31b3], Ybr[id31b3], Yor[id31b3], Yst1[id31b3], beta.or+beta.gor, beta.sc[1] + beta.gsc[1],  beta.orsc[1])  )
    
    Hg0[id31b3]=exp( xbg0 )*( H0(Yst1[id31b3]) + Int1(Yst1[id31b3], Yor[id31b3], beta.sc[1]) +
                                Int01(Yor[id31b3], Ybr[id31b3], Yor[id31b3], Yst1[id31b3], beta.or,  beta.sc[1],  beta.orsc[1])  )
  }
  # ## 4)
  id31b4=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr>=Yst2 & Ybr<Yst3)
  if(length(id31b4)>0){
    Hg1[id31b4]=exp( xbg1 )*( H0(Yst1[id31b4]) + Int1(Yst1[id31b4], Yor[id31b4], beta.sc[1] + beta.gsc[1]) +
                                Int01(Yor[id31b4], Yst2[id31b4], Yor[id31b4], Yst1[id31b4], beta.or+beta.gor, beta.sc[1] + beta.gsc[1], beta.orsc[1]) +
                                Int02(Yst2[id31b4], Ybr[id31b4], Yor[id31b4], Yst2[id31b4], beta.or+beta.gor, beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
    
    Hg0[id31b4]=exp( xbg0 )*( H0(Yst1[id31b4]) + Int1(Yst1[id31b4], Yor[id31b4], beta.sc[1]) +
                                Int01(Yor[id31b4], Yst2[id31b4], Yor[id31b4], Yst1[id31b4], beta.or, beta.sc[1], beta.orsc[1]) +
                                Int02(Yst2[id31b4], Ybr[id31b4], Yor[id31b4], Yst2[id31b4], beta.or, beta.sc[2], beta.orsc[2]) )
  }
  ## 5)                          
  id31b5=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst1 & Yor<Yst2 & Ybr>=Yst3)
  if(length(id31b5)>0){
    Hg1[id31b5]=exp( xbg1 )*( H0(Yst1[id31b5]) + Int1(Yst1[id31b5], Yor[id31b5], beta.sc[1] + beta.gsc[1]) +
                                Int01(Yor[id31b5], Yst2[id31b5], Yor[id31b5], Yst1[id31b5], beta.or+beta.gor, beta.sc[1] + beta.gsc[1], beta.orsc[1]) +
                                Int02(Yst2[id31b5], Yst3[id31b5],Yor[id31b5], Yst2[id31b5],  beta.or+beta.gor, beta.sc[2] + beta.gsc[2], beta.orsc[2]) + 
                                Int03(Yst3[id31b5], Ybr[id31b5], Yor[id31b5], Yst3[id31b5], beta.or+beta.gor, beta.sc[3] + beta.gsc[3], beta.orsc[3]) )
    
    Hg0[id31b5]=exp( xbg0 )*( H0(Yst1[id31b5]) + Int1(Yst1[id31b5], Yor[id31b5], beta.sc[1]) +
                                Int01(Yor[id31b5], Yst2[id31b5],Yor[id31b5], Yst1[id31b5], beta.or, beta.sc[1], beta.orsc[1]) +
                                Int02(Yst2[id31b5], Yst3[id31b5],Yor[id31b5], Yst2[id31b5],  beta.or, beta.sc[2], beta.orsc[2]) + 
                                Int03(Yst3[id31b5], Ybr[id31b5], Yor[id31b5], Yst3[id31b5], beta.or, beta.sc[3], beta.orsc[3]) )
  }
  
  # c)
  ## 1)
  id31c1=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst2 & Yor<Yst3 & Ybr<Yst1)
  Hg1[id31c1]=exp( xbg1 )*( H0(Ybr[id21c1]) )
  Hg0[id31c1]=exp( xbg0 )*( H0(Ybr[id21c1]) )
  
  ## 2)
  id31c2=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst2 & Yor<Yst3  & Ybr>=Yst1 & Ybr<Yst2)
  if(length(id31c2)>0){
  Hg1[id31c2]=exp( xbg1 )*( H0(Yst1[id31c2]) + Int1(Yst1[id31c2], Ybr[id31c2], beta.sc[1] + beta.gsc[1])  )
  Hg0[id31c2]=exp( xbg0 )*( H0(Yst1[id31c2]) + Int1(Yst1[id31c2], Ybr[id31c2], beta.sc[1])  )
  }
  ## 3)
  id31c3=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst2  & Yor<Yst3  & Ybr>=Yst2 & Ybr<Yor)
  if(length(id31c3)>0){
  Hg1[id31c3]=exp( xbg1 )*( H0(Yst1[id31c3]) + Int1(Yst1[id31c3], Yst2[id31c3], beta.sc[1] + beta.gsc[1]) +
                              Int2(Yst2[id31c3], Ybr[id31c3], beta.sc[2] + beta.gsc[2]) )
  Hg0[id31c3]=exp( xbg0 )*( H0(Yst1[id31c3]) + Int1(Yst1[id31c3], Yst2[id31c3], beta.sc[1]) +
                              Int2(Yst2[id31c3], Ybr[id31c3], beta.sc[2]) )
  }
  ## 4)
  id31c4=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst2  & Yor<Yst3  & Ybr>=Yor & Ybr<Yst3)
  if(length(id31c4)>0){
    Hg1[id31c4]=exp( xbg1 )*( H0(Yst1[id31c4]) + Int1(Yst1[id31c4], Yst2[id31c4], beta.sc[1] + beta.gsc[1]) +
                                Int2(Yst2[id31c4], Yor[id31c4], beta.sc[2] + beta.gsc[2]) +
                                Int02(Yor[id31c4], Ybr[id31c4], Yor[id31c4], Yst2[id31c4], beta.or+beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) )
    
    Hg0[id31c4]=exp( xbg0 )*( H0(Yst1[id31c4]) + Int1(Yst1[id31c4], Yst2[id31c4], beta.sc[1]) +
                                Int2(Yst2[id31c4], Yor[id31c4], beta.sc[2]) +
                                Int02(Yor[id31c4], Ybr[id31c4], Yor[id31c4], Yst2[id31c4], beta.or,  beta.sc[2], beta.orsc[2]) )
  }
  
  ## 5)
  id31c5=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst2  & Yor<Yst3  & Ybr>=Yst3)
  if(length(id31c5)>0){
    Hg1[id31c5]=exp( xbg1 )*( H0(Yst1[id31c5]) + Int1(Yst1[id31c5], Yst2[id31c5], beta.sc[1] + beta.gsc[1]) +
                                Int2(Yst2[id31c5], Yor[id31c5], beta.sc[2] + beta.gsc[2]) +
                                Int02(Yor[id31c5], Yst3[id31c5], Yor[id31c5], Yst2[id31c5], beta.or+beta.gor,  beta.sc[2] + beta.gsc[2], beta.orsc[2]) +
                                Int03(Yst3[id31c5], Ybr[id31c5], Yor[id31c5], Yst3[id31c5], beta.or+beta.gor,  beta.sc[3] + beta.gsc[3], beta.orsc[3]) )
    
    Hg0[id31c5]=exp( xbg0 )*( H0(Yst1[id31c5]) + Int1(Yst1[id31c5], Yst2[id31c5], beta.sc[1]) +
                                Int2(Yst2[id31c5], Yor[id31c5], beta.sc[2]) +
                                Int02(Yor[id31c5], Yst3[id31c5], Yor[id31c5], Yst2[id31c5], beta.or,  beta.sc[2], beta.orsc[2]) +
                                Int03(Yst3[id31c5], Ybr[id31c5], Yor[id31c5], Yst3[id31c5],  beta.or,  beta.sc[3], beta.orsc[3]) )
  }
  # d) Yor>=Yst3
  ## 1)
  id31d1=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst3 & Ybr<Yst1)
  Hg1[id31d1]=exp( xbg1 )*( H0(Ybr[id31d1]) )
  Hg0[id31d1]=exp( xbg0 )*( H0(Ybr[id31d1]) )
  
  ## 2)
  id31d2=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst3 & Ybr>=Yst1 & Ybr<Yst2)
  if(length(id31c2)>0){
  Hg1[id31d2]=exp( xbg1 )*( H0(Yst1[id31d2]) + Int1(Yst1[id31d2], Ybr[id31d2], beta.sc[1] + beta.gsc[1])  )
  Hg0[id31d2]=exp( xbg0 )*( H0(Yst1[id31d2]) + Int1(Yst1[id31d2], Ybr[id31d2], beta.sc[1])  )
  }
  ## 3)
  id31d3=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst3 & Ybr>=Yst2 & Ybr<Yst3)
  if(length(id31d3)>0){
  Hg1[id31d3]=exp( xbg1 )*( H0(Yst1[id31d3]) + Int1(Yst1[id31d3], Yst2[id31d3], beta.sc[1] + beta.gsc[1]) +
                              Int2(Yst2[id31d3], Ybr[id31d3], beta.sc[2] + beta.gsc[2]) )
  Hg0[id31d3]=exp( xbg0 )*( H0(Yst1[id31d3]) + Int1(Yst1[id31d3], Yst2[id31d3], beta.sc[1]) +
                              Int2(Yst2[id31d3], Ybr[id31d3], beta.sc[2]) )
  }
  
  ## 4)
  id311d4=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst3 & Ybr>=Yst3 & Ybr<Yor)
  if(length(id311d4)>0){
    Hg1[id311d4]=exp( xbg1 )*( H0(Yst1[id311d4]) + Int1(Yst1[id311d4], Yst2[id311d4], beta.sc[1] + beta.gsc[1]) +
                                 Int2(Yst2[id311d4], Yst3[id311d4], beta.sc[2] + beta.gsc[2]) +
                                 Int03(Yst3[id311d4], Ybr[id311d4], Yor[id311d4], Yst3[id311d4], beta.or+beta.gor,  beta.sc[3] + beta.gsc[3], beta.orsc[3]) )
    
    Hg0[id311d4]=exp( xbg0 )*( H0(Yst1[id311d4]) + Int1(Yst1[id311d4], Yst2[id311d4], beta.sc[1]) +
                                 Int2(Yst2[id311d4], Yst3[id311d4], beta.sc[2] ) +
                                 Int03(Yst3[id311d4], Ybr[id311d4],Yor[id311d4], Yst3[id311d4],  beta.or,  beta.sc[3], beta.orsc[3]) )
  }
  ## 5)
  id311d5=which(nbSC==3 & ORI==1 & BRI==1 & Yor>=Yst3 & Ybr>=Yor)
  if(length(id311d5)>0){
    Hg1[id311d5]=exp( xbg1 )*( H0(Yst1[id311d5]) + Int1(Yst1[id311d5], Yst2[id311d5], beta.sc[1] + beta.gsc[1]) +
                                 Int2(Yst2[id311d5], Yst3[id311d5], beta.sc[2] + beta.gsc[2]) +
                                 Int03(Yst3[id311d5], Yor[id311d5],Yor[id311d5], Yst3[id311d5], beta.or+beta.gor,  beta.sc[3] + beta.gsc[3], beta.orsc[3]) +
                                 Int03(Yor[id311d5], Ybr[id311d5], Yor[id311d5], Yst3[id311d5], beta.or+beta.gor,  beta.sc[3] + beta.gsc[3], beta.orsc[3]) )
    
    Hg0[id311d5]=exp( xbg0 )*( H0(Yst1[id311d5]) + Int1(Yst1[id311d5], Yst2[id311d5], beta.sc[1]) +
                                 Int2(Yst2[id311d5], Yst3[id311d5], beta.sc[2]) +
                                 Int03(Yst3[id311d5], Yor[id311d5], Yor[id311d5], Yst3[id311d5], beta.or,  beta.sc[3], beta.orsc[3]) +
                                 Int03(Yor[id311d5], Ybr[id311d5],Yor[id311d5], Yst3[id311d5],  beta.or,  beta.sc[3], beta.orsc[3]) )
  }
  
  H=(1-cp)*Hg0 + cp*Hg1
  
  return( list(H=H) )
}


### combine loghazard and cumHazard function
marg_fun<- function(mpars, Y, Yor, Yst1, Yst2, Yst3,  Ybr, TDf_or, TDf_sc, cp){
  
 logh=loghaz(mpars, Y, Yor, Yst1, Yst2, Yst3,  Ybr, TDf_or, TDf_sc, cp)$logh
 H= cumHaz(mpars, Y, Yor, Yst1, Yst2, Yst3,  Ybr, TDf_or, TDf_sc, cp)$H
  
return( list(logh=logh, H=H) )
}











