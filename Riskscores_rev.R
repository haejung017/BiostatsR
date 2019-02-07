#Make Risk Scores matrix for different models
RiskScores=function(data, modelobject, modelNumber, n=65){
  b1=data
  RS=matrix(ncol=n,nrow=nrow(b1),byrow=TRUE) 
  
  if(modelNumber==1){
    bG=modelobject$par[7]
    print(bG)
    cov=as.matrix(b1[,33])
    
    XB=cov*bG
    for(i in 1:nrow(RS)){
      RS[i,]= rep(XB[i,],n)
    }
    
  }else if(modelNumber==2){
    bG=modelobject$par[7]
    bS=modelobject$par[8]
    bO=modelobject$par[9]
    print(c(bG,bS,bO))
    cov=as.matrix(b1[,33:35])
    XB=cbind(cov[,1]*bG,cov[,2]*bS,cov[,3]*bO)
    for(i in 1:nrow(RS)){
      RS[i,]= rep(sum(XB[i,]),n)
    }
    
  }else if(modelNumber==3){
    bG=modelobject$par[7]
    bS=modelobject$par[8]
    bO=modelobject$par[9]
    etaS=exp(modelobject$par[12])
    cov=as.matrix(b1[,33:35])
    t=seq(from=0, to=(n-1), by=1)+16
    
    ov=as.matrix(ifelse(b1$time- b1$ov.censortime>0,b1$ov.censortime,99))
    st1=as.matrix(b1$st1); st1[is.na(st1)]=99
    
    xb=function(i){
    xbeta=cov[i,1]*bG+
    ifelse(t>=st1[i,],1,0)*(exp(-etaS*(t-st1[i,])))*bS+
    ifelse(t>=ov[i,],1,0)*bO
    return(xbeta)
    }
    for(i in 1:nrow(RS)){
      RS[i,]= xb(i)
    }
    
  }else if(modelNumber==4){
    bG=modelobject$par[3]
    bS=modelobject$par[4]
    bO=modelobject$par[5]
    print(c(bG,bS,bO))
    cov=as.matrix(b1[,33:35])
    XB=cbind(cov[,1]*bG,cov[,2]*bS,cov[,3]*bO)
    for(i in 1:nrow(RS)){
      RS[i,]= rep(sum(XB[i,]),n)
    }
    
  }else if(modelNumber==5){
    bG=modelobject$par[3]
    bS=modelobject$par[4]
    bO=modelobject$par[5]
    etaS=exp(modelobject$par[6])
    cov=as.matrix(b1[,33:35])
    
    t=seq(from=0, to=(n-1), by=1)+16
    
    ov=as.matrix(ifelse(b1$time- b1$ov.censortime>0,b1$ov.censortime,99))
    st1=as.matrix(b1$st1); st1[is.na(st1)]=99
    
    xb=function(i){
      xbeta=cov[i,1]*bG+
        ifelse(t>=st1[i,],1,0)*(exp(-etaS*(t-st1[i,])))*bS+
        ifelse(t>=ov[i,],1,0)*bO
      return(xbeta)
    }
    for(i in 1:nrow(RS)){
      RS[i,]= xb(i)
    }
  }else if(modelNumber==6){
    
    data2=Data_preparation(b1)
    TVCstatus = cbind(data2$v1+16,data2$v2) 
    parms=tran_parm(modelobject$par)
    cest=parms[[1]]
    phiv=parms[[2]]
    for(i in 1:nrow(RS)){
      RS[i,]= ComputeRS_fullmodel(TVCstatus[i,1:5],TVCstatus[i,6:10],b1$mgene[i],cest,phiv, n=n)
    }
  }else if(modelNumber==7){
    data2=Data_preparation(b1)
    TVCstatus = cbind(data2$v1+16,data2$v2) 
    parms=tran_parm2(modelobject$par)
    cest=parms[[1]]
    phiv=parms[[2]]
    for(i in 1:nrow(RS)){
      RS[i,]= ComputeRS_fullmodel(TVCstatus[i,1:5],TVCstatus[i,6:10],b1$mgene[i],cest,phiv, n=n)
    }
  }else if(modelNumber==8){
    bG=modelobject$par[3]
    print(bG)
    cov=as.matrix(b1[,33])
    XB=cov*bG
    for(i in 1:nrow(RS)){
      RS[i,]= rep(XB[i,],n)
    }
  }
  colnames(RS)<-round(seq(from=0, to=(n-1), by=1)+16,2)
  return(RS)
}
ZED <- function(t,u,phi){
  exp(-(t-u)*phi)
}
#Computing RiskScores for the full model(Gene,3MS,1BO)
ComputeRS_fullmodel=function(v1,v2,G,cest,phiv,n=65){
t=seq(from=0, to=(n-1), by=1)+16
u=t
sc1 = ifelse(length(v1[which(v2==1)])==1,v1[which(v2==1)],99)
sc2 = ifelse(length(v1[which(v2==2)])==1,v1[which(v2==2)],99)
sc3 = ifelse(length(v1[which(v2==3)])==1,v1[which(v2==3)],99)
ov = ifelse(length(v1[which(v2==4)])==1,v1[which(v2==4)],99)
br = ifelse(length(v1[which(v2==5)])==1,v1[which(v2==5)],99)
zsc1 = ifelse(u>=sc1 & u<sc2,1,0)
zsc2 = ifelse(u>=sc2 & u<sc3,1,0)
zsc3 = ifelse(u>=sc3 ,1,0)
zor = ifelse(u<ov,0,1)
zbr = ifelse(u<br,0,1)


zsc1td = ifelse(zsc1==1,ZED(t=u, u=sc1, phiv[1]),0)
zsc2td = ifelse(zsc2==1,ZED(t=u, u=sc2, phiv[2]),0)
zsc3td = ifelse(zsc3==1,ZED(t=u, u=sc3, phiv[3]),0)
zortd = ifelse(zor==1,ZED(t=u, u=ov, phiv[4]),0)

RS=G*cest[1]+(cest[2])*zsc1*zsc1td+(cest[3])*zsc2*zsc2td+(cest[4])*zsc3*zsc3td+cest[5]*zor*zortd 
return(RS)
}
#transform parameters for full model(Gene,3MS,1BO): competing risk
tran_parm=function(theta){
  # vbeta for breast  
  beta.gen1 = theta[7]
  beta.sc1_1 = theta[8]
  beta.sc2_1 = theta[9]
  beta.sc3_1 = theta[10]
  beta.or1 = theta[11]
  phiv <- c(exp(theta[c(14,15,16)]),0)
  fp <- exp(theta[(length(theta)-1):length(theta)])
  cest <- c(beta.gen1,#1
            beta.sc1_1,
            beta.sc2_1,
            beta.sc3_1,
            beta.or1)
  return(list(cest,phiv))
}
#transform parameters for full model(Gene,3MS,1BO): Non-competing risk
tran_parm2=function(theta){
  # vbeta for breast  
  beta.gen1 = theta[3]
  beta.sc1_1 = theta[4]
  beta.sc2_1 = theta[5]
  beta.sc3_1 = theta[6]
  beta.or1 = theta[7]
  phiv <- c(exp(theta[c(8,9,10)]),0)
  fp <- exp(theta[length(theta)])
  cest <- c(beta.gen1,#1
            beta.sc1_1,
            beta.sc2_1,
            beta.sc3_1,
            beta.or1)
  return(list(cest,phiv))
}
#Function to get TVC times and indicators
Data_preparation <- function(data_raw){
  data <- data_raw
  agemin <- 16
  data = data[data$currentage>=agemin,]
  time0 = data$time-agemin
  
  ip <- which(data$proband==1)
  
  #setting up indicator variable for the time dependent variables
  br_time <- data$br.censortime-agemin
  ov_time <- data$ov.censortime-agemin
  sc1_time <- data$st1 - agemin
  sc2_time <- data$st2 - agemin
  sc3_time <- data$st3 - agemin
  delta_BR <- ifelse(time0-br_time > 0,1,0)
  delta_OR <- ifelse(time0-ov_time > 0,1,0)
  
  time_vec <- cbind(sc1_time,sc2_time,sc3_time,ov_time*delta_OR,br_time*delta_BR)
  
  clean_time_vec <- function(time_vec_r){
    time_vec_r[is.na(time_vec_r)] <- 0  
    time_vec_r[time_vec_r==0] <- NA 
    return(time_vec_r)
  }
  
  #sort time vector, make indicator vector for event order
  time_vec2 <- NULL
  for(i in 1:nrow(time_vec)){
    tve <- sort(clean_time_vec(time_vec[i,]), na.last=TRUE)
    time_vec2 <- rbind(time_vec2,tve)
  }
  
  indicator_vec <- NULL
  for(i in 1:nrow(time_vec)){
    tve1 <- clean_time_vec(time_vec[i,])
    tve2 <- order(tve1, na.last=TRUE)
    tve2[is.na(sort(tve1, na.last=TRUE))] <- NA
    indicator_vec <- rbind(indicator_vec,tve2)
  }
  
  cagep <- data$currentage[ip]-agemin
  timep <- data$time[ip]-agemin
  statusp <- data$affect[ip]  # proband disease status at the study entry
  
  
  #setting up indicator variable for the removal status at the time of study entry
  br_time_p <- (data$br.censortimep-agemin)[ip]
  ov_time_p <- (data$ov.censortimep-agemin)[ip]
  sc1_time_p <- (data$st1p - agemin)[ip]
  sc2_time_p <- (data$st2p - agemin)[ip]
  sc3_time_p <- (data$st3p - agemin)[ip]
  
  sc1_time_pn <- (data$st1 - agemin)[ip]
  sc2_time_pn <- (data$st2 - agemin)[ip]
  sc3_time_pn <- (data$st3 - agemin)[ip]
  
  delta_BR_p <- ifelse(br_time_p < cagep, 1, 0) # mastectomy status at study entry
  delta_OR_p <- ifelse(ov_time_p < cagep, 1, 0) # oophorectomy status at study entry
  delta_SC_p <- data$screenp[ip] 
  
  
  # subsetting the probands according to his affection status at study entry
  br_time_p0 <- br_time_p[statusp==0]
  ov_time_p0 <- ov_time_p[statusp==0]
  sc1_time_p0 <- sc1_time_p[statusp==0]
  sc2_time_p0 <- sc2_time_p[statusp==0]
  sc3_time_p0 <- sc3_time_p[statusp==0]
  cagep0 <- cagep[statusp==0]
  timep0 <- timep[statusp==0]
  
  delta_BR_p0 <- delta_BR_p[statusp==0]
  delta_OR_p0 <- delta_OR_p[statusp==0]
  #delta_SC_p0 <- delta_SC_p[statusp==0]
  
  time_vec_p0 <- cbind(sc1_time_p0,sc2_time_p0,sc3_time_p0,ov_time_p0*delta_OR_p0,br_time_p0*delta_BR_p0)
  
  #sort time vector
  time_vec2_p0 <- NULL
  for(i in 1:nrow(time_vec_p0)){
    tve <- sort(clean_time_vec(time_vec_p0[i,]), na.last=TRUE)
    time_vec2_p0 <- rbind(time_vec2_p0,tve)
  }
  indicator_vec_p0 <- NULL
  for(i in 1:nrow(time_vec_p0)){
    tve1 <- clean_time_vec(time_vec_p0[i,])
    tve2 <- order(tve1, na.last=TRUE)
    tve2[is.na(sort(tve1, na.last=TRUE))] <- NA
    indicator_vec_p0 <- rbind(indicator_vec_p0,tve2)
  }
  #
  # br_time_p1 <- br_time_p[statusp==1]
  # ov_time_p1 <- ov_time_p[statusp==1]
  # sc1_time_p1 <- sc1_time_pn[statusp==1]
  # sc2_time_p1 <- sc2_time_pn[statusp==1]
  # sc3_time_p1 <- sc3_time_pn[statusp==1]
  # cagep1 <- cagep[statusp==1]
  # timep1 <- timep[statusp==1]
  # 
  # delta_BR_p1 <- ifelse(timep1-br_time_p1 > 0,1,0)
  # delta_OR_p1 <- ifelse(timep1-ov_time_p1 > 0,1,0)
  
  br_time_p1 <- br_time_p[statusp==1 |statusp==2]
  ov_time_p1 <- ov_time_p[statusp==1 |statusp==2]
  sc1_time_p1 <- sc1_time_p[statusp==1 |statusp==2]#
  sc2_time_p1 <- sc2_time_p[statusp==1 |statusp==2]
  sc3_time_p1 <- sc3_time_p[statusp==1 |statusp==2]
  cagep1 <- cagep[statusp==1 |statusp==2]
  timep1 <- timep[statusp==1 |statusp==2]
  
  delta_BR_p1 <- ifelse(timep1-br_time_p1 > 0,1,0)
  delta_OR_p1 <- ifelse(timep1-ov_time_p1 > 0,1,0)
  #
  # br_time_p1 <- br_time_p[statusp==1 |statusp==2]
  # ov_time_p1 <- ov_time_p[statusp==1 |statusp==2]
  # sc1_time_p1 <- sc1_time_p[statusp==1 |statusp==2]
  # sc2_time_p1 <- sc2_time_p[statusp==1 |statusp==2]
  # sc3_time_p1 <- sc3_time_p[statusp==1 |statusp==2]
  # cagep1 <- cagep[statusp==1 |statusp==2]
  # timep1 <- timep[statusp==1 |statusp==2]
  # 
  # delta_BR_p1 <- delta_BR_p[statusp==1 |statusp==2]
  # delta_OR_p1 <- delta_OR_p[statusp==1 |statusp==2]
  #delta_SC_p1 <- delta_SC_p[statusp==1 |statusp==2]
  
  time_vec_p1 <- cbind(sc1_time_p1,sc2_time_p1,sc3_time_p1,ov_time_p1*delta_OR_p1,br_time_p1*delta_BR_p1)
  
  #sort time vector
  time_vec2_p1 <- NULL
  for(i in 1:nrow(time_vec_p1)){
    tve <- sort(clean_time_vec(time_vec_p1[i,]), na.last=TRUE)
    time_vec2_p1 <- rbind(time_vec2_p1,tve)
  }
  indicator_vec_p1 <- NULL
  for(i in 1:nrow(time_vec_p1)){
    tve1 <- clean_time_vec(time_vec_p1[i,])
    tve2 <- order(tve1, na.last=TRUE)
    tve2[is.na(sort(tve1, na.last=TRUE))] <- NA
    indicator_vec_p1 <- rbind(indicator_vec_p1,tve2)
  }
  
  br_time_p2 <- br_time_p[statusp==2]
  ov_time_p2 <- ov_time_p[statusp==2]
  sc1_time_p2 <- sc1_time_pn[statusp==2]
  sc2_time_p2 <- sc2_time_pn[statusp==2]
  sc3_time_p2 <- sc3_time_pn[statusp==2]
  cagep2 <- cagep[statusp==2]
  timep2 <- timep[statusp==2]
  
  delta_BR_p2 <- ifelse(timep2-br_time_p2 > 0,1,0)
  delta_OR_p2 <- ifelse(timep2-ov_time_p2 > 0,1,0)
  
  # br_time_p1 <- br_time_p1[statusp==1 |statusp==2]
  # ov_time_p1 <- ov_time_p1[statusp==1 |statusp==2]
  # sc1_time_p1 <- sc1_time_p1[statusp==1 |statusp==2]
  # sc2_time_p1 <- sc2_time_p1[statusp==1 |statusp==2]
  # sc3_time_p1 <- sc3_time_p1[statusp==1 |statusp==2]
  # cagep1 <- cagep[statusp==1 |statusp==2]
  # timep1 <- timep[statusp==1 |statusp==2]
  # 
  # delta_BR_p1 <- ifelse(timep1-br_time_p1 > 0,1,0)
  # delta_OR_p1 <- ifelse(timep1-ov_time_p1 > 0,1,0)
  #
  # br_time_p1 <- br_time_p[statusp==1 |statusp==2]
  # ov_time_p1 <- ov_time_p[statusp==1 |statusp==2]
  # sc1_time_p1 <- sc1_time_p[statusp==1 |statusp==2]
  # sc2_time_p1 <- sc2_time_p[statusp==1 |statusp==2]
  # sc3_time_p1 <- sc3_time_p[statusp==1 |statusp==2]
  # cagep1 <- cagep[statusp==1 |statusp==2]
  # timep1 <- timep[statusp==1 |statusp==2]
  # 
  # delta_BR_p1 <- delta_BR_p[statusp==1 |statusp==2]
  # delta_OR_p1 <- delta_OR_p[statusp==1 |statusp==2]
  #delta_SC_p1 <- delta_SC_p[statusp==1 |statusp==2]
  
  time_vec_p2 <- cbind(sc1_time_p2,sc2_time_p2,sc3_time_p2,ov_time_p2*delta_OR_p2,br_time_p2*delta_BR_p2)
  
  #sort time vector
  time_vec2_p2 <- NULL
  for(i in 1:nrow(time_vec_p2)){
    tve <- sort(clean_time_vec(time_vec_p2[i,]), na.last=TRUE)
    time_vec2_p2 <- rbind(time_vec2_p2,tve)
  }
  indicator_vec_p2 <- NULL
  for(i in 1:nrow(time_vec_p2)){
    tve1 <- clean_time_vec(time_vec_p2[i,])
    tve2 <- order(tve1, na.last=TRUE)
    tve2[is.na(sort(tve1, na.last=TRUE))] <- NA
    indicator_vec_p2 <- rbind(indicator_vec_p2,tve2)
  }
  
  return(list(v1=time_vec2,v2=indicator_vec,
              v1p0=time_vec2_p0, v2p0=indicator_vec_p0,
              v1p1=time_vec2_p1, v2p1=indicator_vec_p1,
              v1p2=time_vec2_p2, v2p2=indicator_vec_p2))   
}
Data_preparation2 <- function(data_raw){
  data <- data_raw
  agemin <- 16
  data = data[data$currentage>=agemin,]
  time0 = data$time-agemin

  #setting up indicator variable for the time dependent variables
  br_time <- data$br.censortime-agemin
  ov_time <- data$ov.censortime-agemin
  sc1_time <- data$st1 - agemin
  sc2_time <- sc3_time <- NA
  delta_BR <- ifelse(time0-br_time > 0,1,0)
  delta_OR <- ifelse(time0-ov_time > 0,1,0)
  
  time_vec <- cbind(sc1_time,sc2_time,sc3_time,ov_time*delta_OR,br_time*delta_BR)
  
  clean_time_vec <- function(time_vec_r){
    time_vec_r[is.na(time_vec_r)] <- 0  
    time_vec_r[time_vec_r==0] <- NA 
    return(time_vec_r)
  }
  
  #sort time vector, make indicator vector for event order
  time_vec2 <- NULL
  for(i in 1:nrow(time_vec)){
    tve <- sort(clean_time_vec(time_vec[i,]), na.last=TRUE)
    time_vec2 <- rbind(time_vec2,tve)
  }
  
  indicator_vec <- NULL
  for(i in 1:nrow(time_vec)){
    tve1 <- clean_time_vec(time_vec[i,])
    tve2 <- order(tve1, na.last=TRUE)
    tve2[is.na(sort(tve1, na.last=TRUE))] <- NA
    indicator_vec <- rbind(indicator_vec,tve2)
  }
  
  return(list(v1=time_vec2,v2=indicator_vec))   
}



##################################################################################################

#Compute penetrance for different models
Penetrance=function(data, modelobject, modelNumber, n=65){
  b1=data
  PEN=matrix(ncol=n,nrow=nrow(b1),byrow=TRUE)
  
  if(modelNumber%in%c(1,2,3,6)){#Competing Risk
    base.parms=exp(modelobject$par[1:6])
    fp=exp(modelobject$par[(length(modelobject$par)-1):(length(modelobject$par))])
    
    PEN=ComputePen_comp(b1,PEN,base.parms,fp,modelobject,modelNumber,n)
    
  }else{#Non-competing Risk
    base.parms=exp(modelobject$par[1:2])
    fp=exp(modelobject$par[length(modelobject$par)])
    
    PEN=ComputePen_noncomp(b1,PEN,base.parms,fp,modelobject,modelNumber,n)
    
  }
  colnames(PEN)<-round(seq(from=0, to=(n-1), by=1)+16,2)
  return(PEN)
}  
  
ComputePen_noncomp=function(b1,PEN,base.parms,fp,modelobject,modelNumber,n){
  t=seq(from=0, to=(n-1), by=1)
  
  lambda1=base.parms[1];rho1=base.parms[2]
  bcumhaz1 = (lambda1*t)^rho1 
  G=b1$M;S=b1$S;O=b1$O
  if(modelNumber==4){
    for(i in 1:nrow(b1)){
      H=bcumhaz1*exp(G[i]*modelobject$par[3]+S[i]*modelobject$par[4]+O[i]*modelobject$par[5])
      PEN[i,]=1-((1+H/fp[1])^(-fp[1])) # 3 TIC, G,S,O
    }
  }else if(modelNumber==5){
    mat=cbind(Data_preparation2(b1)[[1]],Data_preparation2(b1)[[2]])
    cest=c(base.parms,modelobject$par[3:4],0,0,modelobject$par[5])
    phiv=c(exp(modelobject$par[6]),0,0,0)
    for(i in 1:nrow(b1)){
      v=mat[i,]
      PENi=vector()
      for(j in 1:length(t)){
        vv1=v[1:5];vv2=v[6:10]
        vv1[is.na(v[1:5])]<-99;vv2[is.na(v[6:10])]<-99
        vv1=matrix(rep(vv1,2),ncol=5,byrow=TRUE);vv2=matrix(rep(vv2,2),ncol=5,byrow=TRUE)
        H=(exp(G[i]*modelobject$par[3]))*CumH_c_ED(v1=vv1,v2=vv2,u=c(t[j],t[j]),cest=cest,affp=FALSE,type=1,phi=phiv)[1,]
        PENi[j]=1-((1-H/fp[1])^(-fp[1]))
      }
      PEN[i,]=PENi
    }
  }else if(modelNumber==7){
    mat=cbind(Data_preparation(b1)[[1]],Data_preparation(b1)[[2]])
    cest=c(base.parms,tran_parm2(modelobject$par)[[1]])
    phiv=tran_parm2(modelobject$par)[[2]]
    for(i in 1:nrow(b1)){
      v=mat[i,]
      PENi=vector()
     for(j in 1:length(t)){
          vv1=v[1:5];vv2=v[6:10]
          vv1[is.na(v[1:5])]<-99;vv2[is.na(v[6:10])]<-99
          vv1=matrix(rep(vv1,2),ncol=5,byrow=TRUE);vv2=matrix(rep(vv2,2),ncol=5,byrow=TRUE)
          H=(exp(G[i]*modelobject$par[3]))*CumH_c_ED(v1=vv1,v2=vv2,u=c(t[j],t[j]),cest=cest,affp=FALSE,type=1,phi=phiv)[1,]
          PENi[j]=1-((1-H/fp[1])^(-fp[1]))
          }
    PEN[i,]=PENi
    }
  }else if(modelNumber==8){
    for(i in 1:nrow(b1)){
      H=bcumhaz1*exp(G[i]*modelobject$par[3]) # 1 TIC, G
      PEN[i,]=1-((1+H/fp[1])^(-fp[1]))
    }
  }
  
  return(PEN) 
}
ComputePen_comp=function(b1,PEN,base.parms,fp,modelobject,modelNumber,n){
  t=seq(from=0, to=(n-1), by=1)
  
  lambda1=base.parms[1];rho1=base.parms[2];lambda2=base.parms[3];rho2=base.parms[4]
  lambda3=base.parms[5];rho3=base.parms[6]
 
  G=b1$M;S=b1$S;O=b1$O
  if(modelNumber==1){
    cov=as.matrix(cbind(G,rep(0,50),rep(0,50)))
    PEN=matrix(ncol=n,nrow=nrow(b1),byrow=TRUE) 
    for (i in 1:nrow(b1)){
      PEN[i,]=PrintPenTIC2(modelobject,cov=cov[i,],type=1,res=n)[,2]
    }
  }else if(modelNumber==2){
    cov=as.matrix(cbind(G,S,O))
    PEN=matrix(ncol=n,nrow=nrow(b1),byrow=TRUE) 
    for (i in 1:nrow(b1)){
      PEN[i,]=PrintPenTIC(modelobject,cov=cov[i,],type=1,res=n)[,2]
    }
  }else if(modelNumber==3){
    data2 <- Data_preparation2(b1) 
    TVCstatus = cbind(data2$v1,data2$v2)
    parms=tran_parm4(modelobject$par)  
    PEN=matrix(ncol=n,nrow=nrow(b1),byrow=TRUE) 
    for (i in 1:nrow(b1)){
      PEN[i,]=PrintPen(parms,eventage=TVCstatus[i,1:5],eventind=TVCstatus[i,6:10],G=G[i],type=1,res=n)[,2]
    }
  }else if(modelNumber==6){
    data2 <- Data_preparation(b1) 
    TVCstatus = cbind(data2$v1,data2$v2)  
    parms=tran_parm3(modelobject$par)  
    PEN=matrix(ncol=n,nrow=nrow(b1),byrow=TRUE) 
    for (i in 1:nrow(b1)){
      PEN[i,]=PrintPen(parms,eventage=TVCstatus[i,1:5],eventind=TVCstatus[i,6:10],G=G[i],type=1,res=n)[,2]
    }
  }
  
  return(PEN) 
}

  
PrintPen=function(parms, eventage, eventind, G, type=1, res=100, data){

  v=c(eventage,eventind)
  cest=parms[[1]]
  phiv=parms[[2]]
  fp=parms[[3]]
  
  
  #estimate
  t=seq(from=0, to=(res-1), by=1)
  mat=cbind(matrix(rep(v,res),ncol=10,nrow=res,byrow=TRUE),t)
  #v=mat[100,]
  
  
  Computation=function(type,G){
    C_F=vector()
    for(i in 1:nrow(mat)){
      Fi=try(proband_penf_ED(v=mat[i,],type=type, affp=FALSE, cest=cest, phi=phiv, frailty=TRUE, fp=fp,G=G),silent = TRUE)
      if(inherits(Fi ,'try-error')){
        Fi=7
      }
      C_F[i]=Fi
    }
    return(C_F)
  }
  Pen1_C=Computation(type,G)
  
  t=t+16
  plotdat1=cbind(t,Pen1_C)
  return(round(plotdat1,3))
  
}
proband_penf_ED <- function(cest, v, affp, type, phi, frailty, fp,G=1){
  if(G==0){
    cest[c(7,12,13)]=0
  }
  
  v1 <- v[1:5]
  v2 <- v[6:10]
  ind=which(v1<v[11])
  v1 =  v1[ind]
  v2 =  v2[ind]
  v1 = c(v1,v[11])
  v1 = v1[!is.na(v1)]
  v2 = v2[!is.na(v2)]
  
  j <- beta_vec1 <- beta_vec2 <- beta_vec3 <- NULL
  v2his <- vector(mode="numeric")
  # v1
  # v2
  # i=2
  beta_sc=0; beta_or=0; lowsc=0; lowor=0; phisc=0;
  for(i in 1:length(v1)){
    #i=2
    
    if(i!=1){
      beta = coeff_c_ED(v2[i-1],type=1, affp, cest);
      if(v2[i-1]==1 | v2[i-1]==2 | v2[i-1]==3 ){
        phisc = phi[v2[i-1]]
        beta_sc = beta;
        lowsc = v1[i-1];
      }else if(v2[i-1]==4){
        beta_or = beta;
        lowor = v2[i-1];
      }else if(v2[i-1]==5 & type==1){
        break
      }
    }else{
      beta = 0;
      beta_sc = beta;
      beta_or = beta;
    }
    
    beta1 <- c(beta_sc,lowsc,beta_or,lowor)
    beta2 <- ifelse(i==1, 0, coeff_c_ED(v2[i-1],type=2, cest=cest, affp=affp) )
    beta3 <- ifelse(i==1, 0, coeff_c_ED(v2[i-1],type=3, cest=cest, affp=affp) )
    
    
    # if(type==1){
    #   inter = ifelse(c(q==1,q==1), c(0,0), coeff_interaction_c_IS(v2[q-1], v2his, bref))
    #   v2his= c(v2his, v2[q-1]); 
    #   intg = c(intg, inter[1])
    #   intw = c(intw, inter[2])
    # }else{
    #   intg <- intw <- 0  
    # }
    low=ifelse(i==1,0,v1[i-1])
    #print(low)
    if(low==0){
      H1=H2=H3=0
      H=c(H1,H2,H3)
    }else{
      vv1=v[1:5];vv2=v[6:10]
      vv1[is.na(v[1:5])]<-99;vv2[is.na(v[6:10])]<-99
      vv1=matrix(rep(vv1,2),ncol=5,byrow=TRUE);vv2=matrix(rep(vv2,2),ncol=5,byrow=TRUE)
      H1=CumH_c_ED(v1=vv1,v2=vv2,u=c(low,low),cest=cest,affp=FALSE,type=1,phi=phi)[1,]
      H2=CumH_c_ED(v1=vv1,v2=vv2,u=c(low,low),cest=cest,affp=FALSE,type=2,phi=phi)[1,]
      H3=CumH_c_ED(v1=vv1,v2=vv2,u=c(low,low),cest=cest,affp=FALSE,type=3,phi=phi)[1,]
      H=c(H1,H2,H3)
    }
    Penpiece = penf.comp_ED(f=integrand_IS_ED, v1=v1, v2=v2, est=cest, affp=affp,beta_vec1=beta1,beta_vec2=beta2,beta_vec3=beta3,
                            type=type,phi=phi,phisc=phisc,lb=low,H=H,frailty=frailty,fp=fp,lower=low,upper=v1[i])
    j = c(j,Penpiece)
  }
  j = sum(j,na.rm = TRUE)
  Penet <- j
  
  # print(round(Penet,3))
  return(Penet)
}
penf.comp_ED <- function(f,v1,v2,est, affp,type,beta_vec1,beta_vec2,beta_vec3,phi,phisc,lb,H,frailty,fp,lower,upper){  
  integrate(f,v1,v2,est=est, affp=affp,beta_vec1=beta_vec1,
            beta_vec2=beta_vec2,beta_vec3=beta_vec3,type=type,H=H,
            phi=phi,phisc=phisc,frailty=frailty,fp=fp,lb=lb,lower=lower,upper=upper)$value
}
integrand_IS_ED <- function(u,v1,v2,est,affp,type,beta_vec1,beta_vec2,beta_vec3,phi,phisc,lb,H,frailty=frailty,fp=c(20,20,20)){
  u <- u
  #u=runif(5,low,v1[i])
  #u
  #  print(u)
  # v1 = v1[!is.na(v1)]
  # v2 = v2[!is.na(v2)]
  sc1 = ifelse(length(v1[which(v2==1)])==1,v1[which(v2==1)],99)
  sc2 = ifelse(length(v1[which(v2==2)])==1,v1[which(v2==2)],99)
  sc3 = ifelse(length(v1[which(v2==3)])==1,v1[which(v2==3)],99)
  ov = ifelse(length(v1[which(v2==4)])==1,v1[which(v2==4)],99)
  br = ifelse(length(v1[which(v2==5)])==1,v1[which(v2==5)],99)
  zsc1 = ifelse(u>=sc1 & u<sc2,1,0)
  zsc2 = ifelse(u>=sc2 & u<sc3,1,0)
  zsc3 = ifelse(u>=sc3 ,1,0)
  zor = ifelse(u<ov,0,1)
  zbr = ifelse(u<br,0,1)
  
  
  zsc1td = ifelse(zsc1==1,ZED(t=u, u=sc1, phi[1]),0)
  zsc2td = ifelse(zsc2==1,ZED(t=u, u=sc2, phi[2]),0)
  zsc3td = ifelse(zsc3==1,ZED(t=u, u=sc3, phi[3]),0)
  zortd = ifelse(zor==1,ZED(t=u, u=ov, phi[4]),0)
  
  if(type==1){
    xbeta = est[7] + (est[8])*zsc1*zsc1td +(est[9])*zsc2*zsc2td+(est[10])*zsc3*zsc3td + 
      est[11]*zor*zortd  
    haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta)
    #print(c(xbeta, est[7]+sum(beta_vec1)))
  }else if(type==2){
    haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(est[12])  
  }
  
  H1 = H[1]-(1-zbr)*exp(est[7])*sapply(u, chaz, lam1=est[1],rho1=est[2],
                                       betasc=beta_vec1[1],betaor=beta_vec1[3], phisc=phisc,phior=phi[4],
                                       lowsc=beta_vec1[2], lowor=beta_vec1[4],low=lb)
  H2 = H[2]-(1-zor)*(est[3]*u)^est[4]*exp(est[12]+ sum(beta_vec2))
  H3 = H[3]-(est[5]*u)^est[6]*exp(est[13]+ sum(beta_vec3) ) 
  
  # H1 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=1)*exp(est[7])
  # H2 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=2)*exp(est[12])
  # H3 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=TRUE,type=3)*exp(est[14])
  
  
  if(frailty==TRUE){
    if(type==1){
      f = haz*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))*(exp(H3))
    }else{
      f = haz*((1-H2/fp[2])^(-fp[2]-1))*((1-H1/fp[1])^(-fp[1]))*(exp(H3))
    }
  } else{
    Hp <- H1+H2+H3
    f = haz*exp(Hp) # h * S   
  }
  return(f)
}
tran_parm3=function(theta){
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  lambda3 = exp(theta[5])
  rho3 = exp(theta[6])
  
  # vbeta for breast  
  beta.gen1 = theta[7]
  beta.sc1_1 = theta[8]
  beta.sc2_1 = theta[9]
  beta.sc3_1 = theta[10]
  beta.or1 = theta[11]
  
  # vbeta for ovarian
  beta.gen2 = theta[12]
  
  # vbeta for death
  beta.gen3 = theta[13]
  
  phiv <- c(exp(theta[c(14,15,16)]),0)
  phi_sc = c(phiv[c(1:3)],rep(0,3)) # 3 rate parms, 3 dummy parms
  phi_or = c(0, 0) # 1 rate parm, 1 dummy parm
  # }else if(tvctype=="CO"){
  #   phiv = c( exp(theta[c(14,15,16,17)]) , theta[c(18,19,20,21)] ) # BO PE version
  #   phi_sc = phiv[c(1:3,5:7)] # 3 rate parms, 3 convergence parms
  #   phi_or = phiv[c(4,8)] # 1 rate parm, 1 convergence parm    
  # }else{
  #   phi_sc=rep(0,6)
  #   phi_or=c(0,0)
  # }
  fp <- exp(theta[(length(theta)-1):length(theta)])
  cest <- c(lambda1,rho1,
            lambda2,rho2,
            lambda3,rho3,
            beta.gen1,#7
            beta.sc1_1,
            beta.sc2_1,
            beta.sc3_1,
            beta.or1,
            beta.gen2,#12
            beta.gen3)
  return(list(cest,phiv,fp))
}
tran_parm4=function(theta){
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  lambda3 = exp(theta[5])
  rho3 = exp(theta[6])
  
  # vbeta for breast  
  beta.gen1 = theta[7]
  beta.sc1_1 = theta[8]
  beta.sc2_1 = 0
  beta.sc3_1 = 0
  beta.or1 = theta[9]
  
  # vbeta for ovarian
  beta.gen2 = theta[10]
  
  # vbeta for death
  beta.gen3 = theta[11]
  
  phiv <- c(exp(theta[12]),0,0,0)
  phi_sc = c(phiv[c(1:3)],rep(0,3)) # 3 rate parms, 3 dummy parms
  phi_or = c(0, 0) # 1 rate parm, 1 dummy parm
  # }else if(tvctype=="CO"){
  #   phiv = c( exp(theta[c(14,15,16,17)]) , theta[c(18,19,20,21)] ) # BO PE version
  #   phi_sc = phiv[c(1:3,5:7)] # 3 rate parms, 3 convergence parms
  #   phi_or = phiv[c(4,8)] # 1 rate parm, 1 convergence parm    
  # }else{
  #   phi_sc=rep(0,6)
  #   phi_or=c(0,0)
  # }
  fp <- exp(theta[(length(theta)-1):length(theta)])
  cest <- c(lambda1,rho1,
            lambda2,rho2,
            lambda3,rho3,
            beta.gen1,#7
            beta.sc1_1,
            beta.sc2_1,
            beta.sc3_1,
            beta.or1,
            beta.gen2,#12
            beta.gen3)
  return(list(cest,phiv,fp))
}


tran_parm5=function(theta){
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  lambda3 = exp(theta[5])
  rho3 = exp(theta[6])
  
  # vbeta for breast  
  beta.gen1 = theta[7]
  beta.sc1_1 = theta[8]
  beta.or1 = theta[9]
  
  # vbeta for ovarian
  beta.gen2 = theta[10]
  
  # vbeta for death
  beta.gen3 = theta[11]
  
  fp <- exp(theta[(length(theta)-1):length(theta)])
  cest <- c(lambda1,rho1,
            lambda2,rho2,
            lambda3,rho3,
            beta.gen1,#7
            beta.sc1_1,
            beta.or1,
            beta.gen2,#12
            beta.gen3)
  return(list(cest,fp))
}
tran_parm6=function(theta){
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  lambda3 = exp(theta[5])
  rho3 = exp(theta[6])
  
  # vbeta for breast  
  beta.gen1 = theta[7]
  beta.sc1_1 = 0
  beta.or1 = 0
  
  # vbeta for ovarian
  beta.gen2 = theta[8]
  
  # vbeta for death
  beta.gen3 = theta[9]
  
  fp <- exp(theta[(length(theta)-1):length(theta)])
  cest <- c(lambda1,rho1,
            lambda2,rho2,
            lambda3,rho3,
            beta.gen1,#7
            beta.sc1_1,
            beta.or1,
            beta.gen2,#12
            beta.gen3)
  return(list(cest,fp))
}
PrintPenTIC=function(modelobject, eventage, eventind, cov, type=1, res=100, data){
  output=modelobject$par
  #RobustVar=RobustVar
  
  #v=c(eventage-16,eventind)
  parms=tran_parm5(output)  
  cest=parms[[1]]
  fp=parms[[2]]
  
  
  #estimate
  t=seq(from=0, to=(res-1), by=1)
  
  Computation=function(t,cest,fp,cov){
    t=t
    Fi=vector()
    for(i in 1:length(t)){
      Fi[i]=penf.comp_TIC(f=integrand_TIC, est=cest, type=1,frailty=TRUE,cov=cov,fp=fp,lower=0,upper=t[i])
    }
    Fi
  }
  
  Pen1_C=Computation(t,cest,fp,cov)
  #Pen1SE_C=Computation_SE(RobustVar,output,mat,G=1,type)
  t=t+16
  plotdat1=cbind(t,Pen1_C)
  return(round(plotdat1,3))
}
penf.comp_TIC <- function(f,est,type,frailty,cov,fp,lower,upper){  
  integrate(f,est=est,type=type,frailty=frailty,cov=cov,fp=fp,lower=lower,upper=upper)$value
}
integrand_TIC <- function(u,est,type,frailty=frailty,cov,fp=c(20,20,20)){
  u <- u
  G=cov[1]
  S=cov[2]
  O=cov[3]
  if(type==1){
    xbeta = est[7]*G + (est[8])*S +(est[9])*O  
    haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta)
    #print(c(xbeta, est[7]+sum(beta_vec1)))
  }else if(type==2){
    haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(est[12])  
  }
  
  H1 = -(est[1]*u)^est[2]*exp(xbeta)
  H2 = -(est[3]*u)^est[4]*exp(est[10]*G)
  H3 = -(est[5]*u)^est[6]*exp(est[11]*G) 
  
  # H1 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=1)*exp(est[7])
  # H2 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=2)*exp(est[12])
  # H3 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=TRUE,type=3)*exp(est[14])
  
  
  if(frailty==TRUE){
    if(type==1){
      f = haz*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))*(exp(H3))
    }else{
      f = haz*((1-H2/fp[2])^(-fp[2]-1))*((1-H1/fp[1])^(-fp[1]))*(exp(H3))
    }
  } else{
    Hp <- H1+H2+H3
    f = haz*exp(Hp) # h * S   
  }
  return(f)
}
PrintPenTIC2=function(modelobject, eventage, eventind, cov, type=1, res=100, data){
  output=modelobject$par
  #RobustVar=RobustVar
  
  #v=c(eventage-16,eventind)
  parms=tran_parm6(output)  
  cest=parms[[1]]
  fp=parms[[2]]
  
  
  #estimate
  t=seq(from=0, to=(res-1), by=1)
  
  Computation=function(t,cest,fp,cov){
    t=t
    Fi=vector()
    for(i in 1:length(t)){
      Fi[i]=penf.comp_TIC(f=integrand_TIC, est=cest, type=1,frailty=TRUE,cov=cov,fp=fp,lower=0,upper=t[i])
    }
    Fi
  }
  
  Pen1_C=Computation(t,cest,fp,cov)
  #Pen1SE_C=Computation_SE(RobustVar,output,mat,G=1,type)
  t=t+16
  plotdat1=cbind(t,Pen1_C)
  return(round(plotdat1,3))
}