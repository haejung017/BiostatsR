PlotPen=function(modelobject, RobustVar, eventage, eventind, type=1, res=100){
  
  output=modelobject$par
  RobustVar=RobustVar
  
  v=c(eventage-16,eventind)
  parms=tran_parm(output)  
  cest=parms[[1]]
  phiv=parms[[2]]
  fp=parms[[3]]
  
  
  #estimate
  t=seq(from=0, to=55, length.out = res)
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
  Pen1_C=Computation(type,1)
  Pen1SE_C=Computation_SE(RobustVar,output,mat,G=1,type)
  Pen1_NC=Computation(type,0)
  Pen1SE_NC=Computation_SE(RobustVar,output,mat,G=0,type)
  
  t=t+16
  plotdat1=cbind(t,Pen1_C)
  plotdat2=cbind(t,Pen1_C+1.96*Pen1SE_C)
  plotdat3=cbind(t,Pen1_C-1.96*Pen1SE_C)
  
  plotdat4=cbind(t,Pen1_NC)
  plotdat5=cbind(t,Pen1_NC+1.96*Pen1SE_NC)
  plotdat6=cbind(t,Pen1_NC-1.96*Pen1SE_NC)
  
  plotdat=cbind(rbind(plotdat1,plotdat2,plotdat3,plotdat4,plotdat5,plotdat6),
                rep(1:6, each=res))#,rep(c(1,2), each=100)),
  #rep(c(1,2), each=200))
  plotdat=as.data.frame(plotdat)
  #plotdat$V3=factor(plotdat$V3, labels=c("PE","ED","CO"))
  colnames(plotdat)=c("t","e","g")#,"cancer")
  
  plotdat[plotdat$e>1 & !is.na(plotdat$e),"e"]<-1
  
  C=plotdat[plotdat$g==1,]
  C_u=plotdat[plotdat$g==2,]
  C_l=plotdat[plotdat$g==3,]
  
  NC=plotdat[plotdat$g==4,]
  NC_u=plotdat[plotdat$g==5,]
  NC_l=plotdat[plotdat$g==6,]
  
  print(c(C_l[nrow(C),"e"],"~",C[nrow(C),"e"],"~",C_u[nrow(C),"e"]))
  print(c(NC_l[nrow(C),"e"],"~",NC[nrow(C),"e"],"~",NC_u[nrow(C),"e"]))
  # 
  # ET0=plotdat[plotdat$t>=0 & plotdat$t<tx & plotdat$g=="ED",]
  # ET1=plotdat[plotdat$t>=tx & plotdat$g=="ED",]
  # 
  # CT0=plotdat[plotdat$t>=0 & plotdat$t<tx & plotdat$g=="CO",]
  # CT1=plotdat[plotdat$t>=tx & plotdat$g=="CO",]
  
  # some variables for ggplot
  y.limits = ifelse(rep(type==1,2),c(0,1),c(0,0.8))
  g.title= ifelse(type==1, "Breast Cancer Penetrance", "Ovarian Cancer Penetrance")
  
  ##output
  ggplot(plotdat, aes(x=t,y=e,group=g))+
    geom_line(data=C,color=2)+
    geom_line(data=C_u,linetype=3,color=2)+
    geom_line(data=C_l,linetype=3,color=2)+
    
    geom_line(data=NC,linetype=4,color=4)+
    geom_line(data=NC_u,linetype=3,color=4)+
    geom_line(data=NC_l,linetype=3,color=4)+
    
    # Screen effect segments
    # geom_line(data=PT0,linetype=1,color=2)+ # No screen segment
    # geom_line(data=PT1,linetype=1,color=2)+ # screen segment
    # 
    # geom_line(data=ET0,linetype=1,color=2)+ # No screen segment
    # geom_line(data=ET1,linetype=1,color=2)+ # screen segment
    # 
    # geom_line(data=CT0,linetype=1,color=2)+ # No screen segment
    # geom_line(data=CT1,linetype=1,color=2)+ # screen segment
    
  # Vertical Line
  # geom_segment(data=PT1, x=tx, y=0, xend=tx, yend=PT1[1,2], linetype=2, color=1)+
  # geom_segment(data=ET1, x=tx, y=0, xend=tx, yend=ET1[1,2], linetype=2, color=1)+
  # geom_segment(data=CT1, x=tx, y=0, xend=tx, yend=CT1[1,2], linetype=2, color=1)+
  # 
  # geom_segment(data=PT1, x=15, y=1.5, xend=70, yend=1.5, linetype=4, color=2)+
  # geom_segment(data=ET1, x=15, y=1.5, xend=70, yend=1.5, linetype=4, color=2)+
  # geom_segment(data=CT1, x=15, y=0.5, xend=70, yend=0.5, linetype=4, color=2)+
  # geom_segment(data=CT1, x=15, y=2, xend=70, yend=2, linetype=4, color=2)+
  # Points
  # geom_point(x=tx, y=0, shape=1,fill="black") + # right before t_x
  # geom_point(data=PT1, x=tx, y=PT1[1,2]) + # right after t_x
  # geom_point(data=ET1, x=tx, y=ET1[1,2]) + # right after t_x
  # geom_point(data=CT1, x=tx, y=CT1[1,2]) + # right after t_x
  
  # # Arrow heads
  # annotate("segment",x=69,xend=70,yend=PT1[nrow(PT1),2],y=PT1[nrow(PT1),2],color="red",arrow = arrow(length = unit(0.2, "cm")))+
  
  # Themes
  scale_x_continuous(name = "Ages in years") +
    scale_y_continuous(name = "Penetrance",limits=y.limits) + # ylab
    ggtitle(g.title) + #title
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size=15), # Centering of the title
          axis.text.x = element_text(face="bold",size=12),  # Change axis ticks labels
          axis.title.x = element_text(size=15,margin = margin(t = 20)),
          axis.title.y = element_text(size=15))
  # facet_wrap(~cancer,scales="free")
  #grid.locator(unit="native") 
  # grid.brackets(200,153,200,67 , lwd=0.5, h=0.03,type=1,col="brown")
  # grid.brackets(200,288,200,202, lwd=0.5, h=0.03,type=1,col="brown")
  # grid.brackets(200,397,200,308, lwd=0.5, h=0.03,type=1,col="brown")
  # grid.brackets(200,424,200,397, lwd=0.4, h=0.03,type=1,col="brown")
}
PrintPen=function(modelobject, RobustVar, eventage, eventind, G, type=1, res=100, data){

  output=modelobject$par
  RobustVar=RobustVar
  
  v=c(eventage,eventind)
  parms=tran_parm(output)  
  cest=parms[[1]]
  phiv=parms[[2]]
  fp=parms[[3]]
  
  
  #estimate
  t=seq(from=0, to=55, length.out = res)
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
  #Pen1SE_C=Computation_SE(RobustVar,output,mat,G=1,type)
  #Pen1_NC=Computation(type,0)
  #Pen1SE_NC=Computation_SE(RobustVar,output,mat,G=0,type)
  
  t=t+16
  plotdat1=cbind(t,Pen1_C)
  return(round(plotdat1,3))
  # plotdat2=cbind(t,Pen1_C+1.96*Pen1SE_C)
  # plotdat3=cbind(t,Pen1_C-1.96*Pen1SE_C)
  # 
  # plotdat4=cbind(t,Pen1_NC)
  # plotdat5=cbind(t,Pen1_NC+1.96*Pen1SE_NC)
  # plotdat6=cbind(t,Pen1_NC-1.96*Pen1SE_NC)
  # 
  # plotdat=cbind(rbind(plotdat1,plotdat2,plotdat3,plotdat4,plotdat5,plotdat6),
  #               rep(1:6, each=res))#,rep(c(1,2), each=100)),
  # #rep(c(1,2), each=200))
  # plotdat=as.data.frame(plotdat)
  # #plotdat$V3=factor(plotdat$V3, labels=c("PE","ED","CO"))
  # colnames(plotdat)=c("t","e","g")#,"cancer")
  # 
  # plotdat[plotdat$e>1 & !is.na(plotdat$e),"e"]<-1
  # 
  # C=plotdat[plotdat$g==1,]
  # C_u=plotdat[plotdat$g==2,]
  # C_l=plotdat[plotdat$g==3,]
  # 
  # NC=plotdat[plotdat$g==4,]
  # NC_u=plotdat[plotdat$g==5,]
  # NC_l=plotdat[plotdat$g==6,]
  # 
  # print(c(C_l[nrow(C),"e"],"~",C[nrow(C),"e"],"~",C_u[nrow(C),"e"]))
  # print(c(NC_l[nrow(C),"e"],"~",NC[nrow(C),"e"],"~",NC_u[nrow(C),"e"]))
  # # 
  # # ET0=plotdat[plotdat$t>=0 & plotdat$t<tx & plotdat$g=="ED",]
  # # ET1=plotdat[plotdat$t>=tx & plotdat$g=="ED",]
  # # 
  # # CT0=plotdat[plotdat$t>=0 & plotdat$t<tx & plotdat$g=="CO",]
  # # CT1=plotdat[plotdat$t>=tx & plotdat$g=="CO",]
  # 
  # # some variables for ggplot
  # y.limits = ifelse(rep(type==1,2),c(0,1),c(0,0.3))
  # g.title= ifelse(type==1, "Breast Cancer Penetrance", "Ovarian Cancer Penetrance")
  # 
  # ##output
  # ggplot(plotdat, aes(x=t,y=e,group=g))+
  #   geom_line(data=C,color=2)+
  #   geom_line(data=C_u,linetype=3,color=2)+
  #   geom_line(data=C_l,linetype=3,color=2)+
  #   
  #   geom_line(data=NC,color=4)+
  #   geom_line(data=NC_u,linetype=3,color=4)+
  #   geom_line(data=NC_l,linetype=3,color=4)+
  #   
  #   # Screen effect segments
  #   # geom_line(data=PT0,linetype=1,color=2)+ # No screen segment
  #   # geom_line(data=PT1,linetype=1,color=2)+ # screen segment
  #   # 
  #   # geom_line(data=ET0,linetype=1,color=2)+ # No screen segment
  #   # geom_line(data=ET1,linetype=1,color=2)+ # screen segment
  #   # 
  #   # geom_line(data=CT0,linetype=1,color=2)+ # No screen segment
  #   # geom_line(data=CT1,linetype=1,color=2)+ # screen segment
  #   
  # # Vertical Line
  # # geom_segment(data=PT1, x=tx, y=0, xend=tx, yend=PT1[1,2], linetype=2, color=1)+
  # # geom_segment(data=ET1, x=tx, y=0, xend=tx, yend=ET1[1,2], linetype=2, color=1)+
  # # geom_segment(data=CT1, x=tx, y=0, xend=tx, yend=CT1[1,2], linetype=2, color=1)+
  # # 
  # # geom_segment(data=PT1, x=15, y=1.5, xend=70, yend=1.5, linetype=4, color=2)+
  # # geom_segment(data=ET1, x=15, y=1.5, xend=70, yend=1.5, linetype=4, color=2)+
  # # geom_segment(data=CT1, x=15, y=0.5, xend=70, yend=0.5, linetype=4, color=2)+
  # # geom_segment(data=CT1, x=15, y=2, xend=70, yend=2, linetype=4, color=2)+
  # # Points
  # # geom_point(x=tx, y=0, shape=1,fill="black") + # right before t_x
  # # geom_point(data=PT1, x=tx, y=PT1[1,2]) + # right after t_x
  # # geom_point(data=ET1, x=tx, y=ET1[1,2]) + # right after t_x
  # # geom_point(data=CT1, x=tx, y=CT1[1,2]) + # right after t_x
  # 
  # # # Arrow heads
  # # annotate("segment",x=69,xend=70,yend=PT1[nrow(PT1),2],y=PT1[nrow(PT1),2],color="red",arrow = arrow(length = unit(0.2, "cm")))+
  # 
  # # Themes
  # scale_x_continuous(name = "Ages in years") +
  #   scale_y_continuous(name = "Penetrance",limits=y.limits) + # ylab
  #   ggtitle(g.title) + #title
  #   theme(plot.title = element_text(hjust = 0.5,size=15), # Centering of the title
  #         axis.text.x = element_text(face="bold",size=12),  # Change axis ticks labels
  #         axis.title.x = element_text(size=15,margin = margin(t = 20)),
  #         axis.title.y = element_text(size=15))
  # 
  # # facet_wrap(~cancer,scales="free")
  # #grid.locator(unit="native") 
  # # grid.brackets(200,153,200,67 , lwd=0.5, h=0.03,type=1,col="brown")
  # # grid.brackets(200,288,200,202, lwd=0.5, h=0.03,type=1,col="brown")
  # # grid.brackets(200,397,200,308, lwd=0.5, h=0.03,type=1,col="brown")
  # # grid.brackets(200,424,200,397, lwd=0.4, h=0.03,type=1,col="brown")
}
Computation_SE=function(Var,output,mat,G,type){
  penSE=vector()
  for(i in 1:nrow(mat)){
    #print(mat[i,])
    FD1=jacobian(proband_penf_ED_jacob, output, v=mat[i,], type=type, affp=TRUE,G=G)
    f1=FD1%*%Var%*%t(FD1)
    f1=try(sqrt(f1),silent = TRUE)
    if(inherits(f1 ,'try-error')){
      f1=NA_real_
    }
    penSE[i]=f1
  }
  return(penSE)
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
proband_penf_ED_jacob <- function(output, v, affp, type,G=1){
  parms=tran_parm(output)
  cest=parms[[1]]
  phi=parms[[2]]
  fp=parms[[3]]
  frailty=TRUE
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
tran_parm=function(theta){
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
ZED <- function(t,u,phi){
  exp(-(t-u)*phi)
}
carrierprobgeno <-function(method="data", data, mode="dominant", q=0.0021){
  carrp <- data$mgene
  id.na<-data$indID[is.na(data$mgene)]
  
  
  mut.ca <- data$relation[data$mgene==1 & !is.na(data$mgene)]
  cfam.id <- data$famID[data$proband==1 & data$mgene==1 & !is.na(data$mgene)]
  nfam.id <- data$famID[data$proband==1 & data$mgene==0 & !is.na(data$mgene)]
  na.fam.id <- data$famID[data$proband==1 & !is.na(data$mgene)]
  i.cfam <- is.element(data$famID,cfam.id)
  i.nfam <- is.element(data$famID,nfam.id)
  
  if(method=="data"){
    for(g in unique(data$relation)){
      for(s in c(0,1)){
        for(a in c(0,1,3,4)){
          carrp[is.na(data$mgene) & data$relation==g & data$gender==s & data$status==a] <- mean(data$mgene[!is.na(data$mgene) & data$relation==g & data$gender==s & data$status==a])
          
        }
      }
    }
  }
  
  # if(method=="data"){
  #   for(g in unique(data$relation)){
  #     for(s in c(0,1)){
  #       for(a in c(0,1,3,4)){
  #         for(te in c("noneremoval","ovaryremoval","breastremoval","bothremoval")){
  #         carrp[is.na(data$mgene) & data$relation==g & data$gender==s & data$status==a & data$type==te] <- mean(data$mgene[!is.na(data$mgene) & data$relation==g & data$gender==s & data$status==a & data$type==te])
  #         }
  #       }  
  #     }
  #   }
  # }
  
  
  else if (method=="mendelian"){
    
    if(is.null(q)) {
      q <- ifelse(mode=="recessive", sqrt(mean(data$mgene[data$generation==0], na.rm=T)),
                  1-sqrt(1-mean(data$mgene[data$generation==0], na.rm=T)) )
      cat("Estimate allele frequency = ", q, "\n")
    }
    else if(q>1 | q<0) stop("The allele frequency (q) should lie between 0 and 1.")
    
    #G1=mutation status of proband
    #G2=mutation status of sib
    #G3=mutation status of child
    #G4=mutation status of parent
    #G5=mutation status of sib's child
    
    if(mode=="dominant"){
      # P1.s1=P(1S0=(AA,Aa,aa) | 0S0=1) = P(parent/offsprint=(AA,Aa,aa) |proband is carrier); 
      # P1.s0=P(1S0=(AA,Aa,aa) | 0S0=0)
      P1.s1 <- c(q/(2-q), (1-q^2)/(2-q), (1-q)^2/(2-q))
      P1.s0 <- c(0, q, 1-q)
      
      # P(0F0=(AA,Aa,aa) | 0S0=1) = P(sib =(AA,Aa,aa) | proband is carrier)
      # P(0F0=(AA,Aa,aa) | 0S0=0)
      P2.s1 <- c((3+2*q-q^2)*q/4/(2-q), (2+3*q-q^2)*(1-q)/2/(2-q), (1-q)^2*(4-q)/4/(2-q) )
      P2.s0 <- c(q^2/4, q*(2-q)/2, (2-q)^2/4 )
      
      # P(2S0=(AA,Aa,aa)|0S0=1) = P(grandparent/grandchild =(AA,Aa,aa) | proband is carrier)
      # this prob is same for niece/nephew, aunt/uncle
      P3.s1 <- c((1+2*q-q^2)*q/2/(2-q), (1+5*q-2*q^2)*(1-q)/2/(2-q), (1-q)^2*(3-q)/2/(2-q) )
      P3.s0 <- c(q^2/2, q*(3-2*q)/2, (1-q)*(2-q)/2 )
      
      GP.AA <- c(q^2+(1-q)*q/2, q^2/2+q/4, q^2/2)             #P(Grandparent/grandchild=AA|0S0=(AA,Aa,aa))
      GP.Aa <- c(q*(1-q)+(1-q)/2, (1-q)*q+1/4, q/2 + q*(1-q)) #P(Grandparent/grandchild=Aa|0S0=(AA,Aa,aa))
      GP.aa <- c((1-q)^2/2, (3-2*q)*(1-q)/4, (1-q)*(2-q)/2)   #P(Grandparent/grandchild=aa|0S0=(AA,Aa,aa))
      
      P4.s1 <- c(sum(GP.AA*P1.s1), sum(GP.Aa*P1.s1), sum(GP.aa*P1.s1)) # P(cousin=(AA,Aa,aa) | 0S0=1)
      P4.s0 <- c(sum(GP.AA*P1.s0), sum(GP.Aa*P1.s0), sum(GP.aa*P1.s0)) # P(cousin=(AA,Aa,aa) | 0S0=0)
      
      p1 <- (1+q-q^2)/(2-q) #p1=sum(P1.s1[1:2]) # parent
      p2 <- (4+5*q-6*q^2+q^3)/(8-4*q) #p2=sum(P2.s1[1:2]) # sib
      p3 <- (1+5*q-5*q^2+q^3)/2/(2-q) #p4=sum(P3.s1[1:2]) # uncle,aunt, niece, nephew, grandparent, grandchild all have same carrier prob
      p4 <- sum(P4.s1[1:2]) # cousin (1F1)
      
      i.rel <- data$relation
      
    }
    else if(mode=="recessive"){
      # recessive
      p1 <- q # P(G3=1|G1=1) = P(G4=1|G1=1)
      p2 <- (1+q)^2/4 # P(G2=1|G1=1)
      p4 <- p1*p2+(1-q)*(1-p2) #P(G5=1|G1=1) = P(G5=1|G2=1)P(G2=1|G1=1) + P(G5=1|G2=0)P(G2=0|G1=1)
      p3 <- (3-2*q-q^2)/4 #P(G2=1|G1=0)
      p0 <- (1+q)^2/4
      p5 <- p1*p3+(1-q)*p0  #P(G5=1|G1=0) = P(G5=1|G2=1)P(G2=1|G1=0) + P(G5=1|G2=0)P(G2=0|G1=0)
    }
    else stop("Unrecognized inheritance mode")
    
    
    # carrp[ is.na(carrp) & (i.rel=="1S0"|i.rel=="0S1")]  <- p1 #parent/offspring
    # carrp[ is.na(carrp) & (i.rel=="0F0")]  <- p2 # sibs
    # carrp[ is.na(carrp) & (i.rel=="0F1" | i.rel=="1F0" | i.rel == "2S0" | i.rel == "0S2")]  <- p3 #gparent/gchild,uncle/aunt,niece/nephew
    # carrp[ is.na(carrp) & (i.rel=="1F1")] <- p4 #cousin
    # carrp[ is.na(carrp) & (i.rel=="NBS" | i.rel=="NSB") ]  <- q*(2-q)
    
    
    carrp[ is.na(carrp) & (i.rel==4|i.rel==3)]  <- p1 #parent/offspring
    carrp[ is.na(carrp) & (i.rel==2)]  <- p2 # sibs
    carrp[ is.na(carrp) & (i.rel==5 | i.rel=="1F0" | i.rel == "2S0" | i.rel == "0S2")]  <- p3 #gparent/gchild,uncle/aunt,niece/nephew
    carrp[ is.na(carrp) & (i.rel=="1F1")] <- p4 #cousin
    carrp[ is.na(carrp) & (i.rel==6 | i.rel==7) ]  <- q*(2-q)
    
    # 
    # j=0
    # imput<-matrix(0,nrow=length(id.na), ncol=3)
    # 
    # for(i in id.na){
    #   i.rel <-data$relation[data$indID==i]
    #   i.fam <- data$famID[data$indID==i]
    #   fam1 <- data[data$famID==i.fam,] # fmily members of the i
    #   
    #   faid <- fam1$fatherID[fam1$indID==i] # father id
    #   moid <- fam1$motherID[fam1$indID==i] # mother id
    #   
    #   chid <- fam1$indID[fam1$fatherID==i | fam1$motherID==i] # children's ids
    #   sbid <- fam1$indID[fam1$fatherID==faid | fam1$motherID==moid] # sib's ids
    #   sbid <- setdiff(sbid, i)
    #   
    #   fag <- ifelse(faid==0, NA, fam1$mgene[fam1$indID==faid]) # father's mutation status
    #   mog <- ifelse(moid==0, NA, fam1$mgene[fam1$indID==moid]) # mother's mutation status
    #   
    #   #if( is.parent(i, faid, moid))  p <- p1
    #   #else if( is.sibling(faid, moid)) p <- p2
    #   #else if( is.niece(i.rel, mut.ca)) p <- p4
    #   #else if( is.grandson(i.rel, mut.ca)) p <- p6
    #   #else if( is.cousin(i.rel, mut.ca)) p <- p8
    #   #else if(i.rel=="NBS") p <- q*(2-q)
    #   
    #   # 1=self, 2=sib, 3=child, 4=parent, 5=sibchild, 6=hus,7=sibspouse
    #   
    #   if(sum(fam1$mgene[is.element(fam1$indID, c(faid,moid,chid))], na.rm=TRUE) > 0) p <- p1
    #   else if(sum(fam1$mgene[is.element(fam1$indID, sbid)], na.rm=TRUE) > 0 ) p <- p2 #P(G2=1|G2=1)
    #   else if(sum(fam1$mgene[ ((i.rel==3) | (i.rel == 5)) & (fam1$generation==1)], na.rm=TRUE) > 0 ) p <- p6 #P(G4=1|G3=1) grandparent is carrier
    #   else if(sum(fam1$mgene[ ((i.rel==3) | (i.rel == 5)) & (fam1$generation==2)], na.rm=TRUE) > 0 ) p <- p4 #P(G3=1|G2=1) aunt or niece is carrier
    #   else if(sum(fam1$mgene[ ((i.rel==3) | (i.rel == 5)) ], na.rm=TRUE) > 0) p <- p8 # cousin is carrier
    #   #else print(paste("no mutation carriers in family",i.fam))
    #   
    #   carrp[data$indID==i] <- p
    
    #} #close for for(i in id.na)
  } # close for else if(method=="mendelian")
  else stop("Unrecognized method")
  
  carrp[is.na(carrp)] <- 0
  data$carrp.geno <- carrp
  return(data)
}
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
impute=function(dat){
  #set.seed(7)
  #set.seed(25)
  set.seed(78)
  for(i in 1:nrow(dat)){
    if(is.na(dat$mgene[i]))dat[i,"mgene"]<-sample(c(0,1),1,prob=c(1-dat$carrp.geno[i],dat$carrp.geno[i]),replace=TRUE)
  }
  return(dat)
}  