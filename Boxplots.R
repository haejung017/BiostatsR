library(ggplot2)
library(latex2exp)
# result #######################################################################################
load_output=function(Fsize){
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/PE_indep"))
  PE_i<<-readRDS("result_mat.rds")
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/PE_medium"))
  PE_m<<-readRDS("result_mat.rds")
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/PE_high"))
  PE_h<<-readRDS("result_mat.rds")
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/ED_indep"))
  ED_i<<-readRDS("result_mat.rds")
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/ED_medium"))
  ED_m<<-readRDS("result_mat.rds")
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/ED_high"))
  ED_h<<-readRDS("result_mat.rds")
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/CO_indep"))
  CO_i<<-readRDS("result_mat.rds")
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/CO_medium"))
  CO_m<<-readRDS("result_mat.rds")
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/CO_high"))
  CO_h<<-readRDS("result_mat.rds")
}
load_variance=function(Fsize,Varmethod){
  #Varmethod choose between "RobustSE_Deltamethod_Variance","HessianSE_MonteCarlo_Variance"
  setwd(paste0("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/??? ??????/",Fsize,"/",Varmethod))
  PE_i_var<<-readRDS("PE_i_var.rds")
  PE_m_var<<-readRDS("PE_m_var.rds")
  PE_h_var<<-readRDS("PE_h_var.rds")
  ED_i_var<<-readRDS("ED_i_var.rds")
  ED_m_var<<-readRDS("ED_m_var.rds")
  ED_h_var<<-readRDS("ED_h_var.rds")
  CO_i_var<<-readRDS("CO_i_var.rds")
  CO_m_var<<-readRDS("CO_m_var.rds")
  CO_h_var<<-readRDS("CO_h_var.rds")  
}
################################################################################################
# Covariate Coefficients
simbox=function(estvec1,estvec2,estvec3,type){
  if(type=="PE"){
    trueparm=c(0.008,2.405,0.007,3.080,0.668,1.952,1.194,10,2.900)
    ori.trueparm = trueparm
    tra.trueparm = trueparm
    tra.trueparm[c(1:4,8:9)] = sapply(c(1:4,8:9), function(s) log(trueparm[s]))
  }else if(type=="ED"){
    trueparm=c(0.008,2.300,0.007,2.932,1.872,1.858,1.224,10,3.240,0.278)
    ori.trueparm = trueparm
    tra.trueparm = trueparm
    tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
  }else{
    trueparm=c(0.008,2.329,0.007,2.906,3.401,2.078,1.566, 10,3.529, 3.530,0.160) 
    ori.trueparm = trueparm
    tra.trueparm = trueparm
    tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
  }
  #estvec=cbind(apply(estvec,2, function(x) Trim(x, trim=0.025,na.rm=TRUE)))
  B=nrow(estvec1)
  #Coefficient diagnostic
  tra.truemat = matrix(rep(tra.trueparm,B), ncol=length(tra.trueparm), nrow=B, byrow = TRUE)
  # sd1=matrix(apply(estvec1,2, sd), ncol=length(tra.trueparm),nrow=B, byrow=TRUE)
  # sd2=matrix(apply(estvec2,2, sd), ncol=length(tra.trueparm),nrow=B, byrow=TRUE)
  # sd3=matrix(apply(estvec3,2, sd), ncol=length(tra.trueparm),nrow=B, byrow=TRUE)
  Bias1=estvec1-tra.truemat
  #Bias1=Bias1/tra.truemat # Relative Bias
  #Bias1=Bias1/sd1 # Normalized Bias
  Bias2=estvec2-tra.truemat
  #Bias2=Bias2/tra.truemat
  #Bias2=Bias2/sd2
  Bias3=estvec3-tra.truemat
  #Bias3=Bias3/tra.truemat
  #Bias3=Bias3/sd3
  
  if(type=="PE"){
    index=c(5,6,7)
    covnames=c('1'=parse(text = TeX('$\\beta_s$')),
               '2'=parse(text = TeX('$\\beta_{g1}$')),
               '3'=parse(text= TeX('$\\beta_{g2}$')))
    }else if(type=="ED"){
    index=c(5,6,7,10)
    covnames=c('1'=parse(text = TeX('$\\beta_s$')),
               '2'=parse(text = TeX('$\\beta_{g1}$')),
               '3'=parse(text= TeX('$\\beta_{g2}$')),
               '4'=parse(text= TeX('log($\\eta$)')))  
    }else if(type=="CO"){ 
      index=c(5,6,7,10,11)
      covnames=c('1'=parse(text = TeX('$\\beta_s$')),
                 '2'=parse(text = TeX('$\\beta_{g1}$')),
                 '3'=parse(text= TeX('$\\beta_{g2}$')),
                 '4'=parse(text= TeX('log($\\eta$)')),
                 '5'=parse(text= TeX('$\\eta_0$')))
    }
  
  boxclear=function(Bias,Boxindex){
    boxdat1=NULL
    for(i in index){
      Bi=cbind(Bias[,i])  
      boxdat1=rbind(boxdat1,Bi)  
    }
    covtype=cbind(rep(1:length(index), each=B)) 
    frailtytype=cbind(rep(Boxindex, B*length(index)))
    boxdat1=cbind(boxdat1,covtype,frailtytype)
    boxdat1
  }
  boxdat1=boxclear(Bias1,1)
  boxdat2=boxclear(Bias2,2)
  boxdat3=boxclear(Bias3,3)
  
  boxdat=as.data.frame(rbind(boxdat1,boxdat2,boxdat3))
  boxdat$V2=as.character(boxdat$V2)
  boxdat$V3=factor(boxdat$V3, labels=c("Low familial dependence","Medium familial dependence","High familial dependence"))
  
  
  a = ggplot(boxdat, aes(V2, V1))+
    geom_boxplot(fill="#4271AE", colour="#1F3552", alpha=0.6) + # box color, line colour, transparency
    geom_hline(yintercept=0, linetype="dashed", color = "red") + # add horizontal line
    scale_x_discrete(name = "Covariate Coefficients",# xlab
                     labels= covnames) + # change the name of the axis ticks
    scale_y_continuous(name = "Bias") + # ylab
    ggtitle(paste0("Boxplot of Bias of the Model Parameters",", ",type, " TVC")) + #title
    theme(plot.title = element_text(hjust = 0.5,size=15), # Centering of the title
          axis.text.x = element_text(face="bold",size=15),  # Change axis ticks labels
          axis.title.x = element_text(size=15,margin = margin(t = 20)),
          axis.title.y = element_text(size=15),
          strip.text.x = element_text(size=15))+     
    facet_grid(. ~ V3) # grouping by third variable
  
  
  return(a)   
}
PEplot=simbox(PE_i,PE_m,PE_h,type="PE")
EDplot=simbox(ED_i,ED_m,ED_h,type="ED")
COplot=simbox(CO_i,CO_m,CO_h,type="CO")

jpeg("boxplotPE.jpg", width=10, height=5, units="in", res=500)
PEplot
dev.off()
jpeg("boxplotED.jpg", width=10, height=5, units="in", res=500)
EDplot
dev.off()
jpeg("boxplotCO.jpg", width=10, height=5, units="in", res=500)
COplot
dev.off()


# Frailty parameters
simbox.frailty=function(estvec1,estvec2,estvec3,type){
  if(type=="PE"){
    trueparm1=c(10,2.900)
    trueparm2=c(3.5,2.900)
    trueparm3=c(1,2.900)
  }else if(type=="ED"){
    trueparm1=c(10,3.240)
    trueparm2=c(3.5,3.240)
    trueparm3=c(1,3.240)
  }else{
    trueparm1=c(10,3.529) 
    trueparm2=c(3.5,3.529) 
    trueparm3=c(1,3.529) 
  }
  tra.trueparm1 = trueparm1
  tra.trueparm2 = trueparm2
  tra.trueparm3 = trueparm3
  tra.trueparm1[1:2] = sapply(1:2, function(s) log(trueparm1[s]))
  tra.trueparm2[1:2] = sapply(1:2, function(s) log(trueparm2[s]))
  tra.trueparm3[1:2] = sapply(1:2, function(s) log(trueparm3[s]))
  #estvec=cbind(apply(estvec,2, function(x) Trim(x, trim=0.025,na.rm=TRUE)))
  estvec1=estvec1[,8:9]
  estvec2=estvec2[,8:9]
  estvec3=estvec3[,8:9]
  B=nrow(estvec1)
  #Coefficient diagnostic
  tra.truemat1 = matrix(rep(tra.trueparm1,B), ncol=length(tra.trueparm1), nrow=B, byrow = TRUE)
  tra.truemat2 = matrix(rep(tra.trueparm2,B), ncol=length(tra.trueparm2), nrow=B, byrow = TRUE)
  tra.truemat3 = matrix(rep(tra.trueparm3,B), ncol=length(tra.trueparm3), nrow=B, byrow = TRUE)
  
  Bias1=estvec1-tra.truemat1
  #Bias1=Bias1/tra.truemat # Relative Bias
  #Bias1=Bias1/sd1 # Normalized Bias
  Bias2=estvec2-tra.truemat2
  #Bias2=Bias2/tra.truemat
  #Bias2=Bias2/sd2
  Bias3=estvec3-tra.truemat3
  #Bias3=Bias3/tra.truemat
  #Bias3=Bias3/sd3
  index=1:2
  covnames=c('1'=parse(text= TeX('log($\\k_1$)')),
             '2'=parse(text= TeX('log($\\k_2$)')))  
  
  boxclear=function(Bias,Boxindex){
    boxdat1=NULL
    for(i in index){
      Bi=cbind(Bias[,i])  
      boxdat1=rbind(boxdat1,Bi)  
    }
    covtype=cbind(rep(1:length(index), each=B)) 
    frailtytype=cbind(rep(Boxindex, B*length(index)))
    boxdat1=cbind(boxdat1,covtype,frailtytype)
    boxdat1
  }
  boxdat1=boxclear(Bias1,1)
  boxdat2=boxclear(Bias2,2)
  boxdat3=boxclear(Bias3,3)
  
  boxdat=as.data.frame(rbind(boxdat1,boxdat2,boxdat3))
  boxdat$V2=as.character(boxdat$V2)
  boxdat$V3=factor(boxdat$V3, labels=c("Low familial dependence","Medium familial dependence","High familial dependence"))
  
  
  a = ggplot(boxdat, aes(V2, V1))+
    geom_boxplot(fill="#4271AE", colour="#1F3552", alpha=0.6) + # box color, line colour, transparency
    geom_hline(yintercept=0, linetype="dashed", color = "red") + # add horizontal line
    scale_x_discrete(name = "Frailty Parameters",# xlab
                     labels= covnames) + # change the name of the axis ticks
    scale_y_continuous(name = "Bias") + # ylab
    ggtitle(paste0("Boxplot of Bias of the Model Parameters",", ",type, " TVC")) + #title
    theme(plot.title = element_text(hjust = 0.5,size=15), # Centering of the title
          axis.text.x = element_text(face="bold",size=15),  # Change axis ticks labels
          axis.title.x = element_text(size=15,margin = margin(t = 20)),
          axis.title.y = element_text(size=15),
          strip.text.x = element_text(size=15))+      
    facet_grid(. ~ V3) # grouping by third variable
  
  
  return(a)   
}
PEplot.frailty=simbox.frailty(PE_i,PE_m,PE_h,type="PE")
EDplot.frailty=simbox.frailty(ED_i,ED_m,ED_h,type="ED")
COplot.frailty=simbox.frailty(CO_i,CO_m,CO_h,type="CO")

jpeg("boxplotPEf.jpg", width=10, height=5, units="in", res=500)
PEplot.frailty
dev.off()
jpeg("boxplotEDf.jpg", width=10, height=5, units="in", res=500)
EDplot.frailty
dev.off()
jpeg("boxplotCOf.jpg", width=10, height=5, units="in", res=500)
COplot.frailty
dev.off()


################################################################################################
# Penetrance boxes
source("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection/BaseCodes/TVCtools.r")
# 1000 families ################################################################################
load_output(1000);load_variance(1000,"RobustSE_Deltamethod_Variance")
################################################################################################
# 500 families #################################################################################
load_output(750);load_variance(750,"RobustSE_Deltamethod_Variance")
################################################################################################
# 250 families #################################################################################
load_output(500);load_variance(500,"RobustSE_Deltamethod_Variance")
################################################################################################
simbox.pen=function(){
  truepenPE=truepenED=truepenCO=list()
  for (f in c(10,3.5,1)){
    trueparmPE=c(0.008,2.405,0.007,3.080,0.668,1.952,1.194,f,2.900)
    trueparmED=c(0.008,2.300,0.007,2.932,1.872,1.858,1.224,f,3.240,0.278)
    trueparmCO=c(0.008,2.329,0.007,2.906,3.401,2.078,1.566,f,3.529, 3.530,0.160) 
    if(f==10){i=1}else if(f==3.5){i=2}else{i=3}
    truepenPE[[i]]=Pen(trueparmPE,"PE",ct=1)
    truepenED[[i]]=Pen(trueparmED,"ED",ct=1)
    truepenCO[[i]]=Pen(trueparmCO,"CO",ct=1)
  }  
  
  func1=function(truepenList, family){
    t= seq(from=35,to=70, length.out = 200)
    if(family=="Low familial dependence"){
      z=1
      }
    else if(family=="Medium familial dependence"){
      z=2
    }else{
        z=3
        }
    penNS=truepenList[[z]][[1]]
    penS=truepenList[[z]][[2]]
    family.g=rep(family,400)
    screen.g=rep(c("NS","S"),each=200)
    piece1=as.data.frame(cbind(t,penNS))
    piece2=as.data.frame(cbind(t,penS))
    colnames(piece1)[2]<-"pen";colnames(piece2)[2]<-"pen"
    testp=rbind(piece1,piece2)
    testp=cbind(testp,family.g,screen.g)
    return(testp)
  }
  func2=function(truepenList, tvcmodel){
    plotDat=func1(truepenList,"Low familial dependence")
    plotDat2=func1(truepenList,"Medium familial dependence")
    plotDat3=func1(truepenList,"High familial dependence")
    tvcmodel=rep(tvcmodel,1200)
    testp=rbind(plotDat,plotDat2,plotDat3)
    testp=cbind(testp,tvcmodel)
    return(testp)
  }
  testp1=func2(truepenPE,"PE")
  testp2=func2(truepenED,"ED")
  testp3=func2(truepenCO,"CO")
  testp=rbind(testp1,testp2,testp3)
  
  NS=testp[testp$screen.g=="NS",]
  S=testp[testp$screen.g=="S",]
  
  PEi=PenASE(PE_i,PE_i_var,"PE",1)
  PEm=PenASE(PE_m,PE_m_var,"PE",1)
  PEh=PenASE(PE_h,PE_h_var,"PE",1)
  EDi=PenASE(ED_i,ED_i_var,"ED",1)
  EDm=PenASE(ED_m,ED_m_var,"ED",1)
  EDh=PenASE(ED_h,ED_h_var,"ED",1)
  COi=PenASE(CO_i,CO_i_var,"CO",1)
  COm=PenASE(CO_m,CO_m_var,"CO",1)
  COh=PenASE(CO_h,CO_h_var,"CO",1)
  
  # > sum(is.na(CO_i_var$RSE[,1]))
  # [1] 25
  # > sum(is.na(CO_m_var$RSE[,1]))
  # [1] 26
  # > sum(is.na(CO_h_var$RSE[,1]))
  # [1] 21
  
  biasbar=function(obj,f,tt,ind){
    if(ind[1]==2){t=c(40,50,60,70)}
    else{t=c(39.5,49.5,59.5,69.5)}
    biasbarS=data.frame(t=t,
                        pen=obj$bias[ind],
                        lower=obj$bias[ind]-1.96*obj$pse1[ind],
                        upper=obj$bias[ind]+1.96*obj$pse1[ind],
                        family.g=f,tvcmodel=tt)
    label=ifelse(f=="indep","Low familial dependence",
                 ifelse(f=="medium","Medium familial dependence",
                        "High familial dependence"))
    biasbarS$family.g=factor(biasbarS$family.g, labels=label)
    biasbarS
  }
  PEib=biasbar(obj=PEi,f="indep",tt="PE",ind=c(2,4,6,8))
  PEmb=biasbar(obj=PEm,f="medium",tt="PE",ind=c(2,4,6,8))
  PEhb=biasbar(obj=PEh,f="high",tt="PE",ind=c(2,4,6,8))
  
  EDib=biasbar(obj=EDi,f="indep",tt="ED",ind=c(2,4,6,8))
  EDmb=biasbar(obj=EDm,f="medium",tt="ED",ind=c(2,4,6,8))
  EDhb=biasbar(obj=EDh,f="high",tt="ED",ind=c(2,4,6,8))
  
  COib=biasbar(obj=COi,f="indep",tt="CO",ind=c(2,4,6,8))
  COmb=biasbar(obj=COm,f="medium",tt="CO",ind=c(2,4,6,8))
  COhb=biasbar(obj=COh,f="high",tt="CO",ind=c(2,4,6,8))
  
  PEib2=biasbar(obj=PEi,f="indep",tt="PE",ind=c(1,3,5,7))
  PEmb2=biasbar(obj=PEm,f="medium",tt="PE",ind=c(1,3,5,7))
  PEhb2=biasbar(obj=PEh,f="high",tt="PE",ind=c(1,3,5,7))
  
  EDib2=biasbar(obj=EDi,f="indep",tt="ED",ind=c(1,3,5,7))
  EDmb2=biasbar(obj=EDm,f="medium",tt="ED",ind=c(1,3,5,7))
  EDhb2=biasbar(obj=EDh,f="high",tt="ED",ind=c(1,3,5,7))
  
  COib2=biasbar(obj=COi,f="indep",tt="CO",ind=c(1,3,5,7))
  COmb2=biasbar(obj=COm,f="medium",tt="CO",ind=c(1,3,5,7))
  COhb2=biasbar(obj=COh,f="high",tt="CO",ind=c(1,3,5,7))
  
  #testp$family.g=factor(testp$family.g, labels=c("Low familial dependence","Medium familial dependence","High familial dependence"))
  
  a= ggplot(testp, aes(x=t, y=pen))+
    geom_line(data=S, color=2)+
    geom_line(data=NS, color=4)+
    
    #geom_line(data=PEib, color=2)+
    geom_point(data=PEib,shape=1)+
    geom_point(data=PEmb,shape=1)+
    geom_point(data=PEhb,shape=1)+
    geom_point(data=EDib,shape=1)+
    geom_point(data=EDmb,shape=1)+
    geom_point(data=EDhb,shape=1)+
    geom_point(data=COib,shape=1)+
    geom_point(data=COmb,shape=1)+
    geom_point(data=COhb,shape=1)+
    
    geom_point(data=PEib2,shape=1)+
    geom_point(data=PEmb2,shape=1)+
    geom_point(data=PEhb2,shape=1)+
    geom_point(data=EDib2,shape=1)+
    geom_point(data=EDmb2,shape=1)+
    geom_point(data=EDhb2,shape=1)+
    geom_point(data=COib2,shape=1)+
    geom_point(data=COmb2,shape=1)+
    geom_point(data=COhb2,shape=1)+
    
    geom_errorbar(data=PEib,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+ 
    geom_errorbar(data=PEmb,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+ 
    geom_errorbar(data=PEhb,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+ 
    geom_errorbar(data=EDib,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+ 
    geom_errorbar(data=EDmb,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+ 
    geom_errorbar(data=EDhb,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+
    geom_errorbar(data=COib,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+ 
    geom_errorbar(data=COmb,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+ 
    geom_errorbar(data=COhb,aes(x=t,ymin=lower,ymax=upper),width=1, color="red")+
    
    geom_errorbar(data=PEib2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+ 
    geom_errorbar(data=PEmb2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+ 
    geom_errorbar(data=PEhb2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+ 
    geom_errorbar(data=EDib2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+ 
    geom_errorbar(data=EDmb2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+ 
    geom_errorbar(data=EDhb2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+
    geom_errorbar(data=COib2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+ 
    geom_errorbar(data=COmb2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+ 
    geom_errorbar(data=COhb2,aes(x=t,ymin=lower,ymax=upper),width=1, color="blue")+
  
    scale_x_continuous(name = "Ages in years") + # change the name of the axis ticks
    scale_y_continuous(name = "Penetrance") + # ylab
    ggtitle(paste0("Carrier Penetrance")) + #title
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size=15), # Centering of the title
          axis.title.x = element_text(size=15,margin = margin(t = 20)),
          axis.title.y = element_text(size=20),
          strip.text.x = element_text(size=15))+   
    facet_grid(tvcmodel~family.g)
return(a)
}


a=simbox.pen()
jpeg("pen1000.jpg", width=10, height=10, units="in", res=250)
a
dev.off()
b=simbox.pen()
jpeg("pen500.jpg", width=10, height=10, units="in", res=250)
b
dev.off()
c=simbox.pen()
jpeg("pen250.jpg", width=10, height=10, units="in", res=250)
c
dev.off()

simbox.cov=function(ind){
  
  if(!(ind%in%c(8,9))){
  f=7
  trueparmPE=c(0.008,2.405,0.007,3.080,0.668,1.952,1.194,f,2.900)
  trueparmED=c(0.008,2.300,0.007,2.932,1.872,1.858,1.224,f,3.240,0.278)
  trueparmCO=c(0.008,2.329,0.007,2.906,3.401,2.078,1.566,f,3.529, 3.530,0.160) 
  TruePE=trueparmPE[ind]
  TrueED=trueparmED[ind]
  TrueCO=trueparmCO[ind]
  }
  
pickcov_i=function(ind,nfam){  
  if(ind%in%c(8,9)){
    TruePE=TrueED=TrueCO=log(10)
  }
  betam=mean(PE_i[,ind])-TruePE
  betase=mean(PE_i_var[[1]][,ind])
  betam2=mean(ED_i[,ind])-TrueED
  betase2=mean(ED_i_var[[1]][,ind])
  betam3=mean(CO_i[,ind],na.rm = TRUE)-TrueCO
  betase3=mean(CO_i_var[[1]][,ind],na.rm = TRUE)
  biasbarS=data.frame(tvc="PE",bias=betam,
                      lower=betam-1.96*betase,
                      upper=betam+1.96*betase,
                      family.g="Low familial dependence",nfam=nfam)
  biasbarS2=data.frame(tvc="ED",
                      bias=betam2,
                      lower=betam2-1.96*betase2,
                      upper=betam2+1.96*betase2,
                      family.g="Low familial dependence",nfam=nfam)
  biasbarS3=data.frame(tvc="CO",
                       bias=betam3,
                       lower=betam3-1.96*betase3,
                       upper=betam3+1.96*betase3,
                       family.g="Low familial dependence",nfam=nfam)
  testdat=rbind(biasbarS,biasbarS2,biasbarS3)
  testdat$nfam=factor(testdat$nfam)
  testdat
}
pickcov_m=function(ind,nfam){  
  if(ind%in%c(8,9)){
    TruePE=TrueED=TrueCO=log(3.5)
  }
  betam=mean(PE_m[,ind])-TruePE
  betase=mean(PE_m_var[[1]][,ind])
  betam2=mean(ED_m[,ind])-TrueED
  betase2=mean(ED_m_var[[1]][,ind])
  betam3=mean(CO_m[,ind],na.rm = TRUE)-TrueCO
  betase3=mean(CO_m_var[[1]][,ind],na.rm = TRUE)
  biasbarS=data.frame(tvc="PE",bias=betam,
                      lower=betam-1.96*betase,
                      upper=betam+1.96*betase,
                      family.g="Medium familial dependence",nfam=nfam)
  biasbarS2=data.frame(tvc="ED",
                       bias=betam2,
                       lower=betam2-1.96*betase2,
                       upper=betam2+1.96*betase2,
                       family.g="Medium familial dependence",nfam=nfam)
  biasbarS3=data.frame(tvc="CO",
                       bias=betam3,
                       lower=betam3-1.96*betase3,
                       upper=betam3+1.96*betase3,
                       family.g="Medium familial dependence",nfam=nfam)
  testdat=rbind(biasbarS,biasbarS2,biasbarS3)
  testdat$nfam=factor(testdat$nfam)
  testdat
}
pickcov_h=function(ind,nfam){  
  if(ind%in%c(8,9)){
    TruePE=TrueED=TrueCO=log(1)
  }
  betam=mean(PE_h[,ind])-TruePE
  betase=mean(PE_h_var[[1]][,ind])
  betam2=mean(ED_h[,ind])-TrueED
  betase2=mean(ED_h_var[[1]][,ind])
  betam3=mean(CO_h[,ind],na.rm = TRUE)-TrueCO
  betase3=mean(CO_h_var[[1]][,ind],na.rm = TRUE)
  biasbarS=data.frame(tvc="PE",bias=betam,
                      lower=betam-1.96*betase,
                      upper=betam+1.96*betase,
                      family.g="High familial dependence",nfam=nfam)
  biasbarS2=data.frame(tvc="ED",
                       bias=betam2,
                       lower=betam2-1.96*betase2,
                       upper=betam2+1.96*betase2,
                       family.g="High familial dependence",nfam=nfam)
  biasbarS3=data.frame(tvc="CO",
                       bias=betam3,
                       lower=betam3-1.96*betase3,
                       upper=betam3+1.96*betase3,
                       family.g="High familial dependence",nfam=nfam)
  testdat=rbind(biasbarS,biasbarS2,biasbarS3)
  testdat$nfam=factor(testdat$nfam)
  testdat
}

  load_output(1000);load_variance(1000,"RobustSE_Deltamethod_Variance")
  d1=pickcov_i(ind,"n=1000")
  d2=pickcov_m(ind,"n=1000")
  d3=pickcov_h(ind,"n=1000")
  load_output(750);load_variance(750,"RobustSE_Deltamethod_Variance")
  d4=pickcov_i(ind,"n=500")
  d5=pickcov_m(ind,"n=500")
  d6=pickcov_h(ind,"n=500")
  load_output(500);load_variance(500,"RobustSE_Deltamethod_Variance")
  d7=pickcov_i(ind,"n=250")
  d8=pickcov_m(ind,"n=250")
  d9=pickcov_h(ind,"n=250")
  
  testdat= rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9)
  
a=ggplot(data=testdat,aes(x=tvc,y=bias))+
    geom_point(shape=1)+
    # geom_point(data=biasbarS2, shape=1)+
    # geom_point(data=biasbarS3, shape=1)+
    geom_errorbar(aes(x=tvc,ymin=lower,ymax=upper),width=0.5, color="blue")+
    # geom_errorbar(data=biasbarS,aes(x=tvc,ymin=lower,ymax=upper),width=0.5, color="blue")+
    # geom_errorbar(data=biasbarS2,aes(x=tvc,ymin=lower,ymax=upper),width=0.5, color="blue")+
    # geom_errorbar(data=biasbarS3,aes(x=tvc,ymin=lower,ymax=upper),width=0.5, color="blue")+
    geom_hline(yintercept = 0)+
    scale_x_discrete(name = "TVC type") + # change the name of the axis ticks
    scale_y_continuous(name = "Bias") + # ylab
    #ggtitle(paste0("")) + #title
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size=15), # Centering of the title
          axis.title.x = element_text(size=15,margin = margin(t = 20)),
          axis.title.y = element_text(size=15),
          strip.text.x = element_text(size=15),
          strip.text.y = element_text(size=15))+     
  facet_grid(nfam ~ family.g)

  return(a)
}

a=simbox.cov(5)
setwd("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection")
jpeg("biasS.jpg", width=10, height=10, units="in", res=250)
a
dev.off()

b=simbox.cov(6)
setwd("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection")
jpeg("biasG1.jpg", width=10, height=10, units="in", res=250)
b
dev.off()

c=simbox.cov(7)
setwd("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection")
jpeg("biasG2.jpg", width=10, height=10, units="in", res=250)
c
dev.off()

d=simbox.cov(8)
setwd("C:/Users/Jay/Desktop/latex_thesis_jay/revise/simulationsection")
jpeg("frailtyk1.jpg", width=10, height=10, units="in", res=250)
d
dev.off()


gridExtra::grid.arrange(a,b,c, nrow=1)


#####
penf <- function(f, X, est, frailty,t,ct){  
  integrate(integrand, X=X, est=est, frailty=frailty, ct=ct,lower=0,upper=t)$value
}
integrand <- function(u, X, est, frailty, ct){
  u <- u
  fp=c(est[8],est[9])
  index=X
  xbeta1 =  index[2]*est[6]
  xbeta2 =  index[2]*est[7]
  
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
  haz2=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
  H1 = - (est[1]*u)^est[2]*exp(xbeta1)
  H2 = - (est[3]*u)^est[4]*exp(xbeta2)
  
  if(ct==1){
    haz=haz1; mo1=-1;mo2=0
  }else{
    haz=haz2; mo1=0;mo2=-1
  }
  
  if(frailty==TRUE){
    f = haz*((1-H1/fp[1])^(-fp[1]+mo1))*((1-H2/fp[2])^(-fp[2]+mo2))
  } else{
    f = haz*exp(H1+H2) # h * S
  }
  return(f)
}
penf.tsc <- function(f, X, est, frailty, tt, t, ct){  
  if(tt=="PE"){cHf=cH_PE;est=c(est,c(0,0))}
  if(tt=="ED"){cHf=cH_ED;est=c(est,0)}
  if(tt=="CO"){cHf=cH_CO}
  integrate(integrand.tsc, X=X, est=est, frailty=frailty, cHf=cHf, ct=ct, lower=20,upper=t)$value
}
integrand.tsc <- function(u, X, est, cHf, frailty, ct){
  u <- u
  fp=c(est[8],est[9])
  index=X
  xbeta1 =  index[2]*est[6]
  xbeta2 =  index[2]*est[7]
  
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)*exp(est[5]*exp(-est[10]*(u-20)) + est[11])
  haz2=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
  H1 = - (sapply(u,cHf,a=20,parm=c(est[1:2],est[5],est[10],est[11]))+(est[1]*20)^est[2])*exp(xbeta1)
  H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
  
  if(ct==1){
    haz=haz1; mo1=-1;mo2=0
  }else{
    haz=haz2; mo1=0;mo2=-1
  }
  
  if(frailty==TRUE){
    f = haz*((1-H1/fp[1])^(-fp[1]+mo1))*((1-H2/fp[2])^(-fp[2]+mo2))
  } else{
    f = haz*exp(H1+H2) # h * S
  }
  return(f)
}
pen4comb=function(X,t,est,tt,ct){
  if(X[1]==0 | t<20){
    pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t,ct=ct)
  }else if(X[1]==1 & t>=20){
    pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20,ct=ct)+
      penf.tsc(integrand.tsc,X=X,est=est,tt=tt,frailty=TRUE,t=t,ct=ct)
  }
  pen
}
Pen=function(trueparm,tt,ct){
  index=rbind(c(0,1),c(1,1))
  t= seq(from=20,to=55, length.out = 200) 
  truepen=list()
  for(i in 1:2){ #iterate for S,G combinations 
    #Calculate true penetrance for given SG 
    truepen[[i]] = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=trueparm,tt=tt,ct=ct))
  }
  return(truepen)
}
PenASE=function(estvec,estvar,tt,ct,si){
  index=rbind(c(0,1),c(1,1))
  ev=estvec
  ev[ev[,8]>=20,8] = 20
  ev[ev[,9]>=20,9] = 20
  if(tt=="PE"){
    ev[,c(1:4,8:9)] <- apply(ev[,c(1:4,8:9)], 2, exp)
  }else{
    ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
  }
  
  
  total=list()
  for(i in 1:2){ #iterate for S,G combinations 
    l=list()
    estpen=NULL
    for(s in 1:nrow(ev)){
      evs = ev[s,]
      if(index[i,1]==0){
        t= c(24.5,34.5,44.5,54.5)
        #t=seq(from=20,to=55, length.out = 50)
      }else{
        t= c(25,35,45,55)
        #t=seq(from=20,to=55, length.out = 50)
      }
      r = try(sapply(t, function(t) pen4comb(X=index[i,],t=t,est=evs,tt=tt,ct=ct)),silent=TRUE)
      if(inherits(r ,'try-error')){
        warning(as.vector(r))
        r=rep(NA_real_,4)
      }
      estpen =rbind(estpen,r)
    }
    print(sum(is.na(estpen[,1])))
    #calculate bias for i
    l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
    Bias = unlist(lapply(l,mean,na.rm=TRUE))
    total[[i]] = cbind(Bias)
  }
  pen40=pen50=pen60=pen70=NULL
  for(i in 1:2){ #i=1,i=2
    pen40=rbind(pen40,total[[i]][1,])
    pen50=rbind(pen50,total[[i]][2,])
    pen60=rbind(pen60,total[[i]][3,])
    pen70=rbind(pen70,total[[i]][4,])
  }
  out=rbind(pen40,pen50,pen60,pen70)  
  
  ASE1=colMeans(estvar$PSE1,na.rm = TRUE)[c(3:4,7:8,11:12,15:16)]
  ASE2=colMeans(estvar$PSE2,na.rm=TRUE)[c(3:4,7:8,11:12,15:16)]
  
  
  
  return(list(bias=out,pse1=ASE1,pse2=ASE2))
}


printresult=function(estvec,estvar,type,dep){
  library(DescTools)
  colmeans=function(vec){
    apply(vec,2,function(x)mean(Trim(x, trim=0.025,na.rm=TRUE)))
  }
  colsds=function(vec){
    apply(vec,2,function(x)sd(Trim(x, trim=0.025,na.rm=TRUE)))
  }
  if(type=="PE"){
    trueparm=c(0.008,2.405,0.007,3.080,0.668,1.952,1.194,dep,2.900)
    Diagnostic_PE= function(estvec, trueparm){
      
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:9)] = sapply(c(1:4,8:9), function(s) log(trueparm[s]))
      
      # Coefficient diagnostic
      tra.truemat = matrix(rep(tra.trueparm,B), ncol=length(tra.trueparm), nrow=B, byrow = TRUE)
      #Mean=colmeans(estvec)
      Mean=colMeans(estvec)
      Bias=Mean-tra.trueparm
      ESE= apply(estvec, 2, sd)
      #ESE= colsds(estvec)
      
      # pentrance diagnostic
      penBiasSE = function(estvec,tp){
        ev=estvec
        ev[,c(1:4,8:9)] <- apply(ev[,c(1:4,8:9)], 2, exp)
        trueparm=tp
        index = rbind(c(0,1),c(1,1))
        
        t= c(25,35,45,55)
        total=list()
        for(i in 1:2){ #iterate for S,G combinations 
          #Calculate true penetrance for given SG 
          truepen = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=trueparm))
          #Calculate estimated penetrnace for given t, SG
          l=list()
          estpen=NULL
          for(s in 1:nrow(ev)){
            evs = ev[s,]
            r = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=evs))
            estpen =rbind(estpen,r)
          }
          
          #calculate bias for i
          l[[1]]=estpen[,1];l[[2]]=estpen[,2]
          #print(l)
          
          Bias = unlist(lapply(l,mean))-truepen
          total[[i]] = cbind(truepen,Bias)
        }
        
        pen40=pen50=pen60=pen70=NULL
        for(i in 1:2){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
          pen40=rbind(pen40,total[[i]][1,])
          pen50=rbind(pen50,total[[i]][2,])
          pen60=rbind(pen60,total[[i]][3,])
          pen70=rbind(pen70,total[[i]][4,])
        }
        out=rbind(pen40,pen50,pen60,pen70)  
        #round(out,4)
      }
      meansd=cbind(tra.trueparm,Bias,ESE) 
      penbias=penBiasSE(estvec=estvec,tp=ori.trueparm)
      #print(round(meansd,4))
      #print(round(penbias,4))
      return(list(meansd=meansd, penbias=penbias))
    }
    Diagnostic_PE_E2= function(estvec, trueparm){
      
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:9)] = sapply(c(1:4,8:9), function(s) log(trueparm[s]))
      
      # Coefficient diagnostic
      tra.truemat = matrix(rep(tra.trueparm,B), ncol=length(tra.trueparm), nrow=B, byrow = TRUE)
      #Mean=colMeans(estvec)
      Mean=colmeans(estvec)
      Bias=Mean-tra.trueparm
      PBias=(Bias/tra.trueparm)*100
      #ESE= apply(estvec, 2, sd)
      ESE= colsds(estvec)
      #RMSE= sqrt( colMeans( (estvec-tra.truemat)^2 ))
      
      # pentrance diagnostic
      penBiasSE = function(estvec,tp){
        ev=estvec
        ev[ev[,8]>=20,8] = 20
        ev[ev[,9]>=20,9] = 20
        ev[,c(1:4,8:9)] <- apply(ev[,c(1:4,8:9)], 2, exp)
        trueparm=tp
        index = as.matrix(expand.grid(c(0,1),c(0,1)))
        #ev=test8$simresult
        
        penf <- function(f, X, est, frailty,  t){  
          integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
        }
        integrand <- function(u, X, est, frailty){
          
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
          H1 = - (est[1]*u)^est[2]*exp(xbeta1)
          H2 = - (est[3]*u)^est[4]*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        penf.tsc <- function(f, X, est, frailty,t){  
          integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
        }
        integrand.tsc <- function(u, X, est, type, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
          H1 = - (sapply(u,cH_PE,a=20,parm=c(est[1:2],est[5]))+(est[1]*20)^est[2])*exp(xbeta1)
          H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        pen4comb=function(X,t,est){
          if(X[1]==0){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
          }else if(X[1]==1){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
              penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
          }
          pen
        }
        #t=45
        t= c(25,35,45,55)
        total=list()
        for(i in 1:4){ #iterate for S,G combinations 
          #Calculate true penetrance for given SG 
          truepen = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=trueparm))
          #Calculate estimated penetrnace for given t, SG
          l=list()
          estpen=NULL
          for(s in 1:nrow(ev)){
            evs = ev[s,]
            r = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=evs))
            estpen =rbind(estpen,r)
          }
          
          #calculate bias for i
          l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
          Bias = unlist(lapply(l,mean))-truepen
          EmpSE = unlist(lapply(l,sd))
          #RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2),mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
          total[[i]] = cbind(truepen,Bias,EmpSE)#,RMSE)
        }
        
        pen40=pen50=pen60=pen70=NULL
        for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
          pen40=rbind(pen40,total[[i]][1,])
          pen50=rbind(pen50,total[[i]][2,])
          pen60=rbind(pen60,total[[i]][3,])
          pen70=rbind(pen70,total[[i]][4,])
        }
        out=rbind(pen40,pen50,pen60,pen70)  
        round(out,4)
      }
      
      #result
      #meansd=cbind(ori.trueparm,tra.trueparm,Mean,Bias,PBias,ESE,RMSE) 
      meansd=cbind(tra.trueparm,Bias,ESE) 
      penbias=penBiasSE(estvec=estvec,tp=ori.trueparm)
      #print(round(meansd,4))
      #print(round(penbias,4))
      return(list(meansd=meansd, penbias=penbias))
    }
    CoverageP=function(estvec,estse,trueparm){
      v1=estvec
      v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:9)] = sapply(c(1:4,8:9), function(s) log(trueparm[s]))
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96 
      sapply(1:length(trueparm), function(s) sum(lo[,s]<=tra.trueparm[s]&up[,s]>=tra.trueparm[s],na.rm = TRUE) )/B
    }
    CoverageP1=function(estvec,estse,trueparm){
      #v1=estvec
      #v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      
      ev=estvec
      ev[ev[,8]>=20,8] = 20
      ev[ev[,9]>=20,9] = 20
      ev[,c(1:4,8:9)] <- apply(ev[,c(1:4,8:9)], 2, exp)
      #ev=test8$simresult
      
      
      penf <- function(f, X, est, frailty,  t){  
        integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
      }
      integrand <- function(u, X, est, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
        H1 = - (est[1]*u)^est[2]*exp(xbeta1)
        H2 = - (est[3]*u)^est[4]*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      penf.tsc <- function(f, X, est, frailty,t){  
        integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
      }
      integrand.tsc <- function(u, X, est, type, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)*exp(est[5])
        H1 = - (sapply(u,cH_PE,a=20,parm=c(est[1:2],est[5]))+(est[1]*20)^est[2])*exp(xbeta1)
        H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      pen4comb=function(X,t,est){
        if(X[1]==0){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
        }else if(X[1]==1){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
            penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
        }
        pen
      }
      #t=45
      
      t= c(25,35,45,55)
      penest=NULL
      for(s in 1:nrow(ev)){
        evs = ev[s,]
        apen=NULL
        for(age in t){
          ap=c(pen4comb(X=c(0,0),t=age,est=evs),
               pen4comb(X=c(1,0),t=age,est=evs),
               pen4comb(X=c(0,1),t=age,est=evs),
               pen4comb(X=c(1,1),t=age,est=evs))
          apen=c(apen,ap)
        }
        penest=rbind(penest,apen)
      }
      
      truepen=NULL
      for(age in t){
        ap=c(pen4comb(X=c(0,0),t=age,est=trueparm),
             pen4comb(X=c(1,0),t=age,est=trueparm),
             pen4comb(X=c(0,1),t=age,est=trueparm),
             pen4comb(X=c(1,1),t=age,est=trueparm))
        truepen=c(truepen,ap)
      }  
      #print(truepen)
      v1=penest
      v2=estse
      
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96
      #print(up)
      #print(lo)
      sapply(1:length(truepen), function(s) sum(lo[,s]<=truepen[s]&up[,s]>=truepen[s],na.rm=TRUE) )/B
    }
    CoverageP2=function(estvec,estse,trueparm){
      #v1=estvec
      #v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      
      ev=estvec
      ev[ev[,8]>=20,8] = 20
      ev[ev[,9]>=20,9] = 20
      
      ev[,c(1:4,8:9)] <- apply(ev[,c(1:4,8:9)], 2, exp)
      #ev=test8$simresult
      
      
      penf <- function(f, X, est, frailty,  t){  
        integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
      }
      integrand <- function(u, X, est, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
        H1 = - (est[1]*u)^est[2]*exp(xbeta1)
        H2 = - (est[3]*u)^est[4]*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      penf.tsc <- function(f, X, est, frailty,t){  
        integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
      }
      integrand.tsc <- function(u, X, est, type, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
        H1 = - (sapply(u,cH_PE,a=20,parm=c(est[1:2],est[5]))+(est[1]*20)^est[2])*exp(xbeta1)
        H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      pen4comb=function(X,t,est){
        if(X[1]==0){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
        }else if(X[1]==1){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
            penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
        }
        pen
      }
      #t=45
      
      t= c(25,35,45,55)
      penest=NULL
      for(s in 1:nrow(ev)){
        evs = ev[s,]
        apen=NULL
        for(age in t){
          ap=c(pen4comb(X=c(0,0),t=age,est=evs),
               pen4comb(X=c(1,0),t=age,est=evs),
               pen4comb(X=c(0,1),t=age,est=evs),
               pen4comb(X=c(1,1),t=age,est=evs))
          apen=c(apen,ap)
        }
        penest=rbind(penest,apen)
      }
      
      truepen=NULL
      for(age in t){
        ap=c(pen4comb(X=c(0,0),t=age,est=trueparm),
             pen4comb(X=c(1,0),t=age,est=trueparm),
             pen4comb(X=c(0,1),t=age,est=trueparm),
             pen4comb(X=c(1,1),t=age,est=trueparm))
        truepen=c(truepen,ap)
      }  
      print(truepen)
      v1=penest
      v2=estse
      
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96
      #print(up)
      #print(lo)
      sapply(1:length(truepen), function(s) sum(lo[,s]<=truepen[s]&up[,s]>=truepen[s],na.rm = TRUE) )/B
    }
    Diagnostic=Diagnostic_PE
    Diagnostic2=Diagnostic_PE_E2
    CoverageP=CoverageP
    CoverageP1=CoverageP1
    CoverageP2=CoverageP2
  }else if(type=="ED"){
    trueparm=c(0.008,2.300,0.007,2.932,1.872,1.858,1.224,dep,3.240,0.278)
    Diagnostic_ED= function(estvec, trueparm){
      
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
      
      # Coefficient diagnostic
      tra.truemat = matrix(rep(tra.trueparm,B), ncol=length(tra.trueparm), nrow=B, byrow = TRUE)
      Mean=colMeans(estvec)
      #Mean=colmeans(estvec)
      Bias=Mean-tra.trueparm
      PBias=(Bias/tra.trueparm)*100
      ESE= apply(estvec, 2, sd)
      #ESE= colsds(estvec)
      #RMSE= sqrt( colMeans( (estvec-tra.truemat)^2 ))
      
      # pentrance diagnostic
      penBiasSE = function(estvec,tp){
        ev=estvec
        ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
        trueparm=tp
        index = as.matrix(expand.grid(c(0,1),c(0,1)))
        #ev=test8$simresult
        
        penf <- function(f, X, est, frailty,  t){  
          integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
        }
        integrand <- function(u, X, est, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
          H1 = - (est[1]*u)^est[2]*exp(xbeta1)
          H2 = - (est[3]*u)^est[4]*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        penf.tsc <- function(f, X, est, frailty,t){  
          integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
        }
        integrand.tsc <- function(u, X, est, type, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)*exp(est[5]*exp(-est[10]*(u-20)))
          H1 = - (sapply(u,cH_ED,a=20,parm=c(est[1:2],est[5],est[10]))+(est[1]*20)^est[2])*exp(xbeta1)
          H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        pen4comb=function(X,t,est){
          if(X[1]==0){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
          }else if(X[1]==1){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
              penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
          }
          pen
        }
        #t=45
        t= c(25,35,45,55)
        total=list()
        for(i in 1:4){ #iterate for S,G combinations 
          #Calculate true penetrance for given SG 
          truepen = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=trueparm))
          #Calculate estimated penetrnace for given t, SG
          l=list()
          estpen=NULL
          for(s in 1:nrow(ev)){
            evs = ev[s,]
            r = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=evs))
            estpen =rbind(estpen,r)
          }
          
          #calculate bias for i
          l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
          Bias = unlist(lapply(l,mean))-truepen
          EmpSE = unlist(lapply(l,sd))
          #RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2) ))#,mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
          total[[i]] = cbind(truepen,Bias,EmpSE)#,RMSE)
        }
        
        pen40=pen50=pen60=pen70=NULL
        for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
          pen40=rbind(pen40,total[[i]][1,])
          pen50=rbind(pen50,total[[i]][2,])
          pen60=rbind(pen60,total[[i]][3,])
          pen70=rbind(pen70,total[[i]][4,])
        }
        out=rbind(pen40,pen50,pen60,pen70)  
        round(out,4)
      }
      
      #result
      #meansd=cbind(ori.trueparm,tra.trueparm,Mean,Bias,PBias,ESE,RMSE)
      meansd=cbind(tra.trueparm,Bias,ESE)
      penbias=penBiasSE(estvec=estvec,tp=ori.trueparm)
      #print(round(meansd,4))
      #print(round(penbias,4))
      return(list(meansd=meansd, penbias=penbias))
    }
    Diagnostic_ED_E2= function(estvec, trueparm){
      
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
      
      # Coefficient diagnostic
      tra.truemat = matrix(rep(tra.trueparm,B), ncol=length(tra.trueparm), nrow=B, byrow = TRUE)
      #Mean=colMeans(estvec)
      Mean=colmeans(estvec)
      Bias=Mean-tra.trueparm
      PBias=(Bias/tra.trueparm)*100
      #ESE= apply(estvec, 2, sd)
      ESE= colsds(estvec)
      #RMSE= sqrt( colMeans( (estvec-tra.truemat)^2 ))
      
      # pentrance diagnostic
      penBiasSE = function(estvec,tp){
        ev=estvec
        ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
        trueparm=tp
        index = as.matrix(expand.grid(c(0,1),c(0,1)))
        #ev=test8$simresult
        
        penf <- function(f, X, est, frailty,  t){  
          integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
        }
        integrand <- function(u, X, est, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
          H1 = - (est[1]*u)^est[2]*exp(xbeta1)
          H2 = - (est[3]*u)^est[4]*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        penf.tsc <- function(f, X, est, frailty,t){  
          integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
        }
        integrand.tsc <- function(u, X, est, type, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
          H1 = - (sapply(u,cH_ED,a=20,parm=c(est[1:2],est[5],est[10]))+(est[1]*20)^est[2])*exp(xbeta1)
          H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        pen4comb=function(X,t,est){
          if(X[1]==0){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
          }else if(X[1]==1){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
              penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
          }
          pen
        }
        #t=45
        t= c(25,35,45,55)
        total=list()
        for(i in 1:4){ #iterate for S,G combinations 
          #Calculate true penetrance for given SG 
          truepen = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=trueparm))
          #Calculate estimated penetrnace for given t, SG
          l=list()
          estpen=NULL
          for(s in 1:nrow(ev)){
            evs = ev[s,]
            r = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=evs))
            estpen =rbind(estpen,r)
          }
          
          #calculate bias for i
          l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
          Bias = unlist(lapply(l,mean))-truepen
          EmpSE = unlist(lapply(l,sd))
          #RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2) ))#,mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
          total[[i]] = cbind(truepen,Bias,EmpSE)#,RMSE)
        }
        
        pen40=pen50=pen60=pen70=NULL
        for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
          pen40=rbind(pen40,total[[i]][1,])
          pen50=rbind(pen50,total[[i]][2,])
          pen60=rbind(pen60,total[[i]][3,])
          pen70=rbind(pen70,total[[i]][4,])
        }
        out=rbind(pen40,pen50,pen60,pen70)  
        round(out,4)
      }
      
      #result
      #meansd=cbind(ori.trueparm,tra.trueparm,Mean,Bias,PBias,ESE,RMSE)
      meansd=cbind(tra.trueparm,Bias,ESE)
      penbias=penBiasSE(estvec=estvec,tp=ori.trueparm)
      #print(round(meansd,4))
      #print(round(penbias,4))
      return(list(meansd=meansd, penbias=penbias))
    }
    CoveragePED=function(estvec,estse,trueparm){
      v1=estvec
      v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96 
      sapply(1:length(trueparm), function(s) sum(lo[,s]<=tra.trueparm[s]&up[,s]>=tra.trueparm[s],na.rm=TRUE))/B
    }
    CoverageP1ED=function(estvec,estse,trueparm){
      #v1=estvec
      #v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      
      ev=estvec
      ev[ev[,8]>=20,8] = 20
      ev[ev[,9]>=20,9] = 20
      ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
      #ev=test8$simresult
      
      
      penf <- function(f, X, est, frailty,  t){  
        integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
      }
      integrand <- function(u, X, est, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
        H1 = - (est[1]*u)^est[2]*exp(xbeta1)
        H2 = - (est[3]*u)^est[4]*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      penf.tsc <- function(f, X, est, frailty,t){  
        integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
      }
      integrand.tsc <- function(u, X, est, type, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)*exp(est[5]*exp(-est[10]*(u-20)))
        H1 = - (sapply(u,cH_ED,a=20,parm=c(est[1:2],est[5],est[10]))+(est[1]*20)^est[2])*exp(xbeta1)
        H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      pen4comb=function(X,t,est){
        if(X[1]==0){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
        }else if(X[1]==1){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
            penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
        }
        pen
      }
      #t=45
      
      t= c(25,35,45,55)
      penest=NULL
      for(s in 1:nrow(ev)){
        evs = ev[s,]
        apen=NULL
        for(age in t){
          ap=c(pen4comb(X=c(0,0),t=age,est=evs),
               pen4comb(X=c(1,0),t=age,est=evs),
               pen4comb(X=c(0,1),t=age,est=evs),
               pen4comb(X=c(1,1),t=age,est=evs))
          apen=c(apen,ap)
        }
        penest=rbind(penest,apen)
      }
      
      truepen=NULL
      for(age in t){
        ap=c(pen4comb(X=c(0,0),t=age,est=trueparm),
             pen4comb(X=c(1,0),t=age,est=trueparm),
             pen4comb(X=c(0,1),t=age,est=trueparm),
             pen4comb(X=c(1,1),t=age,est=trueparm))
        truepen=c(truepen,ap)
      }  
      #print(truepen)
      v1=penest
      v2=estse
      
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96
      #print(up)
      #print(lo)
      sapply(1:length(truepen), function(s) sum(lo[,s]<=truepen[s]&up[,s]>=truepen[s], na.rm = TRUE) )/B
    }
    CoverageP2ED=function(estvec,estse,trueparm){
      #v1=estvec
      #v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      
      ev=estvec
      ev[ev[,8]>=20,8] = 20
      ev[ev[,9]>=20,9] = 20
      ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
      #ev=test8$simresult
      
      
      penf <- function(f, X, est, frailty,  t){  
        integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
      }
      integrand <- function(u, X, est, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
        H1 = - (est[1]*u)^est[2]*exp(xbeta1)
        H2 = - (est[3]*u)^est[4]*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      penf.tsc <- function(f, X, est, frailty,t){  
        integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
      }
      integrand.tsc <- function(u, X, est, type, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
        H1 = - (sapply(u,cH_ED,a=20,parm=c(est[1:2],est[5],est[10]))+(est[1]*20)^est[2])*exp(xbeta1)
        H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      pen4comb=function(X,t,est){
        if(X[1]==0){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
        }else if(X[1]==1){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
            penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
        }
        pen
      }
      #t=45
      
      t= c(25,35,45,55)
      penest=NULL
      for(s in 1:nrow(ev)){
        evs = ev[s,]
        apen=NULL
        for(age in t){
          ap=c(pen4comb(X=c(0,0),t=age,est=evs),
               pen4comb(X=c(1,0),t=age,est=evs),
               pen4comb(X=c(0,1),t=age,est=evs),
               pen4comb(X=c(1,1),t=age,est=evs))
          apen=c(apen,ap)
        }
        penest=rbind(penest,apen)
      }
      
      truepen=NULL
      for(age in t){
        ap=c(pen4comb(X=c(0,0),t=age,est=trueparm),
             pen4comb(X=c(1,0),t=age,est=trueparm),
             pen4comb(X=c(0,1),t=age,est=trueparm),
             pen4comb(X=c(1,1),t=age,est=trueparm))
        truepen=c(truepen,ap)
      }  
      #print(truepen)
      v1=penest
      v2=estse
      
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96
      #print(up)
      #print(lo)
      sapply(1:length(truepen), function(s) sum(lo[,s]<=truepen[s]&up[,s]>=truepen[s],na.rm = TRUE) )/B
    }
    Diagnostic=Diagnostic_ED
    Diagnostic2=Diagnostic_ED_E2
    CoverageP=CoveragePED
    CoverageP1=CoverageP1ED
    CoverageP2=CoverageP2ED
  }else if(type=="CO"){
    trueparm=c(0.008,2.329,0.007,2.906,3.401,2.078,1.566, dep,3.529, 3.530,0.160)  
    Diagnostic_CO= function(estvec, trueparm){
      
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
      
      # Coefficient diagnostic
      tra.truemat = matrix(rep(tra.trueparm,B), ncol=length(tra.trueparm), nrow=B, byrow = TRUE)
      #Mean=colMeans(estvec)
      Mean=colmeans(estvec)
      Bias=Mean-tra.trueparm
      PBias=(Bias/tra.trueparm)*100
      #ESE= apply(estvec, 2, sd)
      ESE= colsds(estvec)
      #RMSE= sqrt( colMeans( (estvec-tra.truemat)^2 ))
      
      # pentrance diagnostic
      penBiasSE = function(estvec,tp){
        ev=estvec
        ev[ev[,8]>=20,8] = 20
        ev[ev[,9]>=20,9] = 20
        ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
        trueparm=tp
        index = as.matrix(expand.grid(c(0,1),c(0,1)))
        #ev=test8$simresult
        
        penf <- function(f, X, est, frailty,  t){  
          integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
        }
        integrand <- function(u, X, est, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
          H1 = - (est[1]*u)^est[2]*exp(xbeta1)
          H2 = - (est[3]*u)^est[4]*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        penf.tsc <- function(f, X, est, frailty,t){  
          integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
        }
        integrand.tsc <- function(u, X, est, type, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)*exp(est[5]*exp(-est[10]*(u-20))+est[11])
          H1 = - (sapply(u,cH_CO,a=20,parm=c(est[1:2],est[5],est[10],est[11]))+(est[1]*20)^est[2])*exp(xbeta1)
          H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        pen4comb=function(X,t,est){
          if(X[1]==0){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
          }else if(X[1]==1){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
              penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
          }
          pen
        }
        #t=45
        t= c(25,35,45,55)
        total=list()
        for(i in 1:4){ #iterate for S,G combinations 
          #Calculate true penetrance for given SG 
          truepen = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=trueparm))
          #Calculate estimated penetrnace for given t, SG
          l=list()
          estpen=NULL
          for(s in 1:nrow(ev)){
            evs = ev[s,]
            #print(s)
            r = try(sapply(t, function(t) pen4comb(X=index[i,],t=t,est=evs)), silent=TRUE)
            if(inherits(r ,'try-error')){
              warning(as.vector(r))
              r=rep(NA_real_, 4)
            }
            estpen =rbind(estpen,r)
          }
          
          #calculate bias for i
          l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
          Bias = unlist(lapply(l,mean, na.rm=TRUE))-truepen
          EmpSE = unlist(lapply(l,sd, na.rm=TRUE))
          #RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2) ))#,mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
          total[[i]] = cbind(truepen,Bias,EmpSE)#,RMSE)
        }
        
        pen40=pen50=pen60=pen70=NULL
        for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
          pen40=rbind(pen40,total[[i]][1,])
          pen50=rbind(pen50,total[[i]][2,])
          pen60=rbind(pen60,total[[i]][3,])
          pen70=rbind(pen70,total[[i]][4,])
        }
        out=rbind(pen40,pen50,pen60,pen70)  
        round(out,4)
      }
      
      #result
      #meansd=cbind(ori.trueparm,tra.trueparm,Mean,Bias,PBias,ESE,RMSE)
      meansd=cbind(tra.trueparm,Bias,ESE)
      penbias=penBiasSE(estvec=estvec,tp=ori.trueparm)
      #print(round(meansd,4))
      #print(round(penbias,4))
      return(list(meansd=meansd, penbias=penbias))
    } 
    Diagnostic_CO_E2= function(estvec, trueparm){
      
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
      
      # Coefficient diagnostic
      tra.truemat = matrix(rep(tra.trueparm,B), ncol=length(tra.trueparm), nrow=B, byrow = TRUE)
      #Mean=colMeans(estvec)
      Mean=colmeans(estvec)
      Bias=Mean-tra.trueparm
      PBias=(Bias/tra.trueparm)*100
      #ESE= apply(estvec, 2, sd)
      ESE= colsds(estvec)
      #RMSE= sqrt( colMeans( (estvec-tra.truemat)^2 ))
      
      # pentrance diagnostic
      penBiasSE = function(estvec,tp){
        ev=estvec
        ev[ev[,8]>=20,8] = 20
        ev[ev[,9]>=20,9] = 20
        ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
        trueparm=tp
        index = as.matrix(expand.grid(c(0,1),c(0,1)))
        #ev=test8$simresult
        
        penf <- function(f, X, est, frailty,  t){  
          integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
        }
        integrand <- function(u, X, est, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
          H1 = - (est[1]*u)^est[2]*exp(xbeta1)
          H2 = - (est[3]*u)^est[4]*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        penf.tsc <- function(f, X, est, frailty,t){  
          integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
        }
        integrand.tsc <- function(u, X, est, type, frailty){
          u <- u
          fp=c(est[8],est[9])
          index=X
          xbeta1 =  index[2]*est[6]
          xbeta2 =  index[2]*est[7]
          
          haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
          H1 = - (sapply(u,cH_CO,a=20,parm=c(est[1:2],est[5],est[10],est[11]))+(est[1]*20)^est[2])*exp(xbeta1)
          H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
          
          if(frailty==TRUE){
            f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
            #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
          } else{
            f = haz1*exp(H1+H2) # h * S
          }
          return(f)
        }
        pen4comb=function(X,t,est){
          if(X[1]==0){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
          }else if(X[1]==1){
            pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
              penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
          }
          pen
        }
        #t=45
        t= c(25,35,45,55)
        total=list()
        for(i in 1:4){ #iterate for S,G combinations 
          #Calculate true penetrance for given SG 
          truepen = sapply(t, function(t) pen4comb(X=index[i,],t=t,est=trueparm))
          #Calculate estimated penetrnace for given t, SG
          l=list()
          estpen=NULL
          for(s in 1:nrow(ev)){
            evs = ev[s,]
            r = try(sapply(t, function(t) pen4comb(X=index[i,],t=t,est=evs)), silent=TRUE)
            if(inherits(r ,'try-error')){
              warning(as.vector(r))
              r=rep(NA_real_, 4)
            }
            estpen =rbind(estpen,r)
          }
          
          #calculate bias for i
          l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
          Bias = unlist(lapply(l,mean, na.rm=TRUE))-truepen
          EmpSE = unlist(lapply(l,sd, na.rm=TRUE))
          #RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2) ))#,mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
          total[[i]] = cbind(truepen,Bias,EmpSE)#,RMSE)
        }
        
        pen40=pen50=pen60=pen70=NULL
        for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
          pen40=rbind(pen40,total[[i]][1,])
          pen50=rbind(pen50,total[[i]][2,])
          pen60=rbind(pen60,total[[i]][3,])
          pen70=rbind(pen70,total[[i]][4,])
        }
        out=rbind(pen40,pen50,pen60,pen70)  
        round(out,4)
      }
      
      #result
      #meansd=cbind(ori.trueparm,tra.trueparm,Mean,Bias,PBias,ESE,RMSE)
      meansd=cbind(tra.trueparm,Bias,ESE)
      penbias=penBiasSE(estvec=estvec,tp=ori.trueparm)
      #print(round(meansd,4))
      #print(round(penbias,4))
      return(list(meansd=meansd, penbias=penbias))
    } 
    CoveragePED=function(estvec,estse,trueparm){
      v1=estvec
      v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      tra.trueparm = trueparm
      tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96 
      sapply(1:length(trueparm), function(s) sum(lo[,s]<=tra.trueparm[s]&up[,s]>=tra.trueparm[s],na.rm=TRUE))/B
    }
    CoverageP1CO=function(estvec,estse,trueparm){
      #v1=estvec
      #v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      
      ev=estvec
      ev[ev[,8]>=20,8] = 20
      ev[ev[,9]>=20,9] = 20
      ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
      #ev=test8$simresult
      
      
      penf <- function(f, X, est, frailty,  t){  
        integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
      }
      integrand <- function(u, X, est, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
        H1 = - (est[1]*u)^est[2]*exp(xbeta1)
        H2 = - (est[3]*u)^est[4]*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      penf.tsc <- function(f, X, est, frailty,t){  
        integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
      }
      integrand.tsc <- function(u, X, est, type, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)*exp(est[5]*exp(-est[10]*(u-20))+est[11])
        H1 = - (sapply(u,cH_CO,a=20,parm=c(est[1:2],est[5],est[10],est[11]))+(est[1]*20)^est[2])*exp(xbeta1)
        H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      pen4comb=function(X,t,est){
        if(X[1]==0){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
        }else if(X[1]==1){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
            penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
        }
        pen
      }
      #t=45
      
      t= c(25,35,45,55)
      penest=NULL
      for(s in 1:nrow(ev)){
        evs = ev[s,]
        apen=NULL
        for(age in t){
          ap=try(c(pen4comb(X=c(0,0),t=age,est=evs),
                   pen4comb(X=c(1,0),t=age,est=evs),
                   pen4comb(X=c(0,1),t=age,est=evs),
                   pen4comb(X=c(1,1),t=age,est=evs)),silent=TRUE)
          if(inherits(ap ,'try-error')){
            warning(as.vector(ap))
            ap=rep(NA_real_, 4)
            #   penSE1=penSE2=rep(NA_real_,16)
          }
          apen=c(apen,ap)
        }
        penest=rbind(penest,apen)
      }
      
      truepen=NULL
      for(age in t){
        ap=c(pen4comb(X=c(0,0),t=age,est=trueparm),
             pen4comb(X=c(1,0),t=age,est=trueparm),
             pen4comb(X=c(0,1),t=age,est=trueparm),
             pen4comb(X=c(1,1),t=age,est=trueparm))
        truepen=c(truepen,ap)
      }  
      #print(truepen)
      v1=penest
      v2=estse
      
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96
      #print(up)
      #print(lo)
      sapply(1:length(truepen), function(s) sum(lo[,s]<=truepen[s]&up[,s]>=truepen[s],na.rm = TRUE) )/B
    }
    CoverageP2CO=function(estvec,estse,trueparm){
      #v1=estvec
      #v2=estse
      trueparm=trueparm
      # Transformation
      B= nrow(estvec)
      ori.trueparm = trueparm
      
      ev=estvec
      ev[ev[,8]>=20,8] = 20
      ev[ev[,9]>=20,9] = 20
      ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
      #ev=test8$simresult
      
      
      penf <- function(f, X, est, frailty,  t){  
        integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
      }
      integrand <- function(u, X, est, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
        H1 = - (est[1]*u)^est[2]*exp(xbeta1)
        H2 = - (est[3]*u)^est[4]*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      penf.tsc <- function(f, X, est, frailty,t){  
        integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=20,upper=t)$value
      }
      integrand.tsc <- function(u, X, est, type, frailty){
        u <- u
        fp=c(est[8],est[9])
        index=X
        xbeta1 =  index[2]*est[6]
        xbeta2 =  index[2]*est[7]
        
        haz1=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)
        H1 = - (sapply(u,cH_CO,a=20,parm=c(est[1:2],est[5],est[10],est[11]))+(est[1]*20)^est[2])*exp(xbeta1)
        H2 = - ((est[3]*u)^est[4]+(est[3]*20)^est[4])*exp(xbeta2)
        
        if(frailty==TRUE){
          f = haz1*((1-H1/fp[1])^(-fp[1]))*((1-H2/fp[2])^(-fp[2]-1))
          #f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
        } else{
          f = haz1*exp(H1+H2) # h * S
        }
        return(f)
      }
      pen4comb=function(X,t,est){
        if(X[1]==0){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=t)
        }else if(X[1]==1){
          pen=penf(integrand,X=X,est=est,frailty=TRUE,t=20)+
            penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
        }
        pen
      }
      #t=45
      
      t= c(25,35,45,55)
      penest=NULL
      for(s in 1:nrow(ev)){
        evs = ev[s,]
        apen=NULL
        for(age in t){
          ap=try(c(pen4comb(X=c(0,0),t=age,est=evs),
                   pen4comb(X=c(1,0),t=age,est=evs),
                   pen4comb(X=c(0,1),t=age,est=evs),
                   pen4comb(X=c(1,1),t=age,est=evs)),silent=TRUE)
          if(inherits(ap ,'try-error')){
            warning(as.vector(ap))
            ap=rep(NA_real_, 4)
            #   penSE1=penSE2=rep(NA_real_,16)
          }
          apen=c(apen,ap)
        }
        penest=rbind(penest,apen)
      }
      
      truepen=NULL
      for(age in t){
        ap=c(pen4comb(X=c(0,0),t=age,est=trueparm),
             pen4comb(X=c(1,0),t=age,est=trueparm),
             pen4comb(X=c(0,1),t=age,est=trueparm),
             pen4comb(X=c(1,1),t=age,est=trueparm))
        truepen=c(truepen,ap)
      }  
      #print(truepen)
      v1=penest
      v2=estse
      
      
      #v1=matrix(rep(tra.trueparm,500),nrow = 500,nco=9,byrow = TRUE)
      up=v1+v2*1.96
      lo=v1-v2*1.96
      #print(up)
      #print(lo)
      sapply(1:length(truepen), function(s) sum(lo[,s]<=truepen[s]&up[,s]>=truepen[s],na.rm = TRUE) )/B
    }
    Diagnostic=Diagnostic_CO
    Diagnostic2=Diagnostic_CO_E2
    CoverageP=CoveragePED
    CoverageP1=CoverageP1CO
    CoverageP2=CoverageP2CO
  }
  
  comb=estvec
  
  RSE=estvar$RSE
  PSE1=estvar$PSE1
  PSE2=estvar$PSE2
  
  A=Diagnostic(comb,trueparm=trueparm)
  B=Diagnostic2(comb,trueparm=trueparm)
  cov=round(cbind(
    #A$meansd,colMeans(RSE, na.rm=TRUE),
    A$meansd,colmeans(RSE),
    CoverageP(comb,RSE,trueparm=trueparm)),4)
  pen1=round(cbind(
    #A$penbias,colMeans(PSE1, na.rm = TRUE),
    A$penbias,colmeans(PSE1),
    CoverageP1(comb,PSE1,trueparm=trueparm)),4)
  pen2=round(cbind(
    #B$penbias,colMeans(PSE2, na.rm=TRUE),
    B$penbias,colmeans(PSE2),
    CoverageP2(comb,PSE2,trueparm=trueparm)),4)
  return(list(cov=cov,pen1=pen1,pen2=pen2))
}
