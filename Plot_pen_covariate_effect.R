estimate <- c(0.009,
              2.407,
              0.006,
              2.381,
              0.015,
              4.065,
              1.649,
              -0.312,
              -0.588,
              1.756,
              -0.883,
              1.23,
              0.113,
              -0.357,
              -1.51)
estimate <- as.vector(brca1.1[,1])
estimate2 <- as.vector(brca1.1int[,1])

plot_proband_penf <- function(v, cest, bref, affp, type){
  v1 <- v[1:5]
  v2 <- v[6:10]
  terminal_time = v[11]
  v1 = v1[v1<terminal_time]
  v1 = c(v1,terminal_time)
  v1 = v1[!is.na(v1)]
  v2 = v2[!is.na(v2)]
  
  j <- beta_vec1 <- beta_vec2 <- beta_vec3 <- inter <- intg <- intw <- NULL
  v2his <- vector(mode="numeric")
  for(q in 1:length(v1)){
    
    beta1 <- ifelse(q==1, 0, coeff_c(v2[q-1],type=1, cest=cest, affp=affp) )
    beta2 <- ifelse(q==1, 0, coeff_c(v2[q-1],type=2, cest=cest, affp=affp) )
    beta3 <- ifelse(q==1, 0, coeff_c(v2[q-1],type=3, cest=cest, affp=affp) )
    beta_vec1 <- c(beta_vec1,beta1)
    beta_vec2 <- c(beta_vec2,beta2)
    beta_vec3 <- c(beta_vec3,beta3)
    
    if(type==1){
      inter = ifelse(c(q==1,q==1), c(0,0), coeff_interaction_c(v2[q-1], v2his, bref))
      v2his= c(v2his, v2[q-1]); 
      intg = c(intg, inter[1])
      intw = c(intw, inter[2])
    }else{
      intg <- intw <- 0  
    }
    
    
    Penpiece = penf.comp(f=integrand, v1=v1, v2=v2, est=cest,bref=bref, affp=affp,beta_vec1=beta_vec1,beta_vec2=beta_vec2,beta_vec3=beta_vec3,
                         intg=intg,intw=intw,type=type,lower=ifelse(q==1,0,v1[q-1]),upper=v1[q])
    j = c(j,Penpiece)
  }
  j = sum(j,na.rm = TRUE)
  Penet <- log(j)
  
  # print(round(Penet,3))
  return(Penet)
}
plot_penf <- function(f,v1,v2,est,bref, affp,type,beta_vec1,beta_vec2,beta_vec3,intg,intw,lower,upper){  
  integrate(f,v1,v2,est=est,bref=bref, affp=affp,beta_vec1=beta_vec1,beta_vec2=beta_vec2,beta_vec3=beta_vec3,type=type,intg=intg,intw=intw,lower=lower,upper=upper)$value
}  
plot_integrand <-function(u,est,G=G,c=c,tdc=NULL,screen_age=NULL,oopho_age=NULL){
  u = u - 16
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(est[7]*G)
  haz2=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(est[12]*G)
  H1 = (est[1]*u)^est[2]*exp(est[7]*G) #cause 1
  H2 = (est[3]*u)^est[4]*exp(est[12]*G) 
  H3 = (est[5]*u)^est[6]*exp(est[14]*G) 
  
  if(c==1){
    h <- haz1
  }else if(c==2){
    h <- haz2
  }
  f <- h*exp(-(H1+H2+H3))
  return(f)
}
plot_penf <- function(f,est,G=G,c=c,tdc,screen_age,oopho_age,lower,upper)  integrate(f,lower,upper,est=est,G=G,c=c,tdc=tdc,screen_age=screen_age,oopho_age=oopho_age)$value
plot_integrand2 <-function(u,est,G=G,c=c,tdc=NULL,screen_age=NULL,oopho_age=NULL){
  u = u - 16
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(est[7]*G)
  haz2=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(est[18]*G)
  H1 = (est[1]*u)^est[2]*exp(est[7]*G) #cause 1
  H2 = (est[3]*u)^est[4]*exp(est[18]*G) 
  H3 = (est[5]*u)^est[6]*exp(est[20]*G) 
  
  if(c==1){
    h <- haz1
  }else if(c==2){
    h <- haz2
  }
  f <- h*exp(-(H1+H2+H3))
  return(f)
}
# 15 , d.scr

#G,3  #S,4   #O,5   #GxS,6  #OxS,7
#estimate <-     c(0.007,2.655,2.624,-3.086,-1.531,4.331,1.317,-1.133,10.874) 
integrand_td_sc <-function(u,est,G=G,c=c,tdc=NULL,screen_age=NULL,oopho_age=NULL){
  u = u - 16
  tdc.age = screen_age-16
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*
    exp(est[7]*G + ifelse(u<tdc.age,0,1)*est[8] + ifelse(u<tdc.age,0,1)*est[10]*G)
  
  H1 = (est[1]*u)^est[2]*
    exp(est[7]*G + ifelse(u<tdc.age,0,1)*est[8] + ifelse(u<tdc.age,0,1)*est[10]*G) #cause 1
  H2 = (est[3]*u)^est[4]*
    exp(est[12]*G + ifelse(u<tdc.age,0,1)*est[13])
  H3 = (est[5]*u)^est[6]*
    exp(est[14]*G + ifelse(u<tdc.age,0,1)*est[15])  
  #cause 1  
  if(c==1){
    h <- haz1 
  }
  f = h*exp(-(H1+H2+H3)) # h * S 
  return(f)
}
integrand_td_oo <-function(u,est,G=G,c=c,tdc=NULL,screen_age=NULL,oopho_age=NULL){
  u = u - 16
  tdc.age = oopho_age-16
  
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*
    exp(est[7]*G + ifelse(u<tdc.age,0,1)*est[9])
  
  H1 = (est[1]*u)^est[2]*
    exp(est[7]*G + ifelse(u<tdc.age,0,1)*est[9]) #cause 1
  H2 = ((est[3]*u)^est[4]*
          exp(est[12]*G + ifelse(u<tdc.age,0,1)*est[13]))*ifelse(u<tdc.age,1,0)
  H3 = (est[5]*u)^est[6]*
    exp(est[14]*G + ifelse(u<tdc.age,0,1)*est[15])  
  if(c==1){
    h <- haz1
  }
  f = h*exp(-(H1+H2+H3)) # h * S 
  return(f)
}
integrand_td_both <-function(u,est,c=c,G=G,tdc=NULL,screen_age=NULL,oopho_age=NULL){
  u = u - 16
  oage = oopho_age-16
  sage = screen_age-16
  haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*
    exp(est[7]*G + ifelse(u<sage,0,1)*est[8] + ifelse(u<oage,0,1)*est[9] + 
          ifelse(u<sage,0,1)*est[10]*G +
          ifelse(u<sage,0,1)*ifelse(u<oage,0,1)*est[11])
  
  
  H1 = (est[1]*u)^est[2]*
    exp(est[7]*G + ifelse(u<sage,0,1)*est[8] + ifelse(u<oage,0,1)*est[9] + 
          ifelse(u<sage,0,1)*est[10]*G + 
          ifelse(u<sage,0,1)*ifelse(u<oage,0,1)*est[11])
  H2 = ((est[3]*u)^est[4]*
          exp(est[12]*G + ifelse(u<sage,0,1)*est[13]))*ifelse(u<oage,1,0)
  H3 = (est[5]*u)^est[6]*
    exp(est[14]*G + ifelse(u<sage,0,1)*est[15])    
  
  if(c==1){
    h <- haz1
  }
  
  f = h*exp(-(H1+H2+H3)) # h * S 
  return(f)
}



Draw.penet3way <- function(estimate,c=1,ceil, tdc=NULL, screen_age=NULL, oopho_age=NULL){
  est <- estimate
  
  C <- NC <- NULL
  for (i in 16:70){
    C[i-15] <- plot_penf(f=plot_integrand,est=est,c=c,G=1,tdc=NULL,screen_age=NULL,oopho_age=NULL,lower=16,upper=i)
    NC[i-15] <- plot_penf(f=plot_integrand,est=est,c=c,G=0,tdc=NULL,screen_age=NULL,oopho_age=NULL,lower=16,upper=i)
  }
  
  if(!is.null(tdc)){
    if(tdc=="oopho"){
      y3 <- y4 <- NULL
      for (i in 16:70){
        y3[i-15] <- penf.comp(f=integrand_td_oo,est=est,c=c,G=1,tdc=tdc,screen_age=screen_age,oopho_age=oopho_age,lower=16,upper=i)
        y4[i-15] <-  penf.comp(f=integrand_td_oo,est=est,c=c,G=0,tdc=tdc,screen_age=screen_age,oopho_age=oopho_age,lower=16,upper=i) 
      }
      plot(16:oopho_age , C[1:(oopho_age-15)], type = "l" ,ylab="Penetrance", xlab="age at onset", ylim=c(0,ceil), xlim=c(20,70), 
           main=paste0("Breast Cancer Penetrance"),col="red",lty=1)
      lines(16:oopho_age , NC[1:(oopho_age-15)], col="blue", lty=1)
      lines(oopho_age:70, C[(oopho_age-15):55], col="red", lty=3)
      lines(oopho_age:70, NC[(oopho_age-15):55], col="blue", lty=3)
      
      lines(oopho_age:70, y3[(oopho_age-15):55], col="red", lty=1)
      lines(oopho_age:70, y4[(oopho_age-15):55], col="blue", lty=1)
      
      abline( v = oopho_age, col = "gray60")
      text(oopho_age,0.7,"oopho",col = "gray60")
      text(67,y3[55]+0.04,paste0(round(y3[55]*100,1),"%"),col="red",cex=1)
      text(67,y4[55]+0.04,paste0(round(y4[55]*100,1),"%"),col="blue",cex=1)
      legend("topleft", c("Carrier", "Noncarrier","Intervention"), lty=c(1,1,1), col=c("red","blue","gray60"),bg="white")   
    } else if(tdc=="screen"){
      y3 <- y4 <- NULL
      
      for (i in 16:70){
        y3[i-15] <- penf.comp(f=integrand_td_sc,est=est,c=c,G=1,tdc=tdc,screen_age=screen_age,lower=16,upper=i)
        y4[i-15] <-  penf.comp(f=integrand_td_sc,est=est,c=c,G=0,tdc=tdc,screen_age=screen_age,lower=16,upper=i) 
      }
      plot(16:screen_age , C[1:(screen_age-15)], type = "l" ,ylab="Penetrance", xlab="age at onset", ylim=c(0,ceil), xlim=c(20,70), 
           main=paste0("Breast Cancer Penetrance"),col="red",lty=1)
      lines(16:screen_age , NC[1:(screen_age-15)], col="blue", lty=1)
      lines(screen_age:70, C[(screen_age-15):55], col="red", lty=3)
      lines(screen_age:70, NC[(screen_age-15):55], col="blue", lty=3)
      
      lines(screen_age:70, y3[(screen_age-15):55], col="red", lty=1)
      lines(screen_age:70, y4[(screen_age-15):55], col="blue", lty=1)
      
      abline( v = screen_age, col = "gray60")
      text(screen_age,0.7,"screen",col = "gray60")
      text(67,y3[55]+0.04,paste0(round(y3[55]*100,1),"%"),col="red",cex=1)
      text(67,y4[55]+0.04,paste0(round(y4[55]*100,1),"%"),col="blue",cex=1)
      legend("topleft", c("Carrier", "Noncarrier","Intervention"), lty=c(1,1,1), col=c("red","blue","gray60"),bg="white")         
    } else if(tdc=="both"){
      y3 <- y4 <- NULL
      
      for (i in 16:70){
        y3[i-15] <- penf.comp(f=integrand_td_both,est=est,c=c,G=1,tdc=tdc,screen_age=screen_age,oopho_age=oopho_age,lower=16,upper=i)
        y4[i-15] <-  penf.comp(f=integrand_td_both,est=est,c=c,G=0,tdc=tdc,screen_age=screen_age,oopho_age=oopho_age,lower=16,upper=i) 
      }
      if(screen_age < oopho_age){
        plot(16:screen_age , C[1:(screen_age-15)], type = "l" ,ylab="Penetrance", xlab="age at onset", ylim=c(0,ceil), xlim=c(20,70), 
             main=paste0("Breast Cancer Penetrance"),col="red",lty=1)
        lines(16:screen_age , NC[1:(screen_age-15)], col="blue", lty=1)
        lines(screen_age:70, C[(screen_age-15):55], col="red", lty=3)
        lines(screen_age:70, NC[(screen_age-15):55], col="blue", lty=3)
        
        lines(screen_age:70, y3[(screen_age-15):55], col="red", lty=1)
        lines(screen_age:70, y4[(screen_age-15):55], col="blue", lty=1)
      }else{
        plot(16:oopho_age , C[1:(oopho_age-15)], type = "l" ,ylab="Penetrance", xlab="age at onset", ylim=c(0,ceil), xlim=c(20,70), 
             main=paste0("Breast Cancer Penetrance"),col="red",lty=1)
        lines(16:oopho_age , NC[1:(oopho_age-15)], col="blue", lty=1)
        lines(oopho_age:70, C[(oopho_age-15):55], col="red", lty=3)
        lines(oopho_age:70, NC[(oopho_age-15):55], col="blue", lty=3)
        
        lines(oopho_age:70, y3[(oopho_age-15):55], col="red", lty=1)
        lines(oopho_age:70, y4[(oopho_age-15):55], col="blue", lty=1)
        
      }
      abline( v = oopho_age, col = "gray60")
      abline( v = screen_age, col = "gray60")
      text(screen_age,0.7,"screen",col = "gray60")
      text(oopho_age,0.6,"oopho",col = "gray60")
      text(67,y3[55]+0.04,paste0(round(y3[55]*100,1),"%"),col="red",cex=1)
      text(67,y4[55]+0.04,paste0(round(y4[55]*100,1),"%"),col="blue",cex=1)
      legend("topleft", c("Carrier", "Noncarrier","Intervention"),  lty=c(1,1,1), col=c("red","blue","gray60"),bg="white")         
    }  # !is.null(tdc)
  }else{
    
    plot(16:70,C, type = "l" ,ylab="Penetrance", xlab="age at onset", ylim=c(0,ceil), xlim=c(20,70), 
         main=paste0(ifelse(c==1,"Breast","Ovarian")," Cancer"),col="red",lty=1)
    
    if(c==2){
      text(70,C[55]+0.01,paste0(round(C[55]*100,1),"%"),col="red",cex=0.8)
      text(70,NC[55]+0.01,paste0(round(NC[55]*100,1),"%"),col="blue",cex=0.8)
    }else{
      text(70,C[55]+0.04,paste0(round(C[55]*100,1),"%"),col="red",cex=0.8)
      text(70,NC[55]+0.04,paste0(round(NC[55]*100,1),"%"),col="blue",cex=0.8)
      
    }
    
    lines(16:70,NC, col="blue", lty=1)
    legend("topleft", c("Carrier", "Noncarrier"), lty=c(1,1), col=c("red","blue"))  
  }
}
Draw.ontop <- function(estimate,c=1,ceil, tdc=NULL, screen_age=NULL, oopho_age=NULL){
  est <- estimate
  
  C <- NC <- NULL
  for (i in 16:70){
    C[i-15] <- plot_penf(f=plot_integrand2,est=est,c=c,G=1,tdc=NULL,screen_age=NULL,oopho_age=NULL,lower=16,upper=i)
    NC[i-15] <- plot_penf(f=plot_integrand2,est=est,c=c,G=0,tdc=NULL,screen_age=NULL,oopho_age=NULL,lower=16,upper=i)
  }
    lines(16:70,C, type = "l",col="red",lty=3)
    
    # if(c==2){
    #   text(70,C[55]+0.01,paste0(round(C[55]*100,1),"%"),col="red",cex=0.8)
    #   text(70,NC[55]+0.01,paste0(round(NC[55]*100,1),"%"),col="blue",cex=0.8)
    # }else{
    #   text(70,C[55]+0.04,paste0(round(C[55]*100,1),"%"),col=32,cex=0.8)
    #   text(70,NC[55]+0.04,paste0(round(NC[55]*100,1),"%"),col=30,cex=0.8)
    #   
    # }
    
    lines(16:70,NC, col="blue", lty=3)
  
}

1.Without any intervention, Left plot is penetrance curve for breast cancer, Right plot is penetrance curve for ovarian cancer

par(mfrow=c(1,2))
Draw.penet3way(estimate,c=1, ceil=1)
Draw.ontop(estimate2,c=1,ceil=1)

Draw.penet3way(estimate,c=2, ceil=0.5)
Draw.ontop(estimate2,c=2,ceil=1)



2. First screening at age 25,30,40,50 followed by Oophorectomy within 5,10,15 years. Dotted lines are the projection without any intervention

par(mfrow=c(2,4))
Draw.penet3way(estimate, ceil=1, tdc="screen", screen_age=25)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=25, oopho_age=30)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=25, oopho_age=35)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=25, oopho_age=40)

Draw.penet3way(estimate, ceil=1, tdc="screen", screen_age=30)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=30, oopho_age=35)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=30, oopho_age=40)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=30, oopho_age=45)

Draw.penet3way(estimate, ceil=1, tdc="screen", screen_age=40)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=40, oopho_age=45)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=40, oopho_age=50)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=40, oopho_age=55)

Draw.penet3way(estimate, ceil=1, tdc="screen", screen_age=50)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=50, oopho_age=55)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=50, oopho_age=65)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=50, oopho_age=70)


3. Oophorectomy at age 30,40,45,50 followed by first screening within 5,10,15 years. Dotted lines are the projection without any intervention

par(mfrow=c(2,4))
Draw.penet3way(estimate, ceil=1, tdc="oopho", oopho_age=30)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=35, oopho_age=30)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=40, oopho_age=30)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=45, oopho_age=30)

Draw.penet3way(estimate, ceil=1, tdc="oopho", oopho_age=40)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=45, oopho_age=40)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=50, oopho_age=40)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=55, oopho_age=40)

Draw.penet3way(estimate, ceil=1, tdc="oopho", oopho_age=45)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=50, oopho_age=45)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=55, oopho_age=45)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=60, oopho_age=45)

Draw.penet3way(estimate, ceil=1, tdc="oopho", oopho_age=50)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=55, oopho_age=50)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=60, oopho_age=50)
Draw.penet3way(estimate, ceil=1, tdc="both", screen_age=65, oopho_age=50)