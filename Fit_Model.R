# FITTING MODEL CODES
FitModelEM <- function( data = data,
                        
                        init.Parms = list( 
                          cause1 = list(base = c(0.008,2.581), time_indep=c(2.084), time_dep=list(main=c(2.364,1.033,1.128,0.411) , interaction=c() )),
                          cause2 = list(base = c(0.006,2.538), time_indep=c(1.437), time_dep=list(main=0.818)),
                          cause3 = list(base = c(0.016,4.193), time_indep=c(-0.094), time_dep=list(main=-1.851)) 
                         ),
                        
                        time.dep.cov = list(time.dep.cov.type = "PE" , #  "ED")
                                                  reccur.type = "IS"  #  "CM") 
                        ),
                        
                        missing.method = "data",
                        base.dist = "Weibull",
                        frailty=FALSE, 
                        weight=FALSE,
                        mutation.prediction=FALSE
){
   
  #initial parameters
  base.parms <- c(init.Parms[[1]]$base, init.Parms[[2]]$base, init.Parms[[3]]$base) # total 6, lambda1,rho1 lambda2,rho2, lambda3,rho3 
    
  vbeta_b <- c( init.Parms[[1]]$time_indep, init.Parms[[1]]$time_dep$main, init.Parms[[1]]$time_dep$interaction ) # breast cancer parameter
  vbeta_o <- c( init.Parms[[2]]$time_indep, init.Parms[[2]]$time_dep$main )                               # ovarian cancer parameter 
  vbeta_d <- c( init.Parms[[3]]$time_indep, init.Parms[[3]]$time_dep$main )                               # death parameter
 
  vbeta <- c(vbeta_b,vbeta_o,vbeta_d)
    
  theta = theta0 = c(log(base.parms),vbeta) # length(theta) = 14

  
  if(time.dep.cov$time.dep.cov.type=="ED"){
    # phi = runif(length(init.Parms[[1]]$time_dep$main), 0,1)
    # theta = c(theta, log(phi))
    phi = c(0.1,0.1,0.1,0.1)
    theta = c(theta, log(phi))
  } else if (time.dep.cov$time.dep.cov.type=="CO"){
    #phi = c(0,2, 0,2, 0,2, 0,0.1)
    phi0 = c(0.415, 0.742 , 0.235, -0.1)
    phi1 = c(0.5 ,0.5, 0.5, 0.5 )
    theta = c(theta, log(phi1), phi0)
  }
    
  if(frailty==TRUE){
      datacopy = data
      fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
      df1 <- aggregate(datacopy$status==1, by=list(datacopy$famID), FUN=sum)[,2]
      df2 <- aggregate(datacopy$status==3, by=list(datacopy$famID), FUN=sum)[,2]
      df3 <- aggregate(datacopy$status==4, by=list(datacopy$famID), FUN=sum)[,2]
      data$df1 <- rep(df1, fsize)
      data$df2 <- rep(df2, fsize)
      data$df3 <- rep(df3, fsize)
      fp <- c(3.35683981, 1.86078790)
    #  theta = theta0 = c(theta, log(fp))
  }
    # data preparation
    data2 <- Data_preparation(data)
    
    # EM
    data.cooked = carrierprobgeno(method = missing.method, data=data)  
    impute=function(dat){
      #set.seed(99)
      for(i in 1:nrow(dat)){
        if(is.na(dat$mgene[i]))dat[i,"mgene"]<-sample(c(0,1),1,prob=c(1-dat$carrp.geno[i],dat$carrp.geno[i]),replace=TRUE)
      }
      return(dat)
    }  
    data.cooked = impute(data.cooked)
    
    est0 <- est <- theta
    # dd <- lval0 <- lval <- 1
    # i <- 0
    # lval <- loglik_Comp_Timedep(theta, data=data, data2=data2, agemin=16)
    # while(dd>0.01){
    #   i <- i+1
    #   est0 <- est
    #   lval0 <- lval
      est1 <- optim(est, loglik_Comp_Timedep, data=data.cooked, data2=data2, 
                    agemin=16, frailty=FALSE, weight=weight,
                    reccur.type = time.dep.cov$reccur.type,
                    time.dep.cov.type = time.dep.cov$time.dep.cov.type,
                    hessian=FALSE, control=list(maxit=100000))
      est2 <- optim(log(fp), loglik_Comp_Timedep_second, data=data.cooked, data2=data2, 
                    agemin=16, frailty=TRUE, weight=weight,
                    reccur.type = time.dep.cov$reccur.type,
                    time.dep.cov.type = time.dep.cov$time.dep.cov.type,stage1=est1$par,
                    hessian=FALSE, control=list(maxit=100000))      
      
    #   lval <- est1$value
    #   dd <- abs(lval0-lval)
    #   est <- est1$par
    #   #dd <- abs(sum(est-est0))
    #   print(c(i, dd, est1$par))
    # }
    # cat("Iterations = ", i, "\n")
   
 
    
    # Results #----------------------------------------
    
    # RobustSE #############################################
    # Var <- try(solve(est1$hessian), TRUE)
    # se <- sqrt(diag(Var))
    # se.exp <- exp(est1$par)*se
    # SE <- c(se.exp[1:6] # 6 for base parms
    #         ,se[7:length(est1$par)]) # 7 for vbeta
    ########################################################  
    
    # AIC ##################################################
    pl=length(est1$par)+length(est2$par)
    AIC <- 2*(pl-(-est2$value))
    Loglik_Value <- -est2$value
    ########################################################
    
    # WALD & Pvalue ########################################
    # wald <- est1$par/SE # add chisquare mixture density for parameters with boundary eg. phi in ED, frailty parameter k
    # 
    # pvalue <- function(wald){
    #   pvalue <- NULL
    #   for (w in 1:length(wald)){
    #     pvalue[w] <- ifelse(wald[w]<0, 2*pnorm(wald[w]), 2*(1-pnorm(wald[w])))
    #   }
    #   return(pvalue)
    # }
    # pval <- pvalue(wald)
    ########################################################
    
    # Penetrance ###########################################
    
    ########################################################
    
    # output ###############################################
    mle=c(est1$par,est2$par)
    
    # fp <- exp(est1$par[length(est1$par)])
    # mle <- round(cbind(as.matrix(c( exp(est1$par[1:6]), 
    #                                     est1$par[7:18],
    #                                 exp(est1$par[19:length(theta)])
    #                                 )), 
    #                      as.matrix(SE),
    #                      as.matrix(pval)),4)
    # colnames(mle) <- c("MLE","SE","P-value")
    # rownames(mle) <- c("Lambda1","Rho1","Lambda2","Rho2","Lambda3","Rho3",
    #                       "Breast.Gene","Breast.Screen1","Breast.Screen2","Breast.Screen3","Breast.oopho",
    #                       #"Breast.S1xG","Breast.S2xG","Breast.S3xG",
    #                       #"Breast.S1xO","Breast.S2xO","Breast.S3xO", # full mode interaction
    #                       "Ovary.Gene","Ovary.Screen1",
    #                       "Death.Gene","Death.Screen1"
    #                       ,"frailty_k" 
    #                     )
      
    
    cat("Printing Model estimates, standard error and P-values... \n")
    #print(mle)
    cat("Resuling model has AIC value of ", AIC , "\n")


    # Mutation prediction diagnostics, prints out stats such as Sensitivity, specificity, AUC, NPV, PPV, ACC    
    if(mutation.prediction==TRUE){
      cat("Computing Mutation covariate prediction based on model estimate...")
      Mut.Pred = mutation_prediction(data= data, parms = est1$par) 
      print(Mut.Pred)
    }else{
      Mut.Pred = NULL
    }
    
    
    return(list( MLE=mle, convergence=c(est1$convergence,est2$convergence), 
                 mutation.prediction = Mut.Pred , 
                 AIC=AIC,
                 Loglik_Value=Loglik_Value
                 ))
}

# FITTING MODEL CODES
FitModelEM_NC <- function( data = data,
                           
                           init.Parms = list( 
                             cause1 = list(base = c(0.008,2.581), time_indep=c(2.084), time_dep=list(main=c(2.364,1.033,1.128,0.411) , interaction=c() ))
                           ),
                           
                           time.dep.cov = list(time.dep.cov.type = "PE" , # "PE", "ED", "CO")
                                               reccur.type = "IS"  #  "CM") 
                           ),
                           
                           missing.method = "data",
                           base.dist = "Weibull",
                           frailty=FALSE, 
                           weight=FALSE,
                           mutation.prediction=FALSE
){
  
  #initial parameters
  base.parms <- c(init.Parms[[1]]$base) # total 2, lambda1,rho1
  vbeta <- c( init.Parms[[1]]$time_indep, init.Parms[[1]]$time_dep$main, init.Parms[[1]]$time_dep$interaction ) # breast cancer parameter                             # death parameter
  
  theta = theta0 = c(log(base.parms),vbeta) # length(theta) = 7
  
  
  # setting up time dependent variables
  if(time.dep.cov$time.dep.cov.type=="ED"){ 
    phi = c(3,2,1)
    theta = c(theta, log(phi))
  }else if (time.dep.cov$time.dep.cov.type=="CO"){
    phi0 = c(0,0,0,0)
    phi1 = c(2,2,2,0.1)
    phi= c(phi0,phi1)
    theta = c(theta, phi0, log(phi1))
  }else{phi=NULL}
  
  
  if(frailty==TRUE){
    datacopy = data
    fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
    df1 <- aggregate(datacopy$status==1, by=list(datacopy$famID), FUN=sum)[,2]
    data$df1 <- rep(df1, fsize)
    fp <- 3.46
    theta = theta0 = c(theta, log(fp))
  }else{fp=NULL}
  
  
  # data preparation
  data2 <- Data_preparation_NC(data)
  
  # EM
  data.cooked = carrierprobgeno_NC(method = missing.method, data=data)  
  
  est0 <- est <- theta
  dd <- lval0 <- lval <- 1#
  i <- 0#
  lval <- loglik_Comp_Timedep_NC(theta, data=data.cooked, data2=data2, agemin=16, frailty=frailty, weight=weight,#
                                 reccur.type = time.dep.cov$reccur.type,#
                                 time.dep.cov.type = time.dep.cov$time.dep.cov.type)#
  plot(0:10,seq(-(lval-100),-lval,length.out=11),col="white")
  points(0,-lval,col="red")
  while(dd>0.01){#
    i <- i+1#
    est0 <- est#
    lval0 <- lval#
    est1 <- optim(est0, loglik_Comp_Timedep_NC, data=data.cooked, data2=data2,
                  agemin=16, frailty=frailty, weight=weight,
                  reccur.type = time.dep.cov$reccur.type,
                  time.dep.cov.type = time.dep.cov$time.dep.cov.type,
                  hessian=FALSE, control=list(maxit=50000))
    lval <- est1$value#
    dd <- abs(lval0-lval)#
    est <- est1$par#
    #dd <- abs(sum(est-est0))#
    points(i,-lval,col="red")
    print(c(i, dd, est1$par, est1$convergence))#
  }#
  cat("Iterations = ", i, "\n")#
  est1 =        optim(est, loglik_Comp_Timedep_NC, data=data.cooked, data2=data2, 
                      agemin=16, frailty=frailty, weight=weight,
                      reccur.type = time.dep.cov$reccur.type,
                      time.dep.cov.type = time.dep.cov$time.dep.cov.type,
                      hessian=TRUE, control=list(maxit=50000))
  
  
  
  # Results 
  print(est1$convergence)
  Var <- try(solve(est1$hessian), TRUE)
  se <- sqrt(diag(Var))
  se.exp <- exp(est1$par)*se
  SE <- c(se.exp[1:2] # 2 for base parms
          ,se[3:length(est1$par)]) # 5 for vbeta
  
  # AIC
  AIC <- 2*(length(est1$par)-(-est1$value))
  Loglik_Value <- -est1$value
  
  # WALD & Pvalue
  wald <- est1$par/SE # add chisquare mixture density for parameters with boundary eg. phi in ED, frailty parameter k
  
  pvalue <- function(wald){
    pvalue <- NULL
    for (w in 1:length(wald)){
      pvalue[w] <- ifelse(wald[w]<0, 2*pnorm(wald[w]), 2*(1-pnorm(wald[w])))
    }
    return(pvalue)
  }
  pval <- pvalue(wald)
  
  # output
  LT = length(theta)
  Lp = length(phi)
  Lf = length(fp)
  
  mle <- round(cbind(as.matrix(c( exp(est1$par[1:2]), 
                                  est1$par[3:(LT-Lp-Lf)],
                                  if(Lp!=0 | Lf!=0){exp(est1$par[(LT-Lp-Lf+1):LT])}else{NULL}
  )), 
  as.matrix(SE),
  as.matrix(pval)),4)
  # colnames(mle) <- c("MLE","SE","P-value")
  # rownames(mle) <- c("Lambda1","Rho1","Lambda2","Rho2","Lambda3","Rho3",
  #                       "Breast.Gene","Breast.Screen1","Breast.Screen2","Breast.Screen3","Breast.oopho",
  #                       #"Breast.S1xG","Breast.S2xG","Breast.S3xG",
  #                       #"Breast.S1xO","Breast.S2xO","Breast.S3xO", # full mode interaction
  #                       "Ovary.Gene","Ovary.Screen1",
  #                       "Death.Gene","Death.Screen1"
  #                       ,"frailty_k" 
  #                     )
  
  
  cat("Printing Model estimates, standard error and P-values... \n")
  print(mle)
  cat("Resuling model has AIC value of ", AIC , "\n")
  
  
  # Mutation prediction diagnostics, prints out stats such as Sensitivity, specificity, AUC, NPV, PPV, ACC    
  if(mutation.prediction==TRUE){
    cat("Computing Mutation covariate prediction based on model estimate...")
    Mut.Pred = mutation_prediction(data= data, parms = est1$par) 
    print(Mut.Pred)
  }else{
    Mut.Pred = NULL
  }
  
  
  return(list( estimation=est1$par, 
               mutation.prediction = Mut.Pred , 
               AIC=AIC))
}


  

## Basic functions
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
carrierprobgeno_NC <-function(method="data", data, mode="dominant", q=0.0021){
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
        for(a in c(0,1)){
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
  br_time_p1 <- br_time_p[statusp==1]
  ov_time_p1 <- ov_time_p[statusp==1]
  sc1_time_p1 <- sc1_time_pn[statusp==1]
  sc2_time_p1 <- sc2_time_pn[statusp==1]
  sc3_time_p1 <- sc3_time_pn[statusp==1]
  cagep1 <- cagep[statusp==1]
  timep1 <- timep[statusp==1]
  
  delta_BR_p1 <- ifelse(timep1-br_time_p1 > 0,1,0)
  delta_OR_p1 <- ifelse(timep1-ov_time_p1 > 0,1,0)
  
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
Data_preparation_NC <- function(data_raw){
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
  
  br_time_p1 <- br_time_p[statusp==1 |statusp==2]
  ov_time_p1 <- ov_time_p[statusp==1 |statusp==2]
  sc1_time_p1 <- sc1_time_pn[statusp==1 |statusp==2]
  sc2_time_p1 <- sc2_time_pn[statusp==1 |statusp==2]
  sc3_time_p1 <- sc3_time_pn[statusp==1 |statusp==2]
  cagep1 <- cagep[statusp==1 |statusp==2]
  timep1 <- timep[statusp==1 |statusp==2]
  
  delta_BR_p1 <- ifelse(timep1-br_time_p1 > 0,1,0)
  delta_OR_p1 <- ifelse(timep1-ov_time_p1 > 0,1,0)
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
  
  return(list(v1=time_vec2,v2=indicator_vec,
              v1p0=time_vec2_p0, v2p0=indicator_vec_p0,
              v1p1=time_vec2_p1, v2p1=indicator_vec_p1))   
}