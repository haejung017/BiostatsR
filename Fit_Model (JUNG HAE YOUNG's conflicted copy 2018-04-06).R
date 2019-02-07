# FITTING MODEL CODES
FitModelEM <- function( data = data,
                        
                        init.Parms = list( 
                          cause1 = list(base = c(0.008,2.053), time_indep=c(1.663), time_dep=list(main=c(1.79,-2.2,1,-1.65) , interaction=c() )),
                          cause2 = list(base = c(0.002,1.643), time_indep=c(2.658), time_dep=list(main=0.064)),
                          cause3 = list(base = c(0.016,4.193), time_indep=c(-0.322), time_dep=list(main=-1.851)) 
                         ),
                        
                        missing.method = "data",
                        base.dist = "Weibull",
                        frailty=FALSE, 
                        weight=FALSE
){
   
  #initial parameters
  base.parms <- c(init.Parms[[1]]$base, init.Parms[[2]]$base, init.Parms[[3]]$base) # total 6, lambda1,rho1 lambda2,rho2, lambda3,rho3 
    

  vbeta_b <- c( init.Parms[[1]]$time_indep, init.Parms[[1]]$time_dep$main, init.Parms[[1]]$time_dep$interaction ) # breast cancer parameter
  vbeta_o <- c( init.Parms[[2]]$time_indep, init.Parms[[2]]$time_dep$main )                               # ovarian cancer parameter 
  vbeta_d <- c( init.Parms[[3]]$time_indep, init.Parms[[3]]$time_dep$main )                               # death parameter
 
  vbeta <- c(vbeta_b,vbeta_o,vbeta_d)
    
  theta = theta0 = c(log(base.parms),vbeta) # length(theta) = 14

  if(frailty==TRUE){
      data[data$status!=0, "status"] <- 1
      fsize <- aggregate(data$status, by=list(data$famID), length)[,2]
      df <- aggregate(data$status, by=list(data$famID), FUN=sum)[,2]
      data$df <- rep(df, fsize)
      fp <- 2.5
      theta = theta0 = c(log(base.parms), vbeta, log(fp))
  }
    
    # data preparation
    data2 <- Data_preparation(data)
    
    # EM
    if(missing.method=="data"){
    data.cooked = carrierprobgeno(data=data)  
    }
    
    
    est0 <- est <- theta
    # dd <- lval0 <- lval <- 1
    # i <- 0
    # lval <- loglik_Comp_Timedep(theta, data=data, data2=data2, agemin=16)
    # while(dd>0.01){
    #   i <- i+1
    #   est0 <- est
    #   lval0 <- lval
      est1 <- optim(est0, loglik_Comp_Timedep, data=data.cooked, data2=data2, agemin=16, frailty=frailty, weight=weight, 
                        hessian=TRUE, control=list(maxit=6000))
    #   lval <- est1$value
    #   dd <- abs(lval0-lval)
    #   est <- est1$par
    #   #dd <- abs(sum(est-est0))    
    #   print(c(i, dd, est1$par))
    # }
    # cat("Iterations = ", i, "\n")
   
 
    
    # Results 
    print(est1$convergence)
    Var <- try(solve(est1$hessian), TRUE)
    se <- sqrt(diag(Var))
    se.exp <- exp(est1$par)*se
    SE <- c(se.exp[1:6] # 6 for base parms
            ,se[7:length(est1$par)]) # 7 for vbeta
    
    # AIC
    AIC <- 2*(length(est1$par)-(-est1$value))
    Loglik_Value <- -est1$value
    
    # WALD & Pvalue
    wald <- est1$par/SE 
    
    pvalue <- function(wald){
      pvalue <- NULL
      for (w in 1:length(wald)){
        pvalue[w] <- ifelse(wald[w]<0, 2*pnorm(wald[w]), 2*(1-pnorm(wald[w])))
      }
      return(pvalue)
    }
    pval <- pvalue(wald)
    
    # output
    if(frailty==TRUE){
      fp <- exp(est1$par[length(est1$par)])
      mle <- round(cbind(as.matrix(c( exp(est1$par[1:6]), 
                                      est1$par[7:length(est1$par)],fp)), 
                         as.matrix(SE),
                         as.matrix(pval)),3)
      colnames(mle) <- c("MLE","SE","P-value")
      rownames(mle) <- c("Lambda1","Rho1","Lambda2","Rho2","Lambda3","Rho3",
                         "Breast.Gene","Breast.Screen","Breast.oopho",
                         "Ovary.Gene","Ovary.Screen",
                         "Death.Gene","Death.Screen",
                         "frailty_k")
    }else{
      mle <- round(cbind(as.matrix(c( exp(est1$par[1:6]), 
                                      est1$par[7:length(est1$par)])), 
                         as.matrix(SE),
                         as.matrix(pval)),3)
      colnames(mle) <- c("MLE","SE","P-value")
      if(length(init.Parms[[1]]$time_dep$interaction)==0){
      rownames(mle) <- c("Lambda1","Rho1","Lambda2","Rho2","Lambda3","Rho3",
                         "Breast.Gene","Breast.Screen1","Breast.Screen2","Breast.Screen3",
                         "Breast.oopho",
                         "Ovary.Gene","Ovary.Screen1",
                         "Death.Gene","Death.Screen1")}else{
      rownames(mle) <- c("Lambda1","Rho1","Lambda2","Rho2","Lambda3","Rho3",
                          "Breast.Gene","Breast.Screen1","Breast.Screen2","Breast.Screen3",
                          "Breast.oopho","Breast.S1xG","Breast.S2xG","Breast.S3xG",
                          "Breast.S1xO","Breast.S2xO","Breast.S3xO",
                          "Ovary.Gene","Ovary.Screen1",
                          "Death.Gene","Death.Screen1")}
      
    }
    print(mle)
    print(AIC)
    return(mle)
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
