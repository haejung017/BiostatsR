library(FamEvent)
write.excel <- function(x,row.names=FALSE,col.names=FALSE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
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
        for(a in c(0,1,2)){
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

############################################ Model1, pop+, time independent covariate, indepen
##############################################################################################
penetrancebias = function(simresult, trueparm){
  ev=simresult
  index = expand.grid(c(0,1),c(0,1))
  xbeta=c(trueparm[3],trueparm[4]) %*% t(index)
  
  t= c(25,35,45,55)
  total=list()
  for(i in 1:4){
  H=(trueparm[1]*t)^trueparm[2]*exp(xbeta[i])#TRUE
  truepen= 1-exp(-H)
  l=lapply(t, function(t) 1-exp(-(ev[,1]*t)^ev[,2]*exp(index[i,1]*ev[,3]+index[i,2]*ev[,4])))
  Bias = c(mean(l[[1]]),mean(l[[2]]),mean(l[[3]]),mean(l[[4]]))-truepen
  EmpSE = c(sd(l[[1]]),sd(l[[2]]),sd(l[[3]]),sd(l[[4]]))
  RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2),mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
  total[[i]] = cbind(truepen,Bias,EmpSE,RMSE)
  }
  
  pen40=NULL;pen50=NULL;pen60=NULL;pen70=NULL
  for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
    pen40=rbind(pen40,total[[i]][1,])
    pen50=rbind(pen50,total[[i]][2,])
    pen60=rbind(pen60,total[[i]][3,])
    pen70=rbind(pen70,total[[i]][4,])
  }
  out=rbind(pen40,pen50,pen60,pen70)  
  return(out)
}
sim1=function(trueparm,mrate=0){
  iter = 0
  lam=trueparm[1];rho=trueparm[2];S=trueparm[3];G=trueparm[4];
  estvec=NULL
    while(iter < 30){
      simdata=simfam(N.fam =150, design="pop+", variation = "none", base.parms=c(lam,rho), vbeta=c(S,G),
                      allelefreq = 0.0021, agemin=16,mrate=mrate)
            k=BC_simulation(data=simdata, initparm=c(lam,rho,S,G),frailty=FALSE)
     estvec = rbind(estvec,k)
       iter = iter+1
        if(iter%%30==0){
            print(iter)
            print(estvec[iter-1,])
            }
    }
    truemat = matrix(rep(trueparm,500), ncol=length(trueparm), nrow=500, byrow = TRUE)
    Bias=colMeans(estvec)-trueparm
    EmpSE=c(sd(estvec[,1]),sd(estvec[,2]),sd(estvec[,3]),sd(estvec[,4]))
    RMSE=sqrt( colMeans( (estvec-truemat)^2 ))
               
    meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
    print(meansd)
    return(list(simresult=estvec, meansd=meansd))
}
##################################################### simulation running & results#############
trueparm=c(0.01,2.5,-1,1.5)
F150_500_mr0   = sim1(trueparm) # mrate=0
F150_500_mr0.5 = sim1(trueparm, mrate=0.5) # mrate=0.5
F150_500_mr0.8 = sim1(trueparm, mrate=0.8) # mrate=0.8
#Result 
write.excel(F150_500_mr0$meansd[,-1])
write.excel(penetrancebias(F150_500_mr0$simresult,trueparm)[,-1]) #penetrance
write.excel(F150_500_mr0.5$meansd[,-1])
write.excel(penetrancebias(F150_500_mr0.5$simresult,trueparm)[,-1]) 
write.excel(F150_500_mr0.8$meansd[,-1])
write.excel(penetrancebias(F150_500_mr0.8$simresult,trueparm)[,-1]) 
###############################################################################################
###############################################################################################


############################################ Model2, pop+, time independent covariate, frailty
##############################################################################################
penetrancebias2 = function(simresult, trueparm){
  ev=simresult
  index = expand.grid(c(0,1),c(0,1))
  xbeta=c(trueparm[3],trueparm[4]) %*% t(index)
  k = trueparm[5]
  
  t= c(25,35,45,55)
  total=list()
  for(i in 1:4){
    H=(trueparm[1]*t)^trueparm[2]*exp(xbeta[i])#TRUE
    truepen= 1-(1+H/k)^(-k)
    
    l=lapply(t, function(t) 1- ( 1+ ((ev[,1]*t)^ev[,2]*exp(index[i,1]*ev[,3]+index[i,2]*ev[,4]))/ev[,5] )^(-ev[,5]))
    
    Bias = c(mean(l[[1]]),mean(l[[2]]),mean(l[[3]]),mean(l[[4]]))-truepen
    EmpSE = c(sd(l[[1]]),sd(l[[2]]),sd(l[[3]]),sd(l[[4]]))
    RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2),mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
    total[[i]] = cbind(truepen,Bias,EmpSE,RMSE)
  }
  
  pen40=NULL;pen50=NULL;pen60=NULL;pen70=NULL
  for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
    pen40=rbind(pen40,total[[i]][1,])
    pen50=rbind(pen50,total[[i]][2,])
    pen60=rbind(pen60,total[[i]][3,])
    pen70=rbind(pen70,total[[i]][4,])
  }
  out=rbind(pen40,pen50,pen60,pen70)  
  return(out)
}
sim2=function(trueparm,mrate=0,nfam){
  iter = 0
  lam=trueparm[1];rho=trueparm[2];S=trueparm[3];G=trueparm[4];k=trueparm[5]
  estvec=NULL
  while(iter < 500){
    simdata=simfam(N.fam =nfam, design="pop+", variation = "frailty", depend=1/k, base.parms=c(lam,rho), vbeta=c(S,G),
                   allelefreq = 0.0021, agemin=16,mrate=mrate)
    simr=tsestfun(data=simdata, initparm=c(lam,rho,S,G,k),frailty=TRUE)
    #simr=FitModelEM_NC_simulation(data=simdata, initparm=c(lam,rho,S,G,k),frailty=TRUE)
    estvec = rbind(estvec,simr)
    iter = iter+1
    if(iter%%30==0){
      print(iter)
      print(estvec[iter-1,])
    }
  }
  truemat = matrix(rep(trueparm,500), ncol=length(trueparm), nrow=500, byrow = TRUE)
  Bias=colMeans(estvec)-trueparm
  EmpSE=c(sd(estvec[,1]),sd(estvec[,2]),sd(estvec[,3]),sd(estvec[,4]),sd(estvec[,5]))
  RMSE=sqrt( colMeans( (estvec-truemat)^2 ))
  
  meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
  print(meansd)
  return(list(simresult=estvec, meansd=meansd))
}
##################################################### simulation running & results#############
######## k = 2.5 ##############################################################################
trueparm=c(0.01,2.5,-1,1.5, 2.5);mrate=0
F150_500_mr0_m2   = sim2(trueparm, nfam=500) # mrate=0
F150_500_mr0.5_m2 = sim2(trueparm, mrate=0.5,nfam=500) # mrate=0.5
F150_500_mr0.8_m2 = sim2(trueparm, mrate=0.8,nfam=500) # mrate=0.8
#Result 
write.excel(F150_500_mr0_m2$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0_m2$simresult,trueparm)[,-1]) #penetrance
write.excel(F150_500_mr0.5_m2$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0.5_m2$simresult,trueparm)[,-1]) 
write.excel(F150_500_mr0.8_m2$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0.8_m2$simresult,trueparm)[,-1]) 
###############################################################################################
######## k = 3.5 ##############################################################################
trueparm=c(0.01,2.5,-1,1.5, 3.5);mrate=0
F150_500_mr0_m2_k3.5   = sim2(trueparm,nfam=500) # mrate=0
F150_500_mr0.5_m2_k3.5 = sim2(trueparm, mrate=0.5,nfam=500) # mrate=0.5
F150_500_mr0.8_m2_k3.5 = sim2(trueparm, mrate=0.8,nfam=500) # mrate=0.8
#Result 
write.excel(F150_500_mr0_m2_k3.5$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0_m2_k3.5$simresult,trueparm)[,-1]) #penetrance
write.excel(F150_500_mr0.5_m2_k3.5$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0.5_m2_k3.5$simresult,trueparm)[,-1]) 
write.excel(F150_500_mr0.8_m2_k3.5$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0.8_m2_k3.5$simresult,trueparm)[,-1]) 
###############################################################################################
######## k = 5 ################################################################################
trueparm=c(0.01,2.5,-1,1.5, 5);mrate=0
F150_500_mr0_m2_k5   = sim2(trueparm,nfam=500) # mrate=0
F150_500_mr0.5_m2_k5 = sim2(trueparm, mrate=0.5,nfam=500) # mrate=0.5
F150_500_mr0.8_m2_k5 = sim2(trueparm, mrate=0.8,nfam=500) # mrate=0.8
#Result 
write.excel(F150_500_mr0_m2_k5$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0_m2_k5$simresult,trueparm)[,-1]) #penetrance
write.excel(F150_500_mr0.5_m2_k5$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0.5_m2_k5$simresult,trueparm)[,-1]) 
write.excel(F150_500_mr0.8_m2_k5$meansd[,-1])
write.excel(penetrancebias2(F150_500_mr0.8_m2_k5$simresult,trueparm)[,-1]) 
###############################################################################################
###############################################################################################



############################ Model3, pop+, competing risk, time independent covariate, no frailty
##############################################################################################
BiasSE= function(estvec, trueparm){
  B= nrow(estvec)
  truemat = matrix(rep(trueparm,B), ncol=length(trueparm), nrow=B, byrow = TRUE)
  Bias=colMeans(estvec)-trueparm
  EmpSE= apply(estvec, 2, sd)
  RMSE= sqrt( colMeans( (estvec-truemat)^2 ))
  meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
  round(meansd,4)
}  #estvec is the simulation result
penBiasSE = function(estvec, trueparm){
  ev=estvec
  index = as.matrix(expand.grid(c(0,1),c(0,1)))
  #k = trueparm[5]
  
  penf <- function(f, X, est,  k,frailty,  t){  
    integrate(integrand, X=X, est=est, k=k, frailty=frailty,lower=0,upper=t)$value
  }
  integrand <- function(u, X, est, k,type, frailty){
    u <- u
    fp=k
    index=X
    xbeta1 = index[1]*est[5] + index[2]*est[6]
    xbeta2 = index[1]*est[7] + index[2]*est[8]
    
    haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
    
    H1 = - (est[1]*u)^est[2]*exp(xbeta1)
    H2 = - (est[3]*u)^est[4]*exp(xbeta2)
    
    if(frailty==TRUE){
      f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
    } else{
      f = haz1*exp(H1+H2) # h * S
    }
    return(f)
  }
  
  t= c(25,35,45,55)
  total=list()
  for(i in 1:4){ #iterate for S,G combinations 
    #Calculate true penetrance for given SG 
    truepen = sapply(t, function(t) penf(integrand, X=index[i,], est=trueparm, k=NULL,frailty=FALSE, t=t))
    #Calculate estimated penetrnace for given t, SG
    l=list()
    estpen=NULL
    for(s in 1:nrow(ev)){
      evs = ev[s,]
      r = sapply(t, function(t) penf(integrand, X=index[i,], est=evs, k=NULL,frailty=FALSE, t=t))
      estpen =rbind(estpen,r)
    }
    
    #calculate bias for i
    l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
    Bias = unlist(lapply(l,mean))-truepen
    EmpSE = unlist(lapply(l,sd))
    RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2),mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
    total[[i]] = cbind(truepen,Bias,EmpSE,RMSE)
  }
  
  pen40=NULL;pen50=NULL;pen60=NULL;pen70=NULL
  for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
    pen40=rbind(pen40,total[[i]][1,])
    pen50=rbind(pen50,total[[i]][2,])
    pen60=rbind(pen60,total[[i]][3,])
    pen70=rbind(pen70,total[[i]][4,])
  }
  out=rbind(pen40,pen50,pen60,pen70)  
  round(out,4)
}
sim3=function(B=30, N=300, trueparm,mrate=0){
  iter = 0
  #lam=trueparm[1];rho=trueparm[2];S=trueparm[3];G=trueparm[4];k=trueparm[5]
  estvec=NULL
  while(iter < B){
    # simdata=simfam.cmp(N.fam =150, design="pop+", variation = "frailty", depend=1/k, base.parms=c(lam,rho), vbeta=c(S,G),
    #                allelefreq = 0.0021, agemin=16,mrate=mrate)
    simdata=simfam.cmp(N.fam = N, trueparm, frailty.dist="none")
    simr = Comp_simulation(data=simdata, initparm=trueparm, frailty=FALSE)
    #simr=FitModelEM_NC_simulation(data=simdata, initparm=c(lam,rho,S,G,k),frailty=TRUE)
    estvec = rbind(estvec,simr)
    iter = iter+1
    #if(iter%%1==0){
      print(iter)
      print(estvec[iter,])
    #}
  }
  truemat = matrix(rep(trueparm,B), ncol=length(trueparm), nrow=B, byrow = TRUE)
  Bias=colMeans(estvec)-trueparm
  EmpSE= apply(estvec, 2, sd)
  RMSE= sqrt( colMeans( (estvec-truemat)^2 ))
  
  meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
  print(round(meansd,4))
  return(list(simresult=estvec, meansd=meansd))
}
##################################################### simulation running & results#############
######## k = 2.5 ##############################################################################
trueparm= c(0.01,2.5, 0.007,3, -1,1.5, -0.5,1 ) ;mrate=0

# mrate=0
F300   = sim3(N=300,trueparm=trueparm) 
F300.2 = sim3(B=40,N=300,trueparm=trueparm) # mrate=
F3_M0.total = rbind(F300$simresult, F300.2$simresult)
BiasSE(F300totla, trueparm)
penBiasSE(F300totla, trueparm)

# mrate=0.5
F3_M5.1    = sim3(B=30,trueparm=trueparm,mrate=0.5) 
F3_M5.2    = sim3(B=40,trueparm=trueparm,mrate=0.5) # mrate=
F3_M5.total = rbind(F3_M5.1$simresult, F3_M5.2$simresult)
BiasSE(F3_M5.total, trueparm)
penBiasSE(F3_M5.total, trueparm)

# mrate=0.8
F3_M8.1    = sim3(B=30,trueparm=trueparm,mrate=0.8) 
F3_M8.2    = sim3(B=40,trueparm=trueparm,mrate=0.8) # mrate=
F3_M8.total = rbind(F3_M8.1$simresult, F3_M8.2$simresult)
BiasSE(F3_M8.total, trueparm)
penBiasSE(F3_M8.total, trueparm)


#save result
write.excel(cbind(BiasSE(F300totla, trueparm),BiasSE(F3_M5.total, trueparm),BiasSE(F3_M8.total, trueparm)))
write.excel(cbind(penBiasSE(F300totla, trueparm),penBiasSE(F3_M5.total, trueparm),penBiasSE(F3_M8.total, trueparm)))
save(F3_M0.total, F3_M5.total, F3_M8.total, file = "simjune22_model3.RData"); 
load("my_data.RData")
###############################################################################################
###############################################################################################


############################ Model4, pop+, competing risk, time independent covariate, frailty
##############################################################################################
tau=function(k){1/(1+2*k)}
BiasSE= function(estvec, trueparm){
  B= nrow(estvec)
  truemat = matrix(rep(trueparm,B), ncol=length(trueparm), nrow=B, byrow = TRUE)
  Bias=colMeans(estvec)-trueparm
  EmpSE= apply(estvec, 2, sd)
  RMSE= sqrt( colMeans( (estvec-truemat)^2 ))
  meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
  round(meansd,4)
}  #estvec is the simulation result
penBiasSE = function(estvec, trueparm){
  ev=estvec
  index = as.matrix(expand.grid(c(0,1),c(0,1)))
  
  
  penf <- function(f, X, est, frailty,  t){  
    integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
  }
  integrand <- function(u, X, est, type, frailty){
    u <- u
    fp=est[9]
    index=X
    xbeta1 = index[1]*est[5] + index[2]*est[6]
    xbeta2 = index[1]*est[7] + index[2]*est[8]
    
    haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
    
    H1 = - (est[1]*u)^est[2]*exp(xbeta1)
    H2 = - (est[3]*u)^est[4]*exp(xbeta2)
    
    if(frailty==TRUE){
     # f = haz1*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
      f = haz1*(1-(H1+H2)/fp[1])^(-fp[1])
    } else{
      f = haz1*exp(H1+H2) # h * S
    }
    return(f)
  }
  
  t= c(25,35,45,55)
  total=list()
  for(i in 1:4){ #iterate for S,G combinations 
    #Calculate true penetrance for given SG 
    truepen = sapply(t, function(t) penf(integrand, X=index[i,], est=trueparm, frailty=TRUE, t=t))
    #Calculate estimated penetrnace for given t, SG
    l=list()
    estpen=NULL
    for(s in 1:nrow(ev)){
      evs = ev[s,]
      r = sapply(t, function(t) penf(integrand, X=index[i,], est=evs, frailty=TRUE, t=t))
      estpen =rbind(estpen,r)
    }
    
    #calculate bias for i
    l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
    Bias = unlist(lapply(l,mean))-truepen
    EmpSE = unlist(lapply(l,sd))
    RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2),mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
    total[[i]] = cbind(truepen,Bias,EmpSE,RMSE)
  }
  
  pen40=NULL;pen50=NULL;pen60=NULL;pen70=NULL
  for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
    pen40=rbind(pen40,total[[i]][1,])
    pen50=rbind(pen50,total[[i]][2,])
    pen60=rbind(pen60,total[[i]][3,])
    pen70=rbind(pen70,total[[i]][4,])
  }
  out=rbind(pen40,pen50,pen60,pen70)  
  round(out,4)
}
sim4=function(B=30, N=300, trueparm, mrate=0){
  iter = 0
  #lam=trueparm[1];rho=trueparm[2];S=trueparm[3];G=trueparm[4];k=trueparm[5]
  estvec=NULL
  while(iter < B){
    # simdata=simfam.cmp(N.fam =150, design="pop+", variation = "frailty", depend=1/k, base.parms=c(lam,rho), vbeta=c(S,G),
    #                allelefreq = 0.0021, agemin=16,mrate=mrate)
    simdata=simfam.cmp(N.fam = N, trueparm, frailty.dist="gamma")
    simr = Comp_simulation(data=simdata, initparm=trueparm, frailty=TRUE)
    #simr=FitModelEM_NC_simulation(data=simdata, initparm=c(lam,rho,S,G,k),frailty=TRUE)
    estvec = rbind(estvec,simr)
    iter = iter+1
    #if(iter%%1==0){
    print(iter)
    print(estvec[iter,])
    #}
  }
  truemat = matrix(rep(trueparm,B), ncol=length(trueparm), nrow=B, byrow = TRUE)
  Bias=colMeans(estvec)-trueparm
  EmpSE= apply(estvec, 2, sd)
  RMSE= sqrt( colMeans( (estvec-truemat)^2 ))
  
  meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
  print(round(meansd,4))
  return(list(simresult=estvec, meansd=meansd))
}
##################################################### simulation running & results#############
######## k = 2.5 ##############################################################################
trueparm= c(0.01,2.5, 0.007,3, -1,1.5, -0.5,1 ,3) ;mrate=0

# mrate=0
M4_F3_G0.1  = sim4(B=30,trueparm=trueparm) 
M4_F3_G0.2  = sim4(B=30,trueparm=trueparm)
M4_F3_G0.3  = sim4(B=30,trueparm=trueparm)#

M4_F3_G0.total = rbind(M4_F3_G0.1$simresult, M4_F3_G0.2$simresult)
BiasSE(M4_F3_G0.total, trueparm) #test=M4_F3_G0.total[M4_F3_G0.total[,9]<50,];BiasSE(test,trueparm);penBiasSE(test,trueparm)
penBiasSE(M4_F3_G0.total, trueparm)

# mrate=0.5
M4_F3_G5.1  = sim4(B=30,trueparm=trueparm,mrate=0.5) 
M4_F3_G5.2  = sim4(B=30,trueparm=trueparm)
M4_F3_G5.3  = sim4(B=30,trueparm=trueparm)#

M4_F3_G5.total = rbind(M4_F3_G5.1$simresult)#, M4_F3_G5.2$simresult)
BiasSE(M4_F3_G5.total, trueparm) #test=M4_F3_G0.total[M4_F3_G0.total[,9]<50,];BiasSE(test,trueparm);penBiasSE(test,trueparm)
penBiasSE(M4_F3_G5.total, trueparm)

# mrate=0.8
F3_M8.1    = sim3(B=30,trueparm=trueparm,mrate=0.8) 
F3_M8.2    = sim3(B=40,trueparm=trueparm,mrate=0.8) # mrate=
F3_M8.total = rbind(F3_M8.1$simresult, F3_M8.2$simresult)
BiasSE(F3_M8.total, trueparm)
penBiasSE(F3_M8.total, trueparm)


#save result
write.excel(cbind(BiasSE(F300totla, trueparm),BiasSE(F3_M5.total, trueparm),BiasSE(F3_M8.total, trueparm)))
write.excel(cbind(penBiasSE(F300totla, trueparm),penBiasSE(F3_M5.total, trueparm),penBiasSE(F3_M8.total, trueparm)))
save(F300.total, F3_M5.total, file = "simjune22.RData"); 
load("my_data.RData")
###############################################################################################
###############################################################################################


############################ Model5, pop+, competing risk, time independent covariate, 2frailty
##############################################################################################
tau=function(k){1/(1+2*k)}
BiasSE= function(estvec, trueparm){
  B= nrow(estvec)
  truemat = matrix(rep(trueparm,B), ncol=length(trueparm), nrow=B, byrow = TRUE)
  Bias=colMeans(estvec)-trueparm
  EmpSE= apply(estvec, 2, sd)
  RMSE= sqrt( colMeans( (estvec-truemat)^2 ))
  meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
  round(meansd,4)
}  #estvec is the simulation result
penBiasSE = function(estvec, trueparm){
  ev=estvec$simresult
  index = as.matrix(expand.grid(c(0,1),c(0,1)))
  
  
  penf <- function(f, X, est, frailty,  t){  
    integrate(integrand, X=X, est=est, frailty=frailty,lower=0,upper=t)$value
  }
  integrand <- function(u, X, est, type, frailty){
    u <- u
    fp=c(est[9],est[10])
    index=X
    xbeta1 = index[1]*est[5] + index[2]*est[6]
    xbeta2 = index[1]*est[7] + index[2]*est[8]
    
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
  
  t= c(25,35,45,55)
  total=list()
  for(i in 1:4){ #iterate for S,G combinations 
    #Calculate true penetrance for given SG 
    truepen = sapply(t, function(t) penf(integrand, X=index[i,], est=trueparm, frailty=TRUE, t=t))
    #Calculate estimated penetrnace for given t, SG
    l=list()
    estpen=NULL
    for(s in 1:nrow(ev)){
      evs = ev[s,]
      r = sapply(t, function(t) penf(integrand, X=index[i,], est=evs, frailty=TRUE, t=t))
      estpen =rbind(estpen,r)
    }
    
    #calculate bias for i
    l[[1]]=estpen[,1];l[[2]]=estpen[,2];l[[3]]=estpen[,3];l[[4]]=estpen[,4]
    Bias = unlist(lapply(l,mean))-truepen
    EmpSE = unlist(lapply(l,sd))
    RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2),mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
    total[[i]] = cbind(truepen,Bias,EmpSE,RMSE)
  }
  
  pen40=NULL;pen50=NULL;pen60=NULL;pen70=NULL
  for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
    pen40=rbind(pen40,total[[i]][1,])
    pen50=rbind(pen50,total[[i]][2,])
    pen60=rbind(pen60,total[[i]][3,])
    pen70=rbind(pen70,total[[i]][4,])
  }
  out=rbind(pen40,pen50,pen60,pen70)  
  round(out,4)
}
sim5=function(B=30, N=300, trueparm, mrate=0,twostage=FALSE){
  iter = 0
  #lam=trueparm[1];rho=trueparm[2];S=trueparm[3];G=trueparm[4];k=trueparm[5]
  estvec=NULL
  while(iter < B){
    # simdata=simfam.cmp(N.fam =150, design="pop+", variation = "frailty", depend=1/k, base.parms=c(lam,rho), vbeta=c(S,G),
    #                allelefreq = 0.0021, agemin=16,mrate=mrate)
    simdata=simfam.cmp(N.fam =500, trueparm, frailty.dist="gamma", mrate=mrate)
    simr = if(twostage==FALSE){
            Comp_simulation2(data=simdata, initparm=trueparm, frailty=TRUE)
           }else{
            Comp_simulation2_2S(data=simdata, initparm=trueparm, frailty=FALSE) 
           }
    
    estvec = rbind(estvec,simr)
    iter = iter+1
    if(iter%%1==0){
    print(iter)
    print(estvec[iter,])
    }
  }
  truemat = matrix(rep(trueparm,B), ncol=length(trueparm), nrow=B, byrow = TRUE)
  Bias=colMeans(estvec)-trueparm
  EmpSE= apply(estvec, 2, sd)
  RMSE= sqrt( colMeans( (estvec-truemat)^2 ))
  
  meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
  print(round(meansd,4))
  return(list(simresult=estvec, meansd=meansd))
}
##################################################### simulation running & results#############
######## k = 2.5 ##############################################################################
trueparm= c(0.01,2.5, 0.01, 2.5, 1,2.5, 0.3,2 ,3.5,2.5) ;mrate=0

# mrate=0
M5_F3_G0.1  = sim5(B=5,trueparm=trueparm)         #trueparm= c(0.01,2.5, 0.007,3, -1,1.5, -0.5,1 ,2.5,5) ;mrate=0
M5_F3_G0.2  = sim5(B=10,trueparm=trueparm) 
M5_F3_G0.3  = sim5(B=15,trueparm=trueparm)
M5_F3_G0.4  = sim5(B=20,trueparm=trueparm)

r1  = sim5(B=20,trueparm=trueparm)         #trueparm= c(0.015,2.5, 0.015,2.5, -1,2.5, -0.5,2 ,3.5,2.5) ;mrate=0
r1.2  = sim5(B=20,trueparm=trueparm)
r1.3  = sim5(B=20,trueparm=trueparm)
r1.4  = sim5(B=20,trueparm=trueparm)
r1.5  = sim5(B=20,trueparm=trueparm)
r=rbind(r1$simresult,r1.2$simresult,r1.3$simresult,r1.4$simresult,r1.5$simresult)
1714  623  527 
r2 = sim5(B=20,trueparm=trueparm) # trueparm= c(0.015,2.5, 0.01,2, -1,2.5, -0.5,2 ,3.5,2.5) ;mrate=0

r3 = sim5(B=20,trueparm=trueparm) #N300,trueparm= c(0.01,2.5, 0.008,3, -1,2.5, -0.5,2 ,3.5,2.5) 
r4 = sim5(N=700,B=20,trueparm=trueparm)#N700,trueparm= c(0.01,2.5, 0.008,3, -1,2.5, -0.5,2 ,3.5,2.5) 

r4 = sim5(N=700,B=20,trueparm=trueparm)#N300,trueparm= c(0.01,2.5, 0.008,3, -1,2.5, -0.5,2 ,3.5,2.5) 

r5 = sim5(B=20,trueparm=trueparm)#N300,trueparm= c(0.015,2.5, 0.015, 2.5, -1,2.5, -0.5,2 ,3.5,2.5) 

r6 = sim5(N=150,B=20,trueparm=trueparm)#N150,trueparm= c(0.015,2.5, 0.015, 2.5, -1,2.5, -0.5,2 ,3.5,2.5) 
                                 #works well family size doesnt matter

r7 = sim5(N=150,B=20,trueparm=trueparm)#N150,trueparm= c(0.015,2.5, 0.015, 2.5, -1,2, -0.5,1.5 ,3.5,2.5)
                                  #works well gene parameter doesnt matter

r8 = sim5(N=150,B=50,trueparm=trueparm)#N150,trueparm= c(0.015,2.5, 0.015, 2.5, -1,2, -0.5,1.5 ,3.5,2.5)
                                        #works well gene parameter doesnt matter

r9 = sim5(N=150,B=50,trueparm=trueparm)


r10 = sim5(N=150,B=50,trueparm=trueparm)#trueparm= c(0.015,2.5, 0.015, 2.5, -1,2, -0.5,1.5 ,1,2.5) 
                                        #high dependence situation 5%penbias

r11 = sim5(N=150,B=50,trueparm=trueparm)#trueparm= c(0.015,2.5, 0.015, 2.5, -1,2.5, -0.5,2 ,1,2.5) 

r12 = sim5(N=150,B=50,trueparm=trueparm)#trueparm= c(0.015,2.5, 0.015, 2.5, -1,2.5, -0.5,2 ,5,2.5) 
                                        # low independence works great <1% penbias

r13 = sim5(N=150,B=50,trueparm=trueparm)#trueparm= c(0.01,2.5, 0.01, 2.5, -1,3.8, -0.5,3.5 ,3.5,2.5)
                                        # lower baseline, higher gene effect <1% penbias

r14 = sim5(N=150,B=50,trueparm=trueparm)#trueparm= c(0.01,2.5, 0.007, 3, -1,3.8, -0.5,3.5 ,3.5,2.5)
                                        # lower ovarian baseline

r15 = sim5(N=400,B=25,trueparm=trueparm)
0.020 2.700 0.015 3.000 1.000 2.000 0.500 1.500 3.500 2.500

r16 = sim5(N=1000,B=10,trueparm=trueparm)#trueparm= c(0.009,2.7, 0.008, 3, 1,5, 0.3,5 ,3.5,2.5)

r17 = sim5(N=500,B=100,trueparm=trueparm)

#change in censorage smaller,
boxp=function(sr,c){
  boxplot(sr$simresult[,c])
  abline(h=trueparm[c])
}

ts1 = sim5(B=20,trueparm=trueparm, twostage = TRUE)

M5_F3_G0.total = rbind(M5_F3_G0.1$simresult, M5_F3_G0.2$simresult, M5_F3_G0.3$simresult,M5_F3_G0.4$simresult)
BiasSE(M5_F3_G0.total, trueparm) #test=M4_F3_G0.total[M4_F3_G0.total[,9]<50,];BiasSE(test,trueparm);penBiasSE(test,trueparm)
penBiasSE(M5_F3_G0.total, trueparm)

# mrate=0.5
M4_F3_G5.1  = sim4(B=30,trueparm=trueparm,mrate=0.5) 
M4_F3_G5.2  = sim4(B=30,trueparm=trueparm)
M4_F3_G5.3  = sim4(B=30,trueparm=trueparm)#

M4_F3_G5.total = rbind(M4_F3_G5.1$simresult)#, M4_F3_G5.2$simresult)
BiasSE(M4_F3_G5.total, trueparm) #test=M4_F3_G0.total[M4_F3_G0.total[,9]<50,];BiasSE(test,trueparm);penBiasSE(test,trueparm)
penBiasSE(M4_F3_G5.total, trueparm)

# mrate=0.8
F3_M8.1    = sim3(B=30,trueparm=trueparm,mrate=0.8) 
F3_M8.2    = sim3(B=40,trueparm=trueparm,mrate=0.8) # mrate=
F3_M8.total = rbind(F3_M8.1$simresult, F3_M8.2$simresult)
BiasSE(F3_M8.total, trueparm)
penBiasSE(F3_M8.total, trueparm)


#save result
write.excel(cbind(BiasSE(F300totla, trueparm),BiasSE(F3_M5.total, trueparm),BiasSE(F3_M8.total, trueparm)))
write.excel(cbind(penBiasSE(F300totla, trueparm),penBiasSE(F3_M5.total, trueparm),penBiasSE(F3_M8.total, trueparm)))
save(F300.total, F3_M5.total, file = "simjune22.RData"); 
load("my_data.RData")
###############################################################################################
###############################################################################################



############################ Model6, pop+, competing risk, time dependent covariate, 2frailty
##############################################################################################
tau=function(k){1/(1+2*k)}
BiasSE= function(estvec, trueparm){
  B= nrow(estvec)
  truemat = matrix(rep(trueparm,B), ncol=length(trueparm), nrow=B, byrow = TRUE)
  Bias=colMeans(estvec)-trueparm
  EmpSE= apply(estvec, 2, sd)
  RMSE= sqrt( colMeans( (estvec-truemat)^2 ))
  meansd=t(rbind(trueparm,Bias,EmpSE,RMSE))
  round(meansd,4)
}  #estvec is the simulation result
penBiasSE = function(){
  #ev=simlist$simresult
  ev=output
  ev[,c(1:4,8:10)] <- apply(ev[,c(1:4,8:10)], 2, exp)
  trueparm=trueparm
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
    integrate(integrand.tsc, X=X, est=est, frailty=frailty,lower=40,upper=t)$value
  }
  integrand.tsc <- function(u, X, est, type, frailty){
    u <- u
    fp=c(est[8],est[9])
    index=X
    xbeta1 =  index[2]*est[6]
    xbeta2 =  index[2]*est[7]
    
    haz1=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)*exp(est[5]*exp(-est[10]*(u-40))+est[11])
    H1 = - sapply(u,cH_CO,a=40,parm=c(est[1:2],est[5],est[10],est[11]))*exp(xbeta1)
    H2 = - (est[3]*u)^est[4]*exp(xbeta2)
    
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
      pen=penf(integrand,X=X,est=est,frailty=TRUE,t=40)+
          penf.tsc(integrand.tsc,X=X,est=est,frailty=TRUE,t=t)
    }
    pen
  }
  #t=45
  t= c(45,55)
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
    l[[1]]=estpen[,1];l[[2]]=estpen[,2]#;l[[3]]=estpen[,3];l[[4]]=estpen[,4]
    Bias = unlist(lapply(l,mean))-truepen
    EmpSE = unlist(lapply(l,sd))
    RMSE = sqrt( c(mean((l[[1]]-truepen[1])^2),mean((l[[2]]-truepen[2])^2) ))#,mean((l[[3]]-truepen[3])^2),mean((l[[4]]-truepen[4])^2)  ))
    total[[i]] = cbind(truepen,Bias,EmpSE,RMSE)
  }
  
  pen60=NULL;pen70=NULL
  for(i in 1:4){ #i=1 NS.NC i=2 S.NC i=3 NS.C i=4 S.C
    pen60=rbind(pen60,total[[i]][1,])
    pen70=rbind(pen70,total[[i]][2,])
  }
  out=rbind(pen60,pen70)  
  round(out,4)
}
sim6=function(B=30, N=300, trueparm, mrate=0, tvctype="PE", estimation="two"){
  iter = 0
  estvec=NULL
  while(iter < B){
     
    realdataestimate1[5]
    realdataestimate1[10]
    trueparm=realdataestimate1
    #trueparm=c(0.008,2.5,0.007,3,3,2,1.5,2.5,1.5,3,0.3)
    #trueparm=c(0.009, 2.6, 0.008, 3.2, 3.2, 2, 1.1, 2.5, 1, 3, 0.1)
    trueparm=c(0.01, 2.5, 0.01, 2.5, 2, 5, 5, 3.5, 2.5, 0.5, -0.5)
    #trueparm=c(0.009276306,2.556138491,0.008440465,3.194603196,3.284644247,3,3,2.798596325,1.509687011,3.194744800,0.119967604)
    round(trueparm,5) 
    0.00928 2.55614 0.00844 3.19460 3.28464 1.94388 0.82047 2.79860 1.50969 3.19474 0.11997)
    trueparm=c(0.008, 2.5, 0.007, 2.7, 3, 2, 1.1, 3.5, 2.5, 0.5, -0.5)
    simdata=simfam.cmp(N.fam = 300, trueparm, frailty.dist="gamma", mrate=0,tvctype="CO")
    simr4= Comp_simulation3(data=simdata, initparm=trueparm, tvctype="CO")#log alpha
    simr8= Comp_simulation3_one(data=simdata, initparm=trueparm, tvctype="CO")#log alpha
   simdata=simfam.cmp()
     backup=simdata_unbiasedbasline
  simdata_unbiasedbasline=simdata  
    table(simdata$status)[2:3]/nrow(simdata)
    table(b1$status)[2:3]/nrow(b1)#0.35,0.07
    table(simdata$status[simdata$proband==1])
    table(simdata$mgene)[2]/nrow(simdata)
    sum(b1$mgene==1)/length(b1$mgene)#0.38

    table(aggregate(simdata$status,by=list(simdata$famID),FUN=length)[,2])

    sum(simdata$TVCstatus)/length(simdata$TVCstatus)
    sum(b1$gender==1)/length(b1$gender)#0.21

  simdata=simdata_unbiasedbasline
    densityplot(simdata$tvc[simdata$TVCstatus==1],xlim=c(10,80))
    densityplot(b1$st1,xlim=c(10,80))

    densityplot(simdata$currentage[simdata$proband==1],xlim=c(10,80))
    densityplot(b1$currentage[b1$proband==1],xlim=c(10,80))

    densityplot(simdata$time[simdata$status==0])
    densityplot(b1$time[b1$status==0])

    densityplot(simdata$time[simdata$status==1],xlim=c(0,100))
    densityplot(b1$time[b1$status==1],xlim=c(0,100))

    densityplot(simdata$time[simdata$status==2])
    densityplot(b1$time[b1$status==2])

    densityplot(simdata$time[simdata$generation==3])
    densityplot(simdata$time[simdata$generation==2])
    densityplot(simdata$time[simdata$generation==1])

    table(aggregate(b1$status,by=list(b1$famID),FUN=length)[,2])

    densityplot(b1$time[b1$generation==3])
    densityplot(b1$time[b1$generation==2])
    densityplot(b1$time[b1$generation==1])

    #0.0093 2.6882 0.0079 3.1522 3.2975 2.1406 1.1128 2.4514 0.9540 3.2502 0.0885
    if(estimation=="two"){
    simr= Comp_simulation3(data=simdata, initparm=trueparm, tvctype="CO")#log alpha
    simr= Comp_simulation3_one(data=simdata, initparm=trueparm, tvctype="CO")#log alpha
    }else{
    }
    # simr= Comp_simulation3_one(data=simdata, initparm=trueparm, tvctype="ED")#log alpha
    # trueparm= c(0.015,2.5, 0.015,2.5, 2,5,5, 3.5,5,  0.5     )
    # simr2 = Comp_simulation3(data=simdata, initparm=trueparm, frailty=TRUE,tvctype="ED")
    # trueparm= c(0.015,2.5, 0.015,2.5, 2,5,5, 3.5,5           )
    # simr3 = Comp_simulation3(data=simdata, initparm=trueparm, frailty=TRUE,tvctype="PE")
    # 
     iter = iter+1

  #y= cbind(trueparm,simr,simr2)
    estvec = rbind(estvec,simr)
    if(iter%%1==0){
      print(iter)
      print(estvec[iter,])
    }
  }
  ori.trueparm = trueparm
  tra.trueparm = trueparm
  tra.trueparm[c(1:4,8:10)] = sapply(c(1:4,8:10), function(s) log(trueparm[s]))
  
  tra.truemat = matrix(rep(tra.trueparm,B), ncol=length(tra.trueparm), nrow=B, byrow = TRUE)
  Bias=colMeans(estvec)-tra.trueparm
  ESE= apply(estvec, 2, sd)
  RMSE= sqrt( colMeans( (estvec-tra.truemat)^2 ))

  meansd=cbind(ori.trueparm,tra.trueparm,Bias,ESE,RMSE)
  print(round(meansd,4))
  return(list(simresult=estvec, summary=meansd, ori.trueparm=ori.trueparm, tra.trueparm=tra.trueparm))
}
##################################################### simulation running & results#############
######## k = 2.5 ##############################################################################

trueparm= c(0.01,2.5, 0.009,3, 1,3,3, 3.5,2.5               ) ;mrate=0#PE&TI
trueparm= c(0.012,2.5, 0.011,3, 4,3,3, 3.5,2.5,  2          ) ;mrate=0#ED
trueparm= c(0.012,2.5, 0.01,3, 4,3,2.5, 3.5,2.5,  2,0.3      ) ;mrate=0#CO

# mrate=0    
test=sim6(B=1,trueparm=trueparm)
test2=sim6(B=1,trueparm=trueparm)
test3=sim6(B=1,trueparm=trueparm)
test4=sim6(B=1,trueparm=trueparm)
test5=sim6(B=10,trueparm=trueparm)
test6=sim6(B=10,trueparm=trueparm)
test7=sim6(B=10,trueparm=trueparm)
test7=sim6(B=10,trueparm=trueparm)
test8=sim6(B=10,trueparm=trueparm)#two stage
test9=sim6(B=10,trueparm=trueparm)#one stage

test10=sim6(B=10,trueparm=trueparm)#one stage
test11=sim6(B=10,trueparm=trueparm)#two stage
test12=sim6(B=10,trueparm=trueparm)#two stagereal
test13=sim6(B=10,trueparm=trueparm)#one
test14=sim6(B=10,trueparm=trueparm)#two

test15=sim6(B=10,trueparm=trueparm)#one
test16=sim6(B=10,trueparm=trueparm)#one
test17=sim6(B=10,trueparm=trueparm)#one

test19=sim6(B=10,trueparm=trueparm)#two

#work? 
#
# trueparm=c(0.009276306,2.556138491,0.008440465,3.194603196,3.284644247,1.943879186,0.820466155,2.798596325,1.509687011,3.194744800,0.119967604)
test20=sim6(B=30,N=500, trueparm=trueparm, estimation="two")#penetranceok
test21=sim6(B=10,N=500, trueparm=trueparm, estimation="one")#penetranceo bad
trueparm=c(0.009276306,2.556138491,0.008440465,3.194603196,3.284644247,1.943879186,0.820466155,2.798596325,1.509687011,3.194744800,0.119967604)
#k1 =5
test22=sim6(B=10,N=500, trueparm=trueparm, estimation="one")#penetranceok

#save result
write.excel(cbind(BiasSE(F300totla, trueparm),BiasSE(F3_M5.total, trueparm),BiasSE(F3_M8.total, trueparm)))
write.excel(cbind(penBiasSE(F300totla, trueparm),penBiasSE(F3_M5.total, trueparm),penBiasSE(F3_M8.total, trueparm)))
save(F300.total, F3_M5.total, file = "simjune22.RData"); 
load("my_data.RData")
###############################################################################################
###############################################################################################

0.009276306 2.556138491 0.008440465 3.194603196 3.284644247 1.943879186 0.820466155 2.798596325 1.509687011
3.194744800 0.119967604


###-----------------------------------------------Functions
# FITTING MODEL CODES
BC_simulation <- function(data = data, initparm=c(0.009,2.4424,1,1.5),missing.method = "data", base.dist = "Weibull", frailty=TRUE){
  
    #initial parameters
    base.parms <- initparm[1:2]
    vbeta <- initparm[3:4]
    theta = theta0 = c(log(base.parms),vbeta) 
    
    
    # setting up time dependent variables
    # if(time.dep.cov$time.dep.cov.type=="ED"){ 
    #   phi = c(3,2,1)
    #   theta = c(theta, log(phi))
    # }else if (time.dep.cov$time.dep.cov.type=="CO"){
    #   phi0 = c(0,0,0,0)
    #   phi1 = c(2,2,2,0.1)
    #   phi= c(phi0,phi1)
    #   theta = c(theta, phi0, log(phi1))
    # }else{phi=NULL}
    
    
    if(frailty==TRUE){
      datacopy = data
      fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
      df1 <- aggregate(datacopy$status==1, by=list(datacopy$famID), FUN=sum)[,2]
      data$df1 <- rep(df1, fsize)
      fp <- initparm[5]
      theta = theta0 = c(theta, log(fp))
    }else{fp=NULL}
    
    # data preparation
    # EM
    data.cooked = carrierprobgeno_NC(method = missing.method, data=data)  
    
    est0 <- est <- theta
    # dd <- lval0 <- lval <- 1#
    # i <- 0#
    # lval <- loglik_Comp_Timedep_NC_simulation(theta, data=data.cooked, agemin=16, frailty=frailty)#
    # plot(0:10,seq(-(lval-100),-lval,length.out=11),col="white")
    # points(0,-lval,col="red")
    # while(dd>0.01){#
    #   i <- i+1#
    #   est0 <- est#
    #   lval0 <- lval#
    # est1 <- optim(est0, loglik_Comp_Timedep_NC_simulation, data=data.cooked,
    #               agemin=16, frailty=frailty,
    #               hessian=FALSE, control=list(maxit=50000))
    #   lval <- est1$value#
    #   dd <- abs(lval0-lval)#
    #   est <- est1$par#
    #   #dd <- abs(sum(est-est0))#
    #   points(i,-lval,col="red")
    #   #print(c(i, dd, est1$par, est1$convergence))#
    # }#
    # cat("Iterations = ", i, "\n")#
    est1 =        optim(est, loglik_BC_simulation, data=data.cooked,  
                        agemin=16, frailty=frailty,
                        hessian=FALSE, control=list(maxit=50000))
    
    
    
    # Results 
    #print(est1$convergence)
    # Var <- try(solve(est1$hessian), TRUE)
    # se <- sqrt(diag(Var))
    # se.exp <- exp(est1$par)*se
    # SE <- c(se.exp[1:2] # 2 for base parms
    #         ,se[3:length(est1$par)]) # 5 for vbeta
    
    # # AIC
    # AIC <- 2*(length(est1$par)-(-est1$value))
    # Loglik_Value <- -est1$value
    
    # WALD & Pvalue
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
  
    # output
    # LT = length(theta)
    # Lf = length(fp)
    if(frailty==FALSE){
    mle <- c( exp(est1$par[1:2]),est1$par[3:4])#,exp(est1$par[5]))
    }else{
    mle <- c( exp(est1$par[1:2]),est1$par[3:4],exp(est1$par[5]))  
    }
    # )), 
    # as.matrix(SE),
    # as.matrix(pval)),4)
    # colnames(mle) <- c("MLE","SE","P-value")
    # rownames(mle) <- c("Lambda1","Rho1","Lambda2","Rho2","Lambda3","Rho3",
    #                       "Breast.Gene","Breast.Screen1","Breast.Screen2","Breast.Screen3","Breast.oopho",
    #                       #"Breast.S1xG","Breast.S2xG","Breast.S3xG",
    #                       #"Breast.S1xO","Breast.S2xO","Breast.S3xO", # full mode interaction
    #                       "Ovary.Gene","Ovary.Screen1",
    #                       "Death.Gene","Death.Screen1"
    #                       ,"frailty_k" 
    #                     )
    
    
    # cat("Printing Model estimates, standard error and P-values... \n")
    # print(mle)
    # cat("Resuling model has AIC value of ", AIC , "\n")
    
    
    # Mutation prediction diagnostics, prints out stats such as Sensitivity, specificity, AUC, NPV, PPV, ACC    
    # if(mutation.prediction==TRUE){
    #   cat("Computing Mutation covariate prediction based on model estimate...")
    #   Mut.Pred = mutation_prediction(data= data, parms = est1$par) 
    #   print(Mut.Pred)
    # }else{
    #   Mut.Pred = NULL
    # }
    # 
    return(mle)
  }
loglik_BC_simulation=function(theta, data, agemin, frailty){
  
  #print(round(c(exp(theta[1:2]),theta[3:length(theta)]),3))
  
  data = data[data$currentage>=agemin,]
  # base parameters
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
   
  
  # vbeta for breast  
  beta.sc = theta[3]
  beta.gen1 = theta[4]

  # frailty parameter
  if(frailty==TRUE){
    fp <- exp(theta[length(theta)])
  }
  gender= data$gender
  
  # Y, delta, carrier prob, proband indicator, sc number
  time0 = data$time-agemin
  status = data$status
  ip = which(data$proband==1)
  ex1 = data$carrp.geno 
  
  #setting up indicator variable for the removal status

  
  #loading up event time(sorted) vectors, event indicator vectors
  # time_vec2 <- data2[[1]]
  # indicator_vec <- data2[[2]]
  # time_vec2_p0 <- data2[[3]]
  # indicator_vec_p0 <- data2[[4]]
  # time_vec2_p1 <- data2[[5]]
  # indicator_vec_p1 <- data2[[6]]
  
  # baseline hazard(weibull), baseline cumulative hazard(weibull)
  bhaz1 = (lambda1^rho1)*rho1*time0^(rho1-1)
  bcumhaz1 = (lambda1*time0)^rho1 
  
  logh1.0 = log(bhaz1) 
  logh1.1 = log(bhaz1) + beta.gen1 
  
  # creating screen,gene,oophorectomy interaction vectors for hazard
  
  # sum of log-hazard
  sum1 = sum((      ex1*(logh1.1 + gender*beta.sc ) +
                      (1-ex1)*(logh1.0 + gender*beta.sc))[status==1], na.rm=TRUE) 
  
  
  # sum of survival until event time for all the individuals in the data

      CH0_bc <- -bcumhaz1*exp(gender*beta.sc)*(1-ex1)
      CH1_bc <- -bcumhaz1*exp(beta.gen1 + gender*beta.sc)*(ex1)
      k <- CH0_bc + CH1_bc
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam1 <- -(CH0_bc+CH1_bc)
    df1 <- data$df1[data$proband==1]
    Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
    sum4.1 <- sum(lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1]), na.rm=T)
    loglik = sum1+sum4.1
  }else{
    sum4.1 = sum(k)
    loglik = sum1+sum4.1
  }
  # Ascertainment correction by design="pop+"
  cagep <- data$currentage[ip]-agemin
  timep <- data$time[ip]-agemin
  # proband disease status at the study entry
  
  # subsetting the probands according to his affection status at study entry
  # cagep0 <- cagep[statusp==0]
  # timep0 <- timep[statusp==0]
  # wtp0 <- wtp[statusp==0]
  # 
  # mat_p0 <- cbind(time_vec2_p0, indicator_vec_p0, cagep0)
  # if(reccur.type=="CM"){
  #   logasc0 <- apply(mat_p0, 1, func_p0, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
  # }else{
  #   if(time.dep.cov.type=="PE"){
  #     #logasc0 <- apply(mat_p0, 1, func_p0_IS, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
  #     testm <- mat_p0
  #     testm[is.na(testm)] <- 99
  #     CH_bc=CumH_c_ED(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=1, phi =c(0,0,0,0))
  #     CH1_bc <- exp(beta.gen1)*CH_bc[,2]
  #     logasc0 <- CH1_bc    
  #   } else{
  #     #logasc0 <- apply(mat_p0, 1, func_p0_IS_ED, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
  #     testm <- mat_p0
  #     testm[is.na(testm)] <- 99
  #     CH_bc=CumH_c_ED(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=1, phi = phiv)
  #     CH1_bc <- exp(beta.gen1)*CH_bc[,2]
  #     logasc0 <- CH1_bc 
  #   }
  # }
  # if(frailty==TRUE){
  #   logS <- log((1-CH1_bc/fp[1])^(-fp[1]))
  # }else {
  #   logS <- logasc0
  # }
  # logasc0 <- sum(wtp0*logS,na.rm=TRUE)
  
  # subsetting the probands according to his affection status at study entry
  CH1_bc <- -exp(beta.gen1+gender[ip]*beta.sc)*(lambda1*cagep)^rho1 
  logasc1 <- CH1_bc
  
  
  if(frailty==TRUE){
    S <-(1-logasc1/fp[1])^(-fp[1])
  }else {
    S <- exp(logasc1)
  }
  logasc1 <- sum((log(1-S)),na.rm=TRUE)
  
  
  slogasc = logasc1
  likelihood  <- loglik - slogasc
  
  #print(c(loglik, slogasc, likelihood))
  return(-likelihood)
  
} 
Comp_simulation <- function(data = data, initparm=c(0.009,2.4424,1,1.5),missing.method = "data", base.dist = "Weibull", frailty=TRUE){
  
  #initial parameters
  base.parms <- initparm[1:4]
  vbeta <- initparm[5:8]
  theta = theta0 = c(log(base.parms),vbeta) 
  
  
  # setting up time dependent variables
  # if(time.dep.cov$time.dep.cov.type=="ED"){ 
  #   phi = c(3,2,1)
  #   theta = c(theta, log(phi))
  # }else if (time.dep.cov$time.dep.cov.type=="CO"){
  #   phi0 = c(0,0,0,0)
  #   phi1 = c(2,2,2,0.1)
  #   phi= c(phi0,phi1)
  #   theta = c(theta, phi0, log(phi1))
  # }else{phi=NULL}
  
  
  if(frailty==TRUE){
    datacopy = data
    fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
    df1 <- aggregate(datacopy$status!=0, by=list(datacopy$famID), FUN=sum)[,2]
    #df2 <- aggregate(datacopy$status==2, by=list(datacopy$famID), FUN=sum)[,2]
    data$df1 <- rep(df1, fsize)
    #data$df2 <- rep(df2, fsize)
    #fp <- c(initparm[9],initparm[10])
    fp <- initparm[9]
    theta = theta0 = c(theta, log(fp))
  }else{fp=NULL}
  
  # data preparation
  # EM
  data.cooked = carrierprobgeno(method = missing.method, data=data)  
  
  est0 <- est <- theta
  # dd <- lval0 <- lval <- 1#
  # i <- 0#
  # lval <- loglik_Comp_Timedep_NC_simulation(theta, data=data.cooked, agemin=16, frailty=frailty)#
  # plot(0:10,seq(-(lval-100),-lval,length.out=11),col="white")
  # points(0,-lval,col="red")
  # while(dd>0.01){#
  #   i <- i+1#
  #   est0 <- est#
  #   lval0 <- lval#
  # est1 <- optim(est0, loglik_Comp_Timedep_NC_simulation, data=data.cooked,
  #               agemin=16, frailty=frailty,
  #               hessian=FALSE, control=list(maxit=50000))
  #   lval <- est1$value#
  #   dd <- abs(lval0-lval)#
  #   est <- est1$par#
  #   #dd <- abs(sum(est-est0))#
  #   points(i,-lval,col="red")
  #   #print(c(i, dd, est1$par, est1$convergence))#
  # }#
  # cat("Iterations = ", i, "\n")#
  est1 =        optim(est, loglik_Comp_simulation, data=data.cooked,  
                      agemin=16, frailty=frailty,
                      hessian=FALSE, control=list(maxit=2000))
  
  
  
  # Results 
  #print(est1$convergence)
  # Var <- try(solve(est1$hessian), TRUE)
  # se <- sqrt(diag(Var))
  # se.exp <- exp(est1$par)*se
  # SE <- c(se.exp[1:2] # 2 for base parms
  #         ,se[3:length(est1$par)]) # 5 for vbeta
  
  # # AIC
  # AIC <- 2*(length(est1$par)-(-est1$value))
  # Loglik_Value <- -est1$value
  
  # WALD & Pvalue
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
  
  # output
  # LT = length(theta)
  # Lf = length(fp)
  if(frailty==FALSE){
    mle <- c( exp(est1$par[1:4]),est1$par[5:8], 10000)#,exp(est1$par[5]))
  }else{
    #mle <- c( exp(est1$par[1:4]),est1$par[5:8],exp(est1$par[5]))  
    mle <- c( exp(est1$par[1:4]) , est1$par[5:8],exp(est1$par[9])) 
  }
  # )), 
  # as.matrix(SE),
  # as.matrix(pval)),4)
  # colnames(mle) <- c("MLE","SE","P-value")
  # rownames(mle) <- c("Lambda1","Rho1","Lambda2","Rho2","Lambda3","Rho3",
  #                       "Breast.Gene","Breast.Screen1","Breast.Screen2","Breast.Screen3","Breast.oopho",
  #                       #"Breast.S1xG","Breast.S2xG","Breast.S3xG",
  #                       #"Breast.S1xO","Breast.S2xO","Breast.S3xO", # full mode interaction
  #                       "Ovary.Gene","Ovary.Screen1",
  #                       "Death.Gene","Death.Screen1"
  #                       ,"frailty_k" 
  #                     )
  
  
  # cat("Printing Model estimates, standard error and P-values... \n")
  # print(mle)
  # cat("Resuling model has AIC value of ", AIC , "\n")
  
  
  # Mutation prediction diagnostics, prints out stats such as Sensitivity, specificity, AUC, NPV, PPV, ACC    
  # if(mutation.prediction==TRUE){
  #   cat("Computing Mutation covariate prediction based on model estimate...")
  #   Mut.Pred = mutation_prediction(data= data, parms = est1$par) 
  #   print(Mut.Pred)
  # }else{
  #   Mut.Pred = NULL
  # }
  # 
  return(mle)
}
loglik_Comp_simulation=function(theta, data, agemin, frailty){
 
  data = data[data$currentage>=agemin,]
  # base parameters
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])

  
  # vbeta for breast  
  beta.sc1_1 = theta[5]
  beta.gen1 = theta[6]

  # vbeta for ovarian
  beta.sc1_2 =  theta[7]
  beta.gen2 = theta[8]


  # time varying covariate type
  # if(time.dep.cov.type=="ED"){
  #   phiv <- exp(theta[c(16,17,18,19)])
  #   phi_sc = phiv[c(1:3)]; phi_sc = c(phi_sc,rep(0,3))
  #   phi_or = phiv[4]; phi_or = c(phi_or, 0)
  #   #phi_or = 0
  # }else if(time.dep.cov.type=="CO"){
  #   phiv = c( exp(theta[c(20,21,22,23)]) , theta[c(16,17,18,19)] )
  #   phi_sc = phiv[c(1:3,5:7)]
  #   phi_or = phiv[c(4,8)]
  #   #phi_or = 0
  # }
  
  # frailty parameter
  if(frailty==TRUE){
    fp <- exp(theta[length(theta)])
  }else{
    fp=10000
  }

  
  # Y, delta, carrier prob, proband indicator, sc number
  time0 = data$time-agemin
  status = data$status
  ip = which(data$proband==1)
  ex1 = data$carrp.geno 
  gender = data$gender

  # baseline hazard(weibull), baseline cumulative hazard(weibull)
  bhaz1 = (lambda1^rho1)*rho1*time0^(rho1-1)
  bcumhaz1 = (lambda1*time0)^rho1 
  bhaz2 = (lambda2^rho2)*rho2*time0^(rho2-1)
  bcumhaz2 = (lambda2*time0)^rho2 

  
  logh1.0 = log(bhaz1) 
  logh1.1 = log(bhaz1) + beta.gen1
  logh2.0 = log(bhaz2) 
  logh2.1 = log(bhaz2) + beta.gen2

  
  # creating screen,gene,oophorectomy interaction vectors for hazard

  # sum of log-hazard
  sum1 = sum((      ex1*(logh1.1 + gender*beta.sc1_1 ) +
                (1-ex1)*(logh1.0 + gender*beta.sc1_1))[status==1], na.rm=TRUE) 
  sum2 = sum((      ex1*(logh2.1 + gender*beta.sc1_2 ) +
                (1-ex1)*(logh2.0 + gender*beta.sc1_2))[status==2], na.rm=TRUE)              
  
  
  # sum of survival until event time for all the individuals in the data
  CH0_bc <- -bcumhaz1*exp(gender*beta.sc1_1)*(1-ex1)
  CH1_bc <- -bcumhaz1*exp(beta.gen1 + gender*beta.sc1_1)*(ex1)
  CH0_ov <- -bcumhaz2*exp(gender*beta.sc1_2)*(1-ex1)
  CH1_ov <- -bcumhaz2*exp(beta.gen2 + gender*beta.sc1_2)*(ex1)
  
  
  k <- CH0_bc + CH1_bc + CH0_ov + CH1_ov
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam1 <- -(CH0_bc+CH1_bc)
    Hfam2 <- -(CH0_ov+CH1_ov)
    Hfam=Hfam1+Hfam2 
    df1 <- data$df1[data$proband==1]
    #df2 <- data$df2[data$proband==1]

    # Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
    # Hfam2 <- aggregate(Hfam2,by=list(data$famID),FUN=sum)[,2]
     Hfam <- aggregate(Hfam,by=list(data$famID),FUN=sum)[,2]
  
       sum4.1 <- sum((lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam)/fp[1])), na.rm=T)
    #sum4.2 <- sum((lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2])), na.rm=T)

    loglik = sum1+sum2 + sum4.1
  }else{
    sum4 <- sum(k,na.rm = TRUE)
    loglik = sum1+sum2+sum4 # numerator in loglikelihood
  }
  
  # Ascertainment correction by design="pop+"
  penf <- function(f,S,est,fp,type,frailty,lower,upper){  
    log(integrate(f,S=S,est=est,fp=fp,type=type,frailty=frailty,lower=lower,upper=upper)$value)
  }
  integrand <- function(u,S,est,fp,type,frailty){
    u <- u
    
    xbeta1 = S*est[5] + est[6]
    xbeta2 = S*est[7] + est[8]
    if(type==1){
      haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
      #print(c(xbeta, est[7]+sum(beta_vec1)))
    }else if(type==2){
      haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)  
    }
    
    H1 = - (est[1]*u)^est[2]*exp(xbeta1)
    H2 = - (est[3]*u)^est[4]*exp(xbeta2)
    
    
    if(frailty==TRUE){
    if(type==1){
      #f = haz*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
      f = haz*(( 1 -  (H1+H2)/fp  ))^(-fp)
    }else{
      #f = haz*((1-H2/fp[2])^(-fp[2]-1))*((1-H1/fp[1])^(-fp[1]))
      f = haz*(( 1 -  (H1+H2)/fp  ))^(-fp)
    }
    } else{
      Hp <- H1+H2
      f = haz*exp(Hp) # h * S
    }
    return(f)
  }
  cagep <- data$currentage[ip]-agemin
  timep <- data$time[ip]-agemin
  statusp <- data$status[ip]  # proband disease status at the study entry
  genderp = gender[ip]
  
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.sc1_2,beta.gen2)
  dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.sc1_2,beta.gen2)
  print(c(dest,fp))
  
  #print(round(c(cest,fp),3))
  cagep1 <- cagep[statusp==1]
  timep1 <- timep[statusp==1]
  genderp1 <- genderp[statusp==1]

  logasc1 = sapply(1:length(cagep1), function(index) penf(integrand,S=genderp1[index],est=cest, fp=fp,frailt=frailty, type=1,lower=0,upper=cagep1[index]))  
  logasc1 <- sum(logasc1,na.rm=TRUE)
  
  # subsetting the probands according to his affection status at study entry
  cagep2 <- cagep[statusp==2]
  timep2 <- timep[statusp==2]
  genderp2 <- genderp[statusp==2]

  logasc2 = sapply(1:length(cagep2), function(index) penf(integrand,S=genderp2[index],est=cest, fp=fp,frailty=frailty, type=2,lower=0,upper=cagep2[index])) 
  logasc2 <- sum(logasc2,na.rm=TRUE)
  
  slogasc = logasc1 + logasc2
  likelihood  <- loglik - slogasc
  
  return(-likelihood)
} 

initparm=trueparm;data=simdata;missing.method="data"
Comp_simulation2 <- function(data = data, initparm=c(0.009,2.4424,1,1.5),missing.method = "data", base.dist = "Weibull", frailty=TRUE){
  
  #initial parameters
  base.parms <- initparm[1:4]
  vbeta <- initparm[5:8]
  theta = theta0 = c(log(base.parms),vbeta) 
  
  
  # setting up time dependent variables
  # if(time.dep.cov$time.dep.cov.type=="ED"){ 
  #   phi = c(3,2,1)
  #   theta = c(theta, log(phi))
  # }else if (time.dep.cov$time.dep.cov.type=="CO"){
  #   phi0 = c(0,0,0,0)
  #   phi1 = c(2,2,2,0.1)
  #   phi= c(phi0,phi1)
  #   theta = c(theta, phi0, log(phi1))
  # }else{phi=NULL}
  
  
  if(frailty==TRUE){
    datacopy = data
    fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
    df1 <- aggregate(datacopy$status==1, by=list(datacopy$famID), FUN=sum)[,2]
    df2 <- aggregate(datacopy$status==2, by=list(datacopy$famID), FUN=sum)[,2]
    data$df1 <- rep(df1, fsize)
    data$df2 <- rep(df2, fsize)
    fp <- c(initparm[9],initparm[10])
    
    theta = theta0 = c(theta, log(fp))
  }else{fp=NULL}
  
  # data preparation
  # EM
  data.cooked = carrierprobgeno(method = missing.method, data=data)  
  
  est0 <- est <- theta
  # dd <- lval0 <- lval <- 1#
  # i <- 0#
  # lval <- loglik_Comp_Timedep_NC_simulation(theta, data=data.cooked, agemin=16, frailty=frailty)#
  # plot(0:10,seq(-(lval-100),-lval,length.out=11),col="white")
  # points(0,-lval,col="red")
  # while(dd>0.01){#
  #   i <- i+1#
  #   est0 <- est#
  #   lval0 <- lval#
  # est1 <- optim(est0, loglik_Comp_Timedep_NC_simulation, data=data.cooked,
  #               agemin=16, frailty=frailty,
  #               hessian=FALSE, control=list(maxit=50000))
  #   lval <- est1$value#
  #   dd <- abs(lval0-lval)#
  #   est <- est1$par#
  #   #dd <- abs(sum(est-est0))#
  #   points(i,-lval,col="red")
  #   #print(c(i, dd, est1$par, est1$convergence))#
  # }#
  # cat("Iterations = ", i, "\n")#
  est1 =        optim(est, loglik_Comp_simulation2, data=data.cooked,  
                      agemin=16, frailty=frailty,
                      hessian=FALSE, control=list(maxit=10000))
  
  
  
  # Results 
  #print(est1$convergence)
  # Var <- try(solve(est1$hessian), TRUE)
  # se <- sqrt(diag(Var))
  # se.exp <- exp(est1$par)*se
  # SE <- c(se.exp[1:2] # 2 for base parms
  #         ,se[3:length(est1$par)]) # 5 for vbeta
  
  # # AIC
  # AIC <- 2*(length(est1$par)-(-est1$value))
  # Loglik_Value <- -est1$value
  
  # WALD & Pvalue
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
  
  # output
  # LT = length(theta)
  # Lf = length(fp)
  if(frailty==FALSE){
    mle <- c( exp(est1$par[1:4]),est1$par[5:8], 77,77)#,exp(est1$par[5]))
  }else{
    #mle <- c( exp(est1$par[1:4]),est1$par[5:8],exp(est1$par[5]))  
    mle <- c( exp(est1$par[1:4]) , est1$par[5:8],exp(c(est1$par[9],est1$par[10]))) 
  }
  # )), 
  # as.matrix(SE),
  # as.matrix(pval)),4)
  # colnames(mle) <- c("MLE","SE","P-value")
  # rownames(mle) <- c("Lambda1","Rho1","Lambda2","Rho2","Lambda3","Rho3",
  #                       "Breast.Gene","Breast.Screen1","Breast.Screen2","Breast.Screen3","Breast.oopho",
  #                       #"Breast.S1xG","Breast.S2xG","Breast.S3xG",
  #                       #"Breast.S1xO","Breast.S2xO","Breast.S3xO", # full mode interaction
  #                       "Ovary.Gene","Ovary.Screen1",
  #                       "Death.Gene","Death.Screen1"
  #                       ,"frailty_k" 
  #                     )
  
  
  # cat("Printing Model estimates, standard error and P-values... \n")
  # print(mle)
  # cat("Resuling model has AIC value of ", AIC , "\n")
  
  
  # Mutation prediction diagnostics, prints out stats such as Sensitivity, specificity, AUC, NPV, PPV, ACC    
  # if(mutation.prediction==TRUE){
  #   cat("Computing Mutation covariate prediction based on model estimate...")
  #   Mut.Pred = mutation_prediction(data= data, parms = est1$par) 
  #   print(Mut.Pred)
  # }else{
  #   Mut.Pred = NULL
  # }
  # 
  return(mle)
}
Comp_simulation3_one <- function(data = data, initparm=c(0.009,2.4424,1,1.5),missing.method = "data", base.dist = "Weibull", tvctype="PE"){
  
  #initial parameters
  base.parms <- initparm[1:4]
  vbeta <- initparm[5:7]
  theta = theta0 = c(log(base.parms),vbeta) 
  
  
  datacopy = data
  fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
  df1 <- aggregate(datacopy$status==1, by=list(datacopy$famID), FUN=sum)[,2]
  df2 <- aggregate(datacopy$status==2, by=list(datacopy$famID), FUN=sum)[,2]
  data$df1 <- rep(df1, fsize)
  data$df2 <- rep(df2, fsize)
  fp <- c(initparm[8:9])
  
   theta = theta0 = c(theta, log(fp))
  
  
  if(tvctype=="ED")theta <- c(theta,log(initparm[10]))
  if(tvctype=="CO")theta <- c(theta,log(initparm[10]),initparm[11])
  # data preparation
  # EM
  data.cooked = carrierprobgeno(method = missing.method, data=data)  
  
  est0 <- est <- theta
 
  est1 =        optim(est, loglik_Comp_simulation3_one, data=data.cooked,
                      agemin=16, frailty=TRUE,tvctype=tvctype,
                       control=list(maxit=100000) )
  # est1 =        nlminb(est,loglik_Comp_simulation3_one,  data=data.cooked,
  #                     agemin=16, frailty=TRUE,tvctype=tvctype,
  #                     hessian=FALSE)#,steptol=1e-5,gradtol=1e-5)

  
  if(tvctype=="TI")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:9]))
  if(tvctype=="PE")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:9]))
  if(tvctype=="ED")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:10]))
  if(tvctype=="CO")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est1$par[8:9]),exp(est1$par[10]),est1$par[11])
  #if(tvctype=="CO")mle = c( est1$par)
  
  return(mle)
}
Comp_simulation3 <- function(data = data, initparm=c(0.009,2.4424,1,1.5),missing.method = "data", base.dist = "Weibull", tvctype="PE"){
  
  #initial parameters
  base.parms <- initparm[1:4]
  vbeta <- initparm[5:7]
  theta = theta0 = c(log(base.parms),vbeta) 
  

    datacopy = data
    fsize <- aggregate(datacopy$status, by=list(datacopy$famID), length)[,2]
    df1 <- aggregate(datacopy$status==1, by=list(datacopy$famID), FUN=sum)[,2]
    df2 <- aggregate(datacopy$status==2, by=list(datacopy$famID), FUN=sum)[,2]
    data$df1 <- rep(df1, fsize)
    data$df2 <- rep(df2, fsize)
    fp <- c(initparm[8:9])
    
   # theta = theta0 = c(theta, log(fp))
  
  
  if(tvctype=="ED")theta <- c(theta,log(initparm[10]))
  if(tvctype=="CO")theta <- c(theta,log(initparm[10]),initparm[11])
  # data preparation
  # EM
  data.cooked = carrierprobgeno(method = missing.method, data=data)  
  
  est0 <- est <- theta

  est1 =        optim(est, loglik_Comp_simulation3, data=data.cooked,  
                      agemin=16, frailty=FALSE,tvctype=tvctype,
                      hessian=FALSE, control=list(maxit=100000))
  est2 =        optim(log(fp), loglik_Comp_simulation3_second, data=data.cooked,  
                      agemin=16, frailty=TRUE,tvctype=tvctype,stage1=est1$par,
                      hessian=FALSE, control=list(maxit=100000))
  print(est2$convergence)


    if(tvctype=="TI")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est2$par[1:2]))
    if(tvctype=="PE")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est2$par[1:2]))
    if(tvctype=="ED")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est2$par[1:2]),exp(est1$par[8]))
    if(tvctype=="CO")mle = c( exp(est1$par[1:4]) , est1$par[5:7],exp(est2$par[1:2]),exp(est1$par[8]),est1$par[9])
    #if(tvctype=="CO")mle = c( est1$par[1:7],est2$par[1:2],est1$par[8:9])

  return(mle)
}
data=data.cooked
loglik_Comp_simulation2=function(theta, data, agemin, frailty){
  
  data = data[data$currentage>=agemin,]
  # base parameters
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  
  
  # vbeta for breast  
  beta.sc1_1 = theta[5]
  beta.gen1 = theta[6]
  
  # vbeta for ovarian
  beta.sc1_2 =  theta[7]
  beta.gen2 = theta[8]
  
  
  # time varying covariate type
  # if(time.dep.cov.type=="ED"){
  #   phiv <- exp(theta[c(16,17,18,19)])
  #   phi_sc = phiv[c(1:3)]; phi_sc = c(phi_sc,rep(0,3))
  #   phi_or = phiv[4]; phi_or = c(phi_or, 0)
  #   #phi_or = 0
  # }else if(time.dep.cov.type=="CO"){
  #   phiv = c( exp(theta[c(20,21,22,23)]) , theta[c(16,17,18,19)] )
  #   phi_sc = phiv[c(1:3,5:7)]
  #   phi_or = phiv[c(4,8)]
  #   #phi_or = 0
  # }
  
  # frailty parameter
  if(frailty==TRUE){
    fp <- exp(theta[c(9,10)])
    
  }else{
    fp=10000
  }
  
  
  # Y, delta, carrier prob, proband indicator, sc number
  time0 = data$time-agemin
  status = data$status
  ip = which(data$proband==1)
  ex1 = data$carrp.geno 
  gender = data$gender
  
  # baseline hazard(weibull), baseline cumulative hazard(weibull)
  bhaz1 = (lambda1^rho1)*rho1*time0^(rho1-1)
  bcumhaz1 = (lambda1*time0)^rho1 
  bhaz2 = (lambda2^rho2)*rho2*time0^(rho2-1)
  bcumhaz2 = (lambda2*time0)^rho2 
  
  
  logh1.0 = log(bhaz1) 
  logh1.1 = log(bhaz1) + beta.gen1
  logh2.0 = log(bhaz2) 
  logh2.1 = log(bhaz2) + beta.gen2
  
  
  # creating screen,gene,oophorectomy interaction vectors for hazard
  
  # sum of log-hazard
  sum1 = sum((      ex1*(logh1.1 + gender*beta.sc1_1 ) +
                      (1-ex1)*(logh1.0 + gender*beta.sc1_1))[status==1], na.rm=TRUE) 
  sum2 = sum((      ex1*(logh2.1 + gender*beta.sc1_2 ) +
                      (1-ex1)*(logh2.0 + gender*beta.sc1_2))[status==2], na.rm=TRUE)              
  
  
  # sum of survival until event time for all the individuals in the data
  CH0_bc <- -bcumhaz1*exp(gender*beta.sc1_1)*(1-ex1)
  CH1_bc <- -bcumhaz1*exp(beta.gen1 + gender*beta.sc1_1)*(ex1)
  CH0_ov <- -bcumhaz2*exp(gender*beta.sc1_2)*(1-ex1)
  CH1_ov <- -bcumhaz2*exp(beta.gen2 + gender*beta.sc1_2)*(ex1)
  
  
  k <- CH0_bc + CH1_bc + CH0_ov + CH1_ov
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam1 <- -(CH0_bc+CH1_bc)
    Hfam2 <- -(CH0_ov+CH1_ov)

    df1 <- data$df1[data$proband==1]
    df2 <- data$df2[data$proband==1]
    
     Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
     Hfam2 <- aggregate(Hfam2,by=list(data$famID),FUN=sum)[,2]
    
    
    sum4.1 <- sum((lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1])), na.rm=T)
    sum4.2 <- sum((lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2])), na.rm=T)
    
    loglik = sum1+sum2 + sum4.1+ sum4.2
  }else{
    sum4 <- sum(k,na.rm = TRUE)
    loglik = sum1+sum2+sum4 # numerator in loglikelihood
  }
  
  # Ascertainment correction by design="pop+"
  penf <- function(f,S,est,fp,type,frailty,lower,upper){  
    log(integrate(f,S=S,est=est,fp=fp,type=type,frailty=frailty,lower=lower,upper=upper)$value)
  }
  integrand <- function(u,S,est,fp,type,frailty){
    u <- u
    
    xbeta1 = S*est[5] + est[6]
    xbeta2 = S*est[7] + est[8]
    if(type==1){
      haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
      #print(c(xbeta, est[7]+sum(beta_vec1)))
    }else if(type==2){
      haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)  
    }
    
    H1 = - (est[1]*u)^est[2]*exp(xbeta1)
    H2 = - (est[3]*u)^est[4]*exp(xbeta2)
    
    
    if(frailty==TRUE){
      if(type==1){
        f = haz*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
        #f = haz*(( 1 -  (H1+H2)/fp  ))^(-fp)
      }else{
        f = haz*((1-H2/fp[2])^(-fp[2]-1))*((1-H1/fp[1])^(-fp[1]))
        #f = haz*(( 1 -  (H1+H2)/fp  ))^(-fp)
      }
    } else{
      Hp <- H1+H2
      f = haz*exp(Hp) # h * S
    }
    return(f)
  }
  cagep <- data$currentage[ip]-agemin
  timep <- data$time[ip]-agemin
  statusp <- data$status[ip]  # proband disease status at the study entry
  genderp = gender[ip]
  
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.sc1_2,beta.gen2)
  dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.sc1_2,beta.gen2)
  print(round(c(dest,fp),5))

  #print(round(c(cest,fp),3))
  cagep1 <- cagep[statusp==1]
  timep1 <- timep[statusp==1]
  genderp1 <- genderp[statusp==1]
  
  logasc1 = sapply(1:length(cagep1), function(index) penf(integrand,S=genderp1[index],est=cest, fp=fp,frailt=frailty, type=1,lower=0,upper=cagep1[index]))  
  logasc1 <- sum(logasc1,na.rm=TRUE)
  
  # subsetting the probands according to his affection status at study entry
  cagep2 <- cagep[statusp==2]
  timep2 <- timep[statusp==2]
  genderp2 <- genderp[statusp==2]
  
  logasc2 = sapply(1:length(cagep2), function(index) penf(integrand,S=genderp2[index],est=cest, fp=fp,frailty=frailty, type=2,lower=0,upper=cagep2[index])) 
  logasc2 <- sum(logasc2,na.rm=TRUE)
  
  slogasc = logasc1 + logasc2
  likelihood  <- loglik - slogasc
  
  return(-likelihood)
} 
loglik_Comp_simulation3_one=function(theta, data, agemin, frailty,tvctype){
  
  #data = data[data$currentage>=agemin,]
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
    eta <- exp(theta[10])
  }else if(tvctype=="CO"){
    eta <- c(exp(theta[10]),theta[11])
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
  #ex1 = data$carrp.geno 
  #gender = data$gender
  mgene=data$mgene
  st1=data$tvc-agemin
  gender=data$TVCstatus
  st1[gender==0]<-time0[gender==0]
  
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
    
    Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
    Hfam2 <- aggregate(Hfam2,by=list(data$famID),FUN=sum)[,2]
    
    
    sum4.1 <- sum((lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1])), na.rm=T)
    sum4.2 <- sum((lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2])), na.rm=T)
    
    loglik = sum1+sum2 + sum4.1+ sum4.2
  }else{
    sum4 <- sum(-k,na.rm = TRUE)
    loglik = sum1+sum2+sum4 # numerator in loglikelihood
  }
  
  # Ascertainment correction by design="pop+"
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.gen2)
  cagep <- data$currentage[ip]-agemin
  statusp <- data$status[ip]  # proband disease status at the study entry
  genderp = gender[ip]
  st1p=st1[ip]
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  
  if(tvctype=="PE"|tvctype=="TI")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp)
  if(tvctype=="ED"|tvctype=="CO")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp,eta)
  print(round(dest,4))
  
  
  # sum of survival until event time for all the individuals in the data
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
    slogasc=sum(log(1- ((1+Hp1/fp[1])^(-fp[1]))*((1+Hp2/fp[2])^(-fp[2]))),na.rm=TRUE)
  }else{
    slogasc=sum(log(1- exp(-Hp1-Hp2)),na.rm=TRUE)  
  }
  
  #print(round(c(cest,fp),3))
  # cagep1.0 <- cagep[statusp==1&genderp==0]
  # cagep1.1 <- cagep[statusp==1&genderp==1]
  # st1p1<-st1p[statusp==1&genderp==1]
  # 
  # if(length(cagep1.0)!=0){
  # logasc1.0 = log(sapply(1:length(cagep1.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep1.0[index])))
  # }else{logasc1.0=0}
  # if(length(cagep1.1)!=0){
  # logasc1.1 = log(sapply(1:length(cagep1.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p1[index]))+
  #                 sapply(1:length(cagep1.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=1,tvc=st1p1[index],tt=tvctype,eta=eta,lower=st1p1[index],upper=cagep1.1[index])))
  # }else{logasc1.1=0}  
  # logasc1 <- sum(logasc1.0,na.rm=TRUE)+sum(logasc1.1,na.rm=TRUE)
  # 
  # # subsetting the probands according to his affection status at study entry
  # cagep2.0 <- cagep[statusp==2&genderp==0]
  # cagep2.1 <- cagep[statusp==2&genderp==1]
  # st1p2<-st1p[statusp==2&genderp==1]
  # 
  # if(length(cagep2.0)!=0){
  # logasc2.0 = log(sapply(1:length(cagep2.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep2.0[index])))
  # }else{logasc2.0=0}
  # if(length(cagep2.1)!=0){
  # logasc2.1 = log(sapply(1:length(cagep2.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p2[index]))+
  #             sapply(1:length(cagep2.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=2,tvc=st1p2[index],tt=tvctype,eta=eta,lower=st1p2[index],upper=cagep2.1[index])))
  # }else{logasc2.1=0}
  # logasc2 <- sum(logasc2.0,na.rm=TRUE)+sum(logasc2.1,na.rm=TRUE)
  
  #slogasc = logasc1 + logasc2
  likelihood  <- loglik - slogasc
  
  return(-likelihood)
}
loglik_Comp_simulation3=function(theta, data, agemin, frailty,tvctype){
    
    data = data[data$currentage>=agemin,]
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
      eta <- exp(theta[8])
    }else if(tvctype=="CO"){
      eta <- c(exp(theta[8]),theta[9])
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
    #ex1 = data$carrp.geno 
    #gender = data$gender
    mgene=data$mgene
    st1=data$tvc-agemin
    gender=data$TVCstatus
    st1[gender==0]<-time0[gender==0]
    
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
      
      Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
      Hfam2 <- aggregate(Hfam2,by=list(data$famID),FUN=sum)[,2]
      
      
      sum4.1 <- sum((lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1])), na.rm=T)
      sum4.2 <- sum((lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2])), na.rm=T)
      
      loglik = sum1+sum2 + sum4.1+ sum4.2
    }else{
      sum4 <- sum(-k,na.rm = TRUE)
      loglik = sum1+sum2+sum4 # numerator in loglikelihood
    }
    
    # Ascertainment correction by design="pop+"
    cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.gen2)
    cagep <- data$currentage[ip]-agemin
    statusp <- data$status[ip]  # proband disease status at the study entry
    genderp = gender[ip]
    st1p=st1[ip]
    # subsetting the probands according to his affection status at study entry
    # subsetting the probands according to his affection status at study entry
    # h1*S
    
    if(tvctype=="PE"|tvctype=="TI")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp)
    if(tvctype=="ED"|tvctype=="CO")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp,eta)
    print(round(dest,4))
    
    
    # sum of survival until event time for all the individuals in the data
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
      slogasc=sum(log(1- ((1+Hp1/fp[1])^(-fp[1]))*((1+Hp2/fp[2])^(-fp[2]))),na.rm=TRUE)
    }else{
      slogasc=sum(log(1- exp(-Hp1-Hp2)),na.rm=TRUE)  
    }
    
    #print(round(c(cest,fp),3))
    # cagep1.0 <- cagep[statusp==1&genderp==0]
    # cagep1.1 <- cagep[statusp==1&genderp==1]
    # st1p1<-st1p[statusp==1&genderp==1]
    # 
    # if(length(cagep1.0)!=0){
    # logasc1.0 = log(sapply(1:length(cagep1.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep1.0[index])))
    # }else{logasc1.0=0}
    # if(length(cagep1.1)!=0){
    # logasc1.1 = log(sapply(1:length(cagep1.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p1[index]))+
    #                 sapply(1:length(cagep1.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=1,tvc=st1p1[index],tt=tvctype,eta=eta,lower=st1p1[index],upper=cagep1.1[index])))
    # }else{logasc1.1=0}  
    # logasc1 <- sum(logasc1.0,na.rm=TRUE)+sum(logasc1.1,na.rm=TRUE)
    # 
    # # subsetting the probands according to his affection status at study entry
    # cagep2.0 <- cagep[statusp==2&genderp==0]
    # cagep2.1 <- cagep[statusp==2&genderp==1]
    # st1p2<-st1p[statusp==2&genderp==1]
    # 
    # if(length(cagep2.0)!=0){
    # logasc2.0 = log(sapply(1:length(cagep2.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep2.0[index])))
    # }else{logasc2.0=0}
    # if(length(cagep2.1)!=0){
    # logasc2.1 = log(sapply(1:length(cagep2.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p2[index]))+
    #             sapply(1:length(cagep2.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=2,tvc=st1p2[index],tt=tvctype,eta=eta,lower=st1p2[index],upper=cagep2.1[index])))
    # }else{logasc2.1=0}
    # logasc2 <- sum(logasc2.0,na.rm=TRUE)+sum(logasc2.1,na.rm=TRUE)
    
    #slogasc = logasc1 + logasc2
    likelihood  <- loglik - slogasc
    
    return(-likelihood)
}
loglik_Comp_simulation3_second=function(theta, data, agemin, frailty,tvctype,stage1){
  
  data = data[data$currentage>=agemin,]
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
    eta <- exp(theta[10])
  }else if(tvctype=="CO"){
    #eta <- c(exp(theta[10]),theta[11])
    eta <- c(exp(stage1[8]),stage1[9])
  }
  
  # frailty parameter
  if(frailty==TRUE){
    fp <- exp(theta[1:2])
  }
  
  
  # Y, delta, carrier prob, proband indicator, sc number
  time0 = data$time-agemin
  status = data$status
  ip = which(data$proband==1)
  #ex1 = data$carrp.geno 
  #gender = data$gender
  mgene=data$mgene
  st1=data$tvc-agemin
  gender=data$TVCstatus
  st1[gender==0]<-time0[gender==0]
  
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
    
    Hfam1 <- aggregate(Hfam1,by=list(data$famID),FUN=sum)[,2]
    Hfam2 <- aggregate(Hfam2,by=list(data$famID),FUN=sum)[,2]
    
    
    sum4.1 <- sum((lfactorial(fp[1]+df1-1)-(df1-1)*log(fp[1])-lfactorial(fp[1]) + (-fp[1]-df1)*log(1+(Hfam1)/fp[1])), na.rm=T)
    sum4.2 <- sum((lfactorial(fp[2]+df2-1)-(df2-1)*log(fp[2])-lfactorial(fp[2]) + (-fp[2]-df2)*log(1+(Hfam2)/fp[2])), na.rm=T)
    
    loglik = sum1+sum2 + sum4.1+ sum4.2
  }else{
    sum4 <- sum(-k,na.rm = TRUE)
    loglik = sum1+sum2+sum4 # numerator in loglikelihood
  }
  
  # Ascertainment correction by design="pop+"
  cest = c(lambda1,rho1,lambda2,rho2,beta.sc1_1,beta.gen1,beta.gen2)
  cagep <- data$currentage[ip]-agemin
  statusp <- data$status[ip]  # proband disease status at the study entry
  genderp = gender[ip]
  st1p=st1[ip]
  # subsetting the probands according to his affection status at study entry
  # subsetting the probands according to his affection status at study entry
  # h1*S
  
  if(tvctype=="PE"|tvctype=="TI")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp)
  if(tvctype=="ED"|tvctype=="CO")dest = c(lambda1*1000,rho1,lambda2*1000,rho2,beta.sc1_1,beta.gen1,beta.gen2,fp,eta)
  print(round(dest,4))
  
  
  # sum of survival until event time for all the individuals in the data
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
    slogasc=sum(log(1- ((1+Hp1/fp[1])^(-fp[1]))*((1+Hp2/fp[2])^(-fp[2]))),na.rm=TRUE)
  }else{
    slogasc=sum(log(1- exp(-Hp1-Hp2)),na.rm=TRUE)  
  }
  
  #print(round(c(cest,fp),3))
  # cagep1.0 <- cagep[statusp==1&genderp==0]
  # cagep1.1 <- cagep[statusp==1&genderp==1]
  # st1p1<-st1p[statusp==1&genderp==1]
  # 
  # if(length(cagep1.0)!=0){
  # logasc1.0 = log(sapply(1:length(cagep1.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep1.0[index])))
  # }else{logasc1.0=0}
  # if(length(cagep1.1)!=0){
  # logasc1.1 = log(sapply(1:length(cagep1.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=1,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p1[index]))+
  #                 sapply(1:length(cagep1.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=1,tvc=st1p1[index],tt=tvctype,eta=eta,lower=st1p1[index],upper=cagep1.1[index])))
  # }else{logasc1.1=0}  
  # logasc1 <- sum(logasc1.0,na.rm=TRUE)+sum(logasc1.1,na.rm=TRUE)
  # 
  # # subsetting the probands according to his affection status at study entry
  # cagep2.0 <- cagep[statusp==2&genderp==0]
  # cagep2.1 <- cagep[statusp==2&genderp==1]
  # st1p2<-st1p[statusp==2&genderp==1]
  # 
  # if(length(cagep2.0)!=0){
  # logasc2.0 = log(sapply(1:length(cagep2.0), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=cagep2.0[index])))
  # }else{logasc2.0=0}
  # if(length(cagep2.1)!=0){
  # logasc2.1 = log(sapply(1:length(cagep2.1), function(index) penf(integrand,S=0,est=cest, fp=fp,frailty=frailty, type=2,tvc=0,tt=tvctype,eta=eta,lower=0,upper=st1p2[index]))+
  #             sapply(1:length(cagep2.1), function(index) penf(integrand,S=1,est=cest, fp=fp,frailty=frailty, type=2,tvc=st1p2[index],tt=tvctype,eta=eta,lower=st1p2[index],upper=cagep2.1[index])))
  # }else{logasc2.1=0}
  # logasc2 <- sum(logasc2.0,na.rm=TRUE)+sum(logasc2.1,na.rm=TRUE)
  
  #slogasc = logasc1 + logasc2
  likelihood  <- loglik - slogasc
  
  return(-likelihood)
}
#------------------------------------------------integration----------------------
library(Rcpp)
sourceCpp("C:/Users/Jay/Desktop/CompRisk/R_CompRisk/simulation/Simulation_Cpp.cpp")


    penf <- function(f,S,est,fp,type,frailty,tvc,tt,eta,lower,upper){  
      integrate(f,S=S,est=est,fp=fp,type=type,frailty=frailty,tvc=tvc,tt=tt,eta=eta,lower=lower,upper=upper)$value
    }
    integrand <- function(u,S,est,fp,type,frailty,tvc,tt,eta){
      u = u
      tvc=tvc
      if(S==0){
      eta=eta  
      xbeta1 =  est[6]
      xbeta2 =  est[7]
      H1 = - (est[1]*u)^est[2]*exp(xbeta1)
      H2 = - (est[3]*u)^est[4]*exp(xbeta2)
      }else{
      if(tt=="ED"){
      xbeta1 =  g(u,tvc,eta[1])*est[5]+est[6]
      xbeta2 =  g(u,tvc,eta[2])*est[7]+est[8]  
      H1 = - exp(est[6])*sapply(u,cH_ED,a=tvc,parm=c(est[1],est[2],est[5],eta[1]))
      H2 = - exp(est[8])*sapply(u,cH_ED,a=tvc,parm=c(est[3],est[4],est[7],eta[2]))
      }else if(tt=="PE"){
      xbeta1 =  S*est[5]+est[6]
      xbeta2 =  S*est[7]+est[8]
      H1 = - exp(xbeta1)*((est[1]*u)^est[2]-(est[1]*tvc)^est[2])
      H2 = - exp(xbeta2)*((est[3]*u)^est[4]-(est[3]*tvc)^est[4])    
      }else if(tt=="CO"){
      xbeta1 =  g(u,tvc,eta[1])*est[5]+est[6]+eta[2]
      xbeta2 =  est[7] 
      H1 = - exp(est[6])*sapply(u,cH_CO,a=tvc,parm=c(est[1],est[2],est[5],eta[1],eta[2]))
      H2 = - (est[3]*u)^est[4]*exp(xbeta2)  
      }
      }
      
      if(type==1){
        haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta1)
      }else if(type==2){
        haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(xbeta2)  
      }
      
      if(frailty==TRUE){
        if(type==1){
          f = haz*((1-H1/fp[1])^(-fp[1]-1))*((1-H2/fp[2])^(-fp[2]))
          #f = haz*(( 1 -  (H1+H2)/fp  ))^(-fp)
        }else{
          f = haz*((1-H2/fp[2])^(-fp[2]-1))*((1-H1/fp[1])^(-fp[1]))
          #f = haz*(( 1 -  (H1+H2)/fp  ))^(-fp)
        }
      } else{
        Hp <- H1+H2
        f = haz*exp(Hp) # h * S
      }
      return(f)
}
#------------------------------------------------integration----------------------
    