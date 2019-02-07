########################################################## somme basic functions
# data.classifier <- function(data){
#   
#   df.comb <- data         
#   brca1 <- df.comb[df.comb$BRCA_FAMILY==1,]
#   brca2 <- df.comb[df.comb$BRCA_FAMILY==2,]
#   brcab <- df.comb[df.comb$BRCA_FAMILY=="Both",]
#   
#   #relation for 1F0
#   brca1[(brca1$relation=="1F0"),"generation"] <- 1
#   brca2[(brca2$relation=="1F0"),"generation"] <- 1
#   brcab[(brcab$relation=="1F0"),"generation"] <- 1
#   
#   b1 <- brca1
#   what <- c(b1$famID[b1$proband==1 & b1$mgene==0],b1$famID[b1$proband==1 & is.na(b1$mgene)])
#   what2 <-c(10000637, 10001232, 10001960 ,10001970 ,10003813 ,21000009 ,22000133 ,22000134, 22000242 ,22000540, 30000036, 30000103,30000237 ,30000277 ,43305355 ,43306006 ,43306512)
#   dfam <- c(what,what2)
#   b1 <<- subset(brca1, !(brca1$famID %in% dfam))
#   
#   #BRCA1 affected and unaffected
#   # b1  <<- brca1
#   #b1.af  <<- brca1[brca1$affect!=0,]
#   
#   #b1.uaf <<- brca1[brca1$affectedFam==0,]
#   
#   #BRCA2 affected and unaffected
#   colnames(brca2)[colnames(brca2)=="mgene"] <- "BRCA1_PERSON_STATUS"
#   colnames(brca2)[colnames(brca2)=="BRCA2_PERSON_STATUS"] <- "mgene"
#   brca2[brca2$mgene==99,"mgene"] <- NA
#   
#   
#   b2 <- brca2
#   what <- c(b2$famID[b2$proband==1 & b2$mgene==0],b2$famID[b2$proband==1 & is.na(b2$mgene)])
#   what2 <-c(10000637, 10001232, 10001960 ,10001970 ,10003813 ,21000009 ,22000133 ,22000134, 22000242 ,22000540, 30000036, 30000103,30000237 ,30000277 ,43305355 ,43306006 ,43306512)
#   dfam <- c(what,what2)
#   b2 <<- subset(brca2, !(brca2$famID %in% dfam))
#   #b2.af  <<- brca2[brca2$affectedFam==1,]
#   #b2.af <<- brca2[brca2$affect!=0,]
#   
#   #bb.af  <<- brcab[brcab$affectedFam==1,]
#   #bb.uaf <<- brcab[brcab$affectedFam==0,]
#   
#   # up to 3rd generation, excluding husband,spouse relationships
#   #  b1.af.3rd <- b1.af[!is.na(b1.af$generation),]
#   #  b1.af.3rd <- b1.af.3rd[(b1.af.3rd$relation != 6 & b1.af.3rd$relation != 7),]
#   #b1.af.3rd <<- b1.af.3rd[!(b1.af.3rd$relation=="1F0"),]
#   #b1.uaf.3rd <- b1.uaf[!is.na(b1.uaf$generation),]
#   
#   #  b2.af.3rd <- b2.af[!is.na(b2.af$generation),]
#   #  b2.af.3rd <- b2.af.3rd[(b2.af.3rd$relation != 6 & b2.af.3rd$relation != 7),]
#   #b2.af.3rd <<- b2.af.3rd[!(b2.af.3rd$relation=="1F0"),]
#   #b2.uaf.3rd <<- b2.uaf[!is.na(b2.uaf$generation),]
#   
#   #bb.af.3rd <- bb.af[!is.na(bb.af$generation),]
#   #bb.uaf.3rd <- bb.uaf[!is.na(bb.uaf$generation),]
# }
# 
# data.classifier(df.comb)

# Basic functions
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
        for(a in c(0,1)){
          carrp[is.na(data$mgene) & data$relation==g & data$gender==s & data$status==a ] <- mean(data$mgene[!is.na(data$mgene) & data$relation==g & data$gender==s & data$status==a ])
          
        }  
      }
    }
  }
  
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





# Basic functions
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
        for(a in c(0,1)){
          carrp[is.na(data$mgene) & data$relation==g & data$gender==s & data$status==a] <- mean(data$mgene[!is.na(data$mgene) & data$relation==g & data$gender==s & data$status==a])
          
        }  
      }
    }
  }
  
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




###### missing data imputed by expected values
