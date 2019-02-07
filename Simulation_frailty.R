# table(simdata$status[simdata$proband==1])
# table(simdata$mgene[simdata$proband==1])

simfam.cmp <- function(){
library(truncnorm)

  iest=c(0.008,2.300,0.007,2.932,1.872,1.858,1.224,3.5,3.240,0.278)

  vbeta=c(iest[5:7],iest[10]) #1 TD

    
  simdat=familyDesign(n=250, affectnum=1,  m.carrier= 1, depend=c(iest[8],iest[9]), 
                      vbeta=vbeta, parms=iest[1:4], variation="frailty", 
                      base.dist="Weibull", frailty.dist="gamma",
                      allelefreq=0.0021, mrate=0, agemin=16,tvctype="ED")
  
  simdat<-data.frame(simdat)
  simdat<-simdat[simdat$time>16,]
  simdat$mut<-simdat$mut1<-simdat$mgene
  names(simdat)[1:5]<-c("famID","indID","gender","moid","faid")
  
  return(simdat)
}



familyDesign <-
  function(n=500, affectnum=0,  m.carrier= 0, dominant.m = TRUE, dominant.s = TRUE, depend=1, base.dist="Weibull", frailty.dist="gamma",
           vbeta= c(-1.126, 2.55, 1.6), parms=c(0.016, 3), variation="none", 
           allelefreq=c(0.002, 0.2), mrate=0.1, age1=c(65,2.5), 
           agemin=14,tvctype="PE") 
  {
    data<-numeric()
    cumind<-0
    i<- 1
    j<- 0
    while (i <= n) {
      j <- j + 1
      dat<-familyStructure(i,cumind=cumind, m.carrier=1, 
                           base.dist="Weibull", frailty.dist=frailty.dist,
                           depend=depend, parms=parms, vbeta=vbeta, dominant.m=TRUE, dominant.s=TRUE,
                           variation=variation, allelefreq=c(0.0021,0.2), mrate=mrate, age1=c(65,1),age2=c(45,1),
                           agemin=16,tvctype=tvctype)
      
      if(is.null(attr(dat, "class"))){

        if(affectnum==3) until<- ifelse(sum(dat[dat[,7]==1,13]) >= 1 & sum(dat[dat[,7]==2,13]) > 1, TRUE, FALSE) 
        # we only keep the families whose proband event status is not 0, currentage < time
        else until <- ifelse( dat[dat[,6]==1,13] > 0, TRUE, FALSE) 
        
        if(!is.null(dim(dat))){
          if(nrow(dat)>0 ){
            if(until){
              data<-rbind(data, dat)
              cumind<-cumind+nrow(dat)
              i<-i+1
            }    
          }
        }
      } # close "is.null(attr(dat, "class"))"
    } # close while
    
    data
  }


familyStructure <-
  function (i, cumind, m.carrier=1, depend=1, parms=c(0.0021,2.9,0.005,3.73), 
            base.dist="Weibull", frailty.dist="gamma",
            vbeta=c(0.39,5.11,-0.69,4.12), dominant.m=TRUE, dominant.s= TRUE, 
            variation="frailty", allelefreq= c(0.02, 0.2), mrate=0, 
            age1=c(65,2.5), age2=c(45,2.5), agemin=14,tvctype="PE")
  {
    tmpdata<-numeric()
    indID<-c(cumind+1,cumind+2)
    motherID<-c(0,0) 
    fatherID<-c(0,0) 
    gender<-c(1,0) ## 0--Female, 1--Male
    proband<-c(0,0) ## 0--Non-proband,  1-- proband
    generation<-c(1,1) ## generation number, 0 in second generation means married from outside
    relation <- c(4,4)  #1=proband, 2=sib, 3=child, 4=parent, 5=sibchild, 6=hasband, 7=sibspouse  
    # SecNum <- 2, we can fix the number of sibs in second generation (eg., 2 sibs)
    
    ## number of sibs in second generation be 2 to 5 
    ## truncated negative binomial (>=2 and <=5) with prob=sibprob=0.4551, and rr=1
    #SecNum <-sample(c(2,3,4,5), 1, replace=TRUE, prob=c(0.4991, 0.2720, 0.1482, 0.0807)) 
    SecNum=2
    ## generate gender by prob=0.5      
    
    tmpgender<-sample(c(1,0), SecNum, replace=TRUE, prob=c(0.5,0.5))
    
    SecNum<-length(tmpgender)
    
    if (SecNum > 0) { ## at least one sample in second generation
      
      tmpproband<-sample(SecNum,1,replace=TRUE)
      
      NumMem<-2*SecNum+cumind+2
      
      for (j in 1:SecNum) {
        
        if (j==tmpproband){ 
          proband <- c(proband, c(1,0))
          relation <- c(relation, c(1,6))
        }
        else {
          proband<-c(proband, c(0,0))
          relation <- c(relation, c(2,7))
        }    
        indID<- c(indID, 2*j+cumind+1, 2*j+cumind+2 )
        fatherID<-c(fatherID, c(cumind+1,0))
        motherID<-c(motherID, c(cumind+2,0))
        
        if (tmpgender[j]==0) gender<-c(gender, c(1,0)) 
        else gender<-c(gender, c(0,1))
        
        generation<-c(generation, c(2,0)) 
        
        #ThiNum<-sample(c(2,3,4,5), 1, replace=TRUE,prob=c(0.4991, 0.2720, 0.1482, 0.0807)) 
        ThiNum=2
        
        for (k in 1:ThiNum) {
          proband<-c(proband,0)
          indID<- c(indID, NumMem+k)
          
          if(j==tmpproband) relation <- c(relation, 3)
          else relation <- c(relation, 5)
          
          if (gender[indID==(2*j+cumind+1)]==0) {
            fatherID<-c(fatherID, 2*j+cumind+2)
            motherID<-c(motherID, 2*j+cumind+1)
          } 
          else {
            fatherID<-c(fatherID, 2*j+cumind+1)
            motherID<-c(motherID, 2*j+cumind+2)
          }
          ## generate gender by prob=0.5
          gender<-c(gender,sample(c(1, 0), 1,replace=TRUE, prob=c(0.5, 0.5)))
          generation<-c(generation, 3)           
        }#Close for for(k in 1:ThiNum)
        NumMem<-NumMem+ThiNum  
      } #Close for  for (j in 1:SecNum)
      
      nIndi <- length(indID)     
      famID<-rep(i, nIndi)    
      ageonset<-rep(0, nIndi)      ## age at onset
      censorage<-rep(0, nIndi)     ## current age
      status<-rep(0, nIndi)        ## disease staus
      tvc<-rep(0, nIndi)           ## time varying covariate
      affected <- rep(0, nIndi)	
      disgene.m <- rep(0, nIndi)	
      disgene.s <- rep(0, nIndi)	
      ParentsG.m <- rep(0, nIndi)  ## parents genotype; it can be 1 to 9; 0 if founders
      pos <- c(1:nIndi) 
      
      
      
      #------------------------------------------------------------------------------
      # Familial correlation will be generated by gamma frailty, log-normal frailty, or secondgene
      
      ## Log-normal frailty (logX~Normal with mean=0, variance=depend)
      if (frailty.dist=="lognormal") alpha <- rnorm(1, mean=0, sd=sqrt(depend))
      ## Gamma Frailty (with mean 1, variance=depend)
      else if (frailty.dist=="gamma") {
        # alpha1 <- log(rgamma(1, shape=depend[1], scale=1/depend[1]))
        # alpha2  <-log(rgamma(1, shape=depend[2], scale=1/depend[2]))}
        alpha1 <- log(rgamma(1, shape=depend[1], scale=1/depend[1]))
        alpha2 <- log(rgamma(1, shape=depend[2], scale=1/depend[2]))}
      else alpha1 =alpha2= 0
      
      #print(alpha)
      #------------------------------------------------------------------------------
      
      #------------------------------------------------------------------------------
      
      ###  Generating time-varying covariate
      # m=20
      # sd=2.5
      # a=(sd/m)^(-1.086)
      # b=m/gamma(1+1/a)
      # tvc=rweibull(nIndi, shape=a, scale=b )
      tvc=runif(length(pos),min=5,max=100)
      #tvc<-rnorm(length(pos), mean=20, sd=2.5)
      #tvc[tvc<=5]=150
      #tvc<-rtruncnorm(length(pos),a=0, b=70, mean=24, sd=10)
      #tvc[sample(1:length(pos),round(length(pos)/1.5))]=150
      
      
      
      ###  Generating the current age
      # First generation with default mean 65 and sd=2.5; Second gen with default mean 45 and sd=2.5; third gene has mean age difference of 20 with their parents and sd=1
      #Probands age generated with truncated normal
      # tvc=vector()
      # tvc[generation==1]<-rnorm(sum(generation==1), mean=age1[1], sd=age1[2]+8)-16
      # tvc[generation==2]<-rnorm(sum(generation==2), mean=age2[1], sd=age2[2]+8)-16
      # tvc[generation==0]<-rnorm(sum(generation==0), mean=age2[1], sd=age2[2]+8)-16
      # tvc[generation==3]<-rnorm(sum(generation==3), mean=20,sd=3)-16
      
      genepos<-pos[generation==1]
      censorage[genepos]<- rnorm(length(genepos), mean=age1[1], sd=age1[2])
      min.page1<- min(censorage[genepos])
      if(min.page1 < agemin) stop("agemin is too large.")
      genepos<-pos[generation==0]
      censorage[genepos]<- rnorm(length(genepos), mean=age2[1], sd=age2[2])
      genepos<-pos[generation==2] 
      for (j in 1:length(genepos)) {
        tmp<-rtruncnorm(1,a=agemin,b=min.page1-14, mean=age2[1]-j, sd=age2[2])
        censorage[genepos[j]]<- tmp
        #min.page2 <- min( c(tmp, censorage[generation==0]) )
        #if(min.page2 < agemin) stop("agemin is too large.")
        sonpos<-pos[fatherID==indID[genepos[j]] | motherID==indID[genepos[j]]]
        min.page2 <- min(censorage[indID==fatherID[sonpos[1]] | indID==motherID[sonpos[1]]])
        #if(min.page2 < agemin) stop("agemin is too large.")
        for (k in 1:length(sonpos)) {
          #censorage[sonpos[k]]<- rtruncnorm(1,a=agemin, b=min.page2, mean=min.page2-20-k, sd=1)
          censorage[sonpos[k]]<- rnorm(1,mean=min.page2-20-k, sd=1)      
        }#Close for for(k in 1:length(sonpos))
      }#Close for for (j in 1:length(genepos))
      
      
      #------------------------------------------------------------------------------
      
      if(length(allelefreq)==1) allelefreq <- c(allelefreq,0) 
      qgene <- allelefreq
      AAq <- qgene^2
      Aaq <- 2*qgene*(1-qgene)
      aaq <- (1-qgene)^2
      G <- cbind(AAq, Aaq, aaq)
      
      #----------------------------------------------------------------------------
      
      
      ## For generating genotypes, we first generate the proband's genotype and based on that, we generate other family members' genotypes
      
      # Generate proband's genotype given its age at onset and gender
      prob.age <- censorage[proband==1] # proband's current age
      prob.sex <- gender[proband==1]
      prob.tvc <- tvc[proband==1]
      probID <- indID[proband==1] # proband's ID
      # generate proband's genotype (always carrier, pop+)
      pGene <- fgene(base.dist, prob.age-agemin, prob.tvc, variation=variation, parms=parms, vbeta=vbeta, alpha1=alpha1,alpha2=alpha2, pg=0, m.carrier=m.carrier, dominant.m=dominant.m, aq = allelefreq,tt=tvctype)
      #------------------------------------------------------------------------------ 
      
      if(variation=="secondgene") ngene <- c(1,2)
      else  ngene <- 1
      
      for(g in ngene){
        #disgene <- rep(0, nIndi)
        if( m.carrier==1 & g==1) {
          if(dominant.m) gg <- 1:2
          else gg <- 1
          prob.G <- sample(gg, 1, replace=TRUE, prob=pGene[1,gg])
        }#Close for if( m.carrier==1 & g==1)
        else prob.G <- sample(c(1,2,3),1, replace=TRUE, prob=pGene[g,])
        #generate proband's parents' genotype based on proband's genotype  
        G1 <- parents.g(prob.G, g=g, allelefreq=allelefreq) 
        
        nsibs <- sum(proband==0 & generation==2)
        sib.G <- kids.g(nsibs, G1) # generate kids' genotype given the proband's genotype
        
        if(g==1){
          # ParentsG.m[generation==2] <- 3*(G1[1]-1) + G1[2] 
          # ParentsG.m = parents' genotype; value can be 1 to 9.
          # let AA==1, Aa==2, aa==3, 
          # 9 possible combinations of parents are 1=(1,1), 2=(1,2), 3=(1,3), 4=(2,1), 5=(2,2), 6=(2,3), 7=(3,1), 8=(3,2), 9=(3,3)
          # disegene.m = geno of marjor gene as 1(AA), 2(Aa), or 3(aa)
          disgene.m[ generation==1 ] <- G1 
          disgene.m[ proband==1 ] <- prob.G
          disgene.m[ proband==0 & generation==2 ] <- sib.G	
        }	
        else{
          # disegene.s = genotype of second gene as 1(AA), 2(Aa), or 3(aa)
          disgene.s[generation==1 ] <- G1
          disgene.s[proband==1] <- prob.G
          disgene.s[proband==0 & generation==2 ] <- sib.G	
        }	 
      }#close for for(g in ngene)
      
      #------------------------------------------------------------------------------
      
      # genotypes of first and second gene for generation=0 which is founder by random selection among AA, Aa, aa 
      disgene.m[generation==0] <- sample(c(1,2,3), sum(generation==0), replace=TRUE, prob=G[1,] )
      
      if(variation=="secondgene"){
        disgene.s[generation==0] <- sample(c(1,2,3), sum(generation==0), replace=TRUE, prob=G[2,] )
      }
      
      # genotypes of third generation
      for(i in indID[generation==3]){
        m.g <- disgene.m[indID==motherID[indID==i]]
        f.g <- disgene.m[indID==fatherID[indID==i]]
        
        disgene.m[indID==i] <- kids.g(1, c(m.g, f.g) )
        #    ParentsG.m[indID==i] <- 3*(m.g-1)+f.g
        
        if(variation=="secondgene"){
          disgene.s[indID==i] <- kids.g(1, c(disgene.s[indID==motherID[indID==i]],
                                             disgene.s[indID==fatherID[indID==i]]))
        }
      }#Close for for(i in indID[generation==3])
      
      ## We dichotomize 3 genotypes into two (carrier or not) depending on the model 
      ## dominant.m=1 if major gene is based on dominant model, 0 for recessive model )
      ## dominant.s=1 if second gene is based on dominant model, 0 for recessive model )
      
      if(dominant.m)  majorgene <- ifelse(disgene.m==3, 0,1)
      else majorgene <- ifelse(disgene.m==1,1,0)
      
      if(variation=="secondgene"){
        if(dominant.s) secondgene <- ifelse(disgene.s==3,0,1)
        else   secondgene <- ifelse(disgene.s==1 , 1, 0)
      }
      else secondgene <- rep(0, nIndi)
      
      #------------------------------------------------------------------------------
      ## simulate the age at on set for the family member: each generation may have different frailty, so simulated seperately.
      
      vbeta1 <- vbeta[1]   ## log relative risk of the gender effect for cause1
      vbeta2 <- vbeta[2]   ## log relative risk of the major gene for cause1
      vbeta3 <- vbeta[3]   ## log relative risk of the gender effect for cause2
      #vbeta4 <- vbeta[4]   ## log relative risk of the major gene for cause2
      
      # gen3 = generation 1,2,3   
      gen3 <- ifelse(generation==2 | generation==0, 2, ifelse(generation==1, 1,3))
      affected <- (proband==1)
      
      xv1beta <-  majorgene*vbeta2+alpha1
      xv2beta <-  majorgene*vbeta3+alpha2
      bz = vbeta1
      if(tvctype=="PE"){
        eta=c(0)
      }else if(tvctype=="ED"){
        eta=c(vbeta[4])  
      }else if(tvctype=="CO"){
        eta=c(vbeta[4],vbeta[5])  
      }
      
      ## generate ageonset for all the members including probands
      genepos <- pos
      uni <- runif(length(genepos), 0,1)
      ageonset[genepos] <- apply(cbind(xv1beta[genepos],xv2beta[genepos], tvc[genepos], uni), 1, inv.surv, base.dist=base.dist, parms=parms, alpha1=alpha1,alpha2=alpha2,bz=bz,tt=tvctype,eta=eta)+agemin
      
      t=ageonset[genepos]-agemin
      tsc=tvc[genepos]
      
      if(tvctype=="PE"){xzb1=bz[1]}
      else if(tvctype=="ED"){xzb1=bz[1]*exp(-eta[1]*(t-tsc))}#;xzb2=bz[2]*exp(-eta[2]*(t-tsc))}
      else if(tvctype=="CO"){xzb1=bz[1]*exp(-eta[1]*(t-tsc))+eta[2]}#;xzb2=bz[2]*exp(-eta[2]*(t-tsc))+eta[4]}
      
      haz1 <- hazards(dist=base.dist, t, parms[1:2])*exp(xv1beta[genepos]+ifelse(t<tsc,0,xzb1))#*alpha1
      haz2 <- hazards(dist=base.dist, t, parms[3:4])*exp(xv2beta[genepos])#*alpha2

      status[genepos] <- ifelse(runif(length(genepos)) < haz1/(haz1+haz2), 1, 2)
      
      currentage<-ifelse(censorage>100, 100, censorage)
      time <- pmin(currentage, ageonset)   
      status <- ifelse(currentage >= ageonset, status, 0)
      tvc <- tvc+agemin
      TVCstatus <- ifelse(time >= tvc, 1, 0)
      genepos <- pos[proband==1]
      TVCstatusp=rep(0, nIndi)
      TVCstatusp[genepos]=ifelse(currentage[genepos]>=tvc[genepos],1,0)
      
      # Generating missing genotypes 
      # mgene is genotype with missing genotype NA
      mgene <-  majorgene 
      mm <- round((length(mgene)-1)*mrate)  # no. of missing 
      mgene[is.element(indID, sample(indID[!proband], mm)) ] <- NA
      
      fsize <- length(famID) # family size
      naff <- sum(status!=0) # number of affected in the family
      tmpdata<-cbind(
        famID, indID, gender, motherID, fatherID, proband, generation,
        majorgene=disgene.m, secondgene=disgene.s, ageonset, currentage, time, status, mgene, relation, fsize, naff,tvc,TVCstatus,TVCstatusp)
    }#Close for if (SecNum > 0)
    #print(tmpdata)
    return(tmpdata)
  }


fgene <-function(base.dist, affage, affsex, variation="secondgene", parms, vbeta, alpha1,alpha2,
                 pg=0, m.carrier=0, dominant.m=TRUE, aq,tt){
  # returns 2x3 matrix , first row for the major gene, 
  #	second row for the second gene
  pAA <- pAa <- paa <- 0
  
  AAq <- Aaq <- aaq <- 0
  AAq[1] <- Pgene(1, pg=pg[1], a.freq=aq[1])
  Aaq[1] <- Pgene(2, pg=pg[1], a.freq=aq[1])
  aaq[1] <- Pgene(3, pg=pg[1], a.freq=aq[1])
  
  av=c(alpha1,alpha2)
  tvc=affsex
  bz=c(vbeta[1])
  if(tt=="PE"){
    eta=c(0)  
  }else if(tt=="ED"){
    eta=c(vbeta[4])
  }else if(tt=="CO"){
    eta=c(vbeta[4],vbeta[5])
  }
  Ft <- 0
  
  pAA <- AAq[1]
  #pAA <- AAq[1]
  if(dominant.m) pAa <-Aaq[1]
  else pAa <- Aaq[1]
  paa <- aaq[1]
  
  #print(Ft)
  return( cbind(pAA, pAa, paa)/c(pAA+pAa+paa))
}

integrand.penf.cmp=function(u,delta, base.dist, parms, xbeta1, xbeta2)
{
  bhaz1 = hazards(base.dist, u, parms=parms[1:2])#*parms[5]
  bhaz2 = hazards(base.dist, u, parms=parms[3:4])#*parms[6]
  
  H1 = cumhaz(base.dist, u, parms=parms[1:2])*exp(xbeta1)#*parms[5] #cause 1
  H2 = cumhaz(base.dist, u, parms=parms[3:4])*exp(xbeta2)#*parms[6] #cause 2
  f= (delta==1)*bhaz1*exp(xbeta1)*exp(-H1-H2) +(delta==2)*bhaz2*exp(xbeta2)*exp(-H1-H2)
  return(f)
}
integrand.penf.cmp2=function(u,tvc,delta, base.dist, parms, xbeta1, xbeta2,bz,tt,eta)
{
  bhaz1 = hazards(base.dist, u, parms=parms[1:2])#*parms[5]
  bhaz2 = hazards(base.dist, u, parms=parms[3:4])#*parms[6]
  
  if(tt=="PE"){
    eta=eta*0;xzb1=bz[1]  
    H1 =exp(xbeta1)*(cumhaz(base.dist, tvc ,parms=parms[1:2])+sapply(u,cH_PE,a=tvc,parm=c(parms[1:2],bz[1])))#*parms[5] #cause 1
    H2 =exp(xbeta2)*(cumhaz(base.dist, tvc ,parms=parms[3:4])+cumhaz(base.dist, u ,parms=parms[3:4])) 
  }else if(tt=="ED"){
   xzb1 =bz[1]*exp(-eta[1]*(u-tvc))#; xzb2=bz[2]*exp(-eta[2]*(u-tvc))
     H1 =exp(xbeta1)*(cumhaz(base.dist, tvc ,parms=parms[1:2])+sapply(u,cH_ED,a=tvc,parm=c(parms[1:2],bz[1],eta[1])))#*parms[5] #cause 1
     H2 =exp(xbeta2)*(cumhaz(base.dist, tvc ,parms=parms[3:4])+cumhaz(base.dist, u ,parms=parms[3:4])) #cause 2
  }else if(tt=="CO"){
    xzb1=bz[1]*exp(-eta[1]*(u-tvc))+eta[2] 
    H1 = exp(xbeta1)*(cumhaz(base.dist, tvc ,parms=parms[1:2])+sapply(u,cH_CO,a=tvc,parm=c(parms[1:2],bz[1],eta[1],eta[2])))#*parms[5] #cause 1
    H2 = exp(xbeta2)*(cumhaz(base.dist, u ,parms=parms[3:4])  +cumhaz(base.dist, tvc ,parms=parms[3:4]))#*parms[6] #cause 2
  }  
  f= (delta==1)*bhaz1*exp(xbeta1+xzb1)*exp(-H1-H2) +(delta==2)*bhaz2*exp(xbeta2)*exp(-H1-H2)
  return(f)
}

penf.cmp=function(t, delta, base.dist, parms, xbeta1,xbeta2)
{
  f=integrate(integrand.penf.cmp, lower=0, upper=t, delta=delta, base.dist=base.dist, parms=parms, xbeta1=xbeta1, xbeta2=xbeta2)$value
  return(f)
}
penf.cmp2=function(t, tvc, delta, base.dist, parms, xbeta1,xbeta2, bz,tt,eta)
{
  f=integrate(integrand.penf.cmp2, lower=tvc, upper=t, tvc=tvc, delta=delta, base.dist=base.dist, parms=parms, xbeta1=xbeta1, xbeta2=xbeta2, bz=bz,tt=tt,eta=eta)$value
  return(f)
}

surv.dist <- function(t, base.dist, parms, xbeta1,xbeta2, alpha1,alpha2, tvc,bz,tt,eta, res){
  if(base.dist=="Weibull") 
    
  ind1<-which(t<tvc)
  ind2<-which(t>=tvc)
  #Int1=cH(A=rep(tvc,length(t)), B=t, parm=c(parms[1],parms[2],bz[1],eta[1]), tt)
  
  if(tt=="PE"){
    Int1=sapply(t[ind2],cH_PE,a=tvc,parm=c(parms[1:2],bz[1]))
  }else if(tt=="ED"){
    Int1=sapply(t[ind2],cH_ED,a=tvc,parm=c(parms[1:2],bz[1],eta[1]))#;Int2=sapply(t[ind2],cH_ED,a=tvc,parm=c(parms[3:4],bz[2],eta[2]))
  }else if(tt=="CO"){
    Int1=sapply(t[ind2],cH_CO,a=tvc,parm=c(parms[1:2],bz[1],eta[1],eta[2]))#;Int2=sapply(t[ind2],cH_CO,a=tvc,parm=c(parms[3:4],bz[2],eta[2],eta[4]))
  }
  chaz1<-rep(NA,length(t))
  chaz2<-rep(NA,length(t))
  
  if(length(ind1)>0){
    chaz1[ind1]= (parms[1]*t[ind1])^parms[2]*exp(xbeta1)#*alpha1 
    chaz2[ind1]= (parms[3]*t[ind1])^parms[4]*exp(xbeta2)#*alpha2
  }
  if(length(ind2)>0){  
    chaz1[ind2]= exp(xbeta1)*((parms[1]*tvc)^parms[2] + Int1)#*alpha1
    chaz2[ind2]= exp(xbeta2)*(parms[3]*t[ind2])^parms[4]#*alpha2
  }
  exp(-chaz1-chaz2)-res
}

surv.dist2 <- function(t, base.dist, parms, xbeta1,xbeta2, alpha1,alpha2, res){
  if(base.dist=="Weibull") 

  chaz1<-rep(NA,length(t))
  chaz2<-rep(NA,length(t))
  chaz1= (parms[1]*t)^parms[2]*exp(xbeta1)*alpha1 
  chaz2= (parms[3]*t)^parms[4]*exp(xbeta2)*alpha2
  exp(-chaz1-chaz2)-res
}

inv.surv <-function(base.dist, parms, alpha1,alpha2,bz,tt,eta, val){
  uniroot(surv.dist, lower=0,upper=10000000, base.dist=base.dist, parms=parms, alpha1=alpha1,alpha2=alpha2, bz=bz,tt=tt,eta=eta,xbeta1=val[1], xbeta2=val[2],tvc=val[3], res=val[4])$root
}

survp.dist <- function(t, base.dist, currentage, parms, xbeta1,xbeta2, alpha1,alpha2,tvc,bz,tt,eta,res){
  A = surv.dist(t=t, base.dist=base.dist, parms=parms, xbeta1=xbeta1,xbeta2=xbeta2, alpha1=alpha1,alpha2=alpha2, tvc=tvc,bz=bz,tt=tt,eta=eta,res=0)
  # haz1 <- hazards(dist=base.dist, t=t, parms[1:2])*exp(xbeta1+ifelse(t<tvc,0,bz[1]*exp(-eta[1]*(t-tvc))))#*alpha1
  # haz2 <- hazards(dist=base.dist, t=t, parms[3:4])*exp(xbeta2)#*alpha2
  # haz=haz1+haz2
  B = surv.dist(t=currentage, base.dist=base.dist, parms=parms, xbeta1=xbeta1, xbeta2=xbeta2, alpha1=alpha1,alpha2=alpha2, tvc=tvc,bz=bz,tt=tt,eta=eta,res=0)
  #B = surv.dist(t=tvc, base.dist=base.dist, parms=parms, xbeta1=xbeta1, xbeta2=xbeta2, alpha1=alpha1,alpha2=alpha2, tvc=100000,bz=bz,tt=tt,eta=eta,res=0)  
  return((A-B)/(1-B)-res)
  #return((haz*A)/(1-B)-res)
} 


inv.survp <- function(val, base.dist, parms, alpha1,alpha2,bz,tt,eta){
  out<-try(uniroot(survp.dist, lower=0,upper=10000000, base.dist=base.dist, parms=parms, 
                   alpha1=alpha1,alpha2=alpha2, bz=bz,tt=tt,eta=eta, xbeta1=val[1],xbeta2=val[2], currentage=val[3], tvc=val[4], res=val[5])$root)
  if(is.null(attr(out,"class"))) return(out)
  else print(c(parms, val))
}

hazards   <- function(dist="Weibull", t, parms){
  if(dist=="Weibull")	 	haz <- (parms[1]^parms[2])*parms[2]*t^(parms[2]-1) 
  haz
}
cumhaz    <- function(dist="Weibull", t, parms){
  if(dist=="Weibull")	 	chaz <- (parms[1]*t)^parms[2] 
  chaz
}

Pgene     <- function(g, pg, a.freq=0.0014){
  #P(g| pg:parents' genotype), a.freq: allele freq.
  qAA <- a.freq^2
  qAa <- 2*a.freq*(1-a.freq)
  qaa <- (1-a.freq)^2
  
  if(length(g)==1) g <- rep(g,length(pg))
  if(length(pg) == 1) pg <- rep(pg, length(g))
  re <- 0
  re[g==1] <- cbind(qAA, 1, 0.5, 0, 0.5, 0.25, 0,0,0,0)[pg[g==1]+1]
  re[g==2] <- cbind(qAa, 0, 0.5, 1, 0.5, 0.5, 0.5, 1, 0.5, 0)[pg[g==2]+1]
  re[g==3] <- cbind(qaa, 0,0,0,0, 0.25,0.5, 0, 0.5,1)[pg[g==3]+1]
  return(re)
}
parents.g <- function(gene, g, allelefreq){
  ## input gene as a vector of (1,2,3) : AA==1, Aa==2, aa==3
  if(length(allelefreq)==1) allelefreq <- c(allelefreq,0)
  qgene <- allelefreq
  AAq <- qgene^2
  Aaq <- 2*qgene*(1-qgene)
  aaq <- (1-qgene)^2
  
  tmp.g <- matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3, 3,1, 3,2, 3,3), ncol=2, byrow=T)
  
  tmp.prob<-c(AAq[g]^2, AAq[g]*Aaq[g], AAq[g]*aaq[g], 
              Aaq[g]*AAq[g], Aaq[g]^2, Aaq[g]*aaq[g], 
              aaq[g]*AAq[g], aaq[g]*Aaq[g], aaq[g]^2)
  
  tmp <- sample(1:9, 1, prob=Pgene(gene, 1:9)*tmp.prob)
  
  return(tmp.g[tmp,])
}
kids.g    <- function(n.kids, p.gene){ 
  ptmp <- sum((p.gene-c(1,0))*c(3,1))
  return(sample(1:3, n.kids, replace=TRUE, prob=Pgene(1:3, ptmp)))
}

######## integrale term in the cum hazard function 

##  if we integrate over one TDC

## under Permanent Exposure (PE) TVC
cH_PE<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  
  bH0 <- function(u){  # Cumulative hazard function
    (lamb*u)^rho
  }
  
  
  return( exp(bz)*(bH0(b) - bH0 (a) ) )
  
}

## under Exponential Decay (ED) TVC
cH_ED<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  phi<-parm[4]
  
  bh0<-function(u){
    (lamb*rho)*(lamb*u)^(rho-1)
  }
  
  integrand <- function(u) { bh0(u)*exp(bz*exp(-(u-a)*phi)) }
  int <- try( integrate( integrand, lower=a, upper=b), silent = TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <-NA_real_
  } else {
    integrated <-int$value
  }
  
  return(integrated)
  
}


## under Cox and Oaks (CO) TVC
cH_CO<-function(a, b, parm){
  
  lamb<-parm[1]
  rho<-parm[2]
  bz<-parm[3]
  phi<-parm[4]
  a0<-parm[5]
  
  bh0<-function(u){
    (lamb*rho)*(lamb*u)^(rho-1)
  }
  
  integrand <- function(u) { bh0( u)*exp(a0+bz*exp(-(u-a)*phi)) }
  int <- try( integrate( integrand, lower=a, upper=b), silent = TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <-NA_real_
  } else {
    integrated <-int$value
  }
  
  return(integrated)
  
}

# combine all these function in one for subject j
cHj<-function(x, parm,TDfunc){
  a=x[1]; b=x[2]
  if(TDfunc=="PE") res<-cH_PE(a, b, parm)
  else if(TDfunc=="ED") res<-cH_ED(a, b, parm)
  else if(TDfunc=="CO") res<-cH_CO(a, b, parm)
  else stop("TDfunc should be PE or ED or CO")
  res
}

#for several subjects :  use lapply function
cH<-function(A, B, parm, TDfunc){
  unlist( lapply(1:length(A), function(i) cHj(x=c(A[i],B[i]), parm=parm, TDfunc=TDfunc) ) )
}

