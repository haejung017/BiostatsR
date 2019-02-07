library(truncnorm)###Two frailty dist
simfam.cmp <- function(N.fam=779, iest, frailty.dist, mrate){
  noerror="error"
  
  # while(!is.null(noerror)){
  simdat=familyDesign(n=N.fam, affectnum=1,  m.carrier= 1, depend=c(iest[8],iest[9]), 
                      vbeta=iest[5:7], parms=iest[1:4], variation="none", 
                      base.dist="Weibull", frailty.dist=frailty.dist,
                      allelefreq=0.0021, mrate=mrate, agemin=16)
  #   noerror <- attr(simdat,"class")
  #   print(noerror)
  # }
  simdat<-data.frame(simdat)
  simdat<-simdat[simdat$time>16,]
  simdat$mut<-simdat$mut1<-simdat$mgene
  simdat$status1<-ifelse(simdat$status==1,1,0)
  simdat$time1 <- simdat$time
  # size<-aggregate(simdat$status,by=list(simdat$famID),length)[,2]
  # df<-aggregate(simdat$status!=0,by=list(simdat$famID),sum)[,2]
  # df1<-aggregate(simdat$status1,by=list(simdat$famID),sum)[,2]
  # simdat$df<-rep(df,size)
  # simdat$df1<-rep(df1,size)
  names(simdat)[1:5]<-c("famID","indID","gender","moid","faid")
  #names(simdat)[11]<-"cage"
  simdat
}

familyDesign <-
  function(n=1000, affectnum=0,  m.carrier= 0, dominant.m = TRUE, dominant.s = TRUE, depend=1, base.dist="Weibull", frailty.dist="gamma",
           vbeta= c(-1.126, 2.55, 1.6), parms=c(0.016, 3), variation="none", 
           allelefreq=c(0.02, 0.2), mrate=0.1, age1=c(65,2.5), 
           agemin=14) 
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
                           variation=variation, allelefreq=c(0.0021,0.2), mrate=mrate, age1=c(65,2.5),age2=c(45,2.5),
                           agemin=16)
      
      if(is.null(attr(dat, "class"))){
        # At least one parent in first gen and two sibs in the second gen should be affected
        if(affectnum==3) until<- ifelse(sum(dat[dat[,7]==1,13]) >= 1 & sum(dat[dat[,7]==2,13]) > 1, TRUE, FALSE) 
        #[,7]=generation, [,13]=status
        else until <- TRUE
        
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
            age1=c(65,2.5), age2=c(45,2.5), agemin=14)
  {
    tmpdata<-numeric()
    indID<-c(cumind+1,cumind+2)
    motherID<-c(0,0) 
    fatherID<-c(0,0) 
    gender<-c(1,0) ## 0--Female, 1--Male
    proband<-c(0,0) ## 0--Non-proband,  1-- proband
    generation<-c(1,1) ## generation number, 0 in second generation means married from outside
    relation <- c(4,4)  #1=proband, 2=sib, 3=child, 4=parent, 5=sibchild, 6=hasband, 7=sibspouse  
    
    ## number of sibs in second generation be 2 to 5 
    ## truncated negative binomial (>=2 and <=5) with prob=sibprob=0.4551, and rr=1
    
    SecNum <-sample(c(2,3,4,5), 1, replace=TRUE, prob=c(0.4991, 0.2720, 0.1482, 0.0807))
    #SecNum <-sample(c(1,2), 1, replace=TRUE, prob=c(0.4991, 0.2720+ 0.1482+0.0807))
    #SecNum=2
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
        
        ThiNum<-sample(c(2,3,4,5), 1, replace=TRUE,prob=c(0.4991, 0.2720, 0.1482, 0.0807)) 
        #ThiNum<-sample(c(1,2), 1, replace=TRUE,prob=c(0.4991, 0.2720+ 0.1482+ 0.0807)) 
        #ThiNum = 2
        
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
          
          gender<-c(gender,sample(c(1, 0), 1,replace=TRUE, prob=c(0.5, 0.5)))
          generation<-c(generation, 3)           
        }#Close for for(k in 1:ThiNum)
        NumMem<-NumMem+ThiNum  
      } #Close for  for (j in 1:SecNum)
      
      nIndi <- length(indID)     
      famID <- rep(i, nIndi)    
      ageonset <- rep(0, nIndi)      ## age at onset
      censorage <- rep(0, nIndi)     ## current age
      status <- rep(0, nIndi)        ## disease staus
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
        alpha1 <- log(rgamma(1, shape=depend[1], scale=1/depend[1]))
        alpha2  <-log(rgamma(1, shape=depend[2], scale=1/depend[2]))}
      else alpha <- c(0,0)
      
      #print(alpha)
      #------------------------------------------------------------------------------
      
      ###  Generating the current age
      # First generation with default mean 65 and sd=2.5; Second gen with default mean 45 and sd=2.5; third gene has mean age difference of 20 with their parents and sd=1
      #Probands age generated with truncated normal
      
      genepos<-pos[generation==1]
      censorage[genepos]<- rnorm(length(genepos), mean=age1[1], sd=age1[2])
      min.page1<- min(censorage[genepos])
      genepos<-pos[generation==0]
      censorage[genepos]<- rnorm(length(genepos), mean=age2[1], sd=age2[2])
      genepos<-pos[generation==2] 
      for (j in 1:length(genepos)) {
        tmp<-rtruncnorm(1,a=agemin,b=min.page1-14, mean=age2[1]-j, sd=age2[2])
        censorage[genepos[j]]<- tmp
        #min.page2 <- min( c(tmp, censorage[generation==0]) )
        sonpos<-pos[fatherID==indID[genepos[j]] | motherID==indID[genepos[j]]]
        min.page2 <- min(censorage[indID==fatherID[sonpos[1]] | indID==motherID[sonpos[1]]])
        for (k in 1:length(sonpos)) {
          #censorage[sonpos[k]]<- rtruncnorm(1,a=0, b=min.page2-14, mean=min.page2-20-k, sd=1)
          censorage[sonpos[k]]<- rnorm(1,mean=min.page2-20-k, sd=1.5)
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
      probID <- indID[proband==1] # proband's ID
      pGene <- fgene(base.dist, prob.age-agemin, prob.sex, variation=variation, parms=parms, vbeta=vbeta, alpha1=alpha1,alpha2=alpha2, pg=0, m.carrier=m.carrier, dominant.m=dominant.m, aq = allelefreq)
      #------------------------------------------------------------------------------ 
      
      if(variation=="secondgene") ngene <- c(1,2)
      else  ngene <- 1
      
      for(g in ngene){
        disgene <- rep(0, nIndi)
        if( m.carrier==1 & g==1) { #pop+
          if(dominant.m) gg <- 1:2
          else gg <- 1
          prob.G <- sample(gg, 1, replace=TRUE, prob=pGene[1,gg])
        }#Close for if( m.carrier==1 & g==1)
        else prob.G <- sample(c(1,2,3),1, replace=TRUE, prob=pGene[g,])
        
        #generate proband's parents' genotype based on proband's genotype  
        G1 <- parents.g(prob.G, g=g, allelefreq=allelefreq) 
        
        sib.G <- kids.g(1, G1) # generate kids' genotype given the proband's genotype
        
        if(g==1){
          ParentsG.m[generation==2] <- 3*(G1[1]-1) + G1[2] 
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

      # genotypes of third generation
      for(i in indID[generation==3]){
        m.g <- disgene.m[indID==motherID[indID==i]]
        f.g <- disgene.m[indID==fatherID[indID==i]]
        
        disgene.m[indID==i] <- kids.g(1, c(m.g, f.g) )
        ParentsG.m[indID==i] <- 3*(m.g-1)+f.g
        
      }#Close for for(i in indID[generation==3])
      
      ## We dichotomize 3 genotypes into two (carrier or not) depending on the model 
      ## dominant.m=1 if major gene is based on dominant model, 0 for recessive model )
      ## dominant.s=1 if second gene is based on dominant model, 0 for recessive model )
      
      if(dominant.m)  majorgene <- ifelse(disgene.m==3, 0,1)
      else majorgene <- ifelse(disgene.m==1,1,0)
      
      if(variation=="secondgene"){
        if(dominant.s) secondgene <- ifelse(disgene.s==3,0,1)
        else   secondgene <- ifelse(disgene.s==1, 1, 0)
      }
      else secondgene <- rep(-1, nIndi)
      
      #------------------------------------------------------------------------------
      ## simulate the age at on set for the family member: each generation may have different frailty, so simulated seperately.
      
      vbeta1 <- vbeta[1]   ## log relative risk of the gender effect for cause1
      vbeta2 <- vbeta[2]   ## log relative risk of the major gene for cause1
      #vbeta3 <- vbeta[3]   ## log relative risk of the gender effect for cause2
      vbeta4 <- vbeta[3]   ## log relative risk of the major gene for cause2
      
      # gen3 = generation 1,2,3   
      gen3 <- ifelse(generation==2 | generation==0, 2, ifelse(generation==1, 1,3))
      affected <- (proband==1)
      
        xv1beta <- gender*vbeta1 + majorgene*vbeta2
        xv2beta <-  majorgene*vbeta4
        
        ## generate ageonset for affected proband
        genepos <- pos[proband==1]
        affage <- censorage[genepos]
        affage.min <- ifelse(affage> agemin, affage-agemin, 0)
        uni <- runif(length(genepos), 0,1)
        ageonset[genepos] <- apply(cbind(xv1beta[genepos],xv2beta[genepos],affage.min, uni), 1, inv.survp, base.dist= base.dist, parms=parms, alpha1=alpha1,alpha2=alpha2)+agemin
        pen1 <- penf.cmp(affage.min, delta=1, base.dist=base.dist, parms, xv1beta[genepos]+alpha1, xv2beta[genepos]+alpha2)
        pen2 <- penf.cmp(affage.min, delta=2, base.dist=base.dist, parms, xv1beta[genepos]+alpha1, xv2beta[genepos]+alpha2)
        
        status[genepos] <- ifelse(runif(length(genepos)) < pen1/(pen1+pen2), 1, 2)
        
        
        ## generate ageonset for non-affecteds
        genepos <- pos[proband==0]
        uni <- runif(length(genepos), 0,1)
        ageonset[genepos] <- apply(cbind(xv1beta[genepos],xv2beta[genepos],uni), 1, inv.surv, base.dist=base.dist, parms=parms, alpha1=alpha1,alpha2=alpha2)+agemin
        haz1 <- hazards(dist=base.dist, ageonset[genepos]-agemin, parms[1:2])*exp(xv1beta[genepos]+alpha1)
        haz2 <- hazards(dist=base.dist, ageonset[genepos]-agemin, parms[3:4])*exp(xv2beta[genepos]+alpha2)
        status[genepos] <- ifelse(runif(length(genepos)) < haz1/(haz1+haz2), 1, 2)
        
        
        
      
      currentage<-ifelse(censorage>100, 100, censorage)
      time <- pmin(currentage, ageonset)   
      status <- ifelse(currentage >= ageonset, status, 0)
      
      
      
      
      # Generating missing genotypes 
      # mgene is genotype with missing genotype NA
      mgene <-  majorgene 
      mm <- round((length(mgene)-1)*mrate)  # no. of missing 
      mgene[is.element(indID, sample(indID[!proband], mm)) ] <- NA
      
      fsize <- length(famID) # family size
      naff <- sum(status!=0) # number of affected in the family
      tmpdata<-cbind(
        famID, indID, gender, motherID, fatherID, proband, generation,
        majorgene=disgene.m, secondgene=disgene.s, ageonset, currentage, time, status, mgene, relation, fsize, naff)
    }#Close for if (SecNum > 0)
    #print(tmpdata)
    return(tmpdata)
    
}

fgene <-function(base.dist, affage, affsex, variation="secondgene", parms, vbeta, alpha1,alpha2,
           pg=0, m.carrier=0, dominant.m=TRUE, aq){
    # returns 2x3 matrix , first row for the major gene, 
    #	second row for the second gene
    pAA <- pAa <- paa <- 0
    
    AAq <- Aaq <- aaq <- 0
    AAq[1] <- Pgene(1, pg=pg[1], a.freq=aq[1])
    Aaq[1] <- Pgene(2, pg=pg[1], a.freq=aq[1])
    aaq[1] <- Pgene(3, pg=pg[1], a.freq=aq[1])
    
    
      Ft <- 0
      for(i in c(0,1)){
        xbeta1 <- affsex*vbeta[1]+i*vbeta[2] + alpha1
        xbeta2 <- i*vbeta[3] + alpha2
        Ft[i+1] <- penf.cmp(affage, delta=1, base.dist=base.dist, parms=parms, xbeta1=xbeta1,xbeta2=xbeta2)
      } 
      pAA <- Ft[2]*AAq[1]
      if(dominant.m) pAa <- Ft[2]*Aaq[1]
      else pAa <- Ft[1]*Aaq[1]
      paa <- Ft[1]*aaq[1]
    
    return( cbind(pAA, pAa, paa)/c(pAA+pAa+paa))
  }

integrand.penf.cmp=function(u,delta, base.dist, parms, xbeta1, xbeta2)
{
  bhaz1 = hazards(base.dist, u, parms=parms[1:2])
  bhaz2 = hazards(base.dist, u, parms=parms[3:4])
  
  H1 = cumhaz(base.dist, u, parms=parms[1:2])*exp(xbeta1) #cause 1
  H2 = cumhaz(base.dist, u, parms=parms[3:4])*exp(xbeta2) #cause 2
  f= (delta==1)*bhaz1*exp(xbeta1)*exp(-H1-H2) +(delta==2)*bhaz2*exp(xbeta2)*exp(-H1-H2)
  return(f)
}

penf.cmp=function(t, delta, base.dist, parms, xbeta1,xbeta2)
{
  f=integrate(integrand.penf.cmp, lower=0, upper=t, delta=delta, base.dist=base.dist, parms=parms, xbeta1=xbeta1, xbeta2=xbeta2)$value
  return(f)
}

surv.dist <- function(t, base.dist, parms, xbeta1,xbeta2, alpha1,alpha2, res){
  if(base.dist=="Weibull") exp(-(parms[1]*t)^parms[2]*exp(xbeta1+alpha1)-(parms[3]*t)^parms[4]*exp(xbeta2+alpha2))-res
}

inv.surv <-
  function(base.dist, parms, alpha1,alpha2, val){
    uniroot(surv.dist, lower=0,upper=10000000, base.dist=base.dist, parms=parms, alpha1=alpha1,alpha2=alpha2, xbeta1=val[1], xbeta2=val[2],res=val[3])$root
  }

survp.dist <- function(t, base.dist, currentage, parms, xbeta1,xbeta2, alpha1,alpha2, res){
  A = surv.dist(t=t, base.dist=base.dist, parms=parms, xbeta1=xbeta1,xbeta2=xbeta2, alpha1=alpha1,alpha2=alpha2, res=0)
  B = surv.dist(t=currentage, base.dist=base.dist, parms=parms, xbeta1=xbeta1, xbeta2=xbeta2, alpha1=alpha1,alpha2=alpha2, res=0)
  return((A-B)/(1-B)-res)  
  
} 

inv.survp <- function(val, base.dist, parms, alpha1,alpha2){
  out<-try(uniroot(survp.dist, lower=0,upper=100000, base.dist=base.dist, parms=parms, 
                   alpha1=alpha1,alpha2=alpha2, xbeta1=val[1],xbeta2=val[2], currentage=val[3], res=val[4])$root)
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
