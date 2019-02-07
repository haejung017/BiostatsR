loglik_Comp_Timedep=function(theta, data, data2, agemin, frailty=FALSE, weight=FALSE, reccur.type, time.dep.cov.type){
  
  print(round(c(exp(theta[1:6]),theta[7:length(theta)]),3))
  
  data = data[data$currentage>=agemin,]
  # base parameters
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
  beta.sc1g_1 = theta[12]
  beta.sc2g_1 = theta[13]
  beta.sc3g_1 = theta[14]
  # beta.sc1o_1 = theta[15]
  # beta.sc2o_1 = theta[16]
  # beta.sc3o_1 = theta[17]
  
  # vbeta for ovarian
  beta.gen2 = theta[15]
  beta.sc1_2 =  theta[16]
  beta.sc2_2 = 0
  beta.sc3_2 = 0
  beta.sc1g_o = 0
  beta.sc2g_o = 0
  beta.sc3g_o = 0
  
  # vbeta for death
  beta.gen3 = theta[17]
  beta.sc1_3 = theta[18] 
  beta.sc2_3 = 0
  beta.sc3_3 = 0
  beta.sc1g_d = 0
  beta.sc2g_d = 0
  beta.sc3g_d = 0
  beta.s1o_d = 0
  beta.s2o_d = 0
  beta.s3o_d = 0
  
  
  # beta.sc1g_1 = 0
  # beta.sc2g_1 = 0
  # beta.sc3g_1 = 0
  beta.sc1o_1 = 0
  beta.sc2o_1 = 0
  beta.sc3o_1 = 0
  
  bref <- c(beta.sc1g_1,beta.sc2g_1,beta.sc3g_1,
            beta.sc1o_1,beta.sc2o_1,beta.sc3o_1)
  
  
  # time varying covariate type
  if(time.dep.cov.type=="ED"){
    phiv <- exp(theta[c(19,20,21)])
    phiv = c(phiv,0) 
    phi_sc = phiv[c(1:3)]
    phi_or = 0
  }
  
  # frailty parameter
  if(frailty==TRUE){
    fp <- exp(theta[length(theta)])
  }
  
  
  # sampling weight 
  if(weight==FALSE){
    wt = rep(1,nrow(data))
  }else{
    wt = 1/data$weight
  } 
  
  # Y, delta, carrier prob, proband indicator, sc number
  time0 = data$time-agemin
  status = data$status
  ip = which(data$proband==1)
  wtp = wt[ip]
  ex1 = data$carrp.geno 
  sc_num = data$screen_number
  sc_num_p = data$screen_number_p
  
  #setting up indicator variable for the removal status
  br_time <- data$br.censortime-agemin
  delta_BR <- ifelse(time0-br_time > 0,1,0)
  
  sc1_time <- data$st1 - agemin
  sc2_time <- data$st2 - agemin
  sc3_time <- data$st3 - agemin
  Z_SC1=  ifelse(is.na(sc1_time), 0,1)
  Z_SC2=  ifelse(is.na(sc2_time), 0,1)
  Z_SC3=  ifelse(is.na(sc3_time), 0,1)
  
  # Variables related to oophorectomy
  ov_time <- data$ov.censortime-agemin
  delta_OR <- ifelse(time0-ov_time > 0,1,0)
  Z_OR <- rep(1,nrow(data))
  
  if(time.dep.cov.type=="ED"){
    
    Z_OR = ifelse( delta_OR==0, 0, ZED(time0, ov_time, phi_or))
    Z_SC1 = ifelse( Z_SC1==0, 0, ZED(time0, sc1_time, phi_sc[1]))
    Z_SC2 = ifelse( Z_SC2==0, 0, ZED(time0, sc2_time, phi_sc[2]))
    Z_SC3 = ifelse( Z_SC3==0, 0, ZED(time0, sc3_time, phi_sc[3]))
  }
  
  #loading up event time(sorted) vectors, event indicator vectors
  time_vec2 <- data2[[1]]
  indicator_vec <- data2[[2]]
  time_vec2_p0 <- data2[[3]]
  indicator_vec_p0 <- data2[[4]]
  time_vec2_p1 <- data2[[5]]
  indicator_vec_p1 <- data2[[6]]
  time_vec2_p2 <- data2[[7]]
  indicator_vec_p2 <- data2[[8]]
  
  # baseline hazard(weibull), baseline cumulative hazard(weibull)
  bhaz1 = (lambda1^rho1)*rho1*time0^(rho1-1)
  bcumhaz1 = (lambda1*time0)^rho1 
  bhaz2 = (lambda2^rho2)*rho2*time0^(rho2-1)
  bcumhaz2 = (lambda2*time0)^rho2 
  bhaz3 = (lambda3^rho3)*rho3*time0^(rho3-1)
  bcumhaz3 = (lambda3*time0)^rho3 
  
  logh1.0 = log(bhaz1) 
  logh1.1 = log(bhaz1) + beta.gen1
  logh2.0 = log(bhaz2) 
  logh2.1 = log(bhaz2) + beta.gen2
  logh3.0 = log(bhaz3) 
  logh3.1 = log(bhaz3) + beta.gen3
  
  # creating screen,gene,oophorectomy interaction vectors for hazard
  sc_vec <-rep(0,nrow(data))
  sc_vec[sc_num==0] <- 0
  sc_vec[sc_num!=0] <- sc_num[sc_num!=0]
  sc_vec_bc  <- vector("numeric", length(sc_num) )
  sc_vec_bc[sc_vec==0] <- 0
  sc_vec_bc[sc_vec==1] <- Z_SC1[sc_vec==1]*beta.sc1_1
  sc_vec_bc[sc_vec==2] <- ifelse(reccur.type=="CM", beta.sc1_1 , 0) + Z_SC2[sc_vec==2]*beta.sc2_1
  sc_vec_bc[sc_vec==3] <- ifelse(reccur.type=="CM", beta.sc1_1 + beta.sc2_1 , 0) + Z_SC3[sc_vec==3]*beta.sc3_1
  sc_vec_ov <- vector("numeric", length(sc_num) )
  sc_vec_ov[sc_vec==0] <- 0
  sc_vec_ov[sc_vec==1] <- beta.sc1_2
  sc_vec_ov[sc_vec==2] <- ifelse(reccur.type=="CM", beta.sc1_2 , 0) + beta.sc1_2
  sc_vec_ov[sc_vec==3] <- ifelse(reccur.type=="CM", beta.sc1_2 + beta.sc2_2 , 0)+ beta.sc1_2
  sc_vec_d <- vector("numeric", length(sc_num) )
  sc_vec_d[sc_vec==0] <- 0
  sc_vec_d[sc_vec==1] <- beta.sc1_3
  sc_vec_d[sc_vec==2] <- ifelse(reccur.type=="CM", beta.sc1_3 , 0) + beta.sc1_3
  sc_vec_d[sc_vec==3] <- ifelse(reccur.type=="CM", beta.sc1_3 + beta.sc2_3 , 0)+ beta.sc1_3
  
  sg_vec_bc <- vector("numeric", length(sc_num) )
  sg_vec_bc[sc_vec==0] <- 0
  sg_vec_bc[sc_vec==1] <- beta.sc1g_1*Z_SC1[sc_vec==1]
  sg_vec_bc[sc_vec==2] <- ifelse(reccur.type=="CM", beta.sc1g_1  , 0) + Z_SC2[sc_vec==2]*beta.sc2g_1
  sg_vec_bc[sc_vec==3] <- ifelse(reccur.type=="CM", beta.sc1g_1 + beta.sc2g_1 , 0) + Z_SC3[sc_vec==3]*beta.sc3g_1
  sg_vec_ov <- vector("numeric", length(sc_num) )
  sg_vec_ov[sc_vec==0] <- 0
  sg_vec_ov[sc_vec==1] <- beta.sc1g_o
  sg_vec_ov[sc_vec==2] <- ifelse(reccur.type=="CM", beta.sc1g_o  , 0) + beta.sc2g_o
  sg_vec_ov[sc_vec==3] <- ifelse(reccur.type=="CM", beta.sc1g_o + beta.sc2g_o , 0) + beta.sc3g_o
  sg_vec_d <- vector("numeric", length(sc_num) )
  sg_vec_d[sc_vec==0] <- 0
  sg_vec_d[sc_vec==1] <- beta.sc1g_d
  sg_vec_d[sc_vec==2] <- ifelse(reccur.type=="CM", beta.sc1g_d  , 0) + beta.sc2g_d
  sg_vec_d[sc_vec==3] <- ifelse(reccur.type=="CM", beta.sc1g_d + beta.sc2g_d , 0) + beta.sc3g_d
  
  os_vec_bc <- vector("numeric", length(sc_num) )
  os_vec_bc[sc_num==1 & delta_OR==1] <- beta.sc1o_1
  os_vec_bc[sc_num==2 & delta_OR==1] <- ifelse(reccur.type=="CM", beta.sc1o_1  , 0) + beta.sc2o_1 # check this 
  os_vec_bc[sc_num==3 & delta_OR==1] <- ifelse(reccur.type=="CM", beta.sc1o_1 + beta.sc2o_1, 0) + beta.sc3o_1 # check this
  os_vec_d <- vector("numeric", length(sc_num) )
  os_vec_d[sc_num==1 & delta_OR==1] <- beta.s1o_d
  os_vec_d[sc_num==2 & delta_OR==1] <- ifelse(reccur.type=="CM", beta.s1o_d  , 0) + beta.s2o_d # check this 
  os_vec_d[sc_num==3 & delta_OR==1] <- ifelse(reccur.type=="CM", beta.s1o_d + beta.s2o_d , 0) + beta.s3o_d # check this
  
  # sum of log-hazard
  sum1 = sum((      wt*ex1*(logh1.1 + delta_OR*Z_OR*beta.or1 +  sc_vec_bc + sg_vec_bc + os_vec_bc ) +
                      wt*(1-ex1)*(logh1.0 + delta_OR*Z_OR*beta.or1 +  sc_vec_bc + os_vec_bc))[status==1], na.rm=TRUE) 
  
  sum2 = sum((      wt*ex1*(logh2.1 + sc_vec_ov + sg_vec_ov ) +
                      wt*(1-ex1)*(logh2.0 + sc_vec_ov ))[status==3], na.rm=TRUE)
  
  sum3 = sum((      wt*ex1*(logh3.1 + sc_vec_d + os_vec_d + sg_vec_d ) +
                      wt*(1-ex1)*(logh3.0 + sc_vec_d + os_vec_d ))[status==4], na.rm=TRUE)
  
  # sum of survival until event time for all the individuals in the data
  cest <- c(lambda1,rho1,
            lambda2,rho2,
            lambda3,rho3,
            beta.gen1,#7
            beta.sc1_1,
            beta.sc2_1,
            beta.sc3_1,
            beta.or1,
            beta.gen2,#12
            beta.sc1_2,
            beta.gen3,#14
            beta.sc1_3)
  
  mat_all <- cbind(time_vec2, indicator_vec, time0, ex1)
  if(reccur.type=="CM"){
    sum4 = apply(mat_all, 1, func_all, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
  }else{
    if(time.dep.cov.type=="PE"){
      #sum4 = apply(mat_all, 1, func_all_IS, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3) 
      testm <- mat_all
      testm[is.na(testm)] <- 99
      CH_bc=CumH_c_IS(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=1)
      CH_ov=CumH_c_IS(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=2)
      CH_d=CumH_c_IS(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=3)
      CH0_bc <- (1-testm[,12])*CH_bc[,1]
      CH1_bc <-    testm[,12]*exp(beta.gen1)*CH_bc[,2]
      CH0_ov <- (1-testm[,12])*CH_ov[,1]
      CH1_ov <-    testm[,12]*exp(beta.gen2)*CH_ov[,2]
      CH0_d <- (1-testm[,12])*CH_d[,1]
      CH1_d <-    testm[,12]*exp(beta.gen3)*CH_d[,2]
      k <- CH0_bc + CH1_bc + CH0_ov + CH1_ov + CH0_d + CH1_d
      sum4=k
    } else {
      # sum4 = apply(mat_all, 1, func_all_IS_ED, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3, phiv=phiv)     
      testm <- mat_all
      testm[is.na(testm)] <-   99
      CH_bc=CumH_c_ED(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=1,phi=phiv)
      CH_ov=CumH_c_ED(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=2,phi=phiv)
      CH_d=CumH_c_ED(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=3,phi=phiv)
      CH0_bc <- (1-testm[,12])*CH_bc[,1]
      CH1_bc <-    testm[,12]*exp(beta.gen1)*CH_bc[,2]
      CH0_ov <- (1-testm[,12])*CH_ov[,1]
      CH1_ov <-    testm[,12]*exp(beta.gen2)*CH_ov[,2]
      CH0_d <- (1-testm[,12])*CH_d[,1]
      CH1_d <-    testm[,12]*exp(beta.gen3)*CH_d[,2]
      k <- CH0_bc + CH1_bc + CH0_ov + CH1_ov + CH0_d + CH1_d
      sum4=k
    }
  }
  
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam <- -sum4
    df <- data$df[data$proband==1]
    Hfam <- aggregate(Hfam,by=list(data$famID),FUN=sum)[,2]
    sum4 <- sum(wtp*(lfactorial(fp+df-1)-(df-1)*log(fp)-lfactorial(fp) + (-fp-df)*log(1+(Hfam)/fp)), na.rm=T)
    loglik = sum1+sum2+sum3+sum4
  }else{
    sum4 <- sum(wt*sum4,na.rm = TRUE)
    loglik = sum1+sum2+sum3+sum4 # numerator in loglikelihood
  }
  
  # Ascertainment correction by design="pop+"
  cagep <- data$currentage[ip]-agemin
  timep <- data$time[ip]-agemin
  statusp <- data$affect[ip]  # proband disease status at the study entry
  
  # subsetting the probands according to his affection status at study entry
  cagep0 <- cagep[statusp==0]
  timep0 <- timep[statusp==0]
  wtp0 <- wtp[statusp==0]
  
  mat_p0 <- cbind(time_vec2_p0, indicator_vec_p0, cagep0)
  if(reccur.type=="CM"){
    logasc0 <- apply(mat_p0, 1, func_p0, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
  }else{
    if(time.dep.cov.type=="PE"){
      #logasc0 <- apply(mat_p0, 1, func_p0_IS, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
      testm <- mat_p0
      testm[is.na(testm)] <- 99
      CH_bc=CumH_c_IS(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=1)
      CH_ov=CumH_c_IS(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=2)
      CH_d=CumH_c_IS(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=3)
      CH1_bc <- exp(beta.gen1)*CH_bc[,2]
      CH1_ov <- exp(beta.gen2)*CH_ov[,2]
      CH1_d <- exp(beta.gen3)*CH_d[,2]
      logasc0 <- CH1_bc + CH1_ov +CH1_d     
    } else{
      #logasc0 <- apply(mat_p0, 1, func_p0_IS_ED, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
      testm <- mat_p0
      testm[is.na(testm)] <- 99
      CH_bc=CumH_c_ED(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=1, phi = phiv)
      CH_ov=CumH_c_ED(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=2, phi = phiv)
      CH_d=CumH_c_ED(v1= testm[,1:5], v2=testm[,6:10], u = testm[,11],cest=cest, bref=bref, affp=FALSE, type=3, phi = phiv)
      CH1_bc <- exp(beta.gen1)*CH_bc[,2]
      CH1_ov <- exp(beta.gen2)*CH_ov[,2]
      CH1_d <- exp(beta.gen3)*CH_d[,2]
      logasc0 <- CH1_bc + CH1_ov +CH1_d    
    }
  }
  if(frailty==TRUE){
    logS <- log((1-logasc0/fp)^(-fp))
  }else {
    logS <- logasc0
  }
  logasc0 <- sum(wtp0*logS,na.rm=TRUE)
  
  # subsetting the probands according to his affection status at study entry
  cagep1 <- cagep[statusp==1]
  timep1 <- timep[statusp==1]
  wtp1 <- wtp[statusp==1]
  
  mat_p1 <- cbind(time_vec2_p1, indicator_vec_p1, cagep1)
  if(reccur.type=="CM"){
    logasc1 <- apply(mat_p1, 1, proband_penf, type=1, affp=TRUE, cest=cest, bref=bref)
  }else{
    if(time.dep.cov.type=="PE"){
      logasc1 <- apply(mat_p1, 1, proband_penf_IS, type=1, affp=TRUE, cest=cest, bref=bref, frailty=frailty, fp=fp)
    } else{
      logasc1 <- apply(mat_p1, 1, proband_penf_ED, type=1, affp=TRUE, cest=cest, bref=bref, phi=phiv, frailty=frailty, fp=fp)
    }
  }
  logasc1 <- sum(wtp1*logasc1,na.rm=TRUE)
  
  # subsetting the probands according to his affection status at study entry
  cagep2 <- cagep[statusp==2]
  timep2 <- timep[statusp==2]
  wtp2 <- wtp[statusp==2]
  
  mat_p2 <- cbind(time_vec2_p2, indicator_vec_p2, cagep2)
  if(reccur.type=="CM"){
    logasc2 <- apply(mat_p2, 1, proband_penf, type=2, affp=TRUE, cest=cest, bref=bref)
  }else{
    if(time.dep.cov.type=="PE"){
      logasc2 <- apply(mat_p2, 1, proband_penf_IS, type=2, affp=TRUE, cest=cest, bref=bref, frailty=frailty, fp=fp)
    } else{
      logasc2 <- apply(mat_p2, 1, proband_penf_ED, type=2, affp=TRUE, cest=cest, bref=bref, phi=phiv, frailty=frailty, fp=fp)
    }
  }
  logasc2 <- sum(wtp2*logasc2,na.rm=TRUE)
  
  slogasc = logasc0 + logasc1 + logasc2
  likelihood  <- loglik - slogasc
  
  #print(c(loglik, slogasc, likelihood))
  return(-likelihood)
  
} 




func_p0 <- function(v, cest, bref, beta.gen1, beta.gen2, beta.gen3){
  v1 = v[1:5]
  v2 = v[6:10]
  cagep0 = v[11]
  v1=v1[!is.na(v1)]
  v2=v2[!is.na(v2)]
  CH_bc <- CumH_proband_c(v1=v1,v2=v2,u=cagep0,cest=cest,bref=bref,affp=FALSE,type=1)
  CH_ov <- CumH_proband_c(v1=v1,v2=v2,u=cagep0,cest=cest,bref=bref,affp=FALSE,type=2)
  CH_d  <- CumH_proband_c(v1=v1,v2=v2,u=cagep0,cest=cest,bref=bref,affp=FALSE,type=3)
  exp(beta.gen1)*CH_bc+exp(beta.gen2)*CH_ov+exp(beta.gen3)*CH_d
}
func_all<- function(v, cest, bref, beta.gen1, beta.gen2, beta.gen3){
  v1 = v[1:5]
  v2 = v[6:10]
  time0  = v[11]
  e1 <- v[12]
  v1 = v1[!is.na(v1)]
  v2 = v2[!is.na(v2)]
  CH_bc <- CumH_c(v1=v1,v2=v2,u=time0,cest=cest,bref=bref,affp=FALSE,type=1)
  CH_ov <- CumH_c(v1=v1,v2=v2,u=time0,cest=cest,bref=bref,affp=FALSE,type=2)
  CH_d  <- CumH_c(v1=v1,v2=v2,u=time0,cest=cest,bref=bref,affp=FALSE,type=3)
  
  CH0_bc <- (1-e1)*CH_bc[1]
  CH1_bc <-    e1*exp(beta.gen1)*CH_bc[2]
  
  CH0_ov <- (1-e1)*CH_ov[1]
  CH1_ov <-    e1*exp(beta.gen2)*CH_ov[2]
  
  CH0_d <- (1-e1)*CH_d[1]
  CH1_d <-    e1*exp(beta.gen3)*CH_d[2]
  
  k <- CH0_bc + CH1_bc + CH0_ov + CH1_ov + CH0_d + CH1_d
  return(k)
}
proband_penf <- function(v, cest, bref, affp, type){
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
penf.comp <- function(f,v1,v2,est,bref, affp,type,beta_vec1,beta_vec2,beta_vec3,intg,intw,phi,lb,lower,upper){  
  integrate(f,v1,v2,est=est,bref=bref, affp=affp,beta_vec1=beta_vec1,
            beta_vec2=beta_vec2,beta_vec3=beta_vec3,type=type,intg=intg,intw=intw,
            phi=phi,lb=lb,lower=lower,upper=upper)$value
}  
integrand <- function(u,v1,v2,est, bref,affp,type,beta_vec1,beta_vec2,beta_vec3,intg,intw){
  u <- u
  #  print(u)
  # v1 = v1[!is.na(v1)]
  # v2 = v2[!is.na(v2)]
  sc1 = ifelse(length(v1[which(v2==1)])==1,v1[which(v2==1)],99)
  sc2 = ifelse(length(v1[which(v2==2)])==1,v1[which(v2==2)],99)
  sc3 = ifelse(length(v1[which(v2==3)])==1,v1[which(v2==3)],99)
  ov = ifelse(length(v1[which(v2==4)])==1,v1[which(v2==4)],99)
  br = ifelse(length(v1[which(v2==5)])==1,v1[which(v2==5)],99)
  
  zsc1 = ifelse(u<sc1,0,1)
  zsc2 = ifelse(u<sc2,0,1)
  zsc3 = ifelse(u<sc3,0,1)
  zor = ifelse(u<ov,0,1)
  zbr = ifelse(u<br,0,1)
  
  if(type==1){
    xbeta = est[7] + (est[8] + bref[1])*zsc1+(est[9] + bref[2])*zsc2+(bref[3] + est[10])*zsc3 + 
      est[11]*zor + bref[4]*zor*zsc1 + bref[5]*zor*zsc2  + bref[6]*zor*zsc3
    haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta)
    #print(c(xbeta, est[7]+sum(beta_vec1)))
  }else if(type==2){
    haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(est[12]+est[13]*zsc1)  
  }
  
  H1 = - (1-zbr)*(est[1]*u)^est[2]*exp(est[7] + sum(beta_vec1) + sum(intg) + sum(intw)) #+ sum(beta_vec1_int_g) + sum(beta_vec1_int_cov) ) #cause 1
  H2 = - (1-zor)*(est[3]*u)^est[4]*exp(est[12]+ sum(beta_vec2))
  H3 = - (est[5]*u)^est[6]*exp(est[14]+ sum(beta_vec3) ) 
  
  
  # H1 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=1)*exp(est[7])
  # H2 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=2)*exp(est[12])
  # H3 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=TRUE,type=3)*exp(est[14])
  
  Hp <- H1+H2+H3
  f = haz*exp(Hp) # h * S 
  
  return(f)
}




# Interval specific hazard
proband_penf_IS <- function(v, cest, bref, affp, type, frailty, fp){
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
    
    beta1 <- ifelse(q==1, 0, coeff_c_IS(v2[q-1],type=1, cest=cest, affp=affp) )
    beta2 <- ifelse(q==1, 0, coeff_c_IS(v2[q-1],type=2, cest=cest, affp=affp) )
    beta3 <- ifelse(q==1, 0, coeff_c_IS(v2[q-1],type=3, cest=cest, affp=affp) )
    beta_vec1 <- c(beta_vec1,beta1)
    beta_vec2 <- c(beta_vec2,beta2)
    beta_vec3 <- c(beta_vec3,beta3)
    
    if(type==1){
      inter = ifelse(c(q==1,q==1), c(0,0), coeff_interaction_c_IS(v2[q-1], v2his, bref))
      v2his= c(v2his, v2[q-1]); 
      intg = c(intg, inter[1])
      intw = c(intw, inter[2])
    }else{
      intg <- intw <- 0  
    }
    
    
    Penpiece = penf.comp_IS(f=integrand_IS, v1=v1, v2=v2, est=cest,bref=bref, affp=affp,beta_vec1=beta_vec1,beta_vec2=beta_vec2,beta_vec3=beta_vec3,
                            intg=intg,intw=intw,type=type,lower=ifelse(q==1,0,v1[q-1]),upper=v1[q])
    j = c(j,Penpiece)
  }
  j = sum(j,na.rm = TRUE)
  Penet <- log(j)
  
  # print(round(Penet,3))
  return(Penet)
}
penf.comp_IS <- function(f,v1,v2,est,bref, affp,type,beta_vec1,beta_vec2,beta_vec3,intg,intw,frailty,fp,lower,upper){  
  integrate(f,v1,v2,est=est,bref=bref, affp=affp,beta_vec1=beta_vec1,
            beta_vec2=beta_vec2,beta_vec3=beta_vec3,type=type,intg=intg,intw=intw,
            lower=lower,upper=upper)$value
}
integrand_IS <- function(u,v1,v2,est, bref,affp,type,beta_vec1,beta_vec2,beta_vec3,intg,intw,frailty=FALSE,fp=1){
  u <- u
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
  zsc1ov = ifelse(u>=sc1,1,0)
  
  if(type==1){
    xbeta = est[7] + (est[8] + bref[1])*zsc1+(est[9] + bref[2])*zsc2+(bref[3] + est[10])*zsc3 + 
      est[11]*zor + bref[4]*zor*zsc1 + bref[5]*zor*zsc2  + bref[6]*zor*zsc3
    haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta)
    #print(c(xbeta, est[7]+sum(beta_vec1)))
  }else if(type==2){
    haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(est[12]+est[13]*zsc1ov)  
  }
  
  H1 = - (1-zbr)*(est[1]*u)^est[2]*exp(est[7] + sum(beta_vec1) + sum(intg) + sum(intw)) #+ sum(beta_vec1_int_g) + sum(beta_vec1_int_cov) ) #cause 1
  H2 = - (1-zor)*(est[3]*u)^est[4]*exp(est[12]+ sum(beta_vec2))
  H3 = - (est[5]*u)^est[6]*exp(est[14]+ sum(beta_vec3) ) 
  
  
  # H1 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=1)*exp(est[7])
  # H2 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=2)*exp(est[12])
  # H3 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=TRUE,type=3)*exp(est[14])
  
  Hp <- H1+H2+H3
  
  if(frailty==TRUE){
    f = haz* ( 1 - (Hp)/fp)^(-fp)
  } else{
    f = haz*exp(Hp) # h * S   
  } 
  
  return(f)
}

# Interval specific hazard + ED
ZED <- function(t,u,phi){
  exp(-(t-u)*phi)
}
proband_penf_ED <- function(v, cest, bref, affp, type, phi, frailty, fp){
  v1 <- v[1:5]
  v2 <- v[6:10]
  terminal_time = v[11]
  v1 = v1[v1<terminal_time]
  v1 = c(v1,terminal_time)
  v1 = v1[!is.na(v1)]
  v2 = v2[!is.na(v2)]
  
  j <- beta_vec1 <- beta_vec2 <- beta_vec3 <- inter <- intg <- intw <- NULL
  v2his <- vector(mode="numeric")
  
  beta_sc=0; beta_or=0; lowsc=0; lowor=0; phisc=0; beta_sg=0;
  for(i in 1:length(v1)){
    
    
    if(i!=1){
      beta = coeff_c_ED(v2[i-1],type=1, affp, cest);
      if(v2[i-1]==1 | v2[i-1]==2 | v2[i-1]==3 ){
        phisc = phi[v2[i-1]]
        beta_sg = bref[v2[i-1]]
        beta_sc = beta;
        lowsc = v2[i-1];
      }else if(v2[i-1]==4){
        beta_or = beta;
        lowor = v2[i-1];
      }
    }else{
      beta = 0;
      beta_sc = beta;
      beta_or = beta;
    }
    
    beta1 <- c(beta_sc,lowsc,beta_or,lowor,beta_sg)
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
    Penpiece = penf.comp_ED(f=integrand_IS_ED, v1=v1, v2=v2, est=cest,bref=bref, affp=affp,beta_vec1=beta1,beta_vec2=beta2,beta_vec3=beta3,
                            intg=intg,intw=intw,type=type,phi=phi,phisc=phisc,lb=low,frailty=frailty,fp=fp,lower=low,upper=v1[i])
    j = c(j,Penpiece)
  }
  j = sum(j,na.rm = TRUE)
  Penet <- log(j)
  
  # print(round(Penet,3))
  return(Penet)
}
penf.comp_ED <- function(f,v1,v2,est,bref, affp,type,beta_vec1,beta_vec2,beta_vec3,intg,intw,phi,phisc,lb,frailty,fp,lower,upper){  
  integrate(f,v1,v2,est=est,bref=bref, affp=affp,beta_vec1=beta_vec1,
            beta_vec2=beta_vec2,beta_vec3=beta_vec3,type=type,intg=intg,intw=intw,
            phi=phi,phisc=phisc,lb=lb,lower=lower,upper=upper)$value
}
integrand_IS_ED <- function(u,v1,v2,est, bref,affp,type,beta_vec1,beta_vec2,beta_vec3,intg,intw,phi,phisc,lb,frailty=FALSE,fp=1){
  u <- u
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
  zsc1ov = ifelse(u>=sc1,1,0)
  
  zsc1td = ifelse(zsc1==1,ZED(t=u, u=sc1, phi[1]),0)
  zsc2td = ifelse(zsc2==1,ZED(t=u, u=sc2, phi[2]),0)
  zsc3td = ifelse(zsc3==1,ZED(t=u, u=sc3, phi[3]),0)
  zortd = ifelse(zor==1,ZED(t=u, u=ov, phi[4]),0)
  
  if(type==1){
    xbeta = est[7] + (est[8] + bref[1])*zsc1*zsc1td +(est[9] + bref[2])*zsc2*zsc2td+(bref[3] + est[10])*zsc3*zsc3td + 
      est[11]*zor*zortd + bref[4]*zor*zsc1 + bref[5]*zor*zsc2  + bref[6]*zor*zsc3
    haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*exp(xbeta)
    #print(c(xbeta, est[7]+sum(beta_vec1)))
  }else if(type==2){
    haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(est[12]+est[13]*zsc1ov)  
  }
  
  H1 = - (1-zbr)*exp(est[7])*sapply(u, chaz, lam1=est[1],rho1=est[2],
                                    betasc=beta_vec1[1],betasg=beta_vec1[5],betaor=beta_vec1[3], phisc=phisc,phior=phi[4],
                                    lowsc=beta_vec1[2], lowor=beta_vec1[4],low=lb)
  H2 = - (1-zor)*(est[3]*u)^est[4]*exp(est[12]+ sum(beta_vec2))
  H3 = - (est[5]*u)^est[6]*exp(est[14]+ sum(beta_vec3) ) 
  
  # H1 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=1)*exp(est[7])
  # H2 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=2)*exp(est[12])
  # H3 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=TRUE,type=3)*exp(est[14])
  
  Hp <- H1+H2+H3
  
  if(frailty==TRUE){
    f = haz* ( 1 - (Hp)/fp)^(-fp)
  } else{
    f = haz*exp(Hp) # h * S   
  }
  return(f)
}
