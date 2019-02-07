loglik_Comp_Timedep=function(theta, data, data2, agemin, frailty=FALSE, weight=FALSE){
  
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
  # beta.sc1g_1 = theta[12]
  # beta.sc2g_1 = theta[13]
  # beta.sc3g_1 = theta[14]
  # beta.sc1o_1 = theta[15]
  # beta.sc2o_1 = theta[16]
  # beta.sc3o_1 = theta[17]
  
  # vbeta for ovarian
  beta.gen2 = theta[12]
  beta.sc1_2 = theta[13]
  beta.sc2_2 = 0
  beta.sc3_2 = 0
  beta.sc1g_o = 0
  beta.sc2g_o = 0
  beta.sc3g_o = 0
  
  # vbeta for death
  beta.gen3 = theta[14]
  beta.sc1_3 = theta[15] 
  beta.sc2_3 = 0
  beta.sc3_3 = 0
  beta.sc1g_d = 0
  beta.sc2g_d = 0
  beta.sc3g_d = 0
  beta.s1o_d = 0
  beta.s2o_d = 0
  beta.s3o_d = 0
  
  #  
  beta.sc1g_1 = 0
  beta.sc2g_1 = 0
  beta.sc3g_1 = 0
  beta.sc1o_1 = 0
  beta.sc2o_1 = 0
  beta.sc3o_1 = 0
  
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
  ov_time <- data$ov.censortime-agemin
  sc1_time <- data$st1 - agemin
  sc2_time <- data$st2 - agemin
  sc3_time <- data$st3 - agemin
  delta_BR <- ifelse(time0-br_time > 0,1,0)
  delta_OR <- ifelse(time0-ov_time > 0,1,0)
  
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
  sc_vec_bc[sc_vec==1] <- beta.sc1_1
  sc_vec_bc[sc_vec==2] <- beta.sc1_1 + beta.sc2_1
  sc_vec_bc[sc_vec==3] <- beta.sc1_1 + beta.sc2_1 + beta.sc3_1
  sc_vec_ov <- NULL
  sc_vec_ov[sc_vec==0] <- 0
  sc_vec_ov[sc_vec==1] <- beta.sc1_2
  sc_vec_ov[sc_vec==2] <- beta.sc1_2 + beta.sc2_2
  sc_vec_ov[sc_vec==3] <- beta.sc1_2 + beta.sc2_2 + beta.sc3_2
  sc_vec_d <- NULL
  sc_vec_d[sc_vec==0] <- 0
  sc_vec_d[sc_vec==1] <- beta.sc1_3
  sc_vec_d[sc_vec==2] <- beta.sc1_3 + beta.sc2_3
  sc_vec_d[sc_vec==3] <- beta.sc1_3 + beta.sc2_3 + beta.sc3_3
  
  sg_vec_bc <- NULL
  sg_vec_bc[sc_vec==0] <- 0
  sg_vec_bc[sc_vec==1] <- beta.sc1g_1
  sg_vec_bc[sc_vec==2] <- beta.sc1g_1 + beta.sc2g_1
  sg_vec_bc[sc_vec==3] <- beta.sc1g_1 + beta.sc2g_1 + beta.sc3g_1
  sg_vec_ov <- NULL
  sg_vec_ov[sc_vec==0] <- 0
  sg_vec_ov[sc_vec==1] <- beta.sc1g_o
  sg_vec_ov[sc_vec==2] <- beta.sc1g_o + beta.sc2g_o
  sg_vec_ov[sc_vec==3] <- beta.sc1g_o + beta.sc2g_o + beta.sc3g_o
  sg_vec_d <- NULL
  sg_vec_d[sc_vec==0] <- 0
  sg_vec_d[sc_vec==1] <- beta.sc1g_d
  sg_vec_d[sc_vec==2] <- beta.sc1g_d + beta.sc2g_d
  sg_vec_d[sc_vec==3] <- beta.sc1g_d + beta.sc2g_d + beta.sc3g_d
  
  os_vec_bc <- vector("numeric", length(sc_num) )
  os_vec_bc[sc_num==1 & delta_OR==1] <- beta.sc1o_1
  os_vec_bc[sc_num==2 & delta_OR==1] <- beta.sc1o_1 + beta.sc2o_1 # check this 
  os_vec_bc[sc_num==3 & delta_OR==1] <- beta.sc1o_1 + beta.sc2o_1 + beta.sc3o_1 # check this
  os_vec_d <- vector("numeric", length(sc_num) )
  os_vec_d[sc_num==1 & delta_OR==1] <- beta.s1o_d
  os_vec_d[sc_num==2 & delta_OR==1] <- beta.s1o_d + beta.s2o_d # check this 
  os_vec_d[sc_num==3 & delta_OR==1] <- beta.s1o_d + beta.s2o_d + beta.s3o_d # check this
  
  # sum of log-hazard
  sum1 = sum((wt*ex1*(logh1.1 + (delta_OR)*(beta.or1)+  sc_vec_bc + sg_vec_bc + os_vec_bc ) +
                wt*(1-ex1)*(logh1.0 + (delta_OR)*beta.or1 +  sc_vec_bc + os_vec_bc))[status==1], na.rm=TRUE) 
  
  sum2 = sum((wt*ex1*(logh2.1 + sc_vec_ov + sg_vec_ov ) +
                wt*(1-ex1)*(logh2.0 + sc_vec_ov ))[status==3], na.rm=TRUE)
  
  sum3 = sum((wt*ex1*(logh3.1 + sc_vec_d + os_vec_d + sg_vec_d ) +
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
  
  bref <- c(beta.sc1g_1,beta.sc2g_1,beta.sc3g_1,
            beta.sc1o_1,beta.sc2o_1,beta.sc3o_1)
  
  sum4 <- NULL # check this: vectorize
  for (i in 1:nrow(data)){
    v1 = time_vec2[i,]
    v2 = indicator_vec[i,]
    v1=v1[!is.na(v1)]
    v2=v2[!is.na(v2)]
    CH_bc <- CumH_c(v1=v1,v2=v2,u=time0[i],cest=cest,bref=bref,affp=FALSE,type=1)
    CH_ov <- CumH_c(v1=v1,v2=v2,u=time0[i],cest=cest,bref=bref,affp=FALSE,type=2)
    CH_d  <- CumH_c(v1=v1,v2=v2,u=time0[i],cest=cest,bref=bref,affp=FALSE,type=3)
    
    CH0_bc <- (1-ex1[i])*CH_bc[1]
    CH1_bc <-    ex1[i]*exp(beta.gen1)*CH_bc[2]
    k_bc <- CH0_bc + CH1_bc
    
    CH0_ov <- (1-ex1[i])*CH_ov[1]
    CH1_ov <-    ex1[i]*exp(beta.gen2)*CH_ov[2]
    k_ov <- CH0_ov + CH1_ov
    
    CH0_d <- (1-ex1[i])*CH_d[1]
    CH1_d <-    ex1[i]*exp(beta.gen3)*CH_d[2]
    k_d <- CH0_d + CH1_d
    
    k <- k_bc + k_ov + k_d
    
    sum4 <- c(sum4,k)
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
  
  #setting up indicator variable for the removal status at the time of study entry
  # br_time_p <- (data$br.censortimep-agemin)[ip]
  # ov_time_p <- (data$ov.censortimep-agemin)[ip]
  # sc_time_p <- (data$st1p - agemin)[ip]
  # 
  # delta_BR_p <- ifelse(br_time_p < cagep, 1, 0) # mastectomy status at study entry
  # delta_OR_p <- ifelse(ov_time_p < cagep, 1, 0) # oophorectomy status at study entry
  # delta_SC_p <- data$screenp[ip] 
  
  # subsetting the probands according to his affection status at study entry
  cagep0 <- cagep[statusp==0]
  timep0 <- timep[statusp==0]
  wtp0 <- wtp[statusp==0]
  
  logasc0 <- NULL
  for (i in 1:length(cagep0)){
    v1 = time_vec2_p0[i,]
    v2 = indicator_vec_p0[i,]
    v1=v1[!is.na(v1)]
    v2=v2[!is.na(v2)]
    CH_bc <- CumH_proband_c(v1=v1,v2=v2,u=cagep0[i],cest=cest,bref=bref,affp=FALSE,type=1)
    CH_ov <- CumH_proband_c(v1=v1,v2=v2,u=cagep0[i],cest=cest,bref=bref,affp=FALSE,type=2)
    CH_d  <- CumH_proband_c(v1=v1,v2=v2,u=cagep0[i],cest=cest,bref=bref,affp=FALSE,type=3)
    k <- exp(beta.gen1)*CH_bc+exp(beta.gen2)*CH_ov+exp(beta.gen3)*CH_d
    logasc0 <- c(logasc0,k)
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
  
  # br_time_p1 <- br_time_p[statusp==1]
  # ov_time_p1 <- ov_time_p[statusp==1]
  # sc_time_p1 <- sc_time_p[statusp==1]
  # #delta_BR_p1 <- delta_BR_p[statusp==1]
  # #delta_OR_p1 <- delta_OR_p[statusp==1]
  # delta_SC_p1 <- delta_SC_p[statusp==1]
  # 
  # delta_BR_p1 <- ifelse(timep1-br_time_p1 > 0,1,0)
  # delta_OR_p1 <- ifelse(timep1-ov_time_p1 > 0,1,0)
  # ov_time_p1= ov_time_p1*delta_OR_p1
  
  
  proband_penf <- function(v1, v2, terminal_time, affp, type){
    
    v1 = v1[v1<terminal_time]
    v1 = c(v1,terminal_time)
    v1 = v1[!is.na(v1)]
    v2 = v2[!is.na(v2)]
    
    j <- beta_vec1 <- beta_vec2 <- beta_vec3 <-  NULL
    for(q in 1:length(v1)){
      
      beta1 <- ifelse(q==1, 0, coeff_c(v2[q-1],type=1, cest=cest, affp=affp) )
      beta2 <- ifelse(q==1, 0, coeff_c(v2[q-1],type=2, cest=cest, affp=affp) )
      beta3 <- ifelse(q==1, 0, coeff_c(v2[q-1],type=3, cest=cest, affp=affp) )
      beta_vec1 <- c(beta_vec1,beta1)
      beta_vec2 <- c(beta_vec2,beta2)
      beta_vec3 <- c(beta_vec3,beta3)
      
      Penpiece = penf.comp(f=integrand, v1=v1, v2=v2, est=cest, affp=affp,beta_vec1=beta_vec1,beta_vec2=beta_vec2,beta_vec3=beta_vec3,
                           type=type,lower=ifelse(q==1,0,v1[q-1]),upper=v1[q])
      j = c(j,Penpiece)
    }
    j = sum(j,na.rm = TRUE)
    Penet <- log(j)
    
    # print(round(Penet,3))
    return(Penet)
  }
  penf.comp <- function(f,v1,v2,est,affp,type,beta_vec1,beta_vec2,beta_vec3,lower,upper){  
    try(integrate(f,v1=v1,v2=v2,est=est,affp=affp,beta_vec1=beta_vec1,beta_vec2=beta_vec2,beta_vec3=beta_vec3,type=type,lower=lower,upper=upper)$value)
  }  
  integrand <- function(u,v1,v2,est,affp,type,beta_vec1,beta_vec2,beta_vec3){
    u <- u
    #print(u)
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
      haz=(est[1]^est[2])*est[2]*(u)^(est[2]-1)*
        exp(est[7] + (est[8] + beta.sc1g_1)*zsc1+(est[9] + beta.sc2g_1)*zsc2+( beta.sc3g_1 + est[10])*zsc3 + 
              est[11]*zor + beta.sc1o_1*zor*zsc1 + beta.sc2o_1*zor*zsc2  + beta.sc3o_1*zor*zsc3 )
    }else if(type==2){
      haz=(est[3]^est[4])*est[4]*(u)^(est[4]-1)*exp(est[12]+est[13]*zsc1)  
    }
    
    H1 = - (1-zbr)*(est[1]*u)^est[2]*exp(est[7] + sum(beta_vec1)) #+ sum(beta_vec1_int_g) + sum(beta_vec1_int_cov) ) #cause 1
    H2 = - (1-zor)*(est[3]*u)^est[4]*exp(est[12]+ sum(beta_vec2))
    H3 = - (est[5]*u)^est[6]*exp(est[14]+ sum(beta_vec3) ) 
    
    
    # H1 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=1)*exp(est[7])
    # H2 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=FALSE,type=2)*exp(est[12])
    # H3 = CumH_proband_c(v1=v1,v2=v2,u=u,cest=est,bref=bref,affp=TRUE,type=3)*exp(est[14])
    
    Hp <- H1+H2+H3
    f = haz*exp(Hp) # h * S 
    
    return(f)
  }
  
  logasc1 <- vector(mode="numeric", length(cagep1))
  for (i in 1:length(cagep1)){
    k <- proband_penf(v1=time_vec2_p1[i,],v2=indicator_vec_p1[i,],
                      terminal_time=cagep1[i], type=1, affp = TRUE)
    logasc1[i] <- k
  }
  logasc1 <- sum(wtp1*logasc1,na.rm=TRUE)
  
  
  # subsetting the probands according to his affection status at study entry
  cagep2 <- cagep[statusp==2]
  timep2 <- timep[statusp==2]
  wtp2 <- wtp[statusp==2]
  # 
  # br_time_p2 <- br_time_p[statusp==2]
  # ov_time_p2 <- ov_time_p[statusp==2]
  # sc_time_p2 <- sc_time_p[statusp==2]
  # delta_BR_p2 <- delta_BR_p[statusp==2]
  # delta_OR_p2 <- delta_OR_p[statusp==2]
  # delta_SC_p2 <- delta_SC_p[statusp==2]
  # 
  # delta_BR_p2 <- ifelse(timep2-br_time_p2 > 0,1,0)
  # #delta_OR_p2 <- ifelse(timep2-ov_time_p1 > 0,1,0)
  # br_time_p2= br_time_p2*delta_BR_p2
  
  logasc2 <- vector(mode="numeric",length(cagep2))
  for (i in 1:length(cagep2)){
    k <- proband_penf(v1=time_vec2_p2[i,],v2=indicator_vec_p2[i,],
                      terminal_time=cagep2[i], type=2, affp = TRUE)
    logasc2[i] <- k
  }
  logasc2 <- sum(wtp2*logasc2,na.rm=TRUE)
  
  slogasc = logasc0 + logasc1 + logasc2
  likelihood  <- loglik - slogasc
  print(c(loglik,slogasc,-likelihood))
  return(-likelihood)
} 
