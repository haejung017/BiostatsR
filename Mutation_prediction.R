sensitivity 0.98790323
specificity 0.02173913
PPV         0.47572816
NPV         0.66666667
ACC         0.47900763

[,1]
sensitivity 0.9395161
specificity 0.1050725
PPV         0.4854167
NPV         0.6590909
ACC         0.5000000


exp_loglik=function(theta, data, data2, agemin, frailty=FALSE, weight=FALSE, ex1){
  
  #print(round(c(exp(theta[1:6]),theta[7:length(theta)]),3))
  
  data = data[data$currentage>=agemin,]
  # base parameters
  lambda1 = theta[1]
  rho1 = theta[2]
  lambda2 = theta[3]
  rho2 = theta[4]
  lambda3 = theta[5]
  rho3 = theta[6]
  
  # vbeta for breast  
  beta.gen1 = theta[7]
  beta.sc1_1 = theta[8]
  beta.sc2_1 = theta[9]
  beta.sc3_1 = theta[10]
  beta.or1 = theta[11]
  beta.sc1g_1 = theta[12]
  beta.sc2g_1 = theta[13]
  beta.sc3g_1 = theta[14]
  beta.sc1o_1 = theta[15]
  beta.sc2o_1 = theta[16]
  beta.sc3o_1 = theta[17]
  
  # vbeta for ovarian
  beta.gen2 = theta[18]
  beta.sc1_2 = theta[19]
  beta.sc2_2 = 0
  beta.sc3_2 = 0
  beta.sc1g_o = 0
  beta.sc2g_o = 0
  beta.sc3g_o = 0
  
  # vbeta for death
  beta.gen3 = theta[20]
  beta.sc1_3 = theta[21] 
  beta.sc2_3 = 0
  beta.sc3_3 = 0
  beta.sc1g_d = 0
  beta.sc2g_d = 0
  beta.sc3g_d = 0
  beta.s1o_d = 0
  beta.s2o_d = 0
  beta.s3o_d = 0
  
  #  
  # beta.sc1g_1 = 0
  # beta.sc2g_1 = 0
  # beta.sc3g_1 = 0
  # beta.sc1o_1 = 0
  # beta.sc2o_1 = 0
  # beta.sc3o_1 = 0
  
  bref <- c(beta.sc1g_1,beta.sc2g_1,beta.sc3g_1,
            beta.sc1o_1,beta.sc2o_1,beta.sc3o_1)
  
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
  ex1 = ex1
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
  sum1 = (wt*ex1*(logh1.1 + (delta_OR)*beta.or1 +  sc_vec_bc + sg_vec_bc + os_vec_bc ) +
                wt*(1-ex1)*(logh1.0 + (delta_OR)*beta.or1 +  sc_vec_bc + os_vec_bc))[status==1] 
  
  sum2 = (wt*ex1*(logh2.1 + sc_vec_ov + sg_vec_ov ) +
                wt*(1-ex1)*(logh2.0 + sc_vec_ov ))[status==3]
  
  sum3 = (wt*ex1*(logh3.1 + sc_vec_d + os_vec_d + sg_vec_d ) +
                wt*(1-ex1)*(logh3.0 + sc_vec_d + os_vec_d ))[status==4]
  
  sumtotal <- rep(0,nrow(data))
  sumtotal[status==1]<- sum1
  sumtotal[status==3]<- sum2
  sumtotal[status==4]<- sum3
  
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
  sum4 = apply(mat_all, 1, func_all, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam <- -sum4
    df <- data$df[data$proband==1]
    Hfam <- aggregate(Hfam,by=list(data$famID),FUN=sum)[,2]
    sum4 <- sum(wtp*(lfactorial(fp+df-1)-(df-1)*log(fp)-lfactorial(fp) + (-fp-df)*log(1+(Hfam)/fp)), na.rm=T)
    loglik = sum1+sum2+sum3+sum4
  }else{
    #sum4 <- sum(wt*sum4,na.rm = TRUE)
    loglik = sumtotal+sum4 # numerator in loglikelihood
  }
  
 
  
  #print(c(loglik, slogasc, likelihood))
  return(exp(loglik))
  
} 
#carrpgeno
IDs <- which(!is.na(b1_comp$mgene)&b1_comp$proband!=1)
theta <- as.vector(brca1.1int[,1])
agemin <- 16  
data <- carrierprobgeno_1(data=b1_comp,method="mendelian")
data2 <- Data_preparation(data)
##
p1 <- exp_loglik(theta=theta, data=data, data2=data2,agemin,ex1=rep(1,nrow(data)))*data$carrp.geno
p0 <- exp_loglik(theta=theta, data=data, data2=data2,agemin,ex1=rep(0,nrow(data)))*(1-data$carrp.geno)

p1[!is.na(data$mgene) & data$mgene==1] <- 1
p1[!is.na(data$mgene) & data$mgene==0] <- 0
p0[!is.na(data$mgene) & data$mgene==1] <- 0
p0[!is.na(data$mgene) & data$mgene==0] <- 1

ex1 <- p1/(p1+p0)
ex1_ori <- ex1
for (id in IDs){
  data <- b1_comp
  data[id,"mgene"] <- NA
  data <- carrierprobgeno_1(data=data,method="mendelian")
  
  p1 <- exp_loglik(theta=theta, data=data, data2=data2,agemin,ex1=rep(1,nrow(data)))*data$carrp.geno
  p0 <- exp_loglik(theta=theta, data=data, data2=data2,agemin,ex1=rep(0,nrow(data)))*(1-data$carrp.geno)
  
  p1[!is.na(data$mgene) & data$mgene==1] <- 1
  p1[!is.na(data$mgene) & data$mgene==0] <- 0
  p0[!is.na(data$mgene) & data$mgene==1] <- 0
  p0[!is.na(data$mgene) & data$mgene==0] <- 1
  cal <- (p1/(p1+p0))[id]
  ex1[id]<- cal
}

diag <- cbind(ex1_ori[IDs],ifelse(ex1[IDs]>=0.2,1,0))
diagstat <- t(c(sum(diag[,1]==diag[,2] & diag[,1]==1)/sum(diag[,1]), #Sensitivity
                sum(diag[,1]==diag[,2] & diag[,1]==0)/sum(diag[,1]==0), #Specificity
                sum(diag[,1]==diag[,2] & diag[,1]==1)/sum(diag[,2]), # Positive predictive value
                sum(diag[,1]==diag[,2] & diag[,1]==0)/sum(diag[,2]==0), # Negative predictive value
                sum(diag[,1]==diag[,2])/length(ex1[IDs]))) # Accuracy
colnames(diagstat) <- c("sensitivity","specificity","PPV","NPV","ACC")
diagstat <- as.data.frame(diagstat)
t(diagstat)
pe_brca1_1 <- mean((abs(ex1_ori[IDs]-ex1[IDs]))^2)


# [,1]
# sensitivity 0.8951613
# specificity 0.3369565
# PPV         0.5481481
# NPV         0.7815126
# ACC         0.6011450
# > mean((abs(ex1_ori[IDs]-ex1[IDs]))^2)
# [1] 0.2164127


exp_loglik=function(theta, data, data2, agemin, frailty=FALSE, weight=FALSE, ex1){
  
  #print(round(c(exp(theta[1:6]),theta[7:length(theta)]),3))
  
  data = data[data$currentage>=agemin,]
  # base parameters
  lambda1 = theta[1]
  rho1 = theta[2]
  lambda2 = theta[3]
  rho2 = theta[4]
  lambda3 = theta[5]
  rho3 = theta[6]
  
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
  
  bref <- c(beta.sc1g_1,beta.sc2g_1,beta.sc3g_1,
            beta.sc1o_1,beta.sc2o_1,beta.sc3o_1)
  
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
  ex1 = ex1
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
  sum1 = (wt*ex1*(logh1.1 + (delta_OR)*beta.or1 +  sc_vec_bc + sg_vec_bc + os_vec_bc ) +
            wt*(1-ex1)*(logh1.0 + (delta_OR)*beta.or1 +  sc_vec_bc + os_vec_bc))[status==1] 
  
  sum2 = (wt*ex1*(logh2.1 + sc_vec_ov + sg_vec_ov ) +
            wt*(1-ex1)*(logh2.0 + sc_vec_ov ))[status==3]
  
  sum3 = (wt*ex1*(logh3.1 + sc_vec_d + os_vec_d + sg_vec_d ) +
            wt*(1-ex1)*(logh3.0 + sc_vec_d + os_vec_d ))[status==4]
  
  sumtotal <- rep(0,nrow(data))
  sumtotal[status==1]<- sum1
  sumtotal[status==3]<- sum2
  sumtotal[status==4]<- sum3
  
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
  sum4 = apply(mat_all, 1, func_all, cest=cest, bref=bref, beta.gen1=beta.gen1, beta.gen2=beta.gen2, beta.gen3=beta.gen3)
  
  if(frailty==TRUE){ # Gamma-frailty
    Hfam <- -sum4
    df <- data$df[data$proband==1]
    Hfam <- aggregate(Hfam,by=list(data$famID),FUN=sum)[,2]
    sum4 <- sum(wtp*(lfactorial(fp+df-1)-(df-1)*log(fp)-lfactorial(fp) + (-fp-df)*log(1+(Hfam)/fp)), na.rm=T)
    loglik = sum1+sum2+sum3+sum4
  }else{
    #sum4 <- sum(wt*sum4,na.rm = TRUE)
    loglik = sumtotal+sum4 # numerator in loglikelihood
  }
  
  
  
  #print(c(loglik, slogasc, likelihood))
  return(exp(loglik))
  
} 
#carrpgeno
IDs <- which(!is.na(b1_comp$mgene)&b1_comp$proband!=1)
theta <- as.vector(brca1.1[,1])
agemin <- 16  
data <- carrierprobgeno_1(data=b1_comp,method="mendelian")
data2 <- Data_preparation(data)
##
p1 <- exp_loglik(theta=theta, data=data, data2=data2,agemin,ex1=rep(1,nrow(data)))*data$carrp.geno
p0 <- exp_loglik(theta=theta, data=data, data2=data2,agemin,ex1=rep(0,nrow(data)))*(1-data$carrp.geno)

p1[!is.na(data$mgene) & data$mgene==1] <- 1
p1[!is.na(data$mgene) & data$mgene==0] <- 0
p0[!is.na(data$mgene) & data$mgene==1] <- 0
p0[!is.na(data$mgene) & data$mgene==0] <- 1

ex1 <- p1/(p1+p0)
ex1_ori <- ex1
for (id in IDs){
  data <- b1_comp
  data[id,"mgene"] <- NA
  data <- carrierprobgeno_1(data=data,method="mendelian")
  
  p1 <- exp_loglik(theta=theta, data=data, data2=data2,agemin,ex1=rep(1,nrow(data)))*data$carrp.geno
  p0 <- exp_loglik(theta=theta, data=data, data2=data2,agemin,ex1=rep(0,nrow(data)))*(1-data$carrp.geno)
  
  p1[!is.na(data$mgene) & data$mgene==1] <- 1
  p1[!is.na(data$mgene) & data$mgene==0] <- 0
  p0[!is.na(data$mgene) & data$mgene==1] <- 0
  p0[!is.na(data$mgene) & data$mgene==0] <- 1
  cal <- (p1/(p1+p0))[id]
  ex1[id]<- cal
}

diag <- cbind(ex1_ori[IDs],ifelse(ex1[IDs]>=0.2,1,0))
diagstat <- t(c(sum(diag[,1]==diag[,2] & diag[,1]==1)/sum(diag[,1]), #Sensitivity
                sum(diag[,1]==diag[,2] & diag[,1]==0)/sum(diag[,1]==0), #Specificity
                sum(diag[,1]==diag[,2] & diag[,1]==1)/sum(diag[,2]), # Positive predictive value
                sum(diag[,1]==diag[,2] & diag[,1]==0)/sum(diag[,2]==0), # Negative predictive value
                sum(diag[,1]==diag[,2])/length(ex1[IDs]))) # Accuracy
colnames(diagstat) <- c("sensitivity","specificity","PPV","NPV","ACC")
diagstat <- as.data.frame(diagstat)
t(diagstat)
pe_brca1_1 <- mean((abs(ex1_ori[IDs]-ex1[IDs]))^2)

[,1]
sensitivity 0.9233871
specificity 0.3079710
PPV         0.5452381
NPV         0.8173077
ACC         0.5992366
mean((abs(ex1_ori[IDs]-ex1[IDs]))^2)
