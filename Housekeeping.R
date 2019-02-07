library(Rcpp)
sourceCpp("C:/Users/Jay/Desktop/CompRisk/R/cumh_ED.cpp")
#sourceCpp("C:/Users/Jay/Desktop/CompRisk/R/cumh_IS.cpp")
source("C:/Users/Jay/Desktop/CompRisk/R/Fit_Model.R")
source("C:/Users/Jay/Desktop/CompRisk/R/Logliklihood.R")
### Raw data
BCraw <- read.csv("C:/Users/jay/Desktop/New Era/July31-2017/July31-2017DM/df6.csv")

### Load data
b1 <- read.csv("C:/Users/jay/Dropbox/Jay/New Era/CompetingRisk Logliklihood/Data/dfcomp_ScreenUpto5_March9.csv")
# load data
b1_NC <- read.csv("C:/Users/jay/Dropbox/Jay/New Era/CompetingRisk Logliklihood/dfbreast_removal_ScreenUpto3_Jan16.csv")


### Fit model Competing Risk
brca1.1frailty6 <- FitModelEM( data = b1,
                        
                      #breast (ovarian, death) ~ mgene (TID) + screen (TD) + oophorectomy (TD)          
                      
                      init.Parms = list( 
                        cause1 = list(base = c(0.008,2.581), time_indep=c(1.960), time_dep=list(main=c(2,-2,1,-1) , interaction=c() )),
                        cause2 = list(base = c(0.008,3.220), time_indep=c(1.331), time_dep=list(main=0)),
                        cause3 = list(base = c(0.016,4.061), time_indep=c(-0.095), time_dep=list(main=0)) 
                      ),
                      
                      time.dep.cov = list(
                                          time.dep.cov.type = "ED" ,    # choose among ("PE","ED","CO")
                                                reccur.type = "IS"      # choose among ("CM","IS")
                                          ),
                      missing.method = "data", # "mendelian"
                      base.dist = "Weibull", # "Piece_Wise", 
                      frailty=TRUE, #
                      weight=FALSE )

brca1.9frailty_NC = FitModelEM_NC( data = b1_NC,#### <- run this line for model fitting ####
                                   #You can change initial parameter Here
                                   init.Parms = list(        #lam1  rho1               gene                  screen1,2,3,oopho         
                                     cause1 = list(base = c(0.008,2.22), time_indep=c(1.4), time_dep=list(main=c(1,1,1,-0.6) , interaction=c() ))
                                   ),
                                   
                                   time.dep.cov = list(time.dep.cov.type = "ED" , #You can change between "PE" and "ED"
                                                       reccur.type = "IS"   #You can change intercal specific(IS) or cumulative effect(CM) 
                                   ),
                                   
                                   missing.method = "data",
                                   base.dist = "Weibull",
                                   frailty=TRUE, # Fit with Gamma frailty? Yes=TRUE   No=FALSE
                                   weight=FALSE,  # Fit with Sampling weight? Yes=TRUE   No=FALSE
                                   mutation.prediction=FALSE) # Mutation status prediction results


7720.937 indep 7581.388 frailty
### Debugging (Comprisk model) Run below lines for variable settings for likelihood function
missing.method = "data";base.dist = "Weibull";frailty=FALSE;weight=FALSE
data=carrierprobgeno(data=b1); data2=Data_preparation(b1); agemin=16;
time.dep.cov = list(time.dep.cov.type = "ED" ,reccur.type = "IS")


base.parms <- c(init.Parms[[1]]$base, init.Parms[[2]]$base, init.Parms[[3]]$base) # total 6, lambda1,rho1 lambda2,rho2, lambda3,rho3 
vbeta_b <- c( init.Parms[[1]]$time_indep, init.Parms[[1]]$time_dep$main )#, init.Parms[[1]]$time_dep$interaction ) # breast cancer parameter
vbeta_o <- c( init.Parms[[2]]$time_indep, init.Parms[[2]]$time_dep$main )                               # ovarian cancer parameter 
vbeta_d <- c( init.Parms[[3]]$time_indep, init.Parms[[3]]$time_dep$main )                               # death parameter
vbeta <- c(vbeta_b,vbeta_o,vbeta_d)
theta = theta0 = c(log(base.parms),vbeta)   

### Debugging (Noncompeting risk model), Run below lines for variable settings for likelihood function
missing.method = "data";base.dist = "Weibull";frailty=FALSE;weight=FALSE
data=carrierprobgeno_NC(data=b1_NC); data2=Data_preparation_NC(b1_NC); agemin=16;
time.dep.cov = list(time.dep.cov.type = "PE" ,reccur.type = "IS")


base.parms <- c(init.Parms[[1]]$base, init.Parms[[2]]$base, init.Parms[[3]]$base) # total 6, lambda1,rho1 lambda2,rho2, lambda3,rho3 
vbeta_b <- c( init.Parms[[1]]$time_indep, init.Parms[[1]]$time_dep$main )#, init.Parms[[1]]$time_dep$interaction ) # breast cancer parameter                           # death parameter
vbeta <- c(vbeta_b)
theta = theta0 = c(log(base.parms),vbeta)   

#Results
plot(seq(0,20,by=0.01),dgamma(seq(0,20,by=0.01), shape=0.1 , scale=1/0.1))#stong dependency
plot(seq(0,5,by=0.01),dgamma(seq(0,5,by=0.01), shape=1.86 , scale=1/1.86))#breast
plot(seq(0,20,by=0.01),dgamma(seq(0,20,by=0.01), shape=2.6 , scale=1/2.6))#ovarian
plot(seq(0,20,by=0.01),dgamma(seq(0,20,by=0.01), shape=30 , scale=1/30))#Death weak dependency?

