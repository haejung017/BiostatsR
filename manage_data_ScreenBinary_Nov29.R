rm(list=ls())

setwd("/Volumes/UNTITLED/projets_de_recherche/2018/Laurent/TVC/Rcode/Rcode_Jan/RcodeOR_SCR_tvc_2Interacts")

library(numDeriv)
library(MASS)

############################### Data management #########################

# load data
df.comb <- read.csv("dfbreast_removal_ScreenBinary_Nov29.csv")

str(df.comb)

## number of missing genotype data
table(df.comb$mgene, useNA = "always")

# missing data imputed by expected values
source("Imputationfunc.R")
#data.classifier(df.comb)

library(survival)
#testset.sim <- list(b1,b2)
#for(k in 1:length(testset.sim)){
  
# p.m1aff <- testset.sim[[k]]
  p.m1aff <- df.comb 

  #Imputing using data
  np.m <- carrierprobgeno(method = "data", data=p.m1aff)
  set.seed(7)
  x <- c(0,1)
  for (i in 1:nrow(np.m)){
    if (is.na(np.m[i,"mgene"])){
      np.m[i,"mgene"] <- sample(x, size=1, prob= c(1-np.m[i,"carrp.geno"], np.m[i,"carrp.geno"]))
    }
  }
  #testset.sim[[k]] <- np.m
#} 


#data <- testset.sim[[1]]
 data <- np.m
  
# check if  missing genotype data
table(data$mgene, useNA = "always")
head(data)
length(unique(data$famID))   # 469

# binary variable that indicates subjects with mastectomy
BR.Ind=ifelse(data$time-data$br.censortime<=0,0,1)

## the number of mastectomy subjects
table(BR.Ind)

# the mastectomy subjects are censored at breast censored time
data$time[BR.Ind==1]=data$br.censortime[BR.Ind==1]
data$status[BR.Ind==1]=0
