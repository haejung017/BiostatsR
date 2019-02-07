
rm(list=ls())

setwd("/Volumes/UNTITLED/projets_de_recherche/2018/Laurent/TVC/Rcode/Rcode_Jan/manage_data")

library(numDeriv)
library(MASS)

############################### Data management #########################

# load data
data <- read.csv("dfbreast_removal_ScreenBinary_Nov29.csv")

length(unique(data$famID))   # 469

str(data)  # 2484 observation

# Check if BRCA1 or BRCA2 ? 
table(data$BRCA_FAMILY)  # only BRCA1 family

###### descriptive statistic

# Event time
Y=data$time   
status=data$status

G=data$mgene         
# check if  missing genotype data
table(G, useNA = "always")  #  1538 missing genotype

#screening
scrI=data$gender
table(scrI)             # 536 screened subjects and 1948 nonscreened

# bilateral mastectomy 
Ybr=data$br.censortime
brI=ifelse(Y-Ybr<=0,0,1)
table(brI)             # 67  mastectomy subjects  (18 in your report)

# bilateral oopharectomy 
Yor=data$ov.censortime
orI=ifelse(Y-Yor<=0,0,1)
table(orI)             # 225  oopharectomy subjects (181 in your report)



#1) Events age mean for censored and breast cancer women
Ym=data.frame(Censored=mean(Y[status==0],na.rm=T), BreastCancer=mean(Y[status==1],na.rm=T), Total=mean(Y,na.rm=T))
Ym

#2) number of carrier and Noncarrier and unknow for censored and breast cancer women
G[is.na(G)]=2     ### for missing genotype       
table(G,status)
GvsBC=rbind(cbind(table(G, status), rowSums(table(G, status)) ), Total=c( colSums(table(G, status)),Total=sum(table(G, status)) )  )

Odds1=round(GvsBC[2,]/GvsBC[1,],digits=2)
OR1=round(Odds1/Odds1[1],digits=2)
tab1=rbind(GvsBC,Odds1,OR1)
tab1



#3) number of screened and Nonscreen for censored and breast cancer women

table(scrI)      # Total = 536 are screened

table(scrI,status)
SCvsBC=rbind(cbind(table(scrI, status), rowSums(table(scrI, status)) ), Total=c( colSums(table(scrI, status)),Total=sum(table(scrI, status)) )  )

Odds2=round(SCvsBC[2,]/SCvsBC[1,],digits=2)
# The OR comparing 1,2,  5 SA to 0
OR2=round(Odds2/Odds2[1],digits=2)
tab2=rbind(SCvsBC,Odds2,OR2)
tab2

#4) number of mastectomy  for censored and breast cancer women

table(brI)      # Total = 67 bilateral mastectomy

table(brI,status)
BRvsBC=rbind(cbind(table(brI, status), rowSums(table(brI, status)) ), Total=c( colSums(table(brI, status)),Total=sum(table(brI, status)) )  )

Odds3=round(BRvsBC[2,]/BRvsBC[1,],digits=2)
OR3=round(Odds3/Odds3[1],digits=2)
tab3=rbind(BRvsBC,Odds3,OR2)
tab3

#5) number of oopharectomy  for censored and breast cancer women

table(orI)      # Total = 225 bilateral oopharectomy

table(orI,status)
BRvsBC=rbind(cbind(table(brI, status), rowSums(table(brI, status)) ), Total=c( colSums(table(brI, status)),Total=sum(table(brI, status)) )  )

Odds3=round(BRvsBC[2,]/BRvsBC[1,],digits=2)
OR3=round(Odds3/Odds3[1],digits=2)
tab3=rbind(BRvsBC,Odds3,OR2)
tab3





# binary variable that indicates subjects with mastectomy
BR.Ind=ifelse(data$time-data$br.censortime<=0,0,1)

## the number of mastectomy subjects
table(BR.Ind)

# the mastectomy subjects are censored at breast censored time
data$time[BR.Ind==1]=data$br.censortime[BR.Ind==1]
data$status[BR.Ind==1]=0
