# joint probability of genotype and phenotype
cprob <- function(theta, data, mut, base.dist, agemin)
{

beta.sex <- theta[3]
beta.gen <- theta[4]

xbeta <- beta.sex*data$gender+beta.gen*mut

y <- ifelse(data$time-agemin<0,0,data$time-agemin)
delta <- data$status
haz<-hazards(base.dist, y, theta)*exp(xbeta)
Haz <-cumhaz(base.dist, y, theta)*exp(xbeta)

if(mut==1) p <- data$carrp.geno
else p <- 1-data$carrp.geno	

return((haz^delta)*exp(-Haz)*p)

}
