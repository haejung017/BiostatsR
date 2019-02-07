# joint probability of genotype and phenotype
cprob.frailty <- function(theta, data, mut, base.dist, frailty.dist, agemin)
{
  
  beta.sex <- theta[3]
  beta.gen <- theta[4]
  kappa = theta[5]
  
  xbeta <- beta.sex*data$gender+beta.gen*mut
  
  y <- ifelse(data$time-agemin<0,0,data$time-agemin)
  delta <- data$status
  haz<-hazards(base.dist, y, theta)*exp(xbeta)
  Haz <-cumhaz(base.dist, y, theta)*exp(xbeta)
  
  S <- laplace(frailty.dist, g=Haz, k=kappa)
  h <- haz*dlaplace(frailty.dist, g=Haz, d=1, k=kappa)
  
  if(mut==1) p <- data$carrp
  else p <- 1-data$carrp	
  
  return((h^delta)*S*p)
  
}
