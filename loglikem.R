loglikem<- function(theta, theta0, data, design, base.dist, agemin, vec=TRUE)
{

etheta <- exp(theta)
etheta0 <- exp(theta0)
beta.sex <- theta[3]
beta.gen <- theta[4]


time0 <- data$time-agemin
status<- data$status
wt <- data$weight

xbeta1 <- beta.sex*data$gender+beta.gen*1
xbeta0 <- beta.sex*data$gender+beta.gen*0

bhaz <- hazards(base.dist, time0, etheta[1:2])
bcumhaz <- cumhaz(base.dist, time0, etheta[1:2])

H1 <- bcumhaz*exp(xbeta1)
H0 <- bcumhaz*exp(xbeta0)
logh1 <- log(bhaz) + xbeta1
logh0 <- log(bhaz) + xbeta0


  p1 <- cprob(c(etheta0[1:2],theta0[3:4]), data=data, mut=1, base.dist=base.dist, agemin=agemin)
  p0 <- cprob(c(etheta0[1:2],theta0[3:4]), data=data, mut=0, base.dist=base.dist, agemin=agemin)
  ex1 <- p1/(p1+p0) #P(x=1|Xp, y)=P(y|x=1)*P(x=1|Xp)/(p1+p0) for EM

#ex1[!is.na(data$mgene)] <- data$mgene[!is.na(data$mgene)]

p1[!is.na(data$mgene) & data$mgene==1] <- 1
p1[!is.na(data$mgene) & data$mgene==0] <- 0
p0[!is.na(data$mgene) & data$mgene==1] <- 0
p0[!is.na(data$mgene) & data$mgene==0] <- 1



  loglik <- wt * (- H1 + status*logh1 ) *p1 + wt * (- H0 + status*logh0 ) * p0
  loglik[data$time<=agemin] <- 0
  
# Ascertainment correction by design

ip <- which(data$proband==1)
cagep <- data$currentage[ip]-agemin
xbeta.p <- beta.sex*data$gender[ip]+beta.gen*data$mgene[ip]
bcumhaz.p <- cumhaz(base.dist, cagep, etheta[1:2])
wt.p <- data$weight[ip]

logasc.p <- wt.p*log(1-exp(-bcumhaz.p*exp(xbeta.p))) 

if(design=="cli" | design=="cli+"){
  
  i.m <- data$generation==1 & data$gender==0
  i.f <- data$generation==1 & data$gender==1
  i.s <- data$generation==2 & data$proband==0 & data$status==1
  
cage.m <- data$currentage[i.m]-agemin
xbeta.m0 <- beta.sex*data$gender[i.m]+beta.gen*0
xbeta.m1 <- beta.sex*data$gender[i.m]+beta.gen*1
bcumhaz.m <- cumhaz(base.dist, cage.m, etheta[1:2])

cage.f <- data$currentage[i.f]-agemin
xbeta.f0 <- beta.sex*data$gender[i.f]+beta.gen*0
xbeta.f1 <- beta.sex*data$gender[i.f]+beta.gen*1
bcumhaz.f <- cumhaz(base.dist, cage.f, etheta[1:2])
                    
cage.s <- data$currentage[i.s]-agemin
xbeta.s0 <- beta.sex*data$gender[i.s]+beta.gen*0
xbeta.s1 <- beta.sex*data$gender[i.s]+beta.gen*1
bcumhaz.s <- cumhaz(base.dist, cage.s, etheta[1:2])
                                        
wt.m <- data$weight[i.m]
wt.f <- data$weight[i.f]
wt.s <- data$weight[i.s]

#loglik <- wt * (- H1 + status*logh1 ) *ex1 + wt * (- H0 + status*logh0 ) * (1-ex1)

logasc.m <-  wt.m*log(1-exp(-bcumhaz.m*exp(xbeta.m1)))*ex1[i.m] + wt.p*log(1-exp(-bcumhaz.m*exp(xbeta.m0)))*(1-ex1[i.m])
logasc.f <-  wt.f*log(1-exp(-bcumhaz.f*exp(xbeta.f1)))*ex1[i.f] + wt.p*log(1-exp(-bcumhaz.f*exp(xbeta.f0)))*(1-ex1[i.f]) 
logasc.s <-  wt.s*log(1-exp(-bcumhaz.s*exp(xbeta.s1)))*ex1[i.s] + wt.p*log(1-exp(-bcumhaz.s*exp(xbeta.s0)))*(1-ex1[i.s])

sumlogasc <- sum(logasc.p, na.rm=TRUE) + sum(logasc.m,na.rm=TRUE) + sum(logasc.f,na.rm=TRUE) + sum(logasc.s,na.rm=TRUE)
loglik[i.m] <- loglik[i.m] - logasc.m
loglik[i.f] <- loglik[i.f] - logasc.f
loglik[i.s] <- loglik[i.s] - logasc.s

}
else sumlogasc <- sum(logasc.p, na.rm=TRUE)

sumloglik<- sum(loglik, na.rm=TRUE)-sumlogasc
loglik[ip] <- loglik[ip] - logasc.p

if(vec) return(-loglik)
else return(-sumloglik)
            
}