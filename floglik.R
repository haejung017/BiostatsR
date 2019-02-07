floglik<- function(theta, data, design, base.dist, frailty.dist, agemin)
{
#data <- data[data$currentage>agemin,]
lambda   <- exp(theta[1])
rho  <- exp(theta[2])
beta.sex <- theta[3]
beta.gen <- theta[4]
kappa <- exp(theta[5])


time0 <- data$time-agemin
status<- data$status
wt <- data$weight
xbeta <- beta.sex*data$gender+beta.gen*data$mgene

bhaz <- hazards(base.dist, time0, c(lambda,rho))
bcumhaz <- cumhaz(base.dist, time0, c(lambda,rho))

H <- bcumhaz*exp(xbeta)
s <- aggregate(H, list(data$famID), sum)[,2]
df <- aggregate(status, list(data$famID), sum)[,2]

logdL <- log( dlaplace(frailty.dist, g=s, d=df, k=kappa) )

logh <- log(bhaz) + xbeta


#loglik <-  wt * (- H + status*logh )

loglik <-  wt*status*logh 

loglik[data$time<=agemin] <- 0

ip <- data$proband==1
wt.p <- wt[ip]

sumiloglik <- sum(loglik, na.rm=T)
sumfloglik <- sum(wt.p*logdL, na.rm=T) 


# Ascertainment correction by design

if(design=="pop" | design=="pop+"){
cagep <- data$currentage[ip]-agemin
xbeta.p <- beta.sex*data$gender[ip]+beta.gen*data$mgene[ip]
bcumhaz.p <- cumhaz(base.dist, cagep,c(lambda,rho))


slogasc <- sum(wt.p*log(1-laplace(frailty.dist, bcumhaz.p*exp(xbeta.p), kappa)  ), na.rm=T)  
}
else{
 stop("Frailty model is only available for POP or POP+ design.")
}

likelihood  <- sumiloglik + sumfloglik - slogasc
return(-likelihood)
}