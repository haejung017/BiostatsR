penplot <- function(base.parms, vbeta, variation="none", base.dist="Weibull", frailty.dist=NULL, depend=1, agemin=20, print=TRUE, ...){

x <- agemin:80
t0 <- x-agemin

if(variation=="frailty" & !any(frailty.dist==c("lognormal",  "gamma"))) 
  stop("Unrecognized frailty distribution; frailty.dist should be either \"gamma\" or \"lognormal\" ")
else if(variation!="frailty" & any(frailty.dist==c("lognormal",  "gamma")) ) stop("Variation should be specified as variation = \"frailty\".")
else if(variation!="frailty" & !is.null(frailty.dist) ) stop("frailty.dist can be specified only with variation=\"frailty\" ")


if(variation=="secondgene") {
  if(length(vbeta)==3) xbeta <- c( vbeta %*% t(expand.grid(c(0,1),c(0,1), c(0,1))) )
  else if(length(vbeta)<3) stop("vbeta should include a second gene effect.")
  else stop("vbeta should be a vector of length 3.")
} 
else {
  if(length(vbeta)==2) xbeta <- c( vbeta %*% t(expand.grid(c(0,1),c(0,1)) ) )
  else stop("vbeta should be a vector of length 2.")
}
s1 <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[1]) # female non-carr
s2 <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[2]) # male non-carr
s3 <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[3]) # female carr
s4 <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[4]) # male carr
s <- list(s1,s2,s3,s4)
if(variation=="secondgene"){
	# second gene carriers
s[[5]] <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[5]) # female non-carr
s[[6]] <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[6]) # male non-carr
s[[7]] <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[7]) # female carr
s[[8]] <- cumhaz(base.dist, t0, base.parms)*exp(xbeta[8]) # male carr	
}

if(variation=="none" | variation=="secondgene") pen <- lapply(s, function(s) 1-exp(-s))
else if(variation=="frailty") pen <- lapply(s, function(dist,s,k) 1-laplace(dist,s,k), dist=frailty.dist, k=depend)
else stop("Unrecognized variation")
 
if(variation=="secondgene"){
par(mfrow=c(1,2))
	plot(x,pen[[5]], ylab="Penetrance", xlab="age at onset", ylim=c(0,1), type="l", lty=2, col="red", main=paste("Penetrances for second gene carriers \n", base.dist, "baseline"), ...)
  lines(x, pen[[6]], lty=2, col="blue")
  lines(x, pen[[7]], lty=1, col="red")
  lines(x, pen[[8]], lty=1, col="blue")
  legend("topleft", c("male carrier", "female carrier", "male noncarrier", "female noncarrier"), bty="n", lty=c(1,1,2,2), col=c("blue","red","blue","red"))

	plot(x,pen[[1]], ylab="Penetrance", xlab="age at onset", ylim=c(0,1), type="l", lty=2, col="red", main=paste("Penetrances for second gene non-carriers \n", base.dist, "baseline"), ...)
  lines(x, pen[[2]], lty=2, col="blue", ...)
  lines(x, pen[[3]], lty=1, col="red", ...)
  lines(x, pen[[4]], lty=1, col="blue", ...)
  legend("topleft", c("male carrier", "female carrier", "male noncarrier", "female noncarrier"), bty="n", lty=c(1,1,2,2), col=c("blue","red","blue","red"))

}
else if(variation=="frailty"){
  par(mfrow=c(1,1))
  plot(x,pen[[1]], ylab="Penetrance", xlab="age at onset", ylim=c(0,1), type="l", lty=2, col="red", main=paste("Penetrance curves \n", base.dist, "baseline and ", frailty.dist,"frailty"), ...)
  lines(x, pen[[2]], lty=2, col="blue", ...)
  lines(x, pen[[3]], lty=1, col="red", ...)
  lines(x, pen[[4]], lty=1, col="blue", ...)
  legend("topleft", c("male carrier", "female carrier", "male noncarrier", "female noncarrier"), bty="n", lty=c(1,1,2,2), col=c("blue","red","blue","red"))
}
else {
  par(mfrow=c(1,1))
  plot(x,pen[[1]], ylab="Penetrance", xlab="age at onset", ylim=c(0,1), type="l", lty=2, col="red", main=paste("Penetrance curves \n", base.dist, "baseline"), ...)
  lines(x, pen[[2]], lty=2, col="blue", ...)
  lines(x, pen[[3]], lty=1, col="red", ...)
  lines(x, pen[[4]], lty=1, col="blue", ...)
  legend("topleft", c("male carrier", "female carrier", "male noncarrier", "female noncarrier"), bty="n", lty=c(1,1,2,2), col=c("blue","red","blue","red"))
}

pen70=c((pen[[4]])[x==70], (pen[[3]])[x==70],(pen[[2]])[x==70],(pen[[1]])[x==70]) 
names(pen70)<-c("male-carrier","famale-carrier","male-noncarr","female-noncarr")

if(variation=="secondgene"){ 
  pen2=c((pen[[8]])[x==70], (pen[[7]])[x==70],(pen[[6]])[x==70],(pen[[5]])[x==70]) 
  pen70 = rbind(pen2, pen70)
  row.names(pen70)<-c("secondgene=1","secondgene=0")
  names(pen) <- c("secondgene=0: male-carrier","secondgene=0: famale-carrier","secondgene=0: male-noncarr","secondgene=0: female-noncarr",
                  "secondgene=1: male-carrier","secondgene=1: famale-carrier","secondgene=1: male-noncarr","secondgene=1: female-noncarr")
  }
else names(pen) <- c("male-carrier","famale-carrier","male-noncarr","female-noncarr")

if(print){
cat("Penetrance by age 70:\n")
print(pen70)
}
invisible(list(pen70=pen70, pen=pen, x.age=x))
}