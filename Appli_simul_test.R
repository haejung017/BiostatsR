setwd("~/Dropbox/2015/Daniel/Simulation/SimulationDaniel/")
library(truncnorm)
source("simul_functions_cmp.R")
source("familyStructure_cmp.R")

# below files are from FamEvent R package
source("parents.g.R")
source("Pgene.R")
source("kids.g.R")
source("hazards.R")
source("cumhaz.R")
source("familydesign.r")

options(width=1000)

iest.cmp <- c(-5.472670754, 0.875468737, -4.688551795, 1.071583616, 0.41, 2.86, -0.72, 1.28, 1.572773928) 

# simulating families
source("bootstraps_functions-rev.R")

set.seed(1)

nfam = 1000 # n=500, 779, 1000

for(logk in log(c(1,2,5,10))){
  cat("\nk =", exp(logk), "\n")
  nsim <- 500
  outnc <- matrix(0,ncol=8,nrow=nsim)
  outnf <- matrix(0,ncol=13,nrow=nsim)
  outc1 <- outc3 <- outc2 <- outc3i <- outc2i <- matrix(0,ncol=14,nrow=nsim) # two-stage
  tval <- c(param.cmp(iest.cmp, k=exp(logk)),0)

    i = 0
while(i < nsim){
  
  simdat<- simfam.cmp(n=nfam, c(iest.cmp[1:8], -logk) )

#  out<-fit.Iterative(simdat, iest.cmp)
#  est.cmp <- optim(iest.cmp, loglikcmp, data=simdat, agemin=14,hessian=hessian,control=list(maxit=4000))
#  out <- fit.Iterative(simdat, iest.cmp)
#  out
  
  
  out<-try(fit7models(simdat, iest.cmp))
  
  if(is.null(attr(out,"class"))){
    i <- i + 1
    outc1[i,]    <- c(out$est.cmp, out$nloglik.cmp)
    outc2[i,]    <- c(out$est.cmp2, out$nloglik.cmp2)
    outc3[i,]    <- c(out$est.cmp3, out$nloglik.cmp3)
    outnf[i,]    <- c(out$est.nofrail, out$nloglik.nofrail)
    outnc[i,]    <- c(out$est.nocmp, out$nloglik.nocmp)
    cat(i," ")
  }
}


    write.table(outc1,  paste0("sim",nfam,"outc1_k=", exp(logk),".txt"))
    write.table(outc2,  paste0("sim",nfam,"outc2_k=", exp(logk),".txt"))
    write.table(outc3,  paste0("sim",nfam,"outc3_k=", exp(logk),".txt"))
    write.table(outnf,  paste0("sim",nfam,"outnf_k=", exp(logk),".txt"))
    write.table(outnc,  paste0("sim",nfam,"outnc_k=", exp(logk),".txt"))
    
    cat("\n True value: ", tval, "\n")
    cat("\n Family numbers: ", n, "\n")
    cat("\n Mean BIAS \n")
    print(apply(outc1, 2, mean) - tval)
    print(apply(outc2, 2, mean) - tval)
    print(apply(outc3, 2, mean) - tval)
    print(apply(outnf, 2, mean) - tval[-13])
    print(apply(outnc, 2, mean) - tval[c(1,2,5,6,9,10,13)])

    cat("\n Median BIAS \n")
    print(apply(outc1, 2, median) - tval)
    print(apply(outc2, 2, median) - tval)
    print(apply(outc3, 2, median) - tval)
    print(apply(outnf, 2, median) - tval[-13])
    print(apply(outnc, 2, median) - tval[c(1,2,5,6,9,10,13)])
    
    cat("\n Empricial SE \n")
    print(apply(outc1, 2, sd))
    print(apply(outc2, 2, sd))
    print(apply(outc3, 2, sd))
    print(apply(outnf, 2, sd))
    print(apply(outnc, 2, sd))
    
}





