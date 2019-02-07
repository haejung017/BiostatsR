summary.simfam <- function(object, digits = max(3, getOption("digits") - 3), ...){

  savedig <- options(digits = digits)
  on.exit(options(savedig))
  design <- attr(object, "design")
  variation <- attr(object, "variation") 
  frailty.dist <- attr(object, "frailty.dist") 
  base.dist <- attr(object, "base.dist") 
  
design.name <- c("pop: population-based study with affected probands",
"pop+: population-based study with affected and mutation carrier probands",
"cli: clinic-based study with affected probands",
"cli+: clinic-based study with affected and mutation carrier probands",
"twostage: two-stage design")[match(design, c("pop","pop+","cli","cli+","twostage"))]

cat("Study design:                          ", design,"\n")
cat("Baseline distribution:                 ", base.dist,"\n")
if(variation == "frailty") 
  cat("Frailty distribution:                  ", attr(object,"frailty.dist"),"\n")
else if(variation=="secondgene") 
  cat("Residucal familial correlation indueced by a second gene variation \n")

k <- ifelse(design=="twostage",2,1)
for(i in 1:k){
  if(design=="twostage"){
	 if(i==1){ 
	 	data <- object[object$weight==1,]
		cat("\n ## High risk families", "\n")
	 }
	 else{
	 	data <- object[object$weight!=1,]
		cat("\n ## Low risk families", "\n")
	 }
  }
else data<- object	 

numfam <- length(data$famID[data$proband==1])
avgnumaffec <- sum(data$status)/numfam
avgnumcarr <- sum(data$mgene,na.rm=T)/numfam
avgfamsize <- mean(data$fsize[data$proband==1])
avgageonset <- mean(data$ageonset[data$status==1])
wt <- unique(data$weight)

ans<-list(num.fam=numfam, avg.num.affected=avgnumaffec, avg.num.carriers=avgnumcarr, avg.family.size=avgfamsize,avgageonset=avgageonset)

if(i==1) ans.high <- ans
nn<-c(
"Number of families:                    ", 
"Average number of affected per family: ", 
"Average number of carriers per family: ", 
"Average family size:                   ",
"Average age of onset for affected:     ")

for(j in 1:5)  cat(nn[j], ans[[j]], "\n")

cat("Sampling weights used:                 ", unique(data$weight))
cat("\n")

if(i==2) ans<-cbind(HighRisk=ans.high, LowRisk=ans)

}
attr(ans,"design") <- design
attr(ans,"base.dist") <- base.dist
if(variation=="frailty") attr(ans,"frailty.dist") <- frailty.dist
else  attr(ans,"frailty.dist") <- NA
attr(ans, "variation") <- variation
class(ans) <- "summary.simfam"
invisible(ans)

}
