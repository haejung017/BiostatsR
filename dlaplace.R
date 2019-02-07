dlaplace <- function(dist, g, d, k){
if(dist=="gamma")  (-1)^d*factorial(k+d-1)*(1+g/k)^(-k)/factorial(k)/k^(d-1)
#else if(dist=="positives") exp(-g^k)
#else if(dist=="power") 1-(1-exp(-g))^k
else if(dist=="lognormal") gh(g,d, k)/sqrt(2*pi)
else stop("Unrecognized frailty distribution")
}