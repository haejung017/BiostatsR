library(Rcpp)
library(pracma)
library(ggplot2)
sourceCpp("cumh_ED.cpp")
source("penplot_funcs.R")
library(gridExtra)
# Data
b1 <- read.csv("BRCA1_comp_final.csv") #BRCA1
#b2 <- read.csv("BRCA2_comp_final.csv") #BRCA2

# modelling result BRCA1
est1 <-readRDS("BRCA1est") # parms estimates
est1var <-readRDS("BRCA1est_robustvar") # parms robust variance 



### Example plots

# eventage = age of TVC occuring, 
# eventind = indicator of TVC, 
          # 1 = first screen
          # 2 = second screen
          # 3 = third screen
          # 4 = BO (Bilateral Oophorectomy)
          # 5 = BM (Bilateral Mastectomy)
# type = 1 (breast), 2 (ovarian)
# res = resolution of the plot

# no screen no BO
PlotPen(est1,est1var,eventage=c(NA,NA,NA,NA,NA),eventind=c(NA,NA,NA,NA,NA),type=1, res=100)
# first screen at age 40
PlotPen(est1,est1var,c(40,NA,NA,NA,NA),c(1,NA,NA,NA,NA),type=1, res=100)
# BO at age 40
PlotPen(est1,est1var,c(40,NA,NA,NA,NA),c(4,NA,NA,NA,NA),type=1, res=100)
# first screen at age 35, second screen at age 40, third screen at age 45
PlotPen(est1,est1var,c(35,40,45,NA,NA),c(1,2,3,NA,NA),type=1, res=100)
# first screen at age 35, BO at age 40, second screen at age 45
PlotPen(est1,est1var,c(35,40,45,NA,NA),c(1,4,2,NA,NA),type=1, res=100)
# BM bilateral mastectomy at age 55
PlotPen(est1,est1var,c(55,NA,NA,NA,NA),c(5,NA,NA,NA,NA),type=1, res=100)
# Screen1 -> BO -> screen2 -> BM
a=PlotPen(est1,est1var,c(35,40,45,50,NA),c(1,4,2,5,NA),type=1, res=300)
b=PlotPen(est1,est1var,c(35,40,45,NA,NA),c(1,4,2,NA,NA),type=1, res=300)


BC=PlotPen(est1,est1var,c(NA,NA,NA,NA,NA),c(NA,NA,NA,NA,NA),type=1, res=100)
OC=PlotPen(est1,est1var,c(NA,NA,NA,NA,NA),c(NA,NA,NA,NA,NA),type=2, res=100)
BCnotitle=BC+ggtitle(NULL)
OCnotitle=OC+ggtitle(NULL)

BCs1=PlotPen(est1,est1var,c(35,NA,NA,NA,NA),c(1,NA,NA,NA,NA),type=1, res=100)
BCs2=PlotPen(est1,est1var,c(40,NA,NA,NA,NA),c(1,NA,NA,NA,NA),type=1, res=100)
BCs3=PlotPen(est1,est1var,c(45,NA,NA,NA,NA),c(1,NA,NA,NA,NA),type=1, res=100)
BCs1=BCs1+ggtitle(NULL)
BCs2=BCs2+ggtitle(NULL)
BCs3=BCs3+ggtitle(NULL)

BCo1=PlotPen(est1,est1var,c(40,NA,NA,NA,NA),c(4,NA,NA,NA,NA),type=1, res=100)
BCo2=PlotPen(est1,est1var,c(45,NA,NA,NA,NA),c(4,NA,NA,NA,NA),type=1, res=100)
BCo3=PlotPen(est1,est1var,c(50,NA,NA,NA,NA),c(4,NA,NA,NA,NA),type=1, res=100)
BCo1=BCo1+ggtitle(NULL)
BCo2=BCo2+ggtitle(NULL)
BCo3=BCo3+ggtitle(NULL)

BCs31=PlotPen(est1,est1var,c(30,35,40,NA,NA),c(1,2,3,NA,NA),type=1, res=100)
BCs32=PlotPen(est1,est1var,c(35,40,45,NA,NA),c(1,2,3,NA,NA),type=1, res=100)

BCs31=BCs31+ggtitle(NULL)
BCs32=BCs32+ggtitle(NULL)

a=PlotPen(est1,est1var,c(30,35,40,45,NA),c(1,2,3,4,NA),type=1, res=100)
b=PlotPen(est1,est1var,c(30,35,40,45,NA),c(1,4,2,3,NA),type=1, res=100)

a=a+ggtitle(NULL)
b=b+ggtitle(NULL)
# OCs1=PlotPen(estimateobj2,c(35,NA,NA,NA,NA),c(1,NA,NA,NA,NA),type=2, res=100)
# 
# PlotPen(estimateobj2,c(35,40,NA,NA,NA),c(1,2,NA,NA,NA),type=1, res=300)
# PlotPen(estimateobj1,c(35,40,NA,NA,NA),c(1,2,NA,NA,NA),type=1, res=100)
# 
# PlotPen(estimateobj2,c(40,NA,NA,NA,NA),c(1,NA,NA,NA,NA),type=2, res=100)

#NONscreened BC,OC
jpeg("BRCA1_BC_OC.jpg", width=7, height=5, units="in", res=500)
grid.arrange(BC,OC,nrow=1)
dev.off()

jpeg("BRCA1_BC_s1.jpg", width=10, height=5, units="in", res=500)
grid.arrange(BCnotitle,BCs1,BCs2,BCs3,nrow=1)
dev.off()

jpeg("BRCA1_BC_o.jpg", width=10, height=5, units="in", res=500)
grid.arrange(BCnotitle,BCo1,BCo2,BCo3,nrow=1)
dev.off()

jpeg("BRCA1_BC_s3.jpg", width=10, height=5, units="in", res=500)
grid.arrange(BCnotitle,BCs31,BCs32,nrow=1)
dev.off()

jpeg("BRCA1.jpg", width=10, height=15, units="in", res=500)

pdf("StackedPlot.pdf",width = 7, height = 20,pointsize=10)
grid.arrange(BCnotitle,OCnotitle,BCs1,BCs3,BCs31,BCs32,a,b,nrow=4,ncol=2)
dev.off()

jpeg("BRCA1_BM.jpg", width=10, height=5, units="in", res=500)
grid.arrange(a,b,ncol=2)
dev.off()
