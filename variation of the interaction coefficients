#### this is a example to generate the Fig2 in my paper.

library(MASS)  ##### this package include the rlm command 
library(robust)
library(reshape)
library(ggplot2)


### load data
setwd("d:/user/yus14/Desktop/R_SweaveWorks/Microbial Community Interaction/git")
load("Abundscale.RData")
load("Parascale.RData")
source("InteractionAnalysis_git.R")


###  calculation of the interaction influnce for  species i
###  take i=1 as example
###

i<-1
derEnv <- rate_change02 (Abundscale, Parascale, i)
derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]


interaction<- interinf02 (Abundscale, Parascale[,colnames(derEnv)], derEnv, i)      

###  interinf_ij is the value    beta_ij\alpha ^ k


###  take j=2 as example, extract the interaction coefficients  beta_12\alpha ^ k
###   the interaction effect from species j=2 on species i=1

j<-2


interinf_ij <- interaction[j,,]



###  chose three parameter as example to show the variation of beta_12\alpha ^ k
###  organic_C, pH,soil_moisture




### In order to make the variation plot, you need to use the original environmental parameter
### load the original environmental parameters,
### here, you need to load your own original paramter


Sample_Info <- read.table("Parameter_Nov2013.txt", header = TRUE, sep = "\t", dec = ",")

EnvPara <- Sample_Info[,colnames(Parascale)]



### extract the interaction coefficients for these three environmental paramters
interinf_organic_C <- as.data.frame(cbind(interinf_ij[2,],EnvPara[1:148,2]))

colnames(interinf_organic_C)<- c("Y","X")

interinf_soil_moisture <- as.data.frame(cbind(interinf_ij[8,],EnvPara[1:148,8]))

colnames(interinf_soil_moisture)<- c("Y","X")

interinf_pH <- as.data.frame(cbind(interinf_ij[9,],EnvPara[1:148,9]))

colnames(interinf_pH)<- c("Y","X")




#### now, you can make the plot for the variation of interaction coefficients 
### example:  the interaction effect from species 2 on species 1, long the pH gradients.


plot(interinf_pH$X,interinf_pH$Y)
