##### this R file is an example to do the interaction analysis on the test data.
##### the species abundance contain 17 taxa
#####  the environmental parameters have 14 parameters, 






###############################  
####         step 0      ######  
###############################  

##### load library, data, function 
#### the interaction calculation need the robust linear regression, therefore, we use the rlm in package MASS
#### both species abundance and environmental data should be data frame in which the rows stand for the samples, and column stand for the species names or environmental parameter names
#### there should be no NA values in the data, and no column which contain too many repeated values (for example, 0), therefore, we need to filter out the data  

library(MASS)  ##### this package include the rlm command 
library(robust)

 
##  load your species abudance  and environmental parameter data


load("taxa_candi.RData")  #### species abundance  part one

load("protist_myxobacteria_preytaxa_realabundances.rdata")  ####  species abundance  part two

Sample_Info <- read.table("Parameter_Nov2013.txt", header = TRUE, sep = "\t", dec = ",")  #### environmental parameters

  

###### load the functions ##### 
######  this file contains all the functions which are needed in the interaction analysis
  source("InteractionAnalysis_git.R")


### filter the original data>>=



###################################################################
#### combine the species abudance from the two part data
taxafull<- cbind(taxa_candi,taxa$Myxomycetes_total)
colnames(taxafull)[42]<-"Myxomycetes_total"



#################################################################
#### extrat the main environmental parametes to the data frame

rownames(Sample_Info)<- Sample_Info$Plot
envapara_tablefull <-Sample_Info

envpara_table <-Sample_Info[,mainpara]


############################################################
## filtering in the the environmental parameters, remove the NA value, remove the columes which contain too many repeated values

Full_Filterenvpara_table <- FilterNA(envpara_table,0)  ####  0 means use the method which can keep the most information when we do the filtering
Full_Filterenvpara_table <- Filterrepeat(Full_Filterenvpara_table,0.33) #### threshold 0.33 means one third of the values are the same value,

###########################################################

#####   the abundance  #######
## filtering in the the environmental parameters, remove the NA value, remove the columes which contain too many repeated values

RelaAbundOrig <- FilterNA(taxafull,0)


###############################
####  common samples after we did the filtering on both species and environmental parameters
chosensample <- intersect(rownames(RelaAbundOrig),rownames(Full_Filterenvpara_table))


##### determin the species abundance, and the environmental parameters data frame 
RelaAbundOrig <-RelaAbundOrig[chosensample,]
Full_Filterenvpara_table <- Full_Filterenvpara_table[chosensample,]






###############################  
####         step 1      ######  
###############################  
##### prepare the data frame for the direct interaction calculation  
#### both species abundance and environmental dataframe should have the same samples 


###### the main environmental parameters we need to consider in our test ####
mainpara <- c( "Inorganic_C", "Organic_C", "CN_ratio", "Fine_Roots_Biomass", "Coarse_Roots_Biomass", "roots_Total_C", "roots_Total_N", "roots_CN_ratio", "soil_moisture", "pH", "NH4", "NO3", "Nmin", "Cmic", "Cmic_Nmic") 


#### the species abundance
RelaAbund <- RelaAbundOrig   
Abundscale <- apply(RelaAbund, 2, function(x) scale(x, center = FALSE, scale = TRUE))  #### rescale the data
rownames(Abundscale) <- rownames(RelaAbund)

nSample<- nrow(RelaAbund)
nSpe <- ncol(RelaAbund)

#### the environmental parameters 
EnvPara <- Full_Filterenvpara_table[rownames(RelaAbund),mainpara[mainpara %in% colnames(Full_Filterenvpara_table)]]
Parascale <- apply(EnvPara, 2, function(x) scale(x, center = FALSE, scale = TRUE))  #### scale the data
rownames(Parascale) <- rownames(EnvPara)

nEnv <- ncol(EnvPara)



###############################  
####         step 2      ######  
###############################  
######  dest the singularity of the original data , and remove the columns which bring the singularity  
RelaAbund<-testsingula(RelaAbund)
Parascale<-testsingula(Parascale)




###############################  
####         step 3      ######  
###############################  
#### define the parameters which are needed in the functions

speciesprecision="high"
envprecision="high"
threshold01<-0.5  #### used to determin the patern in each environmental parameter
threshold<-0.2  ### error level used in the robust test

type<-3  #### use the work flow method to find the global interaction values
count<-10 #### used to calculate the peakvalue of the distribution 
factor02<-2  #### used to compare the difference between the positive and negative part
  
  


###############################  
####         step 4      ######  
###############################  


##### the workflow to calculate the interacion level values are discribed in our paper
### (1) calculate the rate of change of species abundance with respect to environmental parameter: use the function      rate_change, rate_change02,
### (2) calculate the interaction values for each species in the sense of environmental parameter: use the function      interinf, interinf02




##### calculate the interaction Matrix and correlation matrix ###########

### use the work flow method to calculate the global interaction matrix
### function interMatrix02 use the functions:   rate_change02, and interinf02, summ_inter_ij_new
InterMatrix_Summ_150Sample_high_spe <- interMatrix02 (Abundscale,Parascale,3,10,2,speciesprecision,envprecision,threshold01) 


corrMatrix_Summ_150Sample <- cor(Abundscale)


###############################  
####         step 5      ######  
############################### 
###### calculate the main parameters ####
out_150Sample<-main_para(Abundscale,Parascale,0,10,2,speciesprecision,envprecision,threshold01)
main_para_barplot(out_150Sample)

  
###############################  
####         step 6      ######  
###############################   
#### chose the six examples #####  
robust_test_rela<-six_test_robust(InterMatrix_Summ_150Sample_high_spe)




###############################  
####         step 7      ######  
############################### 

#### robust test on the six examples ###
for(k in 1:6){
  species_i<- robust_test_rela[k,1]
  species_j<- robust_test_rela[k,2]
  
    ##### remove one species #####
    permut_type<-"species"
  
    permut_spe<-interMatrix_ij_permut02 (Abundscale,Parascale,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
    
     nam1<-paste("permut_spe", rownames(robust_test_rela)[k],sep="_")  
    assign(nam1, permut_spe);
    
    ###### remove one environmental parameters ######
      permut_type<-"para"
  
    permut_para<-interMatrix_ij_permut02 (Abundscale,Parascale,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
    nam1<-paste("permut_para", rownames(robust_test_rela)[k],sep="_")  
    assign(nam1, permut_para);
  
  
    ###### replace of samples #####
     permut_type<-"samples"
    times<-500
    size<-15
    permut_sample<- interMatrix_ij_permut(Y,X,species_i,species_j,times,size,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
    
    nam1<-paste("permut_sample", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
    assign(nam1, permut_sample);
  
  
  ###### error robust on para ####
  times<- 500
  robust_para<- interMatrix_ij_env_robust(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_para", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_para);
  
  
  ###### error robus on species ####
  
  times<- 500
  robust_spe<- interMatrix_ij_spe_robust(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_spe", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_spe);
  
  
  
}

