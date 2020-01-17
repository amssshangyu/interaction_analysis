##### this R file is an example to do the interaction analysis on the test data.
##### Version 1.0
##### Copyright: Dr. Shang, Yu,  
##### Leibniz-Institut DSMZ-Deutsche Sammlung von Mikroorganismen und Zellkulturen GmbH
##### Inhoffenstraße 7 B
##### 38124 Braunschweig
##### Germany
##### Tel. ++49-531-2616-111
##### Fax ++49-531-2616-418
##### www.dsmz.de
##### email: yus14@dsmz.de

#######################################################################################################
#######################################################################################################

##### load library, data, function 

#######################################################################################################
#######################################################################################################


#### the interaction calculation need the robust linear regression, therefore, we use the rlm in package MASS
#### both species abundance and environmental data should be data frame in which the rows stand for the samples and column stand for the species names or environmental parameter names
#### there should be no NA values in the data, and no column which contains too many repeated values (for example, 0), therefore, we need to filter out the data  

library(MASS)  ##### this package include the rlm command 
library(robust)
library(reshape)
library(ggplot2)

setwd("d:/user/yus14/Desktop/R_SweaveWorks/Microbial Community Interaction/git")





#------------------------------------------------------------------------------------------------------------


#######################################################################################################
#######################################################################################################

################    PartI,   the code used for our original data ananlysis    #########################

#######################################################################################################
#######################################################################################################
###############################  
####         step 1     ######  
###############################  


 
##  load your species abudance  and environmental parameter data
##### the species abundance contain 17 taxa
#####  the environmental parameters have 14 parameters, 


load("Abundscale.RData")  #### species abundance  

load("Parascale.rdata")  ####  environmental parameters


###### load the functions ##### 
######  this file contains all the functions which are needed in the interaction analysis
source("InteractionAnalysis_git.R")


###############################  
####         step 2      ######  
###############################  
######  dest the singularity of the original data , and remove the columns which bring the singularity  
RelaAbund<-testsingula(Abundscale)
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


##### the workflow to calculate the interaction level values are described in our paper
### (1) calculate the rate of change of species abundance with respect to the environmental parameter: use the function      rate_change, rate_change02,
### (2) calculate the interaction values for each species in the sense of environmental parameter: use the function      interinf, interinf02

### !!!! I suggest using rate_change02, interinf02 to do the analysis.  
### !!!! This is because rate_change, interinf construct the linear regreesin modle by taking into account
### !!!! both the environmental parameter and all the species. 
### !!!! Although this is the final full version of the model, it is too complicate for rlm to do the analysis, 
### !!!! it can bring the default error from rlm calculation.
### !!!! rate_change02, interinf02 only use the environmental parameter in the linear regression to simplify the calculation.



##### calculate the interaction Matrix and correlation matrix ###########

### use the work flow method to calculate the global interaction matrix
### the function interMatrix02 use the functions:   rate_change02, and interinf02, summ_inter_ij_new
### this step take long time because the data contain 17 species, 13 parameters, 150 samples, 
### we also provide a smaller dataset which need much short computation time for reader, please see the following parts

InterMatrix_Summ_150Sample_high_spe <- interMatrix02 (RelaAbund,Parascale,3,10,2,speciesprecision,envprecision,threshold01) 


corrMatrix_Summ_150Sample <- cor(Abundscale)



  
###############################  
####         step 5      ######  
###############################   
#### chose the six examples for the robust test::  strong positive, median positive, low positive,####
####  strong negative, median negative, low negative.  #####  
#### the results (row and column) are the index of the two species, 

robust_test_rela<-six_test_robust(InterMatrix_Summ_150Sample_high_spe)  



###############################  
####         step 6      ######  
############################### 

#### robust test on the six examples ###

#### k indicates the index of the six examples
for(k in 1:6){
  species_i<- robust_test_rela[k,1] ##### species i
  species_j<- robust_test_rela[k,2] ##### species j
  
    # ##### remove one species #####
    # permut_type<-"species"
    # 
    # permut_spe<-interMatrix_ij_permut02 (Abundscale,Parascale,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
    # 
    #  nam1<-paste("permut_spe", rownames(robust_test_rela)[k],sep="_")  
    # assign(nam1, permut_spe);
    # 
    # ###### remove one environmental parameters ######
    #   permut_type<-"para"
    # 
    # permut_para<-interMatrix_ij_permut02 (Abundscale,Parascale,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
    # nam1<-paste("permut_para", rownames(robust_test_rela)[k],sep="_")  
    # assign(nam1, permut_para);
    # 
  
    ###### remove  samples #####
     permut_type<-"samples"
    #times<-500
    times<-5 ##### this is only for test of the code
    size<-15
    permut_sample<- interMatrix_ij_permut(RelaAbund,Parascale,species_i,species_j,times,size,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
    
    nam1<-paste("permut_sample", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
    assign(nam1, permut_sample);
  
  
  ###### error robust on para ####
    #times<-500
    times<-5 ##### this is only for test of the code
  robust_para<- interMatrix_ij_env_robust(RelaAbund,Parascale,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_para", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_para);
  
  
  ###### error robus on species ####
  
  #times<-500
  times<-5 ##### this is only for test of the code
  robust_spe<- interMatrix_ij_spe_robust(RelaAbund,Parascale,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_spe", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_spe);
 

}


#######################################################################################################
#######################################################################################################




#-----------------------------------------------------------------------------------------------------------------




#######################################################################################################
#######################################################################################################

################    PartII,   the code used for simplified data               #########################

#######################################################################################################
#######################################################################################################


###############################  
####         step 1     ######  
###############################  



##  load your species abudance  and environmental parameter data
##### the species abundance contain 17 taxa
#####  the environmental parameters have 14 parameters, 


load("Abundscale_short.RData")  #### species abundance  

load("Parascale_short.rdata")  ####  environmental parameters


###### load the functions ##### 
######  this file contains all the functions which are needed in the interaction analysis
source("InteractionAnalysis_git.R")


###############################  
####         step 2      ######  
###############################  
######  dest the singularity of the original data , and remove the columns which bring the singularity  
RelaAbund<-testsingula(Abundscale_short)
Parascale<-testsingula(Parascale_short)




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


##### the workflow to calculate the interaction level values are described in our paper
### (1) calculate the rate of change of species abundance with respect to the environmental parameter: use the function      rate_change, rate_change02,
### (2) calculate the interaction values for each species in the sense of environmental parameter: use the function      interinf, interinf02



### !!!! I suggest using rate_change02, interinf02 to do the analysis.  
### !!!! This is because rate_change, interinf construct the linear regreesin modle by taking into account
### !!!! both the environmental parameter and all the species. 
### !!!! Although this is the final full version of the model, it is too complicate for rlm to do the analysis, 
### !!!! it can bring the default error from rlm calculation.
### !!!! rate_change02, interinf02 only use the environmental parameter in the linear regression to simplify the calculation.





##### calculate the interaction Matrix and correlation matrix ###########

### use the work flow method to calculate the global interaction matrix
### the function interMatrix02 use the functions:   rate_change02, and interinf02, summ_inter_ij_new
### this step take long time because the data contain 17 species, 13 parameters, 150 samples, 
### we also provide a smaller dataset which need much short computation time for reader, please see the following parts

InterMatrix_Summ_high_spe <- interMatrix02 (RelaAbund,Parascale,type,count,factor02,speciesprecision,envprecision,threshold01) 


corrMatrix_Summ <- cor(RelaAbund)




###############################  
####         step 5      ######  
###############################   
#### chose the six examples for the robust test::  strong positive, median positive, low positive,####
####  strong negative, median negative, low negative.  #####  
#### the results (row and column) are the index of the two species, 

robust_test_rela<-six_test_robust(InterMatrix_Summ_high_spe)  



###############################  
####         step 6      ######  
############################### 

#### robust test on the six examples ###

#### k indicates the index of the six examples
for(k in 1:6){
  species_i<- robust_test_rela[k,1] ##### species i
  species_j<- robust_test_rela[k,2] ##### species j
  
  # ##### remove one species #####
  # permut_type<-"species"
  # 
  # permut_spe<-interMatrix_ij_permut02 (Abundscale,Parascale,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
  # 
  # nam1<-paste("permut_spe", rownames(robust_test_rela)[k],sep="_")  
  # assign(nam1, permut_spe);
  # 
  # ###### remove one environmental parameters ######
  # permut_type<-"para"
  # 
  # permut_para<-interMatrix_ij_permut02 (Abundscale,Parascale,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
  # nam1<-paste("permut_para", rownames(robust_test_rela)[k],sep="_")  
  # assign(nam1, permut_para);
  # 
  
  ###### remove  samples #####
  permut_type<-"samples"
  #times<-500
  times<-5 ##### this is only for test of the code
  size<-15
  permut_sample<- interMatrix_ij_permut(RelaAbund,Parascale,species_i,species_j,times,size,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01)
  
  nam1<-paste("permut_sample", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, permut_sample);
  
  
  ###### error robust on para ####
  #times<-500
  times<-5 ##### this is only for test of the code
  robust_para<- interMatrix_ij_env_robust(RelaAbund,Parascale,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_para", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_para);
  
  
  ###### error robus on species ####
  
  #times<-500
  times<-5 ##### this is only for test of the code
  robust_spe<- interMatrix_ij_spe_robust(RelaAbund,Parascale,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01)
  nam1<-paste("robust_spe", rownames(robust_test_rela)[k],"err",threshold,sep="_")  
  assign(nam1, robust_spe);
  
  
}


#######################################################################################################
#######################################################################################################






#------------------------------------------------------------------------------------------------------------







#######################################################################################################
#######################################################################################################

################                  PartII,   Additional tools                  #########################

#######################################################################################################
#######################################################################################################



##############################################################################
##### The follwing code is not used in our manuscript, 
##### we provided here for the possible interest for the readers
###############################################################################


###### figure out the dominant parameters: Some environmental parameters play  the dominant role in determining the interaction level
###### 
###### We calculate (1) the frequency of parameter determining the interaction value; 
###### (2) frequency of  parameter that has the lowest negative value; 
###### (3) frequency of parameters that has the largest positive values.  

out_Sample<-main_para(RelaAbund,Parascale,0,10,2,speciesprecision,envprecision,threshold01)
main_para_barplot(out_Sample)
