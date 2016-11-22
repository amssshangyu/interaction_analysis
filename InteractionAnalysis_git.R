###### this R file contains all the functions which are needed in the interaction analysis


#####  part A --  filtering the data:  FilterNA, FilterZero, Filterrepeat
#####  part B --  calculate derivative and original interaction values  :  rate_change, rate_change02, rate_change_low, interinf, interinf02, interinf03, interinf_low
#####  part C --  calculate the interaction matirx:  interMatrix, interMatrix02, interMatrix_norm, interMatrix_norm_summ, interMatrix_low, interMatrix_ij, interMatrix_ij_low
#####  part D --  statistic summary:  summ_inter_ij, summ_inter_ij_new, summ_interaction, anascatter
#####  part E --  main environmental parameters:  main_para, main_para_barplot
#####  part F --  robust test:  interMatrix_ij_permut, interMatrix_ij_env_spe_robust, interMatrix_ij_env_robust, interMatrix_ij_spe_robust, interMatrix_ij_permut02  
#####  part G --  test singularity: testsingula
#####  part H --  extract six examples from the interaction matrix: six_test_robust 
#####  part I --  make plot of change of interaciton value across the samples:  reshapedata, interinfPlot




####################################################################
######                       part A                            ####
###### filtering the data with the minimum missing information##### 
###################################################################
  FilterNA <- function (data,type){
    #### data is data frame of species abundance or environmental parameter 
    #### type = 0,  remove the NA, with the minimum missing information
    #### type = 1, remove the rows which contain NA
    #### type = other value, remove the columns which contain NA
    
    if (type==0){
      
      while (anyNA(data)){
        
        NAcol_env<- apply(data, 2, function(x) length(which(is.na(x))))
        NArow_env<- apply(data, 1, function(x) length(which(is.na(x))))
        
        max_NAcol_env <- max(NAcol_env,na.rm=TRUE)
        max_NArow_env <-max(NArow_env,na.rm=TRUE)
        
        if (max_NAcol_env>=max_NArow_env){
          data <-data[,-which(NAcol_env==max_NAcol_env)]}  else{
            data <-data[-which(NArow_env==max_NArow_env),]
          }
        
      }
    } else if (type==1){
      NArow<- apply(data, 1, function(x) length(which(is.na(x))))
      
      data <- data[which(NArow==0),] 
      
    }else{
      NAcol<- apply(data, 2, function(x) length(which(is.na(x))))
      data <- data[which(NAcol==0),]  
    }
    
    return(data)
  }


FilterZero <- function (data,type){
  #### data is data frame of species abundance or environmental parameter 
  #### type = 0,  remove the 0, with the minimum missing information
  #### type = 1, remove the rows which contain 0
  #### type = other value, remove the columns which contain 0
  
  if (type==0){
    
    while (any(data==0.0)){
      
      Zerocol_env<- apply(data, 2, function(x) length(which(x==0.0)))
      Zerorow_env<- apply(data, 1, function(x) length(which(x==0.0)))
      
      max_Zerocol_env <- max(Zerocol_env,na.rm=TRUE)
      max_Zerorow_env <-max(Zerorow_env,na.rm=TRUE)
      
      if (max_Zerocol_env>=max_Zerorow_env){
        data <-data[,-which(Zerocol_env==max_Zerocol_env)]}  else{
          data <-data[-which(Zerorow_env==max_Zerorow_env),]
        }
      
    }
  } else if (type==1){
    Zerorow<- apply(data, 1, function(x) length(which(x==0.0)))
    
    data <- data[which(Zerorow==0.0),] 
    
  }else{
    Zerocol<- apply(data, 2, function(x) length(which(x==0.0)))
    data <- data[which(Zerocol==0.0),]  
  }
  
  return(data)
}




Filterrepeat <- function(data,threshold){
  #### remove the columns which contain too many repeat values
  N <- floor(threshold*nrow(data))  #### the threshold for the munber of repeat values
  reptimes <- vector()
  for (j in 1: ncol(data)){
    reptimes <- c(reptimes, max(table(data[,j])))
  }
  
  data <- data[,which(reptimes< N)]
  return(data)
}


  
  
  
#################################################################
######                       part B                          ####
#### calculate the first order derivative and the interaction ####
#################################################################  
  
#### this  part contains all the functions which are used to calculate the rate of change of species abundance 
#### and the original 4 dimensional interaction values. The calculation of the rate of change needs the species abundance and environmental parameters data.
#### the interaction need additionally the results of calculation of the rate of change as the input data.
  
  rate_change <- function(Y, X, i)
  {
    #### This function will calculate the first derivative of species i abundance with respect to
    #### all the environmental parameters in all the samples. Depending on the data structure, this function
    #### will automatically chose the different specision method. The returned results are the 2 dimensional 
    #### dataframe: row is the samples, column is the environmental parameter
    
    
    ####   Y is the abundance, data frame 
    ####   X is the environmental parameters, data frame
    ####   i is the index of species i
    
    nsample <- nrow(Y)
    nSpe <-ncol(Y)
    nEnv <- ncol(X)
    derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
    
    for (n in 1:(nsample -1)){
      #### change of  abundance #####
      
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
      
      deltaSpe_i <- deltaSpe[,i]
      
      ####  change of environmental parameters ###
      
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
      
      if (nsample>=(nSpe+nEnv)){
        regdata <- data.frame(cbind( deltaEnv,deltaSpe[,-i]))
        ###### linear regression on   deltaSpe ~ deltaEnv+deltaSpe,  the influence coming from both 
        ### environmental parameters and species abundance are considered
        
        fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
        fit <- lm (fmal,  regdata)    
        derEnv[n,]<- fit$coefficients[1:nEnv]
      } else if (nsample>=nEnv){
        regdata <- data.frame( deltaEnv)
        ###### linear regression on   deltaSpe ~ deltaEnv, only the influence from 
        ####   environmental parameter is considered
        fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
        fit <- rlm (fmal,  regdata)  
        derEnv[n,]<- fit$coefficients[1:nEnv]
      }else{
        #### low precision, using two point to calculate the derivative, then choose the median 
        temp<-apply(deltaEnv,2,function(x) deltaSpe_i/x)
        temp[which(!is.finite(temp))] <- NA
        
        derEnv[n,] <- apply(temp,2,function(x) mean(x,na.rm=TRUE))
        
      } 
    }
    derEnv <- data.frame(derEnv)
    rownames(derEnv) <-rownames(Y)[1:(nsample-1)]
    colnames(derEnv) <-colnames(X)
    return(derEnv) 
  }


rate_change02 <- function(Y, X, i)
{ #### This function will calculate the first derivative of species i abundance with respect to
  #### all the environmental parameters in all the  samples. Depending on the data structure, this function
  #### only consider the influence of environmental parameters. The returned results are the 2 dimensional 
  #### data frame: row is the samples, column is the environmental parameter
  
  
  
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   i is the species i
  
  #####   this function will calculate the first derivative of species i abundace with respect to all the environmental parameters in the sample n  ######
  nsample <- nrow(Y)
  nSpe <-ncol(Y)
  nEnv <- ncol(X)
  derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
  
  #for (n in 79:80){
  for (n in 1:(nsample -1)){
    #### change of  abundance #####
    #assign("last.warning", NULL, envir = baseenv())
    
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
    
    deltaSpe_i <- deltaSpe[,i]
    
    ####  change of environmental parameters ###
    
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
    

      regdata <- data.frame( deltaEnv)
    
    ##### use robust linear regression  #### 
    #######  lm is not robust at all, the rlm need to test the sigularity of the data 
      ###### linear regression on   deltaSpe ~ deltaEnv
      fmal<- as.formula (paste("deltaSpe_i ~ -1 + ",paste(colnames(regdata),collapse="+")))
     
      fit <- rlm (fmal, regdata, maxit = 200) 
      ####  200 is  the running steps to have a convergence estimation, it can be changed based on the user´s test and experience
      
  
      derEnv[n,]<- fit$coefficients[1:nEnv]
    
  }
  derEnv <- data.frame(derEnv)
  rownames(derEnv) <-rownames(Y)[1:(nsample-1)]
  colnames(derEnv) <-colnames(X)
  return(derEnv) 
}



rate_change_low <- function(Y, X, i)
{ #### This function will calculate the first derivative of species i abundance with respect to
  #### all the environmental parameters in all the  samples. Depending on the data structure, this function
  #### use the low precision method. The returned results are the 2 dimensional 
  #### data frame: row is the samples, column is the environmental parameter
  
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   i is the species i
  
  #####   this function will calculate the first derivative of species i abundace with respect to all the environmental parameters in the sample n  ######
  nsample <- nrow(Y)
  nSpe <-ncol(Y)
  nEnv <- ncol(X)
  derEnv <- matrix(, nrow = nsample-1, ncol = nEnv)
  
  for (n in 1:(nsample -1)){
    #### change of  abundance #####
    
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,]))))
    
    deltaSpe_i <- deltaSpe[,i]
    
    ####  change of environmental parameters ###
    
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,]),1,function(x) x-as.matrix(X[n,]))))
    
   
      #### low precision, using two point to calculate the derivative, then choose the median 
      temp<-apply(deltaEnv,2,function(x) deltaSpe_i/x)
      temp[which(!is.finite(temp))] <- NA
      
      derEnv[n,] <- apply(temp,2,function(x) mean(x,na.rm=TRUE))
      
  }
  derEnv <- data.frame(derEnv)
  rownames(derEnv) <-rownames(Y)[1:(nsample-1)]
  colnames(derEnv) <-colnames(X)
  return(derEnv) 
}




interinf <- function(Y, X, DY, i){
  #### This function will calculate the interaction level of species i abundance 
  #### in all the  samples. Depending on the data structure, this function
  #### will automatically choose the different specific methods. The returned results are the 3-dimensional array
  
  
  
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   DY is the derivative of Y with respect to all the environmental parameters, the rate of change of species abundance
  ####   i is the species i
  
  
  #####   this function will calculate the interaction influence  on species i   #####
  
  
  nsample <- nrow(DY)
  nEnv <- ncol(DY)
  nSpe <- ncol(Y)
  
  if(nsample >=(nSpe+nEnv) ){
    
    #### the model consider the influence from both the species abundacne and environmental parameters,use the Taylor expansion
    
    
    #### array: array[1] is the environmental parameter and the species, 
    ####        array[2] is the environmental parameters, 
    ####        array[3] is the samples
    interaction <- array (0, dim=c(nEnv+nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y),colnames(DY)),colnames(DY),rownames(Y[1:(nsample-1),])))
    # Spe_i <- Y[,i]    ### abundance of species i
    
    for (n in 1: (nsample-1)){
      
      ##### change of the DY ####
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      
      #### change of  abundance #####
      
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      
      ####  change of environmental parameters ###
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      
      regdata <- data.frame(cbind(deltaSpe, deltaEnv))
      ###### linear regression on   deltaSpe ~ deltaEnv
      fmal<- as.formula (paste("deltaDY ~ -1 + ",paste(colnames(regdata),collapse="+")))
      
      fit <- rlm (fmal,  regdata[1:(nsample-1),])
      
      interaction [,,n] <- fit$coefficients/Y[n,i]
      
      if(Y[n,i]==0){
        interaction [,,n]<- NA
      } 
      
      missspecies <- which(Y[n,]==0,useNames = TRUE)
      interaction [missspecies,,n] <-NA
    }
    
  }else if (nsample >= nSpe){
    
    #### the model consider the influence from obly the species abundacne, use the Taylor expansion
    
    
    #### array: array[1] is the species, 
    ####        array[2] is the environmental parameters, 
    ####        array[3] is the samples
    
    interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
    
    for (n in 1: (nsample-1)){
      
      ##### change of the DY ####
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      
      #### change of  abundance #####
      
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      ####  change of environmental parameters ###
      
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      
      regdata <- data.frame(cbind(deltaSpe))
      ###### linear regression on   deltaSpe ~ deltaEnv
      fmal<- as.formula (paste("deltaDY ~ -1 + ",paste(colnames(regdata),collapse="+")))
      
      fit <- lm (fmal,  regdata[1:(nsample-1),])
      
      interaction [,,n] <- fit$coefficients/Y[n,i]
      
      if(Y[n,i]==0){
        interaction [,,n]<- NA
      } 
      
      missspecies <- which(Y[n,]==0,useNames = TRUE)
      interaction [missspecies,,n] <-NA
    }
  }else{
    #### the model use the low precision method
    
    
    #### array: array[1] is the species, 
    ####        array[2] is the environmental parameters, 
    ####        array[3] is the samples
    
    interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
    
    for (n in 1: (nsample-1)){
      
      ##### change of the DY ####
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      
      #### change of  abundance #####
      
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      ####  change of environmental parameters ###
      
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      
      
      temp<-apply(deltaSpe[1:(nsample-1),],c(1,2),function(x) 1/(x*nsample))
      temp[which(!is.finite(temp))] <- 0
      interaction [,,n] <-t(temp) %*% deltaDY    ###### this is equivalent to calculate the mean values
      
      if(Y[n,i]==0){
        interaction [,,n]<- NA
      } 
      
      missspecies <- which(Y[n,]==0,useNames = TRUE)
      interaction [missspecies,,n] <-NA  
    }
  }
  return(interaction)
}






interinf02 <- function(Y, X, DY, i){
  #### This function will calculate the interaction level of species i abundance 
  #### in all the  samples. Depending on the data structure, 
  #### The returned results are the 3 dimensional array
  #### array: array[1] is the species, 
  ####        array[2] is the environmental parameters, 
  ####        array[3] is the samples
  
  #### the model consider the influence from only the species abundance, use the Taylor expansion
  
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   DY is the derivative of Y with respect to all the environmental parameters
  ####   i is the species i
  
    
  
  nsample <- nrow(DY)
  nEnv <- ncol(DY)
  nSpe <- ncol(Y)
  

    
    interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
    
    for (n in 1: (nsample-1)){
      
      ##### change of the DY ####
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      ####  change of environmental parameters ###
      
      #deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      
      regdata <- data.frame(cbind(deltaSpe))
      for(k in 1:ncol(deltaDY)){
        deltaDY_k <- deltaDY[,k]
        
        
        ###### linear regression on   deltaSpe ~ deltaEnv
        fmal<- as.formula (paste("deltaDY_k ~ -1 + ",paste(colnames(regdata),collapse="+")))
        
        fit <- rlm (fmal,  regdata[1:(nsample-1),], maxit=200)
        
        interaction [,k,n] <- fit$coefficients/Y[n,i]
      }
      #### change of  abundance #####
      
   
      
      if(Y[n,i]==0){
        interaction [,,n]<- NA
      } 
      
      missspecies <- which(Y[n,]==0,useNames = TRUE)
      interaction [missspecies,,n] <-NA
    }
  
  return(interaction)
}



interinf03 <- function(Y, X, DY, i){
  #### This function will calculate the interaction level of species i abundance 
  #### in all the  samples. Depending on the data structure, 
  #### The returned results are the 3 dimensional array
  #### array: array[1] is the species, 
  ####        array[2] is the environmental parameters, 
  ####        array[3] is the samples
  
  #### the model consider the influence from only the species abundance, use the Taylor expansion
  ####  but do not divide by the abundance of species i 
  
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   DY is the derivative of Y with respect to all the environmental parameters
  ####   i is the species i
  
  
  
  
  nsample <- nrow(DY)
  nEnv <- ncol(DY)
  nSpe <- ncol(Y)
  
  
  
  interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
  
  for (n in 1: (nsample-1)){
    
    ##### change of the DY ####
    deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
    
    #### change of  abundance #####
    
    deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
    ####  change of environmental parameters ###
    
    deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
    
    regdata <- data.frame(cbind(deltaSpe))
    ###### linear regression on   deltaSpe ~ deltaEnv
    fmal<- as.formula (paste("deltaDY ~ -1 + ",paste(colnames(regdata),collapse="+")))
    
    fit <- lm (fmal,  regdata[1:(nsample-1),])
    
    interaction [,,n] <- fit$coefficients
    
    if(Y[n,i]==0){
      interaction [,,n]<- NA
    } 
    
    missspecies <- which(Y[n,]==0,useNames = TRUE)
    interaction [missspecies,,n] <-NA
  }
  
  return(interaction)
}



interinf_low <- function(Y, X, DY, i){
  #### This function will calculate the interaction level of species i abundance in all the  samples. 
  #### The returned results are the 3 dimensional array
  #### array: array[1] is the species, 
  ####        array[2] is the environmental parameters, 
  ####        array[3] is the samples
  
  #### the model uses the low precision methods
  
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####   DY is the derivative of Y with respect to all the environmental parameters
  ####   i is the species i
  
  
  
  nsample <- nrow(DY)
  nEnv <- ncol(DY)
  nSpe <- ncol(Y)
  
  
  
    interaction <- array (0, dim=c(nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y)),colnames(DY),rownames(Y[1:(nsample-1),])))
    
    for (n in 1: (nsample-1)){
      
      ##### change of the DY ####
      deltaDY <- t(data.frame(apply(as.matrix(DY[-n,]),1,function(x) x-as.matrix(DY[n,]))))
      
      #### change of  abundance #####
      
      deltaSpe <- t(data.frame(apply(as.matrix(Y[-n,]),1,function(x) x-as.matrix(Y[n,])))) 
      ####  change of environmental parameters ###
      
      deltaEnv <- t(data.frame(apply(as.matrix(X[-n,colnames(DY)]),1,function(x) x-as.matrix(X[n,colnames(DY)]))))
      
      
      temp<-apply(deltaSpe[1:(nsample-1),],c(1,2),function(x) 1/(x*nsample))
      temp[which(!is.finite(temp))] <- 0
      interaction [,,n] <-t(temp) %*% deltaDY    ###### this is equivalent to calculate the mean values
      
      if(Y[n,i]==0){
        interaction [,,n]<- NA
      } 
      
      missspecies <- which(Y[n,]==0,useNames = TRUE)
      interaction [missspecies,,n] <-NA  
    }
  
  return(interaction)
}






#################################################################
######                       part C                          ####
####  calculate the global interaction matrix
#################################################################


interMatrix <- function(Y,X,type){
  ##### this function is used to calculate the global interaction matrix between species
  ####  there are three methods to fix the global interaction value
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  ####   Y is the abundance 
  ####   X is the environmental parameter
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  results <-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  for (species_i in 1:nSpe){
    for (species_j in 1:nSpe){
      derEnv <- rate_change (Y, X, species_i)
      derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
      interaction <- interinf (Y, X[,colnames(derEnv)], derEnv, species_i)     
      interinf_ij <- interaction[species_j,,]
      if (type == 0){
        results[species_i,species_j] <- mean(interinf_ij,na.rm = TRUE)
      }
      if (type == 1){
        results[species_i,species_j] <-median(interinf_ij,na.rm = TRUE)
      }
      if (type == 2){
        paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
        results[species_i,species_j]<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
      }
      
    }
    
  }
  return(results) 
}


interMatrix02 <- function(Y,X,type,count,factor02,speciesprecision,envprecision,threshold01){
  ##### this function is used to calculate the global interaction matrix between species
  ##### 
  ####   Y is the abundance 
  ####   X is the environmental parameter
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  #### threshold01: is the parameter used in function summ_inter_ij_new
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  #results <-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  for (species_i in 1:nSpe){
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, species_i)
    }else {
      derEnv <- rate_change_low (Y, X, species_i)
    }
    
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
    interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)
      
    }
    
    for (species_j in 1:nSpe){   
          
      interinf_ij <- interaction[species_j,,]
      if (type == 0){
        results[species_i,species_j] <- mean(interinf_ij,na.rm = TRUE)
      }
      if (type == 1){
        results[species_i,species_j] <-median(interinf_ij,na.rm = TRUE)
      }
      if (type == 2){
        paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
        results[species_i,species_j]<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
      }
      if (type == 3){
       
        Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
        number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
        
        if (number_type_no>5){
          #### if more than 5 parameter have no clear pattern, the global interaction value will be set as NA, 
          ####this threshold can be changed based on the user's data. In our test, we have 13 environmental parameters, therefore we give 5 as the threshold
          results[species_i,species_j] <- NA
        }else{
          #       x1<-Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")]
          #       x1<-as.matrix(x1)
          Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
          Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
          
          posipara<-length(which(Summ_interinf_ij$intertype=="+"))
          negapara<-length(which(Summ_interinf_ij$intertype=="-"))
          
          if((posipara==0)|(negapara>0.8*nEnv)){
            #### if there is no positive part, or the number of negative part is 80% of the total numner,  negative part is dominant 
            results[species_i,species_j] <-  Negamedian
          }else if ((posipara==0)|(posipara>0.8*nEnv)){
            #### if there is no negative part, or the number of positive part is 80% of the total numner,  positive part is dominant 
            results[species_i,species_j] <- Posimedian
          }else if (Posimedian>(factor02*abs(Negamedian))){
            #### compare the median value between the negative and positive part, if positive partis is factir02 times larger than negative part,  positive part is dominant 
            results[species_i,species_j]<- Posimedian
          }else if ((factor02*Posimedian)<(abs(Negamedian))){
            #### compare the median value between the negative and positive part, if negative partis is factir02 times larger than positive part,  negative part is dominant 
            results[species_i,species_j] <- Negamedian
          }else {
            results[species_i,species_j] <- 0
          }           
        }         
      }
    }    
  }
  
  
  return(results) 
}



interMatrix_norm <- function(Y,X,speciesprecision,envprecision){
  ##### this function is used to calculate the global interaction matrix between species
  ##### use the original norm definition
  ####   Y is the abundance 
  ####   X is the environmental parameter
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  #results <-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  for (species_i in 1:nSpe){
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, species_i)
    }else {
      derEnv <- rate_change_low (Y, X, species_i)
    }
    
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
      interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)
      
    }
    
    for (species_j in 1:nSpe){   
      
      interinf_ij <- interaction[species_j,,]
      
      temp<-interinf_ij
      temp[is.na(temp)] <- 0
      results[species_i,species_j] <- norm(temp)
              
      }
            
    }    
      
  return(results) 
}


interMatrix_norm_summ <- function(Y,X,count,speciesprecision,envprecision){
  ##### this function is used to calculate the global interaction matrix between species
  ##### use the norm definition, but exclude the extreme values
  ####   Y is the abundance 
  ####   X is the environmental parameter
  
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  #results <-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  for (species_i in 1:nSpe){
    
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, species_i)
    }else {
      derEnv <- rate_change_low (Y, X, species_i)
    }
    
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
      interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)
      }
    
    for (species_j in 1:nSpe){   
      
      interinf_ij <- interaction[species_j,,]
      
      temp<-c()
      for (para in 1:nrow(interinf_ij)){
        a <-interinf_ij[para,]
        
        counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])]
        temp<- c(temp,norm(peakvalue,type="2"))
        }
      results[species_i,species_j]<- norm(temp,type="2")
    }    
  }    

return(results)

}




interMatrix_low <- function(Y,X,type,count,factor02,threshold01){
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  #results <-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  results <-  as.data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  for (species_i in 1:nSpe){
    derEnv <- rate_change_low (Y, X, species_i)
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i) 
    for (species_j in 1:nSpe){   
      
      interinf_ij <- interaction[species_j,,]
      if (type == 0){
        results[species_i,species_j] <- mean(interinf_ij,na.rm = TRUE)
      }
      if (type == 1){
        results[species_i,species_j] <-median(interinf_ij,na.rm = TRUE)
      }
      if (type == 2){
        paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
        results[species_i,species_j]<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
      }
      if (type == 3){
        #factor02<-5
        Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
        number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
        
        if (number_type_no>5){
          results[species_i,species_j] <- NA
        }else{
          #       x1<-Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")]
          #       x1<-as.matrix(x1)
          Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
          Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
          
          posipara<-length(which(Summ_interinf_ij$intertype=="+"))
          negapara<-length(which(Summ_interinf_ij$intertype=="-"))
          
          if(posipara==0){
            results[species_i,species_j] <-  Negamedian
          }else if (negapara==0){
            results[species_i,species_j] <- Posimedian
          }else if (Posimedian>(factor02*abs(Negamedian))){
            results[species_i,species_j]<- Posimedian
          }else if ((factor02*Posimedian)<(abs(Negamedian))){
            results[species_i,species_j] <- Negamedian
          }else {
            results[species_i,species_j] <- 0
          }           
        }         
      }
    }    
  }
  
  
  return(results) 
}


interMatrix_ij <- function(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01){
  ##### this function is used to calculate the global interaction value beta_ij, one element in the interaction matrix
  ##### 
  ####   Y is the abundance 
  ####   X is the environmental parameter
  #### species_i is the index of species i
  #### species_j is the index of species j
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, integer number: this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  #### threshold01: is the parameter used in function summ_inter_ij_new
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  #results <-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  
  if (envprecision=="high"){
    derEnv <- rate_change02 (Y, X, species_i)
  }else {
    derEnv <- rate_change_low (Y, X, species_i)
  }
  
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  
  if (speciesprecision=="high"){
    interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
  }else{
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)    
  }
    
      interinf_ij <- interaction[species_j,,]
      if (type == 0){
        results <- mean(interinf_ij,na.rm = TRUE)
      }
      if (type == 1){
        results <-median(interinf_ij,na.rm = TRUE)
      }
      if (type == 2){
        paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
        results<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
      }
  
      if (type==3){
        
        Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
        number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
        
        if (number_type_no>5){
          results<- NA
        }else{
          #       x1<-Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")]
          #       x1<-as.matrix(x1)
          Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
          Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
          
          posipara<-length(which(Summ_interinf_ij$intertype=="+"))
          negapara<-length(which(Summ_interinf_ij$intertype=="-"))
          
          if(posipara==0){
            results <-  Negamedian
          }else if (negapara==0){
            results <- Posimedian
          }else if (Posimedian>(factor02*abs(Negamedian))){
            results<- Posimedian
          }else if ((factor02*Posimedian)<(abs(Negamedian))){
            results <- Negamedian
          }else {
            results <- 0
          }
          
          
        }
        
        
      }
  
  
  return(results) 
}



interMatrix_ij_low <- function(Y,X,species_i,species_j,type,count,factor02,threshold01){
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  #results <-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(Y), colnames(Y))),stringsAsFactors=F)
  
  derEnv <- rate_change_low (Y, X, species_i)
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)     
  interinf_ij <- interaction[species_j,,]
  if (type == 0){
    results <- mean(interinf_ij,na.rm = TRUE)
  }
  if (type == 1){
    results <-median(interinf_ij,na.rm = TRUE)
  }
  if (type == 2){
    paratrans <- t(X)[rownames(interinf_ij),colnames(interinf_ij)]
    results<- median(colSums(as.matrix(paratrans)*as.matrix(interinf_ij)),na.rm = TRUE)
  }
  
  if (type==3){
    
    Summ_interinf_ij<- summ_inter_ij_new(interinf_ij,count,0,threshold01)
    number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
    
    if (number_type_no>5){
      results<- NA
    }else{
      #       x1<-Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")]
      #       x1<-as.matrix(x1)
      Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
      Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
      
      posipara<-length(which(Summ_interinf_ij$intertype=="+"))
      negapara<-length(which(Summ_interinf_ij$intertype=="-"))
      
      if(posipara==0){
        results <-  Negamedian
      }else if (negapara==0){
        results <- Posimedian
      }else if (Posimedian>(factor02*abs(Negamedian))){
        results<- Posimedian
      }else if ((factor02*Posimedian)<(abs(Negamedian))){
        results <- Negamedian
      }else {
        results <- 0
      }
      
      
    }
    
    
  }
  
  
  return(results) 
}









################################################################# 
######                       part D                          ####
#### the statistic summary results of interaction results  ######
#################################################################
##### After calculating the original 4-dimensional interaction level value, we need to do some statistic summary to determine the global interaction matrix,
##### based on the pattern, we can decide whether the global interaction is positive,  negative, or no clear pattern


summ_inter_ij <- function(Y,X,species_i,species_j,count,type,speciesprecision,envprecision){
  ##### this function is used to calculate the statistic summary of the original interaction value \beta^k_ij\alpha 
  ##### 
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belong to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  
  if (envprecision=="high"){
    derEnv <- rate_change02 (Y, X, species_i)
  }else {
    derEnv <- rate_change_low (Y, X, species_i)
  }
  
  derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
  
  if (speciesprecision=="high"){
    interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, species_i) 
  }else{
    interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, species_i)
    
  }
  
  
  
  interinf_ij <- interaction[species_j,,]
  
  
  b<- as.character(c(1:count))
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  
  namelist <- c("intertype","min","1stQu","median","mean","3rdQu","max", "peakrange", "peakvalue", "variance",  "numpositive","numnegative",b)
  
  if (type==0){
    #########  acrossing the samples ############
    
    sta_inter <-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
    
    
    for (para in 1:nrow(interinf_ij)){
      a <-interinf_ij[para,]
      temp <- as.numeric(summary(a))[1:6]
      variance <- var(a,na.rm = TRUE) 
      #       Pear_cor <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a)$estimate
      #       Pear_p <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a)$p.value
      #       Spear_cor <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a,method="spearman")$estimate
      #       Spear_p <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a,method="spearman")$p.value
      nposi <- sum(  a > 0,na.rm = TRUE )
      nnega <- nsample-nposi-sum(is.na(a))
      n_non_NA<-nsample- sum(is.na(a))
      if (nposi/nsample>0.8) {
        #####  the contribution from this environmental parameter is positive
        intertype <- "+"    
      }else if (nnega/nsample>0.8){
        #####  the contribution from this environmental parameter is negative
        intertype <-"-"
      }else {
        #####  the contribution from this environmental parameter has no clear pattern
        intertype <-"no"
      }
      
      if (n_non_NA== 1 | n_non_NA== 0){
        sta_inter[para,1:6]
      }else{
        counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
        
        #sta_inter[para,]<-c(temp,peakrange,peakvalue,variance,Pear_cor,Pear_p,Spear_cor,Spear_p,nposi,nnega,counttable)
        sta_inter[para,c(-8,-1)]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
        sta_inter[para,1]<-intertype
        sta_inter[para,8]<-peakrange
      }
      
    }
    
    
  }else {
    
    
    #########  acrossing the environmental parameters #####                  
    sta_inter <-  data.frame(matrix(vector(), ncol(interinf_ij), length(namelist),dimnames=list(colnames(interinf_ij), namelist)))
    
    
    for (n in 1:nrow(interinf_ij)){
      a <-interinf_ij[,n]
      temp <- as.numeric(summary(a))[1:6]
      variance <- var(a,na.rm = TRUE) 
      nposi <- sum(  a > 0 ,na.rm = TRUE )
      nnega <- nsample-2-nposi-sum(is.na(a))
      n_non_NA<-nsample- sum(is.na(a))
      if (n_non_NA== 1 | n_non_NA== 0){
        sta_inter[para,1:6]
      }else{
        counttable<-table(cut(interinf_ij[,n],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
        
        sta_inter[n,-8]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
        sta_inter[n,8]<-peakrange
      }
      
      
    }
    
    
  }
  
  
  return(sta_inter)
}





summ_inter_ij_new <- function(interinf_ij,count,type,threshold01){
  ##### this function is used to calculate the statistic summary of the original interaction value \beta^k_ij\alpha 
  ##### 
  ####  interinf_ij is the interaction data frame of  beta_ij

  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belong to the types of "posi", "nega" or "no"

  
  nsample <- ncol(interinf_ij)
  #nEnv <- ncol(X)
  #nSpe <- ncol(Y)
  
  b<- as.character(c(1:count))
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  
  #namelist <- c("min","1stQu","median","mean","3rdQu","max", "peakrange", "peakvalue", "variance",  "PearsonCorr_ParaInter", "PearsonPvalue_ParaInter","SpearmanCorr_ParaInter","SpearmanPvalue_ParaInter",  "numpositive","numnegative",b)
  namelist <- c("intertype","min","1stQu","median","mean","3rdQu","max", "peakrange", "peakvalue", "variance",  "numpositive","numnegative",b)
  
  if (type==0){
    #########  acrossing the samples ############
    
    sta_inter <-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
    
    
    for (para in 1:nrow(interinf_ij)){
      a <-interinf_ij[para,]
      temp <- as.numeric(summary(a))[1:6]
      variance <- var(a,na.rm = TRUE) 
#       Pear_cor <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a)$estimate
#       Pear_p <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a)$p.value
#       Spear_cor <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a,method="spearman")$estimate
#       Spear_p <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a,method="spearman")$p.value
      nposi <- sum(  a > 0,na.rm = TRUE )
      nnega <- nsample-nposi-sum(is.na(a))
      n_non_NA<-nsample- sum(is.na(a))
      if (nposi/nsample>threshold01) {
        #####  the contribution from this environmental parameter is positive
        intertype <- "+"
      }else if (nnega/nsample>threshold01){
        #####  the contribution from this environmental parameter is negative
        intertype <-"-"
      }else {
        #####  the contribution from this environmental parameter has no clear pattern
        intertype <-"no"
      }

     if (n_non_NA== 1 | n_non_NA== 0){
       sta_inter[para,1:6]
     }else{
       counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
       peakposition <- which.max(counttable)
       peakrange <- names(peakposition)[1]
       peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
       peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
       
       #sta_inter[para,]<-c(temp,peakrange,peakvalue,variance,Pear_cor,Pear_p,Spear_cor,Spear_p,nposi,nnega,counttable)
       sta_inter[para,c(-8,-1)]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
       sta_inter[para,1]<-intertype
       sta_inter[para,8]<-peakrange
     }

    }
      
    
  }else {
    
    
    #########  acrossing the environmental parameters #####                  
    sta_inter <-  data.frame(matrix(vector(), ncol(interinf_ij), length(namelist),dimnames=list(colnames(interinf_ij), namelist)))
    
    
    for (n in 1:nrow(interinf_ij)){
      a <-interinf_ij[,n]
      temp <- as.numeric(summary(a))[1:6]
      variance <- var(a,na.rm = TRUE) 
      nposi <- sum(  a > 0 ,na.rm = TRUE )
      nnega <- nsample-2-nposi-sum(is.na(a))
      n_non_NA<-nsample- sum(is.na(a))
      if (n_non_NA== 1 | n_non_NA== 0){
        sta_inter[para,1:6]
      }else{
        counttable<-table(cut(interinf_ij[,n],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
        peakposition <- which.max(counttable)
        peakrange <- names(peakposition)[1]
        peakrangevalue <- unique(as.numeric(unlist(strsplit(gsub("[(]|]", "", peakrange), ","))))
        peakvalue <- median(a[which(a>peakrangevalue[1]&a<=peakrangevalue[2])])
        
        sta_inter[n,-8]<-c(temp,peakvalue,variance,nposi,nnega,as.numeric(counttable))
        sta_inter[n,8]<-peakrange
      }
      
      
    }
    
    
  }
  
  return(sta_inter)
}





summ_interaction <-function(Y,X,count,type){
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  
  
  #b<- as.character(c(1:count))
  b<- t(as.matrix(as.character(c(1:(count-1)))))
  b<- apply(b,1,function(x) paste("count",x,sep=""))
  
  namelist <- c("min","1stQu","median","mean","3rdQu","max","variance", "PearsonCorr", "PearsonPvalue","SpearmanCorr","SpearmanPvalue",  "numpositive","numnegative",b)
  
  
  if (type==0){
    #########  acrossing the samples ############
    sta_inter <-  array(NA, dim=c(nSpe,nSpe, nEnv, length(namelist)),dimnames=list(colnames(Y),colnames(Y),colnames(X),namelist))
    
    
    #interaction <- array (0, dim=c(nEnv+nSpe,nEnv,nsample-1),dimnames=list(c(colnames(Y),colnames(DY)),colnames(DY),rownames(Y[1:(nsample-1),])))
    for(species_i in 1:nSpe){
      for(species_j in 1:nSpe){
        derEnv <- rate_change (Y, X, species_i)
        derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
        
        interaction <- interinf (Y, X[,colnames(derEnv)], derEnv, species_i)     
        
        interinf_ij <- interaction[species_j,,]
        
        for (para in 1:nrow(interinf_ij)){
          a <-interinf_ij[para,]
          temp <- as.numeric(summary(a))[1:6]
          variance <- var(a,na.rm = TRUE) 
          Pear_cor <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a)$estimate
          Pear_p <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a)$p.value
          Spear_cor <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a,method="spearman")$estimate
          Spear_p <- cor.test(X[colnames(interinf_ij),rownames(interinf_ij)[para]],a,method="spearman")$p.value
          nposi <- sum(  a > 0 ,na.rm = TRUE)
          nnega <- nsample-2-nposi-sum(is.na(a))
          counttable<-table(cut(interinf_ij[para,],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
          sta_inter[species_i,species_j,rownames(interinf_ij)[para],]<-c(temp,variance,Pear_cor,Pear_p,Spear_cor,Spear_p,nposi,nnega,counttable)
        }
        
      }
    }
    
  }else {
    
    
    #########  acrossing the environmental parameters ##### 
    
    #sta_inter <-  data.frame(matrix(vector(), ncol(interinf_ij), length(namelist),dimnames=list(colnames(interinf_ij), namelist)),stringsAsFactors=F)
    sta_inter <-  array(NA, dim=c(nSpe,nSpe, nsample-2, length(namelist)),dimnames=list(colnames(Y),colnames(Y),rownames(Y)[1:(nsample-2)],namelist))
    for(species_i in 1:nSpe){
      for(species_j in 1:nSpe){  
        
        derEnv <- rate_change (Y, X, species_i)
        derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
        
        interaction <- interinf (Y, X[,colnames(derEnv)], derEnv, species_i)     
        
        interinf_ij <- interaction[species_j,,]
        
        for (n in 1:nrow(interinf_ij)){
          a <-interinf_ij[,n]
          temp <- as.numeric(summary(a))[1:6]
          variance <- var(a,na.rm = TRUE) 
          nposi <- sum(  a > 0,na.rm = TRUE )
          nnega <- nsample-2-nposi-sum(is.na(a))
          counttable<-table(cut(interinf_ij[,n],breaks=seq(min(a,na.rm = TRUE),max(a,na.rm = TRUE),length=count)))
          sta_inter[species_i,species_j,colnames(interinf_ij)[para],]<-c(temp,variance,nposi,nnega,counttable)
        }
        
      }
    }
    
  }
  
  return(sta_inter)
  
  
}






anascatter <-function(Y,X,interaction_ij,para,count){
  #####  the interaction is the interaction strength between a pair of species 2 dimensional matrix , rownames is the enviornmental para, colname is ths sample###
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  
  position <- which(rownames(interaction_ij)==para)
  
  window <-5
  N<- floor(nsample/window)
  
  
  b<- t(as.matrix(as.character(c(1:(nsample-1-window)))))
  b<- apply(b,1,function(x) paste("window",x,sep=""))
  
  namelist <- c("extreme",b)
  difference<-  data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
  paradiffer <- data.frame(matrix(vector(), nrow(interinf_ij), length(namelist),dimnames=list(rownames(interinf_ij), namelist)))
  
  
  for (n_para in (1:nrow(interaction_ij))){
    #paraorig<-X[colnames(interaction_ij),colnames(X)==para]
    paraorig <-X[colnames(interaction_ij),n_para]
    paraorder<- order(paraorig)
    #cor(interinf_ij[position,],X[1:(nsample-2),colnames(X)==para])
    
    counttable<-table(cut(paraorig,breaks=seq(min(paraorig,na.rm = TRUE),max(paraorig,na.rm = TRUE),length=count))) #### check the data extrem properties 
    
    
    diff_temp <- vector(mode="numeric", length=0)
    paradiff_temp <- vector(mode="numeric", length=0)
    if (var(counttable)>20){
      breakpoint <- seq(min(paraorig,na.rm = TRUE),max(paraorig,na.rm = TRUE),length=count)
      pp <- which(counttable == max(counttable,na.rm = TRUE))
      
      paraorig <- paraorig[which(paraorig >breakpoint[pp] & paraorig <= breakpoint[pp+1])]
      paraorder<- order(paraorig)
    }
    
    for (k in 1:(nsample-1-window)){
      #if (k!=N){
      #section_k <- interaction_ij[position,paraorder[((k-1)*window+1):(k*window)]]
      section_k <- interaction_ij[position,paraorder[k:(k+window-1)]]
      diff_temp<-c(diff_temp,max(section_k,na.rm = TRUE)-min(section_k,na.rm = TRUE))
      
      
      Parasection_k <- paraorig[paraorder[k:(k+window-1)]]
      paradiff_tem<-c(paradiff_temp,max(Parasection_k,na.rm = TRUE)-min(Parasection_k,na.rm = TRUE))
      #}
      #      else{
      #       section_k <- interaction_ij[position,paraorder[((k-1)*window+1):(nsample-2)]]
      #       diff <-c(diff,max(section_k)-min(section_k))
      #     }
      
    }
    
    difference[n_para,] <- c(max(counttable),diff_temp)
    paradiffer[n_para,] <- c(max(counttable),paradiff_temp)  
  }
  return(difference)
  
  #return(cbind(difference,paradiffer))
}



  


####################################################################
######                       part E                             ####
####  calculate the main environmental  paramters
####################################################################
main_para<-function(Y,X,type,count,factor02,speciesprecision,envprecision,threshold01){
  
  ##### this function is used to find the main parameters
  ##### Y species abundance
  ####  X environmental parameters
  
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belong to the types of "posi", "nega" or "no"
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  
  #InterMatrix_Summ<-  data.frame(matrix(vector(), nSpe, nSpe,dimnames=list(colnames(RelaAbund), colnames(RelaAbund))),stringsAsFactors=F)
  
  sig_para <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
  
  posi_domi_para <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
  nega_domi_para <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
  
  for (i in 1:nSpe){
    if (envprecision=="high"){
      derEnv <- rate_change02 (Y, X, i)
    }else {
      derEnv <- rate_change_low (Y, X, i)
    }
    
    derEnv <- derEnv[, which(apply(derEnv, 2, function(x) !any(is.na(x))))]
    
    if (speciesprecision=="high"){
      interaction <- interinf02 (Y, X[,colnames(derEnv)], derEnv, i) 
    }else{
      interaction <- interinf_low (Y, X[,colnames(derEnv)], derEnv, i)
      
    }
    
    for(j in 1:nSpe){
      interinf_ij <- interaction[species_j,,]
      #Summ_interinf_ij <- summ_inter_ij(Y,X,i,j,count,0)
      Summ_interinf_ij <- summ_inter_ij_new(interinf_ij,count,0,threshold01)
      
      number_type_no <- length((which(Summ_interinf_ij$intertype=="no")))
      
      temp <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
      temp_domi_posi <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
      temp_domi_nega <- data.frame(matrix(vector(), 1, nEnv,dimnames=list("relation", colnames(Parascale))),stringsAsFactors=F)
      
      if (number_type_no>5){
        #InterMatrix_Summ_test[i,j]<- NA
      }else{
        #       x1<-Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")]
        #       x1<-as.matrix(x1)
        #Posimedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
        #Negamedian <- median(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
        
        posipara<-rownames(Summ_interinf_ij)[which(Summ_interinf_ij$intertype=="+")]
        negapara<-rownames(Summ_interinf_ij)[which(Summ_interinf_ij$intertype=="-")]
        
        maxi_posi <- max(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")])
        mini_nega <- min(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")])
        
        domi_posi <- posipara[which(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="+"),c("peakvalue")]>0.7*maxi_posi)]
        domi_nega <- negapara[which(Summ_interinf_ij[which(Summ_interinf_ij$intertype=="-"),c("peakvalue")]<0.7*mini_nega)]
        
        temp_domi_posi[1,domi_posi]<-1
        temp_domi_nega[1,domi_nega]<-1
        
        
        if ((length(posipara)==nEnv)| (Posimedian>(factor02*abs(Negamedian))) ) {
          
          #InterMatrix_Summ[i,j] <- Posimedian
          
          temp[1,posipara]<-1
        } else if ((length(negapara)==nEnv)|((factor02*Posimedian)<(abs(Negamedian)))){
          #InterMatrix_Summ[i,j] <- Negamedian
          
          temp[1,negapara]<--1
          
        }else {
          #InterMatrix_Summ[i,j] <- 0
        }    
        
      }
      
          
      rownames(temp)<- paste("interinf",i,j,sep="_")  
      sig_para<- rbind(sig_para, temp)
      rownames(temp_domi_posi)<- paste("interinf",i,j,sep="_")  
      posi_domi_para<- rbind(posi_domi_para, temp_domi_posi)
      rownames(temp_domi_nega)<- paste("interinf",i,j,sep="_")  
      nega_domi_para<- rbind(nega_domi_para, temp_domi_nega)
      
    }
  }
  
  
out<-list(posi_domi_para,nega_domi_para,sig_para)
  return(out)
  
}



main_para_barplot<- function(out){
  
  posi_domi_para<-out[[1]]
  nega_domi_para<-out[[2]]
  sig_para<-out[[3]]
  ###### make bar plot #### 
  p1<-apply(posi_domi_para, 2, function(x) sum(x,na.rm=TRUE) )
  p2<-apply(nega_domi_para, 2, function(x) sum(x,na.rm=TRUE) )
  p3<-apply(sig_para, 2, function(x) sum(!is.na(x)) )
  
  infor_domi_para<-as.data.frame(rbind(p1,p2,p3))
  
  rownames(infor_domi_para)<- c("positive","negative","determining interaction")
  
  
  infor_domi_para$statinfor <- rownames(infor_domi_para)
  sss<- melt(infor_domi_para,id="statinfor")
  colnames(sss)[1]<-c("dominant_parameters")
  
  
  term<-ggplot(sss, aes(variable, value))+
    geom_bar(stat = "identity",aes(fill=dominant_parameters),position = "dodge",width=0.5)+
    coord_flip()+xlab("environmental parameters")+ylab("number of frequency")+
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.5), angle = 0))+
    theme(axis.text.x = element_text(size = rel(1.5), angle = 0))+
    theme(axis.text.y = element_text(size = rel(1.5), angle = 0))+
    theme(legend.title = element_text(size = rel(1.0), angle = 0))+
    theme(legend.text = element_text(size = rel(1.0), angle = 0))+
    theme(legend.position = "bottom")
  
  pdf(paste("Bar_Plot_statistic_Dominant_Parameters.pdf",sep=""),height = 10, width = 10)
  print(term)
  dev.off()
  return(1)
}






#################################################################
######                       part F                          ####
#####            robust test                        #############
#################################################################




interMatrix_ij_permut <- function(Y,X,species_i,species_j,times,size,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01){
  
  
  ###### randomly remove  size  species for times, save the 2 dimension interaction table into the data frame #####
  #####  permut_type is "samples",  "species", "para"  
  ##### 
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  ####  time, how many times in  the random test
  #### size,  control the number of movement of species or parameters
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  #### permutation type:  "species", "para", "samples"
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belongs to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  species_i_name<-colnames(Y)[species_i]
  species_j_name<-colnames(Y)[species_j]
  
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  
  #interinf_ij_permut <- vector(mode="numeric", length=0)
  testpermut <- c(1:times)
  removenames<-vector()
  for(k in 1:size){
    removenames<- c(removenames,paste(permut_type,k,"remove",sep="_"))
  }
  
  interinf_ij_permut <-  data.frame(matrix(vector(),  times+1,size+1,dimnames=list( c(testpermut,"original"),c("values",removenames))),stringsAsFactors=F)
  
  
  
  if(permut_type=="species"){
    testspecies <- colnames(Y)[c(-species_i,-species_j)]
  
    for (n in 1:times){
      kickoutspe<- sample(testspecies,size,replace = FALSE)
      
      Y_n <- Y[,-which(colnames(Y) %in% kickoutspe)]
      
      newposition_i<- which(colnames(Y_n)==species_i_name)
      newposition_j<- which(colnames(Y_n)==species_j_name)
      
      
      interinf_ij_n <-interMatrix_ij(Y_n,X,newposition_i,newposition_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      
      #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij_n )
      interinf_ij_permut[n,1]<-interinf_ij_n 
      interinf_ij_permut[n,-1]<- kickoutspe
      
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[times+1,1]<-interinf_ij
  }else if (permut_type=="samples"){
    samplesnames <- rownames(Y)
    
    for (n in 1:times){
      kickoutsample<- sample(samplesnames,size,replace = FALSE)
      
      Y_n <- Y[-which(rownames(Y) %in% kickoutsample),]
      X_n<- X[-which(rownames(Y) %in% kickoutsample),]
      interinf_ij_n <-interMatrix_ij(Y_n,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij_n ) 
      interinf_ij_permut[n,1]<-interinf_ij_n 
      interinf_ij_permut[n,-1]<- kickoutsample
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[times+1,1]<-interinf_ij
  }else{
    paranames <- colnames(X)
    
    for (n in 1:times){
      kickoutpara<- sample(paranames,size,replace = FALSE)
      
      
      X_n<- X[,-which(rownames(X) %in% kickoutpara)]
      interinf_ij_n <-interMatrix_ij(Y,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij_n ) 
      interinf_ij_permut[n,1]<-interinf_ij_n 
      interinf_ij_permut[n,-1]<- kickoutsample
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[times+1,1]<-interinf_ij
    
  }
  
  return(interinf_ij_permut)
  
}



######### add the random error to the original data ####
interMatrix_ij_env_spe_robust <- function(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01){
  ##### this function is used to calculate the robust test on species abundance and environmental parameter,  add the random error to the original data
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  #### times, how many times of the random test
  #### threshold, error bar level
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belongs to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  ### threshold means 5% or 10% or etc. 
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  N_para<-nsample*nEnv
  N_spe<-nsample*nSpe
  coef <- 1+threshold
  
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  testpermut <- c(1:times)
  
  
  interinf_ij_robust <-  data.frame(matrix(vector(),  times+1,1,dimnames=list( c(as.vector(testpermut),"original"),"values")),stringsAsFactors=F)
  
  for (n in 1:times){
    randerror_para<- matrix((1+runif(N_para,-1,1)*threshold),nsample,nEnv)
    X_n <- X*randerror_para
    
    randerror_spe<- matrix((1+runif(N_spe,-1,1)*threshold),nsample,nSpe)
    Y_n <- Y*randerror_spe
    
    
    interinf_ij_n <-interMatrix_ij(Y_n,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
    
    interinf_ij_robust[n,1]<-interinf_ij_n 
    
  }
  
  interinf_ij_robust[times+1,1]<-interinf_ij
  
  return(interinf_ij_robust)
}


interMatrix_ij_env_robust <- function(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01){
  
  ##### this function is used to calculate the robust test on environmental parameter,  add the random error to the original data
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  #### times, how many times of the random test
  #### threshold, error bar level
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belongs to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  
  ### threshold means 5% or 10% or etc. 
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  N<-nsample*nEnv
  coef <- 1+threshold
  
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  testpermut <- c(1:times)
  
  
  interinf_ij_robust <-  data.frame(matrix(vector(),  times+1,1,dimnames=list( c(as.vector(testpermut),"original"),"values")),stringsAsFactors=F)
  
  for (n in 1:times){
    randerror<- matrix((1+runif(N,-1,1)*threshold),nsample,nEnv)
    X_n <- X*randerror
    
    interinf_ij_n <-interMatrix_ij(Y,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
    
    interinf_ij_robust[n,1]<-interinf_ij_n 
    
  }
  
  interinf_ij_robust[times+1,1]<-interinf_ij
  
  return(interinf_ij_robust)
}

interMatrix_ij_spe_robust <- function(Y,X,species_i,species_j,times,threshold,type,count,factor02,speciesprecision,envprecision,threshold01){
  
  ### threshold means 5% or 10% or etc. 
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  N<-nsample*nSpe
  coef <- 1+threshold
  
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  testpermut <- c(1:times)
  
  
  interinf_ij_robust <-  data.frame(matrix(vector(),  times+1,1,dimnames=list( c(as.vector(testpermut),"original"),"values")),stringsAsFactors=F)
  
  for (n in 1:times){
    randerror<- matrix((1+runif(N,-1,1)*threshold),nsample,nSpe)
    Y_n <- Y*randerror
    
    interinf_ij_n <-interMatrix_ij(Y_n,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
    
    interinf_ij_robust[n,1]<-interinf_ij_n 
    
  }
  
  interinf_ij_robust[times+1,1]<-interinf_ij
  
  return(interinf_ij_robust)
}



interMatrix_ij_permut02 <- function(Y,X,species_i,species_j,type,count,factor02,permut_type,speciesprecision,envprecision,threshold01){
  ##### randomly remove one species  or parameter , save the 2 dimension interaction table into the array
  ####   Y is the abundance 
  ####   X is the environmental parameter
  ####  species_i, index of species i
  ####  species_j, index of species j
  #### type = 0, use the mean value 
  #### type = 1, use the median value
  #### type = 2, use the weighted mean value
  #### type = 3, use the work flow method to give the statistic estimation
  #### count, this is used to calculate the histogram distribution in count segments based on the minimum and maximum value 
  #### factor02, a factor used to compare the negative and positive parts
  
  #### threshold01: the threshold used to compare positive and negative parts, 
  #### whether the parameter belong to the types of "posi", "nega" or "no"
  
  #### speciesprecision, or envprecision == "high",  use rate_change02 and interinf02
  #### speciesprecision, or envprecision == "low" . use rate_chang_low and interinf_low
  
  ######  remove one  species  or one parameter #####
  #####  permut_type  is species  or para  #####
  nsample <- nrow(Y)
  nEnv <- ncol(X)
  nSpe <- ncol(Y)
  species_i_name<-colnames(Y)[species_i]
  species_j_name<-colnames(Y)[species_j]
  interinf_ij <-  interMatrix_ij(Y,X,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
  
  if(permut_type=="species"){
    testspecies <- colnames(Y)[c(-species_i,-species_j)]
        
    interinf_ij_permut <-  data.frame(matrix(vector(),  length(testspecies)+1,1,dimnames=list( c(testspecies,"original"),c("values"))),stringsAsFactors=F)
    
    for (n in 1:length(testspecies)){
      
      Y_n <- Y[,-(which(colnames(Y)==testspecies[n]))]
      newposition_i<- which(colnames(Y_n)==species_i_name)
      newposition_j<- which(colnames(Y_n)==species_j_name)
      
      interinf_ij_n <-interMatrix_ij(Y_n,X,newposition_i,newposition_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      interinf_ij_permut[n,1]<-interinf_ij_n 
      
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[length(testspecies)+1,1]<-interinf_ij
  }else{
    testpara <- colnames(X)
    
    interinf_ij_permut <-  data.frame(matrix(vector(),  length(testpara)+1,1,dimnames=list( c(testpara,"original"),c("values"))),stringsAsFactors=F)
    
    for (n in 1:length(testpara)){
      
      X_n <- X[,-(which(colnames(X)==testpara[n]))]
      
      interinf_ij_n <-interMatrix_ij(Y,X_n,species_i,species_j,type,count,factor02,speciesprecision,envprecision,threshold01)
      interinf_ij_permut[n,1]<-interinf_ij_n 
      
    }
    #interinf_ij_permut<-c(interinf_ij_permut,interinf_ij)
    interinf_ij_permut[length(testpara)+1,1]<-interinf_ij
  }
  
  return(interinf_ij_permut)
  
}



 
####################################################################
######                       part G                          ####
######## test data singularity
####################################################################
testsingula<-function(data){
  #### test whether the data has singlarities or not, if it has, remove the colums which bring the singularities
  
  N1<- ncol(data)
  N2<-qr(data)$rank
  names<- colnames(data)
  N3<-N1-N2
  
  if (N1==N2){
    results<-data
  }else{
    temp<- data
    while(qr(temp)$rank< ncol(temp)){
      remove<- sample(names,N3,replace = FALSE)
      temp<- data[,-which(names==remove)]
      
    }
    
    results<-temp
  }
  return(results)
}





####################################################################
######                       part H                             ####
###### extrac the six examples from the interaction matrix
####################################################################
 
six_test_robust <- function(data){
  max_positive<-max(data[which(data>0,arr.ind=TRUE)],na.rm=TRUE)
  median_positive<-median(data[which(data>0,arr.ind=TRUE)],na.rm=TRUE)
  min_positive<-min(data[which(data>0,arr.ind=TRUE)],na.rm=TRUE)
  
  min_negative<-min(data[which(data<0,arr.ind=TRUE)],na.rm=TRUE)
  median_negative<-median(data[which(data<0,arr.ind=TRUE)],na.rm=TRUE)
  max_negative<-max(data[which(data<0,arr.ind=TRUE)],na.rm=TRUE)
  
  Strong_posi<- which(data == max_positive, arr.ind =TRUE, useNames = TRUE)
  
  max_median_posi<-min(data[which(data>=median_positive,arr.ind=TRUE)],na.rm=TRUE)
  Median_posi <-which(data == max_median_posi, arr.ind =TRUE, useNames = TRUE)
  
  Weak_posi <-which(data == min_positive, arr.ind =TRUE, useNames = TRUE)
  
  Strong_nega<- which(data == min_negative, arr.ind =TRUE, useNames = TRUE)
  
  max_median_nega<-min(data[which(data>=median_negative,arr.ind=TRUE)],na.rm=TRUE)
  Median_nega <-which(data == max_median_nega, arr.ind =TRUE, useNames = TRUE)
  Weak_nega <-which(data == max_negative, arr.ind =TRUE, useNames = TRUE)
  
  six_type<- c("strong_posi","median_posi","weak_posi","strong_nega","median_nega","weak_nega")
  
  six_rela <-  data.frame(matrix(vector(),  6,2,dimnames=list( six_type,c("row","column"))),stringsAsFactors=F)

  
  six_rela[1,]<-Strong_posi
  six_rela[2,]<-Median_posi
  six_rela[3,]<-Weak_posi
  six_rela[4,]<-Strong_nega
  six_rela[5,]<-Median_nega
  six_rela[6,]<-Weak_nega
  
  return(six_rela)
  
}







####################################################################
######                       part I                             ####
########                make some plots                    #########
####################################################################


#### reshape the data for ggplot #####
reshapedata<- function(interactiondata,Sample_Info){
  
  
  
  ############# do the necessary reshaping of the data frame
  ## also: add the values of the environmental parameter itself
  x <- as.data.frame(t(interactiondata))
  
  needed_colnames <- colnames(x)
  # make sure that it is clear that the column names show the interaction value
  colnames(x) <- paste(colnames(x), "_IAV", sep = "")
  x$plot <- rownames(x)
  
  # add the values of the environmental parameters per plot
  
  x <- merge(x, Sample_Info[, c("Plot", needed_colnames)], by.x = "plot", by.y = "Plot", sort = FALSE)
  
  xx <- melt(x, id = c("plot", needed_colnames))
  colnames(xx)[colnames(xx) == "value"] <- "value_IAV"
  colnames(xx)[colnames(xx) == "variable"] <- "envpara_IAV"
  
  
  xxx <- melt(xx, id = c("plot", "envpara_IAV", "value_IAV"))
  
  
  xxx$check <- gsub("_IAV", "", xxx$envpara_IAV)
  xxx$checktrue <- xxx$check == xxx$variable
  xxx <- droplevels(subset(xxx, checktrue == TRUE))
  return(xxx)
  
  
}  



interinfPlot <- function(xxx,i,j,specieslist){
  ############################# Make each interaction plot separately
  ######################## This code adds the information on who is interacting on whom
  # remove any row with NA
  xxx <- xxx[which(apply(xxx, 1, function(x) !any(is.na(x)))), ]
  species_i<-specieslist[i]
  species_j<-specieslist[j]
  
  
  species_j_on_species_i <- xxx
  
  information <-paste("Interaction of", species_j, "on", species_i, sep=" ")
  species_j_on_species_i$comp <- rep(information,  nrow(species_j_on_species_i))
  
  
  ########### the code below makes the figure, and writes into a separate object
  species_j_on_species_i_P <- ggplot(species_j_on_species_i, aes(value, value_IAV)) +
    theme_bw() +
    facet_wrap( ~ envpara_IAV, scales = "free", nrow = 1) +  # compare to: scales = "free_x"
    geom_hline(yintercept = 0, color = "green", size = 1) +
    geom_point() +
    ggtitle(species_j_on_species_i[1, "comp"])
  
  filenames <-paste(information,".pdf",sep="")
  setwd("./figures")
  # This code compiles all the figures into one plot.
  pdf(filenames, height = 15, width = 25)
  grid.arrange( species_j_on_species_i_P, ncol = 1)
  
  dev.off()
  setwd("..")
  
  return(species_j_on_species_i_P)
}

