#########################
# This R file only contains some library functions.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

library(mice)
library(Rcpp)

sourceCpp("NAP.cpp")

### Calculate mean difference test statistic

mean_diff <- function(scores.a, scores.b, method, direction)
{
  if(method == "marker")
  {
    scores.a <- as.vector(na.omit(scores.a))
    scores.b <- as.vector(na.omit(scores.b))
    
    if((length(scores.a) * length(scores.b)) == 0)
      return(Inf)
  }
  
  result <- abs(mean(scores.b) - mean(scores.a)) 
  return(result)
}

### Calculate NAP

NAP_calculation <- function(scores.a, scores.b, method, direction)
{
  if(method == "marker")
  {
    scores.a <- as.vector(na.omit(scores.a))
    scores.b <- as.vector(na.omit(scores.b))
    
    if((length(scores.a) * length(scores.b)) == 0)
      return(Inf)
  }
  
  ### Use C++ function for fast comparisons
  NAP.observed <- NAP_cpp(scores.a, scores.b)

  ### Make NAP symmetric
  if (NAP.observed < 0.5) 
  {
    NAP.observed <- 1 - NAP.observed
  }
  
  return(NAP.observed)
}

### Calculate test statistic (either MD or NAP)

ESM_calc <- function(ESM, scores.a, scores.b, method, direction)
{
  if(ESM == "MD")
  {
    stat <- mean_diff(scores.a = scores.a, scores.b = scores.b, method = method, direction = direction)
  } else 
  if(ESM == "NAP")
  {
    stat <- NAP_calculation(scores.a = scores.a, scores.b = scores.b, method = method, direction = direction)
  } else
  {
    stop(paste("Unexpected ESM:", ESM))
  }
  
  return(stat)
}

### Add leading and lagging vectors as covariates for univariate data

Add_covariates <- function(data, model)
{
  if(!(model %in% c("mvn.3", "mvn.6")))
  {
    data <- data[,1:2]
    nrw <- nrow(data)
    data["Lead"] <- c(data[-1, 2], NA)
    data["Lag"] <- c(NA, data[-nrw, 2])
  }
  
  return(data)
}

### Generate multiple imputations from dataset

MI_handler <- function(data, nMI, model)
{
  data <- Add_covariates(data, model)
  data <- data[,2:4]
  
  mi <- mice(data, m = nMI, method = "norm", remove_collinear = FALSE, printFlag = FALSE)
  return(mi)
}

### Perform randomization test and calculate p-value

Compute_RT <- function(
  data_list,   # Simulated dataset(s)
  design,      # SCE design type
  model,       # Data model
  ESM,         # Test statistic
  method,      # Missing data handling method
  alfa,        # Level of significance
  direction,   # Direction of test statistic (Only used if randomization test is one-sided)
  limit_phase, # Minimum number of measurements in a phase
  nMBD,        # Number of participants in MBD
  nMC,         # Number of randomizations in Monte Carlo randomization test
  nMI,         # Number of imputations in multiple imputation
  ABAB_idx     # Phase change indices for all possible ABAB randomizations
)
{ 
  ### Randomization test for RBD
  
  if (design == "RBD") 
  { 
    ### Calculate observed test statistic
    
    data <- data_list[[1]]
    
    if(method == "MI")
    {
      ### Generate multiple imputations
      
      mi <- MI_handler(data, nMI, model)
      index.a <- data[, 1] == "A"
      index.b <- data[, 1] == "B"
      
      completes <- matrix(nrow = nrow(data), ncol = nMI)
      imps <- numeric(nMI)
      
      ### Calculate average statistic from multiple imputations
      
      for(i in 1:nMI)
      {
        temp <- complete(mi, i)[,1]
        observed.a <- temp[index.a]
        observed.b <- temp[index.b]
        completes[,i] <- temp
        
        imps[i] <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
      }
      ESM_obs <- mean(imps)
      
    } else
    {
      ### Calculate observed test statistic for complete data or randomized marker method
      
      observed.a <- data[, 2][data[, 1] == "A"]
      observed.b <- data[, 2][data[, 1] == "B"]
      
      ESM_obs <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
    }
    
    if(is.infinite(ESM_obs))
      stop("Observed test statistic invalid")
    
    ### Generate Monte Carlo randomizations and calculate statistic for each randomization
    
    observed <- data[, 2]
    MT <- length(observed)
    RD <- numeric(nMC)
    ab <- c("A", "B")
    ba <- c("B", "A")
    
    for (i in 1:(nMC-1))
    {
      assignment <- character(MT) 
      rand <- runif(MT/2)
      
      ### Randomly generate a Monte Carlo randomization
      
      for(j in 1:(MT/2))
      {
        index <- if(rand[j] < 0.5) ab else ba
        assignment[(j*2-1):(j*2)] <- index
      }
      
      if(method == "MI")
      {
        ### Calculate average statistic from multiple imputations
        
        index.a <- assignment == "A"
        index.b <- assignment == "B"
        
        imps <- numeric(nMI)
        for(k in 1:nMI)
        {
          ascores <- completes[,k][index.a]
          bscores <- completes[,k][index.b]
          
          imps[k] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
        }
        RD[i] <- mean(imps)
        
      } else
      {
        ### Calculate statistic for complete data or randomized marker method
        
        ascores <- observed[assignment == "A"] 
        bscores <- observed[assignment == "B"]
        
        RD[i] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
      }
      
    }
    
    ### Add observed statistic to Monte Carlo distribution
    RD[nMC] <- ESM_obs
  }
  
  ### Randomization test for ABAB
  
  if(design == "ABAB")
  {
    ### Calculate observed test statistic
    
    data <- data_list[[1]]
    
    if(method == "MI")
    {
      ### Generate multiple imputations
      
      mi <- MI_handler(data, nMI, model)
      index.a <- data[, 1] == "A"
      index.b <- data[, 1] == "B"
      
      completes <- matrix(nrow = nrow(data), ncol = nMI)
      imps <- numeric(nMI)
      
      ### Calculate average statistic from multiple imputations
      
      for(i in 1:nMI)
      {
        temp <- complete(mi, i)[,1]
        observed.a <- temp[index.a]
        observed.b <- temp[index.b]
        completes[,i] <- temp
        
        imps[i] <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
      }
      ESM_obs <- mean(imps)
      
    } else
    {
      ### Calculate observed test statistic for complete data or randomized marker method
      
      observed.a <- data[, 2][data[, 1] == "A"]
      observed.b <- data[, 2][data[, 1] == "B"]
      
      ESM_obs <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
    }
    
    if(is.infinite(ESM_obs))
      stop("Observed test statistic invalid")
    
    ### Generate Monte Carlo randomizations 
    
    RD <- numeric(nMC)
    observed <- data[, 2]
    MT <- length(observed)
    
    ### Check if list of all possible ABAB randomizations has been passed down from power calculation function
    ### If yes, use that list and randomly select Monte Carlo randomizations
    ### If not, calculate list of all possible randomizations and randomly select Monte Carlo randomizations
    
    if(ncol(ABAB_idx) == 3)
    {
      quantity <- nrow(ABAB_idx)
      selection <- sample(1:quantity,nMC-1,replace=TRUE) # Randomly select Monte Carlo randomizations
      
      index.a1 <- ABAB_idx[,1]
      index.b1 <- ABAB_idx[,2]
      index.a2 <- ABAB_idx[,3]
    } else
    {
      quantity<-choose(MT-4*limit_phase+3,3)
      selection<-sample(1:quantity,nMC-1,replace=TRUE) # Randomly select Monte Carlo randomizations
      index1<-1:(MT-4*limit_phase+1)
      index2<-rev(cumsum(index1))
      
      index.a1<-numeric()
      for(it in 1:length(index1))
      {
        index.a1<-c(index.a1,(rep((limit_phase+index1[it]-1),index2[it])))
      }
      
      index.b1<-numeric()
      for(itr in index1)
      {
        for(it in (itr-1):(MT-4*limit_phase))
        {
          index.b1<-c(index.b1,rep((2*limit_phase+it),(MT-4*limit_phase+1-it)))
        }
      }
      
      indexa2<-numeric()
      for(it in 1:length(index1))
      {
        indexa2<-c(indexa2,(index1[it]:length(index1)))
      }
      index.a2<-numeric()
      for(it in 1:length(indexa2))
      {
        index.a2<-c(index.a2,(4*limit_phase-limit_phase-1+(indexa2[it]:length(index1))))
      }
    }
    
    ### Calculate statistic for each randomization
    
    for(it in 1:(nMC-1))
    {
      if(method == "MI")
      {
        ### Calculate average statistic from multiple imputations
        
        index.a <- c(1:(index.a1[selection[it]]), (1+index.b1[selection[it]]):index.a2[selection[it]])
        index.b <- c((1+index.a1[selection[it]]):index.b1[selection[it]], (1+index.a2[selection[it]]):MT)
        
        imps <- numeric(nMI)
        for(k in 1:nMI)
        {
          ascores <- completes[,k][index.a]
          bscores <- completes[,k][index.b]
          
          imps[k] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
        }
        RD[it] <- mean(imps)
        
      } else
      {
        ### Calculate statistic for complete data or randomized marker method
        
        ascores <- c(observed[1:(index.a1[selection[it]])], observed[(1+index.b1[selection[it]]):index.a2[selection[it]]])
        bscores <- c(observed[(1+index.a1[selection[it]]):index.b1[selection[it]]], observed[(1+index.a2[selection[it]]):MT])
        
        RD[it] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
      }
    }
    
    ### Add observed statistic to Monte Carlo distribution
    RD[nMC] <- ESM_obs
  }
  
  ### Randomization test for MBD
  
  if(design == "MBD") 
  { 
    ### Calculate observed test statistic
    
    MT <- nrow(data_list[[1]])
    diff <- numeric(nMBD)
    imps <- numeric(nMI)
    
    if(method == "MI")
    {
      ### Impute data for each participant separately
      
      completes_list <- list()
      
      for(p in 1:nMBD)
      {
        data <- data_list[[p]]
        
        ### Generate multiple imputations
        
        mi <- MI_handler(data, nMI, model)
        index.a <- data[, 1] == "A"
        index.b <- data[, 1] == "B"
        completes <- matrix(nrow = MT, ncol = nMI)
        
        ### Calculate average statistic from multiple imputations
        
        for(i in 1:nMI)
        {
          temp <- complete(mi, i)[,1]
          observed.a <- temp[index.a]
          observed.b <- temp[index.b]
          completes[,i] <- temp
          
          imps[i] <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
        }
        
        completes_list[[p]] <- completes
        diff[p] <- mean(imps)
      }
      
      ESM_obs <- mean(diff)
      
    } else
    {
      ### Calculate observed test statistic for complete data or randomized marker method
      
      for(p in 1:nMBD)
      {
        data <- data_list[[p]]
        
        observed.a <- data[, 2][data[, 1] == "A"]
        observed.b <- data[, 2][data[, 1] == "B"]
        
        diff[p] <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
      }
      
      ESM_obs <- mean(diff[is.finite(diff)])
    }
    
    if(is.nan(ESM_obs))
      stop("Observed test statistic invalid")
    
    ### Generate Monte Carlo randomizations and calculate statistic for each randomization
    
    RD <- numeric(nMC)
    startpoints <- (limit_phase + 1):(MT - limit_phase + 1)
    
    for (i in 1:(nMC-1))
    {
      ### Randomly select start points for each participant
      combstartpts <- sample(startpoints, nMBD) 
      
      if(method == "MI")
      {
        ### Calculate average statistic from multiple imputations for each participant
        
        for(p in 1:nMBD)
        {
          index.a <- 1:(combstartpts[p] - 1)
          index.b <- combstartpts[p]:MT
          completes <- completes_list[[p]]
          
          for(k in 1:nMI)
          {
            ascores <- completes[,k][index.a]
            bscores <- completes[,k][index.b]
            
            imps[k] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
          }
          diff[p] <- mean(imps)
        }
        
        RD[i] <- mean(diff)
        
      } else
      {
        ### Calculate statistic for each participant (for complete data or randomized marker method)
        
        for(p in 1:nMBD)
        {
          observed <- data_list[[p]][,2]
          
          ascores <- observed[1:(combstartpts[p] - 1)] 
          bscores <- observed[combstartpts[p]:MT]
          
          diff[p] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
        }
        
        RD[i] <- mean(diff[is.finite(diff)])
      }
      
      if(is.nan(RD[i]))
        RD[i] <- Inf
    }
    
    ### Add observed statistic to Monte Carlo distribution
    RD[nMC] <- ESM_obs
  }

  ### Calculate p-value from randomization distribution and determine if H0 can be rejected
  
  test <- RD >= ESM_obs
  pvalue <- sum(test) / nMC  
  reject <- if(pvalue <= alfa) 1 else 0
  
  output <- c(ESM_obs, pvalue, reject)
  return(output)
}

