#########################
# This R file only contains some library functions.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

library(imputeTS)
library(data.table)
library(mice)
library(Rcpp)
sourceCpp("NAP.cpp")

### Calculate mean difference test statistic

mean_diff <- function(scores.a, scores.b, method, direction)
{
  if(method == "marker" || method == "TS")
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
  if(method == "marker" || method == "TS")
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

### Impute missing data using time series method
### Return data with missing data in place if model fitting fails

TS_handler <- function(data, num_reps)
{
  for(i in 1:num_reps)
  {
    data[, i*2] <- tryCatch(
      na.kalman(data[, i*2], model = "auto.arima", smooth = TRUE), 
      error = function(e) { 
        print("ARIMA failed!")
        print(e)
        return(data[, i*2])
      }
    )
  }
  
  return(data)
}

### Generate multiple imputations from dataset

MI_handler <- function(data, num_MI)
{
  data <- data[,1:2]
  
  data["Lead"] <- shift(data[,2], type = "lead")
  data["Lag"] <- shift(data[,2], type = "lag")
  data <- data[,2:4]
  
  mi <- mice(data, m = num_MI, method = c(rep("norm", 3)), remove_collinear = FALSE, printFlag = FALSE)
  return(mi)
}

### Perform randomization test and calculate p-value

ESM.random.RT <- function(
  data,        # Simulated dataset
  design,      # SCE design type
  number,      # Number of randomizations in Monte Carlo randomization test
  ESM,         # Test statistic
  method,      # Missing data handling method
  limit_phase, # Minimum number of measurements in a phase
  reps_MBD,    # Number of participants in MBD
  num_MI,      # Number of imputations in multiple imputation
  direction,   # Direction of test statistic (Only used if randomization test is one-sided)
  alfa,        # Level of significance
  ABAB_idx     # Phase change indices for all possible ABAB randomizations
)
{ 
  ### Impute missing data using time series method
  
  if(method == "TS")
  {
    data <- TS_handler(data = data, num_reps = (if(design == "MBD") reps_MBD else 1)) 
  }
  
  ### Randomization test for RBD
  
  if (design == "RBD") 
  { 
    ### Calulate observed test statistic
    
    if(method == "MI")
    {
      ### Generate multiple imputations
      
      mi <- MI_handler(data, num_MI)
      index.a <- data[, 1] == "A"
      index.b <- data[, 1] == "B"
      
      completes <- matrix(nrow = nrow(data), ncol = num_MI)
      imps <- numeric(num_MI)
      
      ### Calculate average statistic from multiple imputations
      
      for(i in 1:num_MI)
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
      ### Calulate observed test statistic for complete data, randomized marker method, or time series method
      
      observed.a <- data[, 2][data[, 1] == "A"]
      observed.b <- data[, 2][data[, 1] == "B"]
      
      ESM_obs <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
    }
    
    if(is.infinite(ESM_obs))
      stop("Observed test statistic invalid")
    
    ### Generate Monte Carlo randomizations and calculate statistic for each randomization
    
    observed <- data[, 2]
    MT <- length(observed)
    RD <- numeric(number)
    ab <- c("A", "B")
    ba <- c("B", "A")
    
    for (i in 1:(number-1))
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
        
        imps <- numeric(num_MI)
        for(k in 1:num_MI)
        {
          ascores <- completes[,k][index.a]
          bscores <- completes[,k][index.b]
          
          imps[k] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
        }
        RD[i] <- mean(imps)
        
      } else
      {
        ### Calulate statistic for complete data, randomized marker method, or time series method
        
        ascores <- observed[assignment == "A"] 
        bscores <- observed[assignment == "B"]
        
        RD[i] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
      }
      
    }
    
    ### Add observed statistic to Monte Carlo distribution
    RD[number] <- ESM_obs
  }
  
  ### Randomization test for ABAB
  
  if(design == "ABAB")
  {
    ### Calulate observed test statistic
    
    if(method == "MI")
    {
      ### Generate multiple imputations
      
      mi <- MI_handler(data, num_MI)
      index.a <- data[, 1] == "A"
      index.b <- data[, 1] == "B"
      
      completes <- matrix(nrow = nrow(data), ncol = num_MI)
      imps <- numeric(num_MI)
      
      ### Calculate average statistic from multiple imputations
      
      for(i in 1:num_MI)
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
      ### Calulate observed test statistic for complete data, randomized marker method, or time series method
      
      observed.a <- data[, 2][data[, 1] == "A"]
      observed.b <- data[, 2][data[, 1] == "B"]
      
      ESM_obs <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
    }
    
    if(is.infinite(ESM_obs))
      stop("Observed test statistic invalid")
    
    ### Generate Monte Carlo randomizations 
    
    RD <- numeric(number)
    observed <- data[, 2]
    MT <- length(observed)
    
    ### Check if list of all possible ABAB randomizations has been passed down from power calculation function
    ### If yes, use that list and randomly select Monte Carlo randomizations
    ### If not, calculate list of all possible randomizations and randomly select Monte Carlo randomizations
    
    if(ncol(ABAB_idx) == 3)
    {
      quantity <- nrow(ABAB_idx)
      selection <- sample(1:quantity,number-1,replace=TRUE) # Randomly select Monte Carlo randomizations
      
      index.a1 <- ABAB_idx[,1]
      index.b1 <- ABAB_idx[,2]
      index.a2 <- ABAB_idx[,3]
    } else
    {
      quantity<-choose(MT-4*limit_phase+3,3)
      selection<-sample(1:quantity,number-1,replace=TRUE) # Randomly select Monte Carlo randomizations
      index1<-1:(MT-4*limit_phase+1)
      index2<-rev(cumsum(index1))
      
      index.a1<-numeric()
      for(it in 1:length(index1)){
        index.a1<-c(index.a1,(rep((limit_phase+index1[it]-1),index2[it])))
      }
      
      index.b1<-numeric()
      for(itr in index1){
        for(it in (itr-1):(MT-4*limit_phase)){
          index.b1<-c(index.b1,rep((2*limit_phase+it),(MT-4*limit_phase+1-it)))
        }
      }
      
      indexa2<-numeric()
      for(it in 1:length(index1)){
        indexa2<-c(indexa2,(index1[it]:length(index1)))
      }
      index.a2<-numeric()
      for(it in 1:length(indexa2)){
        index.a2<-c(index.a2,(4*limit_phase-limit_phase-1+(indexa2[it]:length(index1))))
      }
    }
    
    ### Calculate statistic for each randomization
    
    for(it in 1:(number-1))
    {
      if(method == "MI")
      {
        ### Calculate average statistic from multiple imputations
        
        index.a <- c(1:(index.a1[selection[it]]), (1+index.b1[selection[it]]):index.a2[selection[it]])
        index.b <- c((1+index.a1[selection[it]]):index.b1[selection[it]], (1+index.a2[selection[it]]):MT)
        
        imps <- numeric(num_MI)
        for(k in 1:num_MI)
        {
          ascores <- completes[,k][index.a]
          bscores <- completes[,k][index.b]
          
          imps[k] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
        }
        RD[it] <- mean(imps)
        
      } else
      {
        ### Calulate statistic for complete data, randomized marker method, or time series method
        
        ascores <- c(observed[1:(index.a1[selection[it]])], observed[(1+index.b1[selection[it]]):index.a2[selection[it]]])
        bscores <- c(observed[(1+index.a1[selection[it]]):index.b1[selection[it]]], observed[(1+index.a2[selection[it]]):MT])
        
        RD[it] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
      }
    }
    
    ### Add observed statistic to Monte Carlo distribution
    RD[number] <- ESM_obs
  }
  
  ### Randomization test for MBD
  
  if(design == "MBD")
  {
    N<-reps_MBD
    MT<-nrow(data)
    diff<-numeric(N)
    
    ### Calulate observed test statistic
    
    if(method == "MI")
    {
      ### Impute data for each participant separately
      
      completes <- list()
      for(it in 1:N)
      {
        ### Generate multiple imputations
        
        mi <- MI_handler(data[,(it*2-1):(it*2)], num_MI)
        index.a <- data[, (it*2-1)] == "A"
        index.b <- data[, (it*2-1)] == "B"
        
        completes[[it]] <- matrix(nrow = nrow(data), ncol = num_MI)
        imps <- numeric(num_MI)
        
        ### Calculate average statistic from multiple imputations
        
        for(i in 1:num_MI)
        {
          temp <- complete(mi, i)[,1]
          observed.a <- temp[index.a]
          observed.b <- temp[index.b]
          completes[[it]][,i] <- temp
          
          imps[i] <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
        }
        diff[it] <- mean(imps)
      }
      
    } else
    {
      ### Calulate observed test statistic for complete data, randomized marker method, or time series method
      
      for(it in 1:N)
      {
        observed.a<-data[,it*2][data[,(it*2)-1]=="A"]
        observed.b<-data[,it*2][data[,(it*2)-1]=="B"]
        
        diff[it] <- ESM_calc(ESM = ESM, scores.a = observed.a, scores.b = observed.b, method = method, direction = direction)
      }
    }
    
    ### Calculate average test statistic from multiple participants
    ESM_obs <- mean(diff[is.finite(diff)])
    
    if(is.nan(ESM_obs))
      stop("Observed test statistic invalid")
    
    startpoints <- (limit_phase + 1):(MT - limit_phase + 1)
    RD <- numeric(number)
    
    ### Generate Monte Carlo randomizations and calculate statistic for each randomization
    
    for(i in 1:(number-1))
    {
      combstartpts <- sample(startpoints, N) # Randomly select start points for each participant
      diff<-numeric(N)
      
      if(method == "MI")
      {
        ### Calculate average statistic from multiple imputations for each participant
        
        for(it in 1:N)
        {
          index.a <- 1:(combstartpts[it]-1)
          index.b <- combstartpts[it]:MT
          
          imps <- numeric(num_MI)
          for(k in 1:num_MI)
          {
            ascores <- completes[[it]][,k][index.a]
            bscores <- completes[[it]][,k][index.b]
            
            imps[k] <- ESM_calc(ESM = ESM, scores.a = ascores, scores.b = bscores, method = method, direction = direction)
          }
          diff[it] <- mean(imps)
        }
        
      } else
      {
        ### Calulate statistic for each participant (for complete data, randomized marker method, or time series method)
        
        for(it in 1:N)
        {
          scores.a <- data[1:(combstartpts[it]-1), it*2]
          scores.b <- data[combstartpts[it]:MT, it*2]
          
          diff[it] <- ESM_calc(ESM = ESM, scores.a = scores.a, scores.b = scores.b, method = method, direction = direction)
        }
      }
      
      ### Calculate average statistic from multiple participants
      RD[i] <- mean(diff[is.finite(diff)])
      
      if(is.nan(RD[i]))
        RD[i] <- Inf
    }
    
    ### Add observed statistic to Monte Carlo distribution
    RD[number] <- ESM_obs
  }
  
  ### Calculate p-value from randomization distribution and determine if H0 can be rejected
  
  test <- RD >= ESM_obs
  pvalue <- sum(test) / number  
  reject <- if(pvalue <= alfa) 1 else 0
  
  output <- c(ESM_obs, pvalue, reject)
  return(output)
}

