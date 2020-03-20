#########################
# This R file only contains some library functions.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

library(MASS)

source("RT.R")

### Generate simulation data

Generate_data <- function(model, N, missprop, misstype, AR, corr) 
{
  data_sim <- switch(model,
    "normal" = rnorm(n = N, mean = 0, sd = 1),
    "uniform" = runif(n = N, min = 0, max = sqrt(12)),
    "AR1" = arima.sim(model = list(ar = AR), n = N),
    "mvn" = mvrnorm(
      n = N, 
      mu = c(0,0,0), 
      Sigma = matrix(c(1, corr, corr, corr, 1, corr, corr, corr, 1), ncol = 3)
    )
  )
  
  data_sim <- as.data.frame(data_sim)
  
  ### Introduce missingness
  
  if(misstype %in% c("censor+", "censor-", "mvn+", "mvn-"))
  {
    top_censor <- (misstype %in% c("censor+", "mvn+"))
    col <- if(misstype %in% c("mvn+", "mvn-")) 2 else 1
    
    sorted <- sort(data_sim[, col], decreasing = top_censor)
    trunc <- sorted[missprop * N]
    miss <- if(top_censor) (data_sim[, col] >= trunc) else (data_sim[, col] <= trunc)
    
    data_sim[miss, 1] <- NA
  }
  
  return(data_sim)
}

### Check if data for any treatment level is completely missing

Check_data_insufficient <- function(data_input)
{
  check <- !is.na(data_input[, 2])
  check.a <- check[data_input[, 1] == "A"]
  check.b <- check[data_input[, 1] == "B"]
  checkall <- sum(check.a) * sum(check.b)
  return(checkall == 0)
}

### Estimate power and type I error rate

Calculate_power_RT <- function(
  design,          # SCE design type
  model,           # Data model
  ESM,             # Test statistic
  ES,              # Effect size
  N,               # Number of measurement
  method,          # Missing data handling method
  missprop = 0,    # Proportion of missing data 
  misstype = "",   # Mechanism used to generate missing data
  alfa,            # Level of significance
  AR,              # Autocorrelation
  corr,            # Correlation between bivariate normal vectors
  direction = "+", # Direction of test statistic (Only used if randomization test is one-sided)
  limit_phase,     # Minimum number of measurements in a phase
  nCP,             # Number of randomizations per simulated dataset (>1 only if calcualting conditional power)
  nMC,             # Number of randomizations in Monte Carlo randomization test
  nMI              # Number of imputations in multiple imputation
)
{
  ### Simulate dataset
  data_obs <- Generate_data(model = model, N = N, missprop = missprop, misstype = misstype, AR = AR, corr = corr)
  
  ABAB_idx <- matrix()
  
  ### Randomly select randomization for RBD
  
  if (design == "RBD")
  {
    ab <- c("A", "B")
    ba <- c("B", "A")
    assignments <- matrix(nrow = nCP, ncol = N)
    
    for(i in 1:nCP)
    {
      assignment <- character(N) 
      rand <- runif(N/2)
      
      for(j in 1:(N/2))
      {
        index <- if(rand[j] < 0.5) ab else ba
        assignment[(j*2-1):(j*2)] <- index
      }
      assignments[i,] <- assignment
    }
  }
  
  ### Randomly select randomization for ABAB design
  
  if(design == "ABAB")
  {
    ### Calculate all possible randomizations for ABAB
    
    quantity<-choose(N-4*limit_phase+3,3)
    selection<-sample(1:quantity,nCP,replace=TRUE)
    index1<-1:(N-4*limit_phase+1)
    index2<-rev(cumsum(index1))
    
    index.a1<-numeric()
    for(it in 1:length(index1))
    {
      index.a1<-c(index.a1,(rep((limit_phase+index1[it]-1),index2[it])))
    }
    
    index.b1<-numeric()
    for(itr in index1)
    {
      for(it in (itr-1):(N-4*limit_phase))
      {
        index.b1<-c(index.b1,rep((2*limit_phase+it),(N-4*limit_phase+1-it)))
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
    
    ### Select assignment from list of all possible randomizations
    
    assignments <- matrix(nrow=nCP, ncol=N)
    for(it in 1:nCP)
    {
      assignments[it, ] <- c(
        rep("A", index.a1[selection[it]]),
        rep("B", index.b1[selection[it]] - index.a1[selection[it]]),
        rep("A", index.a2[selection[it]] - index.b1[selection[it]]),
        rep("B", N - index.a2[selection[it]])
      )
    }
    
    ### Save all possible randomizations for ABAB
    ABAB_idx <- matrix(c(index.a1, index.b1, index.a2), ncol = 3, byrow = FALSE)
  }
  
  ### Calculate p-value and whether H0 is rejected for simulated dataset
  
  count_rejections <- numeric()
  pvalues <- numeric()
  
  for(i in 1:nCP)
  {
    ### Add effect size
    
    data_shift <- data_obs
    labels <- assignments[i,] 
    data_shift[labels == "B", 1] <- data_obs[labels == "B", 1] + ES 
    
    data_input <- cbind(as.data.frame(labels, stringsAsFactors = FALSE), data_shift)
    
    #print(data_input)
    
    if(method != "full" && Check_data_insufficient(data_input = data_input))
    {
      ### Reject H0 is data entirely missing for a treatment level
      
      count_rejections[i] <- 0
      pvalues[i] <- NA
      #print("missing statistic!")
    } else
    {
      ### Dispatch data to randomization test function
      
      output <- Compute_RT(
        data = data_input, 
        design = design, 
        model = model,
        ESM = ESM, 
        method = method,
        alfa = alfa,
        direction = direction, 
        limit_phase = limit_phase,
        nMC = nMC, 
        nMI = nMI,
        ABAB_idx = ABAB_idx
      )
      
      count_rejections[i] <- output[3]
      pvalues[i] <- output[2]
    }
  }
  
  ### Return results (H0 rejection or not)
  
  conditional_power <- mean(count_rejections) * 100
  return(conditional_power)
}
