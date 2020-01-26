#########################
# This R file only contains some library functions.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

source("RT.R")

### Generate simulation data

Generate_data <- function(model, AR, N, method, missing) 
{
  if(model == "normal")
  {
    data_sim <- rnorm(n=N, mean=0, sd=1)
  }
  
  if(model == "uniform")
  {
    data_sim <- runif(n=N, min = 0, max = sqrt(12))
  }
  
  if(model == "AR1")
  {
    data_sim <- arima.sim(model=list(ar=AR),n=N)
  }
  
  ### Introduce missingness
  
  if(method != "full")
  {
    miss <- sample(N, missing*N)
    data_sim[miss] <- NA
  }
  
  return(data_sim)
}

### Check if data for any treatment level is completely missing

Check_data_insufficient <- function(data_input, num_reps)
{
  checkall <- numeric(num_reps)
  for(i in 1:num_reps)
  {
    check <- !is.na(data_input[, (i*2)])
    check.a <- check[data_input[, (i*2-1)] == "A"]
    check.b <- check[data_input[, (i*2-1)] == "B"]
    checkall[i] <- sum(check.a) * sum(check.b)
  }
  
  return(sum(checkall) == 0)
}

### Estimate power and type I error rate

Calculate_conditional_power_random <- function(
  design,        # SCE design type
  model,         # Data model
  number,        # Number of randomizations per simulated dataset (>1 only if calcualting conditional power)
  ESM,           # Test statistic
  limit_phase,   # Minimum number of measurements in a phase
  reps_MBD,      # Number of participants in MBD
  num_MI,        # Number of imputations in multiple imputation
  direction="+", # Direction of test statistic (Only used if randomization test is one-sided)
  alfa,          # Level of significance
  N,             # Number of measurement
  ES,            # Effect size
  AR,            # Autocorrelation
  number_MC,     # Number of randomizations in Monte Carlo randomization test
  method,        # Missing data handling method
  missing=0      # Proportion of missing data 
)
{
  ### Simulate dataset
  
  num_reps <- if(design == "MBD") reps_MBD else 1
  data_obs <- numeric(num_reps*N)
  
  for(i in 1:num_reps)
  {
    data_obs[((i-1)*N+1):(i*N)] <- Generate_data(model=model, AR=AR, N=N, method=method, missing=missing)
  }
  
  ABAB_idx <- matrix()
  
  ### Randomly select randomization for RBD
  
  if (design == "RBD")
  {
    MT <- N
    ab <- c("A", "B")
    ba <- c("B", "A")
    assignments <- matrix(nrow = number, ncol = MT)
    
    for(i in 1:number)
    {
      assignment <- character(MT) 
      rand <- runif(MT/2)
      
      for(j in 1:(MT/2))
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
    
    MT<-N
    quantity<-choose(MT-4*limit_phase+3,3)
    selection<-sample(1:quantity,number,replace=TRUE)
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
    
    ### Select assignment from list of all possible randomizations
    
    assignments <- matrix(nrow=number, ncol=MT)
    for(it in 1:number)
    {
      assignments[it, ] <- c(
        rep("A", index.a1[selection[it]]),
        rep("B", index.b1[selection[it]] - index.a1[selection[it]]),
        rep("A", index.a2[selection[it]] - index.b1[selection[it]]),
        rep("B", MT - index.a2[selection[it]])
      )
    }
    
    ### Save all possible randomizations for ABAB
    ABAB_idx <- matrix(c(index.a1, index.b1, index.a2), ncol = 3, byrow = FALSE)
  }
  
  ### Randomly select randomization for MBD
  
  if(design == "MBD")
  {
    MT <- N
    startpoints <- (limit_phase + 1):(MT - limit_phase + 1)
    assignments <- matrix(nrow=number, ncol=MT*reps_MBD)
    
    for(i in 1:number)
    {
      combstartpts <- sample(startpoints, reps_MBD)
      ass <- numeric()
      
      for(it in 1:reps_MBD)
      {
        ass <- c(ass, rep("A", combstartpts[it]-1), rep("B", MT - combstartpts[it] + 1))
      }
      
      assignments[i, ] <- ass
    }
  }
  
  ### Calculate p-value and whether H0 is rejected for simulated dataset
  
  count_rejections <- numeric()
  pvalues <- numeric()
  
  for(i in 1:number)
  {
    ### Add effect size
    
    data_shift <- data_obs
    labels_obs <- assignments[i,] 
    data_shift[labels_obs == "B"] <- data_obs[labels_obs == "B"] + ES 
    
    if(design == "MBD")
    {
      data_input <- data.frame(1:N)
      for(j in 1:reps_MBD)
      {
        data_input[paste("Labs", j, sep = "")] <- labels_obs[((j - 1) * N + 1):(j * N)]
        data_input[paste("Vals", j, sep = "")] <- data_shift[((j - 1) * N + 1):(j * N)]
      }
      data_input <- data_input[, -1]
    } else
    {
      data_input <- data.frame(Labs = labels_obs, Vals = data_shift)
    }
    
    
    if(method != "full" && Check_data_insufficient(data_input = data_input, num_reps = num_reps))
    {
      ### Reject H0 is data entirely missing for a treatment level
      
      count_rejections[i] <- 0
      pvalues[i] <- NA
      #print("missing statistic!")
    } else
    {
      ### Dispatch data to randomization test function
      
      output <- ESM.random.RT(
        data = data_input, 
        design = design, 
        number = number_MC, 
        ESM = ESM, 
        method = method,
        limit_phase = limit_phase,
        reps_MBD = reps_MBD,
        num_MI = num_MI,
        direction = direction, 
        alfa = alfa,
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

