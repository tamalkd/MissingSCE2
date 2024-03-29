#########################
# This R file only contains some library functions.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

source("RT.R")

### Generate unique seed for every set of simulation conditions

Generate_seed <- function(design, model, ESM, ES, N, method, missprop, misstype)
{
  designs <- c("RBD", "ABAB", "MBD")
  models <- c("AR.6", "normal", "uniform", "mvn.3", "mvn.6", "AR.3")
  ESMs <- c("MD", "NAP")
  ESs <- c(0, 1, 2)
  Ns <- c(20, 30, 40)
  methods <- c("full", "marker", "MI")
  missprops <- c(0, 0.1, 0.3, 0.5)
  misstypes <- c("none", "trunc+", "trunc-", "mvn+", "mvn-", "mcar")
  
  idx_design <- match(design, designs)
  idx_model <- match(model, models)
  idx_ESM <- match(ESM, ESMs)
  idx_ES <- match(ES, ESs)
  idx_N <- match(N, Ns)
  idx_method <- match(method, methods)
  idx_missprop <- match(missprop, missprops)
  idx_misstype <- match(misstype, misstypes)
  
  seed <- idx_design + idx_model * 10 + idx_ESM * 100 + idx_ES * 1e3 + 
    idx_N * 1e4 + idx_method * 1e5 + idx_missprop * 1e6 + idx_misstype * 1e7
  
  return(seed)
}

### Generate simulation data

Generate_data <- function(model, N) 
{
  data_sim <- switch(model,
    "normal" = rnorm(n = N, mean = 0, sd = 1),
    "uniform" = runif(n = N, min = 0, max = sqrt(12)),
    "AR.3" = arima.sim(model = list(ar = 0.3), n = N),
    "AR.6" = arima.sim(model = list(ar = 0.6), n = N),
    "mvn.3" = matrix(rnorm(n = 3 * N, mean = 0, sd = 1), ncol = 3),
    "mvn.6" = matrix(rnorm(n = 3 * N, mean = 0, sd = 1), ncol = 3),
    stop("Unknown model: Cannot generate data!")
  )
  
  data_sim <- as.data.frame(data_sim)
  return(data_sim)
}

### Introduce missingness

Add_missing <- function(data, N, missprop, misstype)
{
  if(misstype %in% c("trunc+", "trunc-", "mvn+", "mvn-")) # MNAR and MAR mechanisms
  {
    top_trunc <- (misstype %in% c("trunc+", "mvn+"))
    col <- if(misstype %in% c("mvn+", "mvn-")) 2 else 1
    
    sorted <- sort(data[, col], decreasing = top_trunc)
    trunc <- sorted[missprop * N]
    miss <- if(top_trunc) (data[, col] >= trunc) else (data[, col] <= trunc)
    
    data[miss, 1] <- NA
  } else
  if(misstype == "mcar") # MCAR mechanism
  {
    miss <- sample(N, missprop * N)
    data[miss, 1] <- NA
  }
  
  return(data)
}

### Adjust correlations of observed data and covariates for multivariate normal models

Adjust_mvn_corr <- function(data, model)
{
  corr = switch(model,
    "mvn.3" = 0.3,
    "mvn.6" = 0.6,
    stop("Unknown model: Cannot calculate correlation!")
  )
  
  norm <- data[, 1]
  norm <- (norm - mean(norm)) / sd(norm)
  coef <- corr / sqrt(1 - (corr ^ 2))
  
  data[, 2] <- (coef * norm + data[, 2]) / sqrt(1 + (coef ^ 2))
  data[, 3] <- (coef * norm + data[, 3]) / sqrt(1 + (coef ^ 2))
  
  return(data)
}

### Check if data for any treatment level is completely missing

Check_data_insufficient <- function(data_list, nReps)
{
  checkall <- numeric(nReps)
  for(i in 1:nReps)
  {
    data_input <- data_list[[i]]
    check <- !is.na(data_input[, 2])
    check.a <- check[data_input[, 1] == "A"]
    check.b <- check[data_input[, 1] == "B"]
    checkall[i] <- sum(check.a) * sum(check.b)
  }
  
  return(sum(checkall) == 0)
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
  misstype = "",   # Missing data mechanism
  alfa,            # Level of significance
  direction = "+", # Direction of test statistic (Only used if randomization test is one-sided)
  limit_phase,     # Minimum number of measurements in a phase
  nCP,             # Number of randomizations per simulated dataset (>1 only if calculating conditional power)
  nMBD,            # Number of participants in MBD
  nMC,             # Number of randomizations in Monte Carlo randomization test
  nMI              # Number of imputations in multiple imputation
)
{
  ### Define matrices for assignments
  
  assignments <- list()
  ABAB_idx <- matrix()
  
  ### Randomly select randomization for RBD
  
  if (design == "RBD")
  {
    ab <- c("A", "B")
    ba <- c("B", "A")
    
    for(i in 1:nCP)
    {
      assignment <- character(N) 
      rand <- runif(N/2)
      
      for(j in 1:(N/2))
      {
        index <- if(rand[j] < 0.5) ab else ba
        assignment[(j*2-1):(j*2)] <- index
      }
      assignments[[i]] <- matrix(assignment, ncol = 1)
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
    
    for(it in 1:nCP)
    {
      assignments[[it]] <- matrix(
        c(
          rep("A", index.a1[selection[it]]),
          rep("B", index.b1[selection[it]] - index.a1[selection[it]]),
          rep("A", index.a2[selection[it]] - index.b1[selection[it]]),
          rep("B", N - index.a2[selection[it]])
        ), 
        ncol = 1
      )
    }
    
    ### Save all possible randomizations for ABAB
    ABAB_idx <- matrix(c(index.a1, index.b1, index.a2), ncol = 3, byrow = FALSE)
  }
  
  ### Randomly select randomization for MBD
  
  if(design == "MBD")
  {
    startpoints <- (limit_phase + 1):(N - limit_phase + 1)
    asgn_matrix <- matrix(nrow = N, ncol = nMBD)
    
    for(i in 1:nCP)
    {
      combstartpts <- sample(startpoints, nMBD)
      for(j in 1:nMBD)
      {
        asgn_matrix[,j] <- c(rep("A", combstartpts[j] - 1), rep("B", N - combstartpts[j] + 1))
      }
      
      assignments[[i]] <- asgn_matrix
    }
  }
  
  ### Simulate dataset
  
  nReps <- if(design == "MBD") nMBD else 1
  data_obs <- list()
  
  for(i in 1:nReps)
  {
    data_obs[[i]] <- Generate_data(model = model, N = N)
  }
  
  ### Calculate p-value and whether H0 is rejected for simulated dataset
  
  count_rejections <- numeric(nCP)
  pvalues <- numeric(nCP)
  
  for(i in 1:nCP)
  {
    labels <- assignments[[i]] 
    data_list <- list()
    
    for(j in 1:nReps)
    {
      ### Add effect size
      
      label_col <- labels[,j]
      data_shift <- data_obs[[j]]
      data_shift[label_col == "B", 1] <- data_shift[label_col == "B", 1] + ES
      
      ### Add missingness and create dataset
      
      if(model %in% c("mvn.3", "mvn.6"))
      {
        data_shift <- Adjust_mvn_corr(data = data_shift, model = model)
      }
      if(method != "full")
      {
        data_shift <- Add_missing(data = data_shift, N = N, missprop = missprop, misstype = misstype)
      }
      data_list[[j]] <- cbind(as.data.frame(label_col, stringsAsFactors = FALSE), data_shift)
    }
    
    ### Calculate randomization test result
    
    if(method != "full" && Check_data_insufficient(data_list = data_list, nReps = nReps))
    {
      ### Reject H0 is data entirely missing for a treatment level
      
      count_rejections[i] <- 0
      pvalues[i] <- NA
      
      #print("missing statistic!")
    } else
    {
      ### Dispatch data to randomization test function
      
      output <- Compute_RT(
        data_list = data_list, 
        design = design, 
        model = model,
        ESM = ESM, 
        method = method,
        alfa = alfa,
        direction = direction, 
        limit_phase = limit_phase,
        nMBD = nMBD,
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
