#####################
# PLEASE READ CAREFULLY BEFORE EXECUTING!
#
# This code would normally be run on supercomputers. The entirety of the computation is expected to take 
# 20000 hours on a single core processor running at 2.5 GHz. 
#
# Please only run one condition at a time if running on a personal computer. This means only one choice 
# among 'designs' (single-case design), 'models' (simulation data model), 'ESMs' (effect size 
# measure/test statistic), 'ESs' (effect size applied), 'Ns' (number of measurements), 'methods' (missing 
# data handling method), and 'missprops' (missing data proportion) should be kept, all others should to 
# be deleted.
#
# We strongly recommend setting parameter 'replications' (number of simulated datasets) to 1000 or lower.
# This would lead to lower accuracy for power/type I error estimate, however otherwise the computation 
# time would be too long.
#
# PLease ensure the following files are in the same folder (or in the working directory):
# 'RT.R' 'power_calc.R' 'NAP.cpp'
#
# Please ensure the following R packages are installed:
# 'foreach' 'parallel' 'doParallel' 'imputeTS' 'data.table' 'mice' 'Rcpp' 
#####################

### All simulation conditions (for reference)

# designs <- c("RBD", "ABAB")                            # SCE design types
# models <- c("AR1", "normal", "uniform", "mvn")         # Data models
# ESMs <- c("MD", "NAP")                                 # Test statistics
# ESs <- c(0, 1, 2)                                      # Effect sizes
# Ns <- c(40, 30, 20)                                    # Number of measurements
# methods <- c("full", "marker", "TS", "MI")             # Missing data handling methods
# missprops <- c(0.1, 0.3, 0.5)                          # Proportion of missing data 
# misstypes <- c("censor+", "censor-", "mvn+", "mvn-")   # Mechanism used to generate missing data

### Simulation conditions to test

designs <- c("RBD", "ABAB")                            # SCE design types
models <- c("AR1", "normal", "uniform", "mvn")         # Data models
ESMs <- c("MD", "NAP")                                 # Test statistics
ESs <- c(0, 1, 2)                                      # Effect sizes
Ns <- c(40, 30, 20)                                    # Number of measurements
methods <- c("full", "marker", "TS", "MI")             # Missing data handling methods
missprops <- c(0.1, 0.3, 0.5)                          # Proportion of missing data
misstypes <- c("censor+", "censor-", "mvn+", "mvn-")   # Mechanism used to generate missing data

### Other parameters

alfa = 0.05                 # Level of significance
AR = 0.6                    # Autocorrelation
corr = 0.6                  # Correlation between bivariate normal vectors
direction = "+"             # Direction of test statistic (Only used if test statistic is one-sided)
limit_phase = 3             # Minimum number of measurements in a phase
nCP = 1                     # Number of assignments per simulated dataset (>1 only if calcualting conditional power)
nMC = 1000                  # Number of randomizations in Monte Carlo randomization test
nMI = 10                    # Number of imputations in multiple imputation
replications = 10000        # Number of simulated datasets

### Run simulations

library(foreach)
library(parallel)
library(doParallel)
library(Rcpp)

cores <- detectCores()
print(cores)

source("RT.R")
source("power.R")

Result_table <- data.frame(
  design = character(),
  model = character(),
  ESM = character(),
  ES = numeric(),
  N = integer(),
  method = character(),
  missprop = numeric(),
  misstype = character(),
  nMC = integer(),
  Reps = integer(),
  timer = numeric(),
  power = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(designs))
{
  for (m in 1:length(Ns)) 
  {
    for (k in 1:length(ESMs)) 
    {
      for (j in 1:length(models)) 
      {
        for (l in 1:length(ESs)) 
        {
          for(n in 1:length(methods))
          {
            if(methods[n]=="full")
            {
              outlist <- data.frame(
                design = designs[i], 
                model = models[j], 
                ESM = ESMs[k], 
                ES = ESs[l], 
                N = Ns[m], 
                method = methods[n],
                missprop = 0,
                misstype = "",
                nMC = nMC,
                Reps = replications,
                timer = NA,
                power = NA,
                stringsAsFactors = FALSE
              )
              Result_table <- rbind(Result_table, outlist)
            } else
            {
              for(o in 1:length(missprops))
              {
                for(p in 1: length(misstypes))
                {
                  if((models[j] == "mvn") && (misstypes[p] %in% c("mvn+", "mvn-"))
                     || (models[j] != "mvn") && !(misstypes[p] %in% c("mvn+", "mvn-")))
                  {
                    outlist <- data.frame(
                      design = designs[i], 
                      model = models[j], 
                      ESM = ESMs[k], 
                      ES = ESs[l], 
                      N = Ns[m], 
                      method = methods[n],
                      missprop = missprops[o],
                      misstype = misstypes[p],
                      nMC = nMC,
                      Reps = replications,
                      timer = NA,
                      power = NA,
                      stringsAsFactors = FALSE
                    )
                    Result_table <- rbind(Result_table, outlist)
                  }
                  
                }
              }
            }
            
          }
        }
      }
    }
  }
}
print(Result_table)

for(rnum in 1:nrow(Result_table))
{
  ### Run simulation with missing data
  
  start <- Sys.time()
  
  set.seed(1000)  # Set random seed to make results exactly reproducible
  row <- Result_table[rnum,]
  
  ### Run parallel
  
  # cl <- makeCluster(max(cores - 3, 1))
  # registerDoParallel(cl)
  # 
  # output <- foreach(
  #   it = 1:replications,
  #   .inorder = FALSE,
  #   .combine = 'c',
  #   .packages = c("MASS", "Rcpp", "data.table", "imputeTS", "mice"),
  #   .noexport = c("NAP_cpp")
  # ) %dopar%
  # {
  #   sourceCpp("NAP.cpp") # Explicitly source C++ script inside loop to satisfy foreach
  # 
  #   result <- Calculate_power_RT(
  #     design = row$design,
  #     model = row$model,
  #     ESM = row$ESM,
  #     ES = row$ES,
  #     N = row$N,
  #     method = row$method,
  #     missprop = row$missprop,
  #     misstype = row$misstype,
  #     alfa = alfa,
  #     AR = AR,
  #     corr = corr,
  #     direction = direction,
  #     limit_phase = limit_phase,
  #     nCP = nCP,
  #     nMC = nMC,
  #     nMI = nMI
  #   )
  # }
  # 
  # stopCluster(cl)
  
  ### Run without parallelization
  
  output <- numeric()
  for(it in 1:replications)
  {
    output[it] <- Calculate_power_RT(
      design = row$design,
      model = row$model,
      ESM = row$ESM,
      ES = row$ES,
      N = row$N,
      method = row$method,
      missprop = row$missprop,
      misstype = row$misstype,
      alfa = alfa,
      AR = AR,
      corr = corr,
      direction = direction,
      limit_phase = limit_phase,
      nCP = nCP,
      nMC = nMC,
      nMI = nMI
    )
  }
  
  power <- mean(output)
  end <- Sys.time()
  
  Result_table[rnum,"power"] <- power
  Result_table[rnum,"timer"] <- as.numeric(difftime(end, start, units = "secs"))
  print(power)
}

### Publish results

print(Result_table)
write.table(Result_table, "Results", append = TRUE, row.names = FALSE, col.names = FALSE)
