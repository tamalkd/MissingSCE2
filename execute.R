#########################
# This code is only intended for executing the simulation on supercomputers.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

### Import simulation conditions

args <- commandArgs(TRUE)

print(args)

design <- as.character(args[1])      # SCE design type
model <- as.character(args[2])       # Data model
ESM <- as.character(args[3])         # Test statistic
ES <- as.integer(args[4])            # Effect size
N <- as.integer(args[5])             # Number of measurement
method <- as.character(args[6])      # Missing data handling method
missprop <- as.numeric(args[7])      # Proportion of missing data 
misstype <- as.character(args[8])    # Missing data mechanism

options(warn=-1)

### Other parameters

alfa <- 0.05          # Level of significance
direction <- "+"      # Direction of test statistic (Only used if test statistic is one-sided)
limit_phase <- 3      # Minimum number of measurements in a phase
nCP <- 1              # Number of assignments per simulated dataset (>1 only if calculating conditional power)
nMBD <- 4             # Number of participants in MBD
nMC <- 1000           # Number of randomizations in Monte Carlo randomization test
nMI <- 10             # Number of imputations in multiple imputation
replications <- 10000 # Number of simulated datasets

### Run simulation

strt <- Sys.time()

source("RT.R")
source("power.R")

seed <- Generate_seed(
  design = design,
  model = model,
  ESM = ESM,
  ES = ES,
  N = N,
  method = method,
  missprop = missprop,
  misstype = misstype
)

set.seed(seed)  # Set unique seed to make results exactly reproducible

output <- numeric(replications)
for(it in 1:replications)
{
  output[it] <- Calculate_power_RT(
    design = design,
    model = model,
    ESM = ESM,
    ES = ES,
    N = N,
    method = method,
    missprop = missprop,
    misstype = misstype,
    alfa = alfa,
    direction = direction,
    limit_phase = limit_phase,
    nCP = nCP,
    nMBD = nMBD,
    nMC = nMC,
    nMI = nMI
  )
}

power <- mean(output)
end <- Sys.time()

Result_table <- data.frame(
  design = design, 
  model = model, 
  ESM = ESM, 
  ES = ES, 
  N = N, 
  method = method,
  missprop = missprop,
  misstype = misstype,
  nMC = nMC,
  Reps = replications,
  timer = as.numeric(difftime(end, strt, units = "secs")),
  power = power,
  stringsAsFactors = FALSE
)

### Save results

print(Result_table)
write.table(Result_table, "Results", append = TRUE, row.names = FALSE, col.names = FALSE)
