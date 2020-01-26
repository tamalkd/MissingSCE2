#########################
# This code is only intended for executing the simulation on supercomputers.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

### Import simulation conditions

args <- commandArgs(TRUE)

print(args)

design <- as.character(args[1])   # SCE design type
model <- as.character(args[2])    # Data model
ESM <- as.character(args[3])      # Test statistic
ES <- as.integer(args[4])         # Effect size
N <- as.integer(args[5])          # Number of measurement
method <- as.character(args[6])   # Missing data handling method
missing <- as.numeric(args[7])    # Proportion of missing data 

options(warn=-1)

### Other parameters

AR = 0.6                    # Autocorrelation
limit_phase = 3             # Minimum number of measurements in a phase
reps_MBD = 4                # Number of participants in MBD
num_MI = 10                 # Number of imputations in multiple imputation
number = 1                  # Number of randomizations per simulated dataset (>1 only if calcualting conditional power)
number_MC = 1000            # Number of randomizations in Monte Carlo randomization test
alfa = 0.05                 # Level of significance
replications = 100000       # Number of simulated datasets
direction = "+"             # Direction of test statistic (Only used if randomization test is one-sided)

### Run simulation

strt <- Sys.time()

source("RT.R")
source("power_calc.R")

set.seed(1000)  # Set random seed to make results exactly reproducible

output <- numeric()
for(it in 1:replications)
{
  output[it] <- Calculate_conditional_power_random(
    design=design, 
    model=model, 
    number=number, 
    ESM=ESM, 
    limit_phase=limit_phase,
    reps_MBD=reps_MBD,
    num_MI=num_MI,
    direction=direction, 
    alfa=alfa, 
    N=N, 
    ES=ES, 
    AR=AR,
    number_MC=number_MC,
    method=method,
    missing=missing
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
  missing = missing,
  Reps = replications,
  nMC = number_MC,
  nCP = number,
  timer = difftime(end, strt, units = "secs"),
  power = power,
  stringsAsFactors = FALSE
)

### Save results

print(Result_table)
write.table(Result_table, "Results", append = TRUE, row.names = FALSE, col.names = FALSE)
