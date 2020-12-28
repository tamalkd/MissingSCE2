#########################
# This code is only intended for generating the simulation conditions for 
# running the simulations on supercomputers.
# For running the simulations on a personal computer, please refer to 'test.R'.
#########################

### All simulation conditions

designs <- c("RBD", "ABAB", "MBD")                                    # SCE design types
models <- c("AR.3", "AR.6", "normal", "uniform", "mvn.3", "mvn.6")    # Data models
ESMs <- c("MD", "NAP")                                                # Test statistics
ESs <- c(0, 1, 2)                                                     # Effect sizes
Ns <- c(40, 30, 20)                                                   # Number of measurements
methods <- c("full", "marker", "MI")                                  # Missing data handling methods
missprops <- c(0.1, 0.3, 0.5)                                         # Proportion of missing data
misstypes <- c("trunc+", "trunc-", "mvn+", "mvn-", "mcar")            # Missing data mechanism

### Generate all possible combinations of simulation conditions

outfile <- "start.csv"

Result_table <- data.frame(
  design = character(),
  model = character(),
  ESM = character(),
  ES = numeric(),
  N = integer(),
  method = character(),
  missprop = numeric(),
  misstype = character(),
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
                misstype = "none"
              )
              Result_table <- rbind(Result_table, outlist)
            } else
            {
              for(o in 1:length(missprops))
              {
                for(p in 1: length(misstypes))
                {
                  if((models[j] %in% c("mvn.3", "mvn.6")) || (misstypes[p] %in% c("trunc+", "trunc-")))
                  {
                    outlist <- data.frame(
                      design = designs[i], 
                      model = models[j], 
                      ESM = ESMs[k], 
                      ES = ESs[l], 
                      N = Ns[m], 
                      method = methods[n],
                      missprop = missprops[o],
                      misstype = misstypes[p]
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

### Publish table containing all possible combinations

print(Result_table)
write.csv(Result_table, outfile, row.names = FALSE)
