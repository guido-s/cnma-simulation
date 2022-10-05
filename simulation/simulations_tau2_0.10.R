library(netmeta)


## Load functions for running the simulation
##
source("R/funcs.R")
source("R/disconnect.R")
source("R/disconnect_additional.R")
source("R/fwdselection.R")
source("R/run_simulation.R")
source("R/sim-funcs.R")


## Scenarios:
## A  = additive effects
## B1 = additivity of A+B mildly violated
## B2 = additivity of C+D mildly violated
## C1 = additivity of A+B strongly violated
## C2 = additivity of C+D strongly violated
##
scenario <- c("A", "B1", "B2", "C1", "C2")


## Simulation settings
##
M <- 1000 # number of generated data sets
##
tau2 <- 0.1    # moderate heterogeneity
equal <- FALSE # equal treatment effects? (i.e., all OR = 1)
p.baseline <- 0.1
min.grpsize <- 50
max.grpsize <- 200

set.seed(1910)
##
seed <- sample(1:10000, length(scenario))
seed


for (i in seq_along(scenario))
  run_simulation(scenario[i], equal, tau2, p.baseline,
                 min.grpsize, max.grpsize, M, seed[i])
