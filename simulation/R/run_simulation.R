run_simulation <- function(scenario, equal, tau2, p.baseline,
                           min.grpsize, max.grpsize,
                           M = 1000, seed) {
  
  ##
  ##
  ## (1) Simulation settings
  ##
  ##
  
  ## Treatment effects
  ##
  if (equal)
    OR.A <- OR.B <- OR.C <- OR.D <- 1
  else {
    OR.A <- 1.40
    OR.B <- 1.20
    OR.C <- 2.30
    OR.D <- 1.50
  }
  ##
  ## Scenarios
  ##
  I.AB <- I.AC <- I.CD <- 1
  ##
  if (scenario == "B1") {
    I.AB <- 1.5
  }
  else if (scenario == "B2") {
    I.CD <- 1.5
  }
  else if (scenario == "C1") {
    I.AB <- 2.0
  }
  else if (scenario == "C2") {
    I.CD <- 2.0
  }
  ##
  ## Effects of multicomponent interventions
  ##
  OR.AB <- OR.A * OR.B * I.AB
  OR.AC <- OR.A * OR.C * I.AC
  OR.CD <- OR.C * OR.D * I.CD
  ##
  ## Treatment labels etc.
  ##
  trts <- c("P", "A", "B", "C", "D", "A+B", "A+C", "C+D")
  reference.group <- "P"
  ##
  ## Vector with treatment effects
  ##
  ORs = c(OR.A, OR.B, OR.C, OR.D, OR.AB, OR.AC, OR.CD)
  names(ORs) <- trts[-1]


  ## Create directory for data and save settings
  ##
  dir <- simpath(scenario, equal, tau2, p.baseline, min.grpsize, max.grpsize)
  dir_create(dir)
  ##
  sm <- "OR"
  ##
  save(scenario, equal, tau2, p.baseline,
       min.grpsize, max.grpsize,
       M, seed,
       sm, ORs, reference.group,
       I.AB, I.AC, I.CD,
       file = paste0(dir, "settings.rda"))
  
  
  ##
  ##
  ## (2) Generate simulation data sets
  ##
  ##
  
  set.seed(seed)
  ##
  simdata <- list()
  ##
  for (i in seq_len(M)) {
    simdata.i <-
      sim.dataset(ORs, tau2, p.baseline, min.grpsize, max.grpsize,
                  reference.group, sm)
    ##
    attr(simdata.i, "scenario") <- scenario
    attr(simdata.i, "equal") <- equal
    attr(simdata.i, "iteration") <- i
    attr(simdata.i, "M") <- M
    ##
    simdata[[i]] <- simdata.i
  }
  ##
  save(simdata, file = paste0(dir, "simdata.rda"))
  
  
  ##
  ##
  ## (3) Run the simulation
  ##
  ##
  
  ## Run simulations
  ##
  lapply(simdata, sim.cnma)
  
  if (FALSE) {
    
    ## Parallel computing
    ##
    library(parallel)
    ##
    mc.cores <- max(1, detectCores()) # number of cores
    
    if (.Platform$OS.type == "unix") { # macOS or Linux
      mclapply(simdata, sim.cnma, mc.cores = mc.cores)
    }
    else { # MS Windows
      cl <- makeCluster(mc.cores)
      ##
      clusterExport(cl = cl,
                    varlist = c("simdata", "simpath",
                                "sim.cnma",
                                "dir_create",
                                "netmeta", "netcomb", "discomb", "modselection",
                                "gs", "netconnection", "combinations", "selcomp",
                                "selint", "eval.ms","netgraph", "sim.networks",
                                "disconnect", "sim.events", "sim.ms.disc",
                                "sim.ms.disconnected"),
                    envir = environment())
      ##
      ## create P initial subsets with parallel computations
      ##
      paral <- parLapply(cl, simdata, sim.cnma)
      ##
      stopCluster(cl)
    }
  }
  
  
  invisible(NULL)
}
