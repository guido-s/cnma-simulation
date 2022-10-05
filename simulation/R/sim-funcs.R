sim.cnma <- function(x) {
  
  scenario <- attr(x, "scenario")
  equal <- attr(x, "equal")
  tau2 <- attr(x, "tau2")
  p.baseline <- attr(x, "p.baseline")
  min.grpsize <- attr(x, "min.grpsize")
  max.grpsize <- attr(x, "max.grpsize")
  iteration <- attr(x, "iteration")
  M <- attr(x, "M")
  ##
  reference.group <- attr(x, "reference.group")
  sm <- attr(x, "sm")
  
  
  ##
  ##
  ## (1) Create directories for simulation scenario
  ##
  ##
  
  ## Directory for simulation results
  ##
  dir <-
    simpath(scenario, equal, tau2,
            p.baseline, min.grpsize, max.grpsize,
            iteration, M)
  dir_create(dir)
  
  
  ##
  ##
  ## (2) Connected network
  ##
  ##
  
  ## Standard network meta-analysis
  ##
  NMA.full <-
    netmeta(TE, seTE, treat1, treat2, studlab, data = x,
            reference.group = reference.group, sm = sm,
            fixed = FALSE)
  
  
  ## Additive model for full connected network ####
  ##
  aCNMA.full <- netcomb(NMA.full)

  ## Model selection for the full network
  ##
  ms.full <- fwdselection(aCNMA.full)
  ##
  save(NMA.full, aCNMA.full, ms.full,
       file = paste0(dir, "full-modsel.rda"))
  
  
  ##
  ##
  ## (3) Disconnected networks
  ##
  ##
  
  ## Simulate disconnected networks
  ##
  ## No unique treatment comparisons with the reference treatment
  ##
  discdata <-
    disconnect_additional(TE, seTE, treat1, treat2, studlab, data = x,
                          minimal.set = reference.group)
  ##
  ## Randomly select a disconnected network from all possible
  ## disconnected networks
  ##
  sel.disc <- sample(seq_len(length(discdata$set)), 1)
  ##
  discdata$data <- discdata$data[[sel.disc]]
  discdata$set <- discdata$set[[sel.disc]]
  discdata$fulldata <- discdata$fulldata[[sel.disc]]
  ##

  aCNMA.disc <- 
    discomb(TE, seTE, treat1, treat2, studlab,
            data = discdata$data, sm = sm)
  ##
  ms.disc <- fwdselection(aCNMA.disc)
  ##  
  save(discdata, aCNMA.disc, ms.disc, sel.disc,
       file = paste0(dir, "disc-modsel.rda"))
  
  
  invisible(NULL)
}
## Function to simulate dataset

sim.dataset <- function(ORs,
                        tau2, p.baseline, min.grpsize, max.grpsize,
                        reference.group, sm) {

  ## Construct network geometry
  ##
  treat1 <- c("P", "P", "P", "P", "P", "P", "P", "A", "A", "A", "B",
              "B", "B", "B", "C", "C", "C", "C", "D", "D", "D", "A+B",
              "A+B", "A+C", "P", "P", "P", "P")
  treat2 <- c("A", "B", "C", "D", "A+B", "A+C", "C+D", "C", "D", "A+C",
              "C", "D", "A+C", "C+D", "D", "A+B", "A+C", "C+D", "A+B",
              "A+C", "C+D", "A+C", "C+D", "C+D", "A", "B", "A+B", "A+C")
  ##
  k <- length(treat1)
  n <- length(unique(c(treat1, treat2)))


  ## studlab (no multi-arm studies)
  ##
  studlab <- 1:k


  ## Add placebo (reference.group) to ORs
  ##
  ORs <- c(1, ORs)
  names(ORs)[1] <- reference.group


  ## Sample size arm 1: discrete uniform
  ##
  n1 <- sample(seq(min.grpsize, max.grpsize), k, replace = TRUE)
  ##
  ## Equal sample size in arm 2
  ##
  n2 <- n1


  ## Simulate event1 and event2
  ##
  event1 <- rep(NA, k)
  event2 <- rep(NA, k)
  ##
  ##i="A"
  for (i in names(ORs)) {
    sel1 <- treat1 == i
    sel2 <- treat2 == i
    ##
    if (sum(sel1) > 0)
      event1[sel1] <- sim.events(n1[sel1], p.baseline, ORs[i])
    ##
    if (sum(sel2) > 0)
      event2[sel2] <- sim.events(n2[sel2], p.baseline, ORs[i], tau2)
  }


  simdata <- data.frame(studlab, event1, n1, treat1, event2, n2, treat2)
  ##
  res <- pairwise(list(treat1, treat2),
                  event = list(event1, event2),
                  n = list(n1, n2),
                  data = simdata,
                  sm = sm, reference.group = reference.group)
  ##
  attr(res, "ORs") <- ORs
  attr(res, "tau2") <- tau2
  attr(res, "p.baseline") <- p.baseline
  attr(res, "min.grpsize") <- min.grpsize
  attr(res, "max.grpsize") <- max.grpsize
  attr(res, "reference.group") <- reference.group

  res
}
## function for the events simulation

sim.events <- function(size, p.baseline, effect, tau2 = 0) {

  ## study-specific log-odds
  lnOR <- rnorm(length(size), log(effect), sqrt(tau2))

  ## study-specific baseline probabilities
  prob = (p.baseline * exp(lnOR)) / (1 - p.baseline * (1 - exp(lnOR)))

  r <- rbinom(length(size), size, prob)

  r
}
