library("netmeta")


keepall <- TRUE
verbose <- TRUE
##
inactive <- NULL


##
##
## (1) Load functions and data set
##
##

source("R/combinations.R")
source("R/fwdselection.R")
source("R/disconnect.R")
source("R/disconnect_additional.R")
source("R/sortdisc.R")
source("R/separate.R")
##
load("AEdata.rda")


##
##
## (2) Analyze the full network
##
##

## Standard network meta-analysis of treatment combinations
##
NMA.full <- netmeta(lnRR, selnRR, t1, t2, id, data = AEdata,
  sm = "RR", common = FALSE, reference.group = "plac")
##
## Additive model for full connected network
##
aCNMA.full <- netcomb(NMA.full, inactive = inactive)
##
## List of treatment combinations, i.e., potential interactions
##
combinations(aCNMA.full) # eleven combinations of two components
combinations(aCNMA.full, n = 3) # no combination of three components
##
## CNMA model selection in the full network
##
ms.full <-
  fwdselection(aCNMA.full, keepall = keepall, verbose = verbose)
##
## Save results for full network
##
save(NMA.full, aCNMA.full, ms.full,
  file = "results/full-modsel.rda")


##
##
## (3) Create disconnected networks
##
##

## Initial disconnected network consisting of minimal set
## (with placebo as reference intervention)
##
minset <- c("plac",
            "amis", "beta", "dola", "scop", "meto+scop",
            "dexa", "drop", "dexa+drop", "gran", "dexa+gran","drop+gran",
            "dexa+trop", "trop",
            "dexa+onda")
##
disc.init <-
  disconnect(lnRR, selnRR, t1, t2, id, data = AEdata,
    main.trts = minset)
##
## Generate additional disconnected networks
##
disc.add <-
  disconnect_additional(lnRR, selnRR, t1, t2, id, data = AEdata,
    minimal.set = minset, verbose = verbose)
##
## Sort disconnect networks by
## - decreasing number of pairwise comparisons
## - decreasing number of studies
## - decreasing number of pairwise comparisons in main subnetwork
##
AEdata.discs <- sortdisc(disc.init, disc.add, verbose = TRUE)
##
## Save datasets of disconnected networks
##
save(AEdata.discs, file = "results/discs-AEdata.rda")


##
##
## (4) Run CNMA model selection in disconnected networks
##
##

n.discs <- length(AEdata.discs$data)
aCNMA.discs <- ms.discs <- vector(mode = "list", length = n.discs)
##
for (i in seq_len(n.discs))
  aCNMA.discs[[i]] <-
    discomb(TE, seTE, treat1, treat2, studlab,
            data = AEdata.discs$data[[i]], sm = "RR",
            inactive = inactive)
warnings()
##
ms.discs <-
  lapply(aCNMA.discs, fwdselection, keepall = keepall, verbose = verbose)
warnings()
##
names(ms.discs) <- paste0("disc", seq_len(n.discs))
##
## Save results of model selection in disconnected networks
##
save(ms.discs, file = paste0("results/discs-modsel.rda"))


##
##
## (5) Separate NMA analyses of main subnetworks
##
##

sep.discs <- lapply(AEdata.discs$data, separate)
##
## Save results of separate analyses
##
save(sep.discs, file = "results/discs-sep.rda")
