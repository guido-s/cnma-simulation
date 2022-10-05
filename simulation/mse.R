source("R/funcs.R")


M <- 1000
equal <- FALSE
p.baseline <- 0.1
min.grpsize <- 50
max.grpsize <- 200


##
##
## (1) Connected network, tau2 = 0.00
##
##

mse_A_0.00 <-
  read_mse_full("A", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B1_0.00 <-
  read_mse_full("B1", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B2_0.00 <-
  read_mse_full("B2", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C1_0.00 <-
  read_mse_full("C1", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C2_0.00 <-
  read_mse_full("C2", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_NMA_0.00 <-
  rbind(mse_A_0.00$NMA,
        mse_B1_0.00$NMA, mse_B2_0.00$NMA,
        mse_C1_0.00$NMA, mse_C2_0.00$NMA)
##
mse_aCNMA_0.00 <-
  rbind(mse_A_0.00$aCNMA,
        mse_B1_0.00$aCNMA, mse_B2_0.00$aCNMA,
        mse_C1_0.00$aCNMA, mse_C2_0.00$aCNMA)
##
mse_iCNMA_0.00 <-
  rbind(mse_A_0.00$iCNMA,
        mse_B1_0.00$iCNMA, mse_B2_0.00$iCNMA,
        mse_C1_0.00$iCNMA, mse_C2_0.00$iCNMA)


##
##
## (2) Connected network, tau2 = 0.01
##
##

mse_A_0.01 <-
  read_mse_full("A", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B1_0.01 <-
  read_mse_full("B1", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B2_0.01 <-
  read_mse_full("B2", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C1_0.01 <-
  read_mse_full("C1", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C2_0.01 <-
  read_mse_full("C2", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_NMA_0.01 <-
  rbind(mse_A_0.01$NMA,
        mse_B1_0.01$NMA, mse_B2_0.01$NMA,
        mse_C1_0.01$NMA, mse_C2_0.01$NMA)
##
mse_aCNMA_0.01 <-
  rbind(mse_A_0.01$aCNMA,
        mse_B1_0.01$aCNMA, mse_B2_0.01$aCNMA,
        mse_C1_0.01$aCNMA, mse_C2_0.01$aCNMA)
##
mse_iCNMA_0.01 <-
  rbind(mse_A_0.01$iCNMA,
        mse_B1_0.01$iCNMA, mse_B2_0.01$iCNMA,
        mse_C1_0.01$iCNMA, mse_C2_0.01$iCNMA)


##
##
## (3) Connected network, tau2 = 0.10
##
##

mse_A_0.10 <-
  read_mse_full("A", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B1_0.10 <-
  read_mse_full("B1", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B2_0.10 <-
  read_mse_full("B2", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C1_0.10 <-
  read_mse_full("C1", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C2_0.10 <-
  read_mse_full("C2", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_NMA_0.10 <-
  rbind(mse_A_0.10$NMA,
        mse_B1_0.10$NMA, mse_B2_0.10$NMA,
        mse_C1_0.10$NMA, mse_C2_0.10$NMA)
##
mse_aCNMA_0.10 <-
  rbind(mse_A_0.10$aCNMA,
        mse_B1_0.10$aCNMA, mse_B2_0.10$aCNMA,
        mse_C1_0.10$aCNMA, mse_C2_0.10$aCNMA)
##
mse_iCNMA_0.10 <-
  rbind(mse_A_0.10$iCNMA,
        mse_B1_0.10$iCNMA, mse_B2_0.10$iCNMA,
        mse_C1_0.10$iCNMA, mse_C2_0.10$iCNMA)


##
##
## (4) Connected network, MSE
##
##

mse.full.data <-
  rbind(mse_NMA_0.00, mse_NMA_0.01, mse_NMA_0.10,
        mse_aCNMA_0.00, mse_aCNMA_0.01, mse_aCNMA_0.10,
        mse_iCNMA_0.00, mse_iCNMA_0.01, mse_iCNMA_0.10)
##
mse.full.long <- tolong(mse.full.data, varname = "MSE")
##
mse.full <- by_full(mse.full.data, "MSE")


##
##
## (5) Disconnected network, tau2 = 0.00
##
##

mse_A_0.00 <-
  read_mse_disc("A", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B1_0.00 <-
  read_mse_disc("B1", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B2_0.00 <-
  read_mse_disc("B2", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C1_0.00 <-
  read_mse_disc("C1", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C2_0.00 <-
  read_mse_disc("C2", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
mse_aCNMA_0.00 <-
  rbind(mse_A_0.00$aCNMA,
        mse_B1_0.00$aCNMA, mse_B2_0.00$aCNMA,
        mse_C1_0.00$aCNMA, mse_C2_0.00$aCNMA)
##
mse_iCNMA_0.00 <-
  rbind(mse_A_0.00$iCNMA,
        mse_B1_0.00$iCNMA, mse_B2_0.00$iCNMA,
        mse_C1_0.00$iCNMA, mse_C2_0.00$iCNMA)


##
##
## (6) Disconnected network, tau2 = 0.01
##
##

mse_A_0.01 <-
  read_mse_disc("A", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B1_0.01 <-
  read_mse_disc("B1", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B2_0.01 <-
  read_mse_disc("B2", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C1_0.01 <-
  read_mse_disc("C1", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C2_0.01 <-
  read_mse_disc("C2", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
mse_aCNMA_0.01 <-
  rbind(mse_A_0.01$aCNMA,
        mse_B1_0.01$aCNMA, mse_B2_0.01$aCNMA,
        mse_C1_0.01$aCNMA, mse_C2_0.01$aCNMA)
##
mse_iCNMA_0.01 <-
  rbind(mse_A_0.01$iCNMA,
        mse_B1_0.01$iCNMA, mse_B2_0.01$iCNMA,
        mse_C1_0.01$iCNMA, mse_C2_0.01$iCNMA)


##
##
## (7) Disconnected network, tau2 = 0.10
##
##

mse_A_0.10 <-
  read_mse_disc("A", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B1_0.10 <-
  read_mse_disc("B1", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_B2_0.10 <-
  read_mse_disc("B2", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C1_0.10 <-
  read_mse_disc("C1", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_C2_0.10 <-
  read_mse_disc("C2", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
mse_aCNMA_0.10 <-
  rbind(mse_A_0.10$aCNMA,
        mse_B1_0.10$aCNMA, mse_B2_0.10$aCNMA,
        mse_C1_0.10$aCNMA, mse_C2_0.10$aCNMA)
##
mse_iCNMA_0.10 <-
  rbind(mse_A_0.10$iCNMA,
        mse_B1_0.10$iCNMA, mse_B2_0.10$iCNMA,
        mse_C1_0.10$iCNMA, mse_C2_0.10$iCNMA)


##
##
## (8) Disconnected network, MSE
##
##

mse.disc.data <-
  rbind(mse_aCNMA_0.00, mse_aCNMA_0.01, mse_aCNMA_0.10,
        mse_iCNMA_0.00, mse_iCNMA_0.01, mse_iCNMA_0.10)
##
mse.disc.long <- tolong(mse.disc.data, varname = "MSE")
##
mse.disc <- by_disc(mse.disc.data, "MSE")


##
##
## (9) Save MSE
##
##

save(mse.full, mse.full.data, mse.full.long,
     mse.disc, mse.disc.data, mse.disc.long,
     file = "results/mse.rda")
