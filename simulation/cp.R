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

cp_A_0.00 <-
  read_cp_full("A", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B1_0.00 <-
  read_cp_full("B1", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B2_0.00 <-
  read_cp_full("B2", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C1_0.00 <-
  read_cp_full("C1", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C2_0.00 <-
  read_cp_full("C2", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_NMA_0.00 <-
  rbind(cp_A_0.00$NMA,
        cp_B1_0.00$NMA, cp_B2_0.00$NMA,
        cp_C1_0.00$NMA, cp_C2_0.00$NMA)
##
cp_aCNMA_0.00 <-
  rbind(cp_A_0.00$aCNMA,
        cp_B1_0.00$aCNMA, cp_B2_0.00$aCNMA,
        cp_C1_0.00$aCNMA, cp_C2_0.00$aCNMA)
##
cp_iCNMA_0.00 <-
  rbind(cp_A_0.00$iCNMA,
        cp_B1_0.00$iCNMA, cp_B2_0.00$iCNMA,
        cp_C1_0.00$iCNMA, cp_C2_0.00$iCNMA)


##
##
## (2) Connected network, tau2 = 0.01
##
##

cp_A_0.01 <-
  read_cp_full("A", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B1_0.01 <-
  read_cp_full("B1", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B2_0.01 <-
  read_cp_full("B2", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C1_0.01 <-
  read_cp_full("C1", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C2_0.01 <-
  read_cp_full("C2", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_NMA_0.01 <-
  rbind(cp_A_0.01$NMA,
        cp_B1_0.01$NMA, cp_B2_0.01$NMA,
        cp_C1_0.01$NMA, cp_C2_0.01$NMA)
##
cp_aCNMA_0.01 <-
  rbind(cp_A_0.01$aCNMA,
        cp_B1_0.01$aCNMA, cp_B2_0.01$aCNMA,
        cp_C1_0.01$aCNMA, cp_C2_0.01$aCNMA)
##
cp_iCNMA_0.01 <-
  rbind(cp_A_0.01$iCNMA,
        cp_B1_0.01$iCNMA, cp_B2_0.01$iCNMA,
        cp_C1_0.01$iCNMA, cp_C2_0.01$iCNMA)


##
##
## (3) Connected network, tau2 = 0.10
##
##

cp_A_0.10 <-
  read_cp_full("A", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B1_0.10 <-
  read_cp_full("B1", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B2_0.10 <-
  read_cp_full("B2", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C1_0.10 <-
  read_cp_full("C1", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C2_0.10 <-
  read_cp_full("C2", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_NMA_0.10 <-
  rbind(cp_A_0.10$NMA,
        cp_B1_0.10$NMA, cp_B2_0.10$NMA,
        cp_C1_0.10$NMA, cp_C2_0.10$NMA)
##
cp_aCNMA_0.10 <-
  rbind(cp_A_0.10$aCNMA,
        cp_B1_0.10$aCNMA, cp_B2_0.10$aCNMA,
        cp_C1_0.10$aCNMA, cp_C2_0.10$aCNMA)
##
cp_iCNMA_0.10 <-
  rbind(cp_A_0.10$iCNMA,
        cp_B1_0.10$iCNMA, cp_B2_0.10$iCNMA,
        cp_C1_0.10$iCNMA, cp_C2_0.10$iCNMA)


##
##
## (4) Connected network, CP
##
##

cp.full.data <-
  rbind(cp_NMA_0.00, cp_NMA_0.01, cp_NMA_0.10,
        cp_aCNMA_0.00, cp_aCNMA_0.01, cp_aCNMA_0.10,
        cp_iCNMA_0.00, cp_iCNMA_0.01, cp_iCNMA_0.10)
##
cp.full.long <- tolong(cp.full.data, varname = "CP")
##
cp.full <- by_full(cp.full.data, "CP")


##
##
## (5) Disconnected network, tau2 = 0.00
##
##

cp_A_0.00 <-
  read_cp_disc("A", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B1_0.00 <-
  read_cp_disc("B1", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B2_0.00 <-
  read_cp_disc("B2", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C1_0.00 <-
  read_cp_disc("C1", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C2_0.00 <-
  read_cp_disc("C2", equal, "0", p.baseline, min.grpsize, max.grpsize, M)
##
cp_aCNMA_0.00 <-
  rbind(cp_A_0.00$aCNMA,
        cp_B1_0.00$aCNMA, cp_B2_0.00$aCNMA,
        cp_C1_0.00$aCNMA, cp_C2_0.00$aCNMA)
##
cp_iCNMA_0.00 <-
  rbind(cp_A_0.00$iCNMA,
        cp_B1_0.00$iCNMA, cp_B2_0.00$iCNMA,
        cp_C1_0.00$iCNMA, cp_C2_0.00$iCNMA)


##
##
## (6) Disconnected network, tau2 = 0.01
##
##

cp_A_0.01 <-
  read_cp_disc("A", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B1_0.01 <-
  read_cp_disc("B1", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B2_0.01 <-
  read_cp_disc("B2", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C1_0.01 <-
  read_cp_disc("C1", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C2_0.01 <-
  read_cp_disc("C2", equal, "0.01", p.baseline, min.grpsize, max.grpsize, M)
##
cp_aCNMA_0.01 <-
  rbind(cp_A_0.01$aCNMA,
        cp_B1_0.01$aCNMA, cp_B2_0.01$aCNMA,
        cp_C1_0.01$aCNMA, cp_C2_0.01$aCNMA)
##
cp_iCNMA_0.01 <-
  rbind(cp_A_0.01$iCNMA,
        cp_B1_0.01$iCNMA, cp_B2_0.01$iCNMA,
        cp_C1_0.01$iCNMA, cp_C2_0.01$iCNMA)


##
##
## (7) Disconnected network, tau2 = 0.10
##
##

cp_A_0.10 <-
  read_cp_disc("A", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B1_0.10 <-
  read_cp_disc("B1", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_B2_0.10 <-
  read_cp_disc("B2", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C1_0.10 <-
  read_cp_disc("C1", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_C2_0.10 <-
  read_cp_disc("C2", equal, "0.1", p.baseline, min.grpsize, max.grpsize, M)
##
cp_aCNMA_0.10 <-
  rbind(cp_A_0.10$aCNMA,
        cp_B1_0.10$aCNMA, cp_B2_0.10$aCNMA,
        cp_C1_0.10$aCNMA, cp_C2_0.10$aCNMA)
##
cp_iCNMA_0.10 <-
  rbind(cp_A_0.10$iCNMA,
        cp_B1_0.10$iCNMA, cp_B2_0.10$iCNMA,
        cp_C1_0.10$iCNMA, cp_C2_0.10$iCNMA)


##
##
## (8) Disconnected network, CP
##
##

cp.disc.data <-
  rbind(cp_aCNMA_0.00, cp_aCNMA_0.01, cp_aCNMA_0.10,
        cp_iCNMA_0.00, cp_iCNMA_0.01, cp_iCNMA_0.10)
##
cp.disc.long <- tolong(cp.disc.data, varname = "CP")
##
cp.disc <- by_disc(cp.disc.data, "CP")


##
##
## (9) Save CP
##
##

save(cp.full, cp.full.data, cp.full.long,
     cp.disc, cp.disc.data, cp.disc.long,
     file = "results/cp.rda")
