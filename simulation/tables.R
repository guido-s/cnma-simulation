library(gmodels)


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

## Scenario A
##
full_A_0.00 <-
  selected_full("A", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B1
##
full_B1_0.00 <-
  selected_full("B1", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B2
##
full_B2_0.00 <-
  selected_full("B2", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C1
##
full_C1_0.00 <-
  selected_full("C1", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C2
##
full_C2_0.00 <-
  selected_full("C2", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
full_0.00 <-
  rbind(full_A_0.00, full_B1_0.00, full_B2_0.00,
        full_C1_0.00, full_C2_0.00)


##
##
## (2) Connected network, tau2 = 0.01
##
##

## Scenario A
##
full_A_0.01 <-
  selected_full("A", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B1
##
full_B1_0.01 <-
  selected_full("B1", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B2
##
full_B2_0.01 <-
  selected_full("B2", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C1
##
full_C1_0.01 <-
  selected_full("C1", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C2
##
full_C2_0.01 <-
  selected_full("C2", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
full_0.01 <-
  rbind(full_A_0.01, full_B1_0.01, full_B2_0.01,
        full_C1_0.01, full_C2_0.01)


##
##
## (3) Connected network, tau2 = 0.10
##
##

## Scenario A
##
full_A_0.10 <-
  selected_full("A", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B1
##
full_B1_0.10 <-
  selected_full("B1", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B2
##
full_B2_0.10 <-
  selected_full("B2", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C1
##
full_C1_0.10 <-
  selected_full("C1", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C2
##
full_C2_0.10 <-
  selected_full("C2", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
full_0.10 <-
  rbind(full_A_0.10, full_B1_0.10, full_B2_0.10,
        full_C1_0.10, full_C2_0.10)


##
##
## (4) Table 3
##
##

sink("results/Table3.txt")
##
## n_diff = sum(pval.Q.diff < 0.05)
##
cat("*** tau2 = 0.00 ***\n")
with(full_0.00, table(scenario, pval.Q.diff < 0.05))
with(full_0.00, table(scenario, interactions))
##
cat("*** tau2 = 0.01 ***\n")
with(full_0.01, table(scenario, pval.Q.diff < 0.05))
with(full_0.01, table(scenario, interactions))
##
cat("*** tau2 = 0.10 ***\n")
with(full_0.10, table(scenario, pval.Q.diff < 0.05))
with(full_0.10, table(scenario, interactions))
##
sink()


##
##
## (5) Disconnected network, tau2 = 0.00
##
##

## Scenario A
##
disc_A_0.00 <-
  selected_disc("A", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B1
##
disc_B1_0.00 <-
  selected_disc("B1", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B2
##
disc_B2_0.00 <-
  selected_disc("B2", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C1
##
disc_C1_0.00 <-
  selected_disc("C1", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C2
##
disc_C2_0.00 <-
  selected_disc("C2", equal, "0", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
disc_0.00 <-
  rbind(disc_A_0.00, disc_B1_0.00, disc_B2_0.00,
        disc_C1_0.00, disc_C2_0.00)


##
##
## (6) Disconnected network, tau2 = 0.01
##
##

## Scenario A
##
disc_A_0.01 <-
  selected_disc("A", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B1
##
disc_B1_0.01 <-
  selected_disc("B1", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B2
##
disc_B2_0.01 <-
  selected_disc("B2", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C1
##
disc_C1_0.01 <-
  selected_disc("C1", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C2
##
disc_C2_0.01 <-
  selected_disc("C2", equal, "0.01", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
disc_0.01 <-
  rbind(disc_A_0.01, disc_B1_0.01, disc_B2_0.01,
        disc_C1_0.01, disc_C2_0.01)


##
##
## (7) Disconnected network, tau2 = 0.10
##
##

## Scenario A
##
disc_A_0.10 <-
  selected_disc("A", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B1
##
disc_B1_0.10 <-
  selected_disc("B1", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario B2
##
disc_B2_0.10 <-
  selected_disc("B2", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C1
##
disc_C1_0.10 <-
  selected_disc("C1", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
## Scenario C2
##
disc_C2_0.10 <-
  selected_disc("C2", equal, "0.1", p.baseline,
                min.grpsize, max.grpsize, M,
                details = TRUE)
##
disc_0.10 <-
  rbind(disc_A_0.10, disc_B1_0.10, disc_B2_0.10,
        disc_C1_0.10, disc_C2_0.10)
##
with(disc_0.10, table(scenario, selected))
with(disc_0.10, table(scenario, interactions))


##
##
## (8) Table 4
##
##

sink("results/Table4.txt")
##
cat("*** tau2 = 0.00 ***\n")
with(disc_0.00, table(scenario, interactions))
##
cat("*** tau2 = 0.01 ***\n")
with(disc_0.01, table(scenario, interactions))
##
cat("*** tau2 = 0.10 ***\n")
with(disc_0.10, table(scenario, interactions))
##
sink()
