library("netmeta")
library("dplyr")

source("R/summary_modselect.R")
source("R/getvar.R")

options(width = 150)

loadit <- TRUE

## Load (C)NMA results for connected network
##
if (loadit)
  load("results/full-modsel.rda")


## Interaction CNMA model for full connected network
##
ms.full$selected1
ms.full$selected2
ms.full$selected3
ms.full$selected4
ms.full$selected # iCNMA3
##
summary_modselect(ms.full)
##
##ms.full$iCNMA3[[ms.full$selected3$id]]


## Load CNMA results for disconnected networks
##
if (loadit)
  load("results/discs-AEdata.rda")
if (loadit)
  load("results/discs-modsel.rda")


##
## Model selection for network with minimal set
##
getinit <- function(x) unique(x$id) == 1
id <- (1:length(AEdata.discs$fulldata))[sapply(AEdata.discs$fulldata, getinit)]
disc.init <- paste0("disc", id)
disc.init
##
## Information on subnetworks
##
dat.init <- AEdata.discs$data[[as.numeric(gsub("disc", "", disc.init))]]
dis.init <-
  discomb(TE, seTE, treat1, treat2, studlab, data = dat.init, sm = "RR")
nc.init <- netconnection(dat.init)
nc.init
## Interventions in subnetworks
with(nc.init, sort(unique(c(treat1[subnet == 1], treat2[subnet == 1]))))
with(nc.init, sort(unique(c(treat1[subnet == 2], treat2[subnet == 2]))))
## Number of interventions, studies and pairwise comparisons in main subnetwork
length(
  with(nc.init, sort(unique(c(treat1[subnet == 2], treat2[subnet == 2])))))
length(unique(nc.init$studlab[nc.init$subnet == 2]))
length(nc.init$studlab[nc.init$subnet == 2])
## Number of interventions, studies and pairwise comparisons in auxiliary
## subnetwork
length(
  with(nc.init, sort(unique(c(treat1[subnet == 1], treat2[subnet == 1])))))
length(unique(nc.init$studlab[nc.init$subnet == 1]))
length(nc.init$studlab[nc.init$subnet == 1])
##
ms.discs[[disc.init]]$selected1
ms.discs[[disc.init]]$selected2
ms.discs[[disc.init]]$selected3
ms.discs[[disc.init]]$selected4
ms.discs[[disc.init]]$selected # iCNMA4
##
summary_modselect(ms.discs[[disc.init]])


## Results of model selection in disconnected CNMA models
##
modsel <- data.frame(
  k = sapply(ms.discs, getvar, var1 = "k"),
  m = sapply(ms.discs, getvar, var1 = "m"),
  s = sapply(ms.discs, getvar, var1 = "subnets"),
  selected = sapply(ms.discs, getvar),
  modsel = "", aCNMA = "",
  Q0 = sapply(ms.discs, getvar, var1 = "aCNMA", var2 = "Q.additive"),
  df0 = sapply(ms.discs, getvar, var1 = "aCNMA", var2 = "df.Q.additive"))
##
## Models with 1 two-way interaction
##
modsel$iCNMA1 <- modsel$int11 <-
  sapply(ms.discs, getvar, var1 = "selected1", var2 = "interaction")
modsel$Q1 <- sapply(ms.discs, getvar, var1 = "selected1", var2 = "Q")
modsel$df1 <- sapply(ms.discs, getvar, var1 = "selected1", var2 = "df")
##
## Models with 2 two-way interaction
##
modsel$int21 <-
  sapply(ms.discs, getvar, var1 = "selected2", var2 = "interaction1")
modsel$int22 <-
  sapply(ms.discs, getvar, var1 = "selected2", var2 = "interaction2")
##
modsel$iCNMA2 <- paste(modsel$int21, modsel$int22, sep = " + ")
modsel$Q2 <- sapply(ms.discs, getvar, var1 = "selected2", var2 = "Q")
modsel$df2 <- sapply(ms.discs, getvar, var1 = "selected2", var2 = "df")
##
## Models with 3 two-way interaction
##
modsel$int31 <-
  sapply(ms.discs, getvar, var1 = "selected3", var2 = "interaction1")
modsel$int32 <-
  sapply(ms.discs, getvar, var1 = "selected3", var2 = "interaction2")
modsel$int33 <-
  sapply(ms.discs, getvar, var1 = "selected3", var2 = "interaction3")
modsel$Q3 <- sapply(ms.discs, getvar, var1 = "selected3", var2 = "Q")
modsel$df3 <- sapply(ms.discs, getvar, var1 = "selected3", var2 = "df")
##
modsel$iCNMA3 <-
  paste(modsel$int31, modsel$int32, modsel$int33, sep = " + ")
##
## Models with 4 two-way interaction
##
modsel$int41 <-
  sapply(ms.discs, getvar, var1 = "selected4", var2 = "interaction1")
modsel$int42 <-
  sapply(ms.discs, getvar, var1 = "selected4", var2 = "interaction2")
modsel$int43 <-
  sapply(ms.discs, getvar, var1 = "selected4", var2 = "interaction3")
modsel$int44 <-
  sapply(ms.discs, getvar, var1 = "selected4", var2 = "interaction4")
##
modsel$iCNMA4 <-
  paste(modsel$int41, modsel$int42, modsel$int43, modsel$int44, sep = " + ")
modsel$Q4 <- sapply(ms.discs, getvar, var1 = "selected4", var2 = "Q")
modsel$df4 <- sapply(ms.discs, getvar, var1 = "selected4", var2 = "df")
##
for (i in seq_len(nrow(modsel)))
  modsel$modsel[i] <- modsel[i, modsel$selected[i]]
##
modsel %>% select(k, m, s, selected, iCNMA1, iCNMA2, iCNMA3, iCNMA4)
##
save(modsel, file = "results/discs-modsel-table.rda")


## Detailed results for iCNMA model selection
##
lapply(ms.discs, summary_modselect)


## Table with selected iCNMA model
##
table(modsel$modsel, modsel$selected)
