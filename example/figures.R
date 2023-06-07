library("netmeta")
library("dplyr")


source("R/getvar.R")


## Load (C)NMA results for (dis)connected networks
##
load("results/full-modsel.rda")
load("results/discs-modsel.rda")
load("results/discs-modsel-table.rda")


##
##
## (1) Network graphics
##
##

## Figure 1 (network plot for connected Cochrane data)
##
pdf("graphics/Figure1.pdf")
##
netgraph(NMA.full, number = TRUE, plastic = FALSE,
    col = "black", seq = "o", srt.labels = "o",
    rotate = 360 * (3 / NMA.full$n))
##
dev.off()


## Figure A5 (network plots for disconnected networks 1-9)
##
load("results/discs-AEdata.rda")
n.discs <- length(AEdata.discs$data)
##
pdf("graphics/FigureA5.pdf", width = 12, height = 9)
##
par(mfrow = c(3, 3))
for (i in seq_len(n.discs)) {
  nc.i <- netconnection(AEdata.discs$data[[i]])
  netgraph(nc.i,
    reference = "plac", srt.label = "o", lwd = 1)
  title(sub = paste0("Disconnected ", i,
                     " (m=", nc.i$m, ", k=", nc.i$k, ")"))
}
##
invisible(dev.off())


##
##
## (2) Forest plots
##
##

## Figure 5
##

f2 <- ms.full$hetstats3[order(ms.full$hetstats3$Q), ]
sel <-
  f2$interaction1 == "onda*scop" | f2$interaction2 == "onda*scop" |
  f2$interaction3 == "onda*scop"
f2 <- f2[sel, ]
head(f2, 21)
ms.full$selected
ms.full$selected3


## int1 <- "dexa*gran"
## int2 <- "dexa*trop"
## int3 <- "onda*scop"
## int123 <- paste(c(int3, int1, int2), collapse = " + ")
## sel2 <-
##   with(ms.full$hetstats3,
##        interaction1 == int1 &
##        interaction2 == int2 &
##        interaction3 == int3)
## ##
## id2 <- ms.full$hetstats3$id[sel2]
## id2
##
cnmas <-
  c(list(NMA.full,                      # standard NMA
         aCNMA.full,                    # additive CNMA
         ms.full$iCNMA.selected),       # full connected NMA
    lapply(ms.discs, getvar, var1 = "iCNMA.selected"))
                                        # iCNMA, disconnected
 ##
names.cnmas <-
  c("Standard NMA", "Additive CNMA",
    paste0("Interaction CNMA (", 
           paste(ms.full$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. ", seq_len(nrow(modsel)), " (",
           modsel$modsel, ")"))
##
nb.iCNMA3 <-
  netbind(cnmas, name = names.cnmas,
          common = FALSE,
          col.study = seq_len(length(names.cnmas)),
          col.square = seq_len(length(names.cnmas)))
##
pdf("graphics/Figure5.pdf", height = 25.5, width = 10)
forest(nb.iCNMA3, subset = c("amis", "apre", "apre+scop", "palo", "ramo"),
  xlim = c(0.05, 10), at = c(0.1, 0.5, 1, 2, 10),
  leftlabs = "Intervention",
  col.by = "black")
dev.off()

## Figure A6 (all results)
##
pdf("graphics/FigureA6.pdf", height = 126, width = 10)
forest(nb.iCNMA3,
  xlim = c(0.05, 10), at = c(0.1, 0.5, 1, 2, 10),
  leftlabs = "Intervention",
  col.by = "black")
dev.off()

## Figure A7 (results for main subnetwork)
##
load("results/discs-sep.rda")
##
cnmas <-
  c(list(NMA.full,    # standard NMA
         aCNMA.full), # full connected NMA
    sep.discs)
##
names.cnmas <-
  c("Standard NMA", "Additive CNMA connected network", 
    paste("Separate NMA disc.", seq_len(length(sep.discs))))
##
sepNMAs <-
  netbind(cnmas, name = names.cnmas, common = FALSE,
    col.study = seq_along(names.cnmas),
    col.square = seq_along(names.cnmas))
##
pdf("graphics/FigureA7.pdf", height = 87.5, width = 7)
forest(sepNMAs, xlim = c(0.01, 10), leftlabs = "Intervention")
dev.off()
