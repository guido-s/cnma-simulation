library(netmeta)


## Load (C)NMA results for (dis)connected networks
##
load("results/full-modsel.rda")
load("results/discs-modsel.rda")


##
##
## (1) Network graphics
##
##

## Figure 4 (network plot for connected Cochrane data)
##
pdf("graphics/Figure4.pdf")
netgraph(NMA.full, number = TRUE, multi = FALSE, plastic = FALSE,
  col = "black", thickness = "se.random",
  seq = "o", srt.labels = "o")
dev.off()


## Figure A5 (network plots for disconnected networks 1-9)
##
load("results/discs-AEdata.rda")
##
pdf("graphics/FigureA5.pdf")
par(mfrow = c(3, 3))
for (i in 1:9) {
  netgraph(netconnection(treat1, treat2, studlab,
             data = AEdata.discs$data[[i]]),
    reference = "plac", srt.label = "o", lwd = 3)
  title(sub = paste("Disconnected", i))
}
invisible(dev.off())


##
##
## (2) Forest plots
##
##

## Figure 5
##

int1 <- "dexa*gran"
int2 <- "dexa*trop"
int3 <- "onda*scop"
int123 <- paste(c(int3, int1, int2), collapse = " + ")
sel2 <-
  with(ms.full$hetstats3,
       interaction1 == int1 &
       interaction2 == int2 &
       interaction3 == int3)
##
id2 <- ms.full$hetstats3$id[sel2]
id2
##
cnmas <-
  list(NMA.full,                      # standard NMA
       aCNMA.full,                    # additive CNMA
       ms.full$iCNMA.selected,        # full connected NMA
       ms.full$iCNMA3[[id2]],         # full connected NMA
       ms.discs$disc1$iCNMA.selected, # iCNMA, disconnected 1
       ms.discs$disc2$iCNMA.selected, # iCNMA, disconnected 2
       ms.discs$disc3$iCNMA.selected, # iCNMA, disconnected 3
       ms.discs$disc4$iCNMA.selected, # iCNMA, disconnected 4
       ms.discs$disc5$iCNMA.selected, # iCNMA, disconnected 5
       ms.discs$disc6$iCNMA.selected, # iCNMA, disconnected 6
       ms.discs$disc7$iCNMA.selected, # iCNMA, disconnected 7
       ms.discs$disc8$iCNMA.selected, # iCNMA, disconnected 8
       ms.discs$disc9$iCNMA.selected) # iCNMA, disconnected 9
 ##
 names.cnmas <-
  c("Standard NMA", "Additive CNMA",
    paste0("Interaction CNMA (",
           paste(ms.full$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA (",
           int123, ")"),
    paste0("Interaction CNMA disc. 1 (",
           paste(ms.discs$disc1$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. 2 (",
           paste(ms.discs$disc2$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. 3 (",
           paste(ms.discs$disc3$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. 4 (",
           paste(ms.discs$disc4$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. 5 (",
           paste(ms.discs$disc5$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. 6 (",
           paste(ms.discs$disc6$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. 7 (",
           paste(ms.discs$disc7$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. 8 (",
           paste(ms.discs$disc8$selected3[, 1:3], collapse = " + "), ")"),
    paste0("Interaction CNMA disc. 9 (",
           paste(ms.discs$disc9$selected3[, 1:3], collapse = " + "), ")"))
##
nb.iCNMA3 <-
  netbind(cnmas, name = names.cnmas, common = FALSE,
          col.study = seq_len(length(names.cnmas)),
          col.square = seq_len(length(names.cnmas)))
##
pdf("graphics/Figure5.pdf", height = 16.5, width = 9)
forest(nb.iCNMA3, subset = c("amis", "apre", "apre+scop","palo", "ramo"),
  xlim = c(0.05, 10), at = c(0.1, 0.5, 1, 2, 10),
  leftlabs = "Intervention",
  col.by = "black")
dev.off()

## Figure A6 (all results)
##
pdf("graphics/FigureA6.pdf", height = 79.5, width = 9)
forest(nb.iCNMA3,
  xlim = c(0.05, 10), at = c(0.1, 0.5, 1, 2, 10),
  leftlabs = "Intervention",
  col.by = "black")
dev.off()

## Figure A7 (results for main subnetwork)
##
load("results/discs-sep.rda")
##
sep.disc1 <- sep.discs[[1]]
sep.disc2 <- sep.discs[[2]]
sep.disc3 <- sep.discs[[3]]
sep.disc4 <- sep.discs[[4]]
sep.disc5 <- sep.discs[[5]]
sep.disc6 <- sep.discs[[6]]
sep.disc7 <- sep.discs[[7]]
sep.disc8 <- sep.discs[[8]]
sep.disc9 <- sep.discs[[9]]
##
cnmas <- list(NMA.full, # standard NMA
  aCNMA.full,           # full connected NMA
  sep.disc1, sep.disc2, sep.disc3, sep.disc4,
  sep.disc5, sep.disc6, sep.disc7, sep.disc8, 
  sep.disc9)
##
names.cnmas <-
  c("Standard NMA", "Additive CNMA connected network", 
    paste("Separate NMA disc.", 1:9))
##
sepNMAs <-
  netbind(cnmas, name = names.cnmas, common = FALSE,
    col.study = seq_along(names.cnmas),
    col.square = seq_along(names.cnmas))
##
pdf("graphics/FigureA7.pdf", height = 54, width = 7)
forest(sepNMAs, xlim = c(0.01, 5), leftlabs = "Intervention")
dev.off()
