library(netmeta)
library(lattice)
##library(grid)
##library(gridBase)
library(dplyr)


source("R/funcs.R")


## Load MSE and CP
##
load("results/mse.rda")
mse.full$model <-
  factor(mse.full$model,
         levels = c("NMA", "iCNMA", "aCNMA"),
         labels = c("NMA", "iCNMA", "aCNMA"))
mse.disc$model <-
  factor(mse.disc$model,
         levels = c("iCNMA", "aCNMA"),
         labels = c("iCNMA", "aCNMA"))
##
load("results/cp.rda")
cp.full$model <-
  factor(cp.full$model,
         levels = c("NMA", "iCNMA", "aCNMA"),
         labels = c("NMA", "iCNMA", "aCNMA"))
cp.disc$model <-
  factor(cp.disc$model,
         levels = c("iCNMA", "aCNMA"),
         labels = c("iCNMA", "aCNMA"))


## Load simulated connected data sets and select one
##
load("results/Scenario=A, equal=0, tau2=0.1, p=0.1, min=50, max=200/simdata.rda")
##
net <- netmeta(simdata[[3]])

  
##
##
## (1) Network plot (Figure 2)
##
##

pdf("graphics/Figure2.pdf")
netgraph(net, number = TRUE, multi = FALSE, plastic = FALSE,
  col = "black", thickness = "number", lwd = 3,
  seq = net$trts)
dev.off()

  
##
##
## (2) Figure 3
##
##

pdf("graphics/Figure3a.pdf", width = 8, height = 4)
xyplot(MSE ~ I(as.numeric(het)) | scenario,
       data = mse.full, groups = model, type = "b",
       layout = c(5, 1),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "moderate")),
         y = list(at = seq(0, 0.1, by = 0.02))),
       ylim = c(-0.005, 0.07),
       xlab = "", ylab = "Average Mean Squared Error",
       col = c("blue", "green", "red"), lty = 3:1, lwd = 2)
dev.off()
##
pdf("graphics/Figure3b.pdf", width = 8, height = 4)
xyplot(CP ~ I(as.numeric(het)) | scenario,
       data = cp.full, groups = model, type = "b",
       layout = c(5, 1),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "moderate"))),
       ylim = c(0.89, 1.01),
       xlab = "Heterogeneity", ylab = "Average Coverage Probability",
       col = c("blue", "green", "red"), lty = 3:1, lwd = 2,
       panel = function(...) {
         panel.abline(h = 0.95)
         panel.abline(h = c(0.936, 0.964), lty = 2)
         panel.xyplot(...) 
       },
       key =
         list(text = list(c("NMA", "Selected CNMA", "Additive CNMA")),
              space = "bottom",
              col = c("blue", "green", "red"), pch = 1, lty = 3:1, lwd = 2,
              points = FALSE, lines = FALSE,
              columns = 3))
dev.off()

mse <- rbind(mse.full, mse.disc)
mse$network <-
  factor(mse$network, c("disconnected", "connected"), c("Disconnected", "Connected"))
##
pdf("graphics/Figure3.pdf", width = 9, height = 6)
xyplot(MSE ~ I(as.numeric(het)) | scenario * network,
       data = mse, groups = model, type = c("b", "g"),
       layout = c(5, 2),
       par.settings = list(strip.background = list(col = "transparent")),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "mod")),
         y = list(relation = "free",
                  alternating = FALSE,
                  limits =
                    list(c(0, 0.40), c(0, 0.40), c(0, 0.40),
                         c(0, 0.40), c(0, 0.40),
                         c(0, 0.07), c(0, 0.07), c(0, 0.07),
                         c(0, 0.07), c(0, 0.07)),
                  at = list(seq(0, 0.4, by = 0.1), "", "", "", "",
                            seq(0, 0.06, by = 0.02), "", "", "", ""))
       ),
       ##ylim = c(-0.01, 0.40),
       xlab = "Heterogeneity", ylab = "Average Mean Squared Error",
       col = c("blue", "green", "red"), lty = 3:1, lwd = 2,
       key =
         list(text = list(c("NMA", "Selected CNMA", "Additive CNMA")),
              space = "bottom",
              col = c("blue", "green", "red"), pch = 1, lty = 3:1, lwd = 2,
              points = FALSE, lines = FALSE,
              columns = 3))
dev.off()


##
##
## (3) Figure 4
##
##

pdf("graphics/Figure4a.pdf", width = 8, height = 4)
xyplot(MSE ~ I(as.numeric(het)) | scenario,
       data = mse.disc, groups = model, type = "b",
       layout = c(5, 1),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "moderate")),
         y = list(at = seq(0, 0.4, by = 0.1))),
       ylim = c(-0.02, 0.42),
       xlab = "Heterogeneity", ylab = "Average Mean Squared Error",
       col = c("green", "red"), lty = 2:1, lwd = 2)
dev.off()
##
pdf("graphics/Figure4b.pdf", width = 8, height = 4)
xyplot(CP ~ I(as.numeric(het)) | scenario,
       data = cp.disc, groups = model, type = "b",
       layout = c(5, 1),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "moderate"))),
       ylim = c(0.83, 1.01),
       xlab = "Heterogeneity", ylab = "Average Coverage Probability",
       col = c("green", "red"), lty = 2:1, lwd = 2,
       panel = function(...) {
         panel.abline(h = 0.95)
         panel.abline(h = c(0.936, 0.964), lty = 2)
         panel.xyplot(...) 
       },
       key =
         list(text = list(c("Selected CNMA", "Additive CNMA")),
              space = "bottom",
              col = c("green", "red"), pch = 1, lty = 2:1, lwd = 2,
              points = FALSE, lines = FALSE,
              columns = 2))
dev.off()

cp <- rbind(cp.full, cp.disc)
cp$network <-
  factor(cp$network, c("disconnected", "connected"), c("Disconnected", "Connected"))
##
pdf("graphics/Figure4.pdf", width = 9, height = 6)
xyplot(CP ~ I(as.numeric(het)) | scenario * network,
       data = cp, groups = model, type = c("b", "g"),
       layout = c(5, 2),
       par.settings = list(strip.background = list(col = "transparent")),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "mod"))),
       ylim = c(0.82, 1.01),
       xlab = "Heterogeneity", ylab = "Average Coverage Probability",
       col = c("blue", "green", "red"), lty = 3:1, lwd = 2,
       panel = function(...) {
         panel.abline(h = 0.95)
         panel.abline(h = c(0.936, 0.964), lty = 2)
         panel.xyplot(...) 
       },
       key =
         list(text = list(c("NMA", "Selected CNMA", "Additive CNMA")),
              space = "bottom",
              col = c("blue", "green", "red"), pch = 1, lty = 3:1, lwd = 2,
              points = FALSE, lines = FALSE,
              columns = 3))
dev.off()


##
##
## (4) Supplement
##
##

## Connected network, MSE, individual effects
##
mse.full.trts <-
  by_full(subset(mse.full.long, treat == "A"), "MSE") %>%
  mutate(treat = "A", .after = het)
##
for (trt.i in unique(mse.full.long$treat)[-1])
  mse.full.trts <-
    rbind(mse.full.trts,
          by_full(subset(mse.full.long, treat == trt.i), "MSE") %>%
          mutate(treat = trt.i, .after = het))
##
scens <- c("A", "B1", "B2", "C1", "C2")
mse.full.trts$scenario <-
  factor(mse.full.trts$scenario,
         levels = scens,
         labels = paste("Scenario:", scens))
##
trts <- c("A", "B", "C", "D", "A+B", "A+C", "C+D")
mse.full.trts$treat <-
    factor(mse.full.trts$treat,
           levels = rev(trts),
           labels = paste0("OR: ", rev(trts), " vs P"))
##
mse.full.trts$model <-
  factor(mse.full.trts$model,
         levels = c("NMA", "iCNMA", "aCNMA"),
         labels = c("NMA", "iCNMA", "aCNMA"))
##
pdf("graphics/FigureA1.pdf", width = 9, height = 15)
xyplot(MSE ~ I(as.numeric(het)) | scenario * treat, 
       data = mse.full.trts, groups = model, type = c("b", "g"),
       layout = c(5, 7),
       par.settings = list(strip.background = list(col = "transparent")),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "mod"))),
       xlab = "Heterogeneity", ylab = "Mean Squared Error",
       col = c("blue", "green", "red"), lty = 3:1, lwd = 2,
       key =
         list(text = list(c("NMA", "Selected CNMA", "Additive CNMA")),
              space = "bottom",
              col = c("blue", "green", "red"), pch = 1, lty = 3:1, lwd = 2,
              points = FALSE, lines = FALSE,
              columns = 3))
dev.off()


## Connected network, CP, individual effects
##
cp.full.trts <-
  by_full(subset(cp.full.long, treat == "A"), "CP") %>%
  mutate(treat = "A", .after = het)
##
for (trt.i in unique(cp.full.long$treat)[-1])
  cp.full.trts <-
    rbind(cp.full.trts,
          by_full(subset(cp.full.long, treat == trt.i), "CP") %>%
          mutate(treat = trt.i, .after = het))
##
scens <- c("A", "B1", "B2", "C1", "C2")
cp.full.trts$scenario <-
  factor(cp.full.trts$scenario,
         levels = scens,
         labels = paste("Scenario:", scens))
##
trts <- c("A", "B", "C", "D", "A+B", "A+C", "C+D")
cp.full.trts$treat <-
    factor(cp.full.trts$treat,
           levels = rev(trts),
           labels = paste0("OR: ", rev(trts), " vs P"))
##
cp.full.trts$model <-
  factor(cp.full.trts$model,
         levels = c("NMA", "iCNMA", "aCNMA"),
         labels = c("NMA", "iCNMA", "aCNMA"))
##
pdf("graphics/FigureA2.pdf", width = 9, height = 15)
xyplot(CP ~ I(as.numeric(het)) | scenario * treat,
       data = cp.full.trts, groups = model, type = c("b", "g"),
       layout = c(5, 7),
       par.settings = list(strip.background = list(col = "transparent")),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "mod"))),
       xlab = "Heterogeneity", ylab = "Coverage Probability",
       col = c("blue", "green", "red"), lty = 3:1, lwd = 2,
       panel = function(...) {
         panel.abline(h = 0.95)
         panel.abline(h = c(0.936, 0.964), lty = 2)
         panel.xyplot(...) 
       },
       key =
         list(text = list(c("NMA", "Selected CNMA", "Additive CNMA")),
              space = "bottom",
              col = c("blue", "green", "red"), pch = 1, lty = 3:1, lwd = 2,
              points = FALSE, lines = FALSE,
              columns = 3))
dev.off()


## Disconnected network, MSE, individual effects
##
mse.disc.trts <-
  by_disc(subset(mse.disc.long, treat == "A"), "MSE") %>%
  mutate(treat = "A", .after = het)
##
for (trt.i in unique(mse.disc.long$treat)[-1])
  mse.disc.trts <-
    rbind(mse.disc.trts,
          by_disc(subset(mse.disc.long, treat == trt.i), "MSE") %>%
          mutate(treat = trt.i, .after = het))
##
scens <- c("A", "B1", "B2", "C1", "C2")
mse.disc.trts$scenario <-
  factor(mse.disc.trts$scenario,
         levels = scens,
         labels = paste("Scenario:", scens))
##
trts <- c("A", "B", "C", "D", "A+B", "A+C", "C+D")
mse.disc.trts$treat <-
    factor(mse.disc.trts$treat,
           levels = rev(trts),
           labels = paste0("OR: ", rev(trts), " vs P"))
##
mse.disc.trts$model <-
  factor(mse.disc.trts$model,
         levels = c("iCNMA", "aCNMA"),
         labels = c("iCNMA", "aCNMA"))
##
pdf("graphics/FigureA3.pdf", width = 9, height = 15)
xyplot(MSE ~ I(as.numeric(het)) | scenario * treat,
       data = mse.disc.trts, groups = model, type = c("b", "g"),
       layout = c(5, 7),
       par.settings = list(strip.background = list(col = "transparent")),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "mod"))),
       xlab = "Heterogeneity", ylab = "Mean Squared Error",
       col = c("green", "red"), lty = 2:1, lwd = 2,
       key =
         list(text = list(c("Selected CNMA", "Additive CNMA")),
              space = "bottom",
              col = c("green", "red"), pch = 1, lty = 2:1, lwd = 2,
              points = FALSE, lines = FALSE,
              columns = 2))
dev.off()


## Disconnected network, CP, individual effects
##
cp.disc.trts <-
  by_disc(subset(cp.disc.long, treat == "A"), "CP") %>%
  mutate(treat = "A", .after = het)
##
for (trt.i in unique(cp.disc.long$treat)[-1])
  cp.disc.trts <-
    rbind(cp.disc.trts,
          by_disc(subset(cp.disc.long, treat == trt.i), "CP") %>%
          mutate(treat = trt.i, .after = het))
##
scens <- c("A", "B1", "B2", "C1", "C2")
cp.disc.trts$scenario <-
  factor(cp.disc.trts$scenario,
         levels = scens,
         labels = paste("Scenario:", scens))
##
trts <- c("A", "B", "C", "D", "A+B", "A+C", "C+D")
cp.disc.trts$treat <-
    factor(cp.disc.trts$treat,
           levels = rev(trts),
           labels = paste0("OR: ", rev(trts), " vs P"))
##
cp.disc.trts$model <-
  factor(cp.disc.trts$model,
         levels = c("iCNMA", "aCNMA"),
         labels = c("iCNMA", "aCNMA"))
##
pdf("graphics/FigureA4.pdf", width = 9, height = 15)
xyplot(CP ~ I(as.numeric(het)) | scenario * treat, 
       data = cp.disc.trts, groups = model, type = c("b", "g"),
       layout = c(5, 7),
       par.settings = list(strip.background = list(col = "transparent")),
       scales = list(
         x = list(at = 1:3, labels = c("no", "low", "mod"))),
       xlab = "Heterogeneity", ylab = "Coverage Probability",
       col = c("green", "red"), lty = 2:1, lwd = 2,
       panel = function(...) {
         panel.abline(h = 0.95)
         panel.abline(h = c(0.936, 0.964), lty = 2)
         panel.xyplot(...) 
       },
       key =
         list(text = list(c("Selected CNMA", "Additive CNMA")),
              space = "bottom",
              col = c("green", "red"), pch = 1, lty = 2:1, lwd = 2,
              points = FALSE, lines = FALSE,
              columns = 2))
dev.off()
