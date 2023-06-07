dir_create <- function(directory,
                       show_warnings = FALSE,
                       recursive = TRUE) {

  if (!file.exists(directory))
    dir.create(directory, show_warnings, recursive)

  invisible(NULL)
}


combinations <- function(x, n = NULL) {
  
  meta:::chkclass(x, "netcomb")
  
  comps <- x$comps
  trts <- x$trts
  ##
  combs <- trts[!(trts %in% comps)]
  ##
  if (!is.null(x$inactive))
    combs <- combs[!(combs %in% x$inactive)]
  
  if (!is.null(n)) {
    if (!meta:::is_wholenumber(n))
      stop("Argument 'n' must be a whole number.")
    sel <- unlist(lapply(strsplit(combs, x$sep.comps,
                                  fixed = TRUE), length)) == n
    combs <- combs[sel]
  }
  
  combs
}


selcomp <- function(x, n = 1, sep = "+") {
  if (!meta:::is_wholenumber(n))
    stop("Argument 'n' must be a whole number.")
  ##
  res <- lapply(strsplit(x, sep, fixed = TRUE),
                gsub, pattern = " ", replacement = "")
  n.comps <- sapply(res, length)
  ##
  ## if (any(n > n.comps))
  ##  stop("Argument 'n' must be smaller equal than ", max(n.comps), ".")
  ##
  res <- sapply(res, "[", n)
  
  res
}


selint <- function(x, n = 1, sep = "*") {
  
  if (!meta:::is_wholenumber(n))
    stop("Argument 'n' must be a whole number.")
  
  res <- lapply(strsplit(x, sep, fixed = TRUE),
                gsub, pattern = " ", replacement = "")
  n.comps <- sapply(res, length)
  ##
  ## if (any(n > n.comps))
  ##  stop("Argument 'n' must be smaller equal than ", max(n.comps), ".")
  ##
  res <- sapply(res, "[", n)
  
  res
}


getint <- function(x) {
  comps <- x$comps
  paste(comps[grepl("*", comps, fixed = TRUE)], collapse = "+")
}


simpath <- function(scenario, equal, tau2,
                    p.baseline, min.grpsize, max.grpsize,
                    iteration, M) {
  
  path <- paste0("results/Scenario=", scenario,
                 ", equal=", as.numeric(equal),
                 ", tau2=", tau2,
                 ", p=", p.baseline,
                 ", min=", min.grpsize,
                 ", max=", max.grpsize,
                 "/")
  ##
  if (!missing(iteration)) {
    len.M <- nchar(M)
    len.it <- nchar(iteration)
    path <- paste0(path, "iteration",
                   paste(rep("0", len.M - len.it), collapse = ""),
                   iteration, "/")
  }
  
  path
}


selected_full <- function(scenario, equal, tau2, p.baseline,
                          min.grpsize, max.grpsize, M, seq = NULL,
                          details = FALSE) {
  
  if (missing(seq))
    seq <- seq(1, M)
  
  res <- data.frame(scenario, equal, tau2, p.baseline,
                    min.grpsize, max.grpsize,
                    run = seq,
                    pval.Q.diff = NA,
                    selected = "", interactions = "")
  ##
  for (i in seq) {
    path <- simpath(scenario, equal, tau2,
                    p.baseline, min.grpsize, max.grpsize,
                    i, M)
    load(paste0(path, "full-modsel.rda"))
    ##
    res$pval.Q.diff[i] <- ms.full$aCNMA$pval.Q.diff
    res$selected[i] <- ms.full$selected
    res$interactions[i] <-
      switch(ms.full$selected,
             aCNMA = "--",
             iCNMA1 = ms.full$selected1$interaction,
             iCNMA2 = paste(ms.full$selected2$interaction1,
                            ms.full$selected2$interaction2,
                            sep = " + "),
             iCNMA3 = paste(ms.full$selected3$interaction1,
                            ms.full$selected3$interaction2,
                            ms.full$selected3$interaction3,
                            sep = " + "))
  }

  res$interactions <-
    factor(res$interactions,
           levels = c("--",
                      "A*B", "A*C", "C*D",
                      "A*B + C*D", "A*B + A*C", "A*C + C*D"))
  ##
  if (details)
    CrossTable(res$selected,
               prop.r = TRUE, prop.c = FALSE, prop.t = FALSE,
               prop.chisq = FALSE,
               digits = 1, format = "SPSS")
  
  res
}


selected_disc <- function(scenario, equal, tau2, p.baseline,
                          min.grpsize, max.grpsize, M, seq = NULL,
                          details = FALSE) {
  
  if (missing(seq))
    seq <- seq(1, M)
  
  res <- data.frame(scenario, equal, tau2, p.baseline,
                    min.grpsize, max.grpsize,
                    run = seq,
                    selected = "", interactions = "")
  ##
  for (i in seq) {
    path <- simpath(scenario, equal, tau2,
                    p.baseline, min.grpsize, max.grpsize,
                    i, M)
    load(paste0(path, "disc-modsel.rda"))
    ##
    res$selected[i] <- ms.disc$selected
    res$interactions[i] <-
      switch(ms.disc$selected,
             aCNMA = "--",
             iCNMA1 = ms.disc$selected1$interaction,
             iCNMA2 = paste(ms.disc$selected2$interaction1,
                            ms.disc$selected2$interaction2,
                            sep = " + "),
             iCNMA3 = paste(ms.disc$selected3$interaction1,
                            ms.disc$selected3$interaction2,
                            ms.disc$selected3$interaction3,
                            sep = " + "))
  }

  res$interactions <-
    factor(res$interactions,
           levels = c("--",
                      "A*B", "A*C", "C*D",
                      "A*B + C*D", "A*B + A*C", "A*C + C*D"))
  ##
  if (details)
    CrossTable(res$selected,
               prop.r = TRUE, prop.c = FALSE, prop.t = FALSE,
               prop.chisq = FALSE,
               digits = 1, format = "SPSS")
  
  res
}


read_mse_full <- function(scenario, equal, tau2, p.baseline,
                          min.grpsize, max.grpsize, M, seq = NULL) {
  
  if (missing(seq))
    seq <- seq(1, M)
  
  path <- simpath(scenario, equal, tau2,
                  p.baseline, min.grpsize, max.grpsize)
  load(paste0(path, "settings.rda"))
  ##
  lnORs <- log(ORs)
  trts <- names(lnORs)
  
  ## Variable names for MSEs and lnORs
  ##
  mses <- paste0("MSE.", gsub("+", "", trts, fixed = TRUE))
  ests <- paste0("lnOR.", gsub("+", "", trts, fixed = TRUE))
  
  ## Sceletions of data sets for results
  ## (variable MSE contains the overall MSE)
  ##
  NMA <- aCNMA <- iCNMA <-
    data.frame(network = "", model = "",
               scenario, equal, tau2, p.baseline,
               min.grpsize, max.grpsize,
               run = seq,
               MSE = NA)
  
  ## Add variables for MSEs and log odds ratio estimates
  ##
  for (j in seq_along(trts)) {
    NMA[[mses[j]]] <- aCNMA[[mses[j]]] <- iCNMA[[mses[j]]] <- NA
  }
  for (j in seq_along(trts)) {
    NMA[[ests[j]]] <- aCNMA[[ests[j]]] <- iCNMA[[ests[j]]] <- NA
  }
  
  ## Calculate overall MSE
  ##
  for (i in seq) {
    path <- simpath(scenario, equal, tau2,
                    p.baseline, min.grpsize, max.grpsize,
                    i, M)
    load(paste0(path, "full-modsel.rda"))
    ##
    if (ms.full$selected == "aCNMA")
      iCNMA.i <- ms.full$aCNMA
    else
      iCNMA.i <- ms.full$iCNMA.selected
    ##
    for (j in seq_along(trts)) {
      ##
      ## Estimated log odds ratios from NMA, aCNMA and selected iCNMA
      ##
      NMA[[ests[j]]][i] <- NMA.full$TE.random[trts[j], "P"]
      aCNMA[[ests[j]]][i] <- aCNMA.full$TE.random[trts[j], "P"]
      iCNMA[[ests[j]]][i] <- iCNMA.i$TE.random[trts[j], "P"]
      ##
      ## Calculate MSEs for individual parameters
      ##
      NMA[[mses[j]]][i] <- (NMA[[ests[j]]][i] - lnORs[trts[j]])^2
      aCNMA[[mses[j]]][i] <- (aCNMA[[ests[j]]][i] - lnORs[trts[j]])^2
      iCNMA[[mses[j]]][i] <- (iCNMA[[ests[j]]][i] - lnORs[trts[j]])^2
    }
    ##
    ## MSEs
    ##
    NMA$MSE[i] <- sum((NMA[i, ests] - lnORs)^2) / length(trts)
    aCNMA$MSE[i] <- sum((aCNMA[i, ests] - lnORs)^2) / length(trts)
    iCNMA$MSE[i] <- sum((iCNMA[i, ests] - lnORs)^2) / length(trts)
  }
  
  
  ## Use same abbreviated names for true log odds ratios
  ##
  names(lnORs) <-
    paste0("lnOR.", gsub("+", "", names(lnORs), fixed = TRUE))
  
  
  NMA$network <- aCNMA$network <- iCNMA$network <- "connected"
  ##
  NMA$model <- "NMA"
  aCNMA$model <- "aCNMA"
  iCNMA$model <- "iCNMA"
  
  
  list(NMA = NMA, aCNMA = aCNMA, iCNMA = iCNMA, lnORs = lnORs)
}


read_mse_disc <- function(scenario, equal, tau2, p.baseline,
                          min.grpsize, max.grpsize, M, seq = NULL) {
  
  if (missing(seq))
    seq <- seq(1, M)
  
  path <- simpath(scenario, equal, tau2,
                  p.baseline, min.grpsize, max.grpsize)
  load(paste0(path, "settings.rda"))
  ##
  lnORs <- log(ORs)
  trts <- names(lnORs)
  
  ## Variable names for MSEs and lnORs
  ##
  mses <- paste0("MSE.", gsub("+", "", trts, fixed = TRUE))
  ests <- paste0("lnOR.", gsub("+", "", trts, fixed = TRUE))
  
  ## Sceletions of data sets for results
  ## (variable MSE contains the overall MSE)
  ##
  aCNMA <- iCNMA <-
    data.frame(network = "", model = "",
               scenario, equal, tau2, p.baseline,
               min.grpsize, max.grpsize,
               run = seq,
               MSE = NA)
  
  ## Add variables for MSEs and log odds ratio estimates
  ##
  for (j in seq_along(trts)) {
    aCNMA[[mses[j]]] <- iCNMA[[mses[j]]] <- NA
  }
  for (j in seq_along(trts)) {
    aCNMA[[ests[j]]] <- iCNMA[[ests[j]]] <- NA
  }
  
  ## Calculate overall MSE
  ##
  for (i in seq) {
    path <- simpath(scenario, equal, tau2,
                    p.baseline, min.grpsize, max.grpsize,
                    i, M)
    load(paste0(path, "disc-modsel.rda"))
    ##
    if (ms.disc$selected == "aCNMA")
      iCNMA.i <- ms.disc$aCNMA
    else
      iCNMA.i <- ms.disc$iCNMA.selected
    ##
    for (j in seq_along(trts)) {
      ##
      ## Estimated log odds ratios from aCNMA and selected iCNMA
      ##
      aCNMA[[ests[j]]][i] <- aCNMA.disc$TE.random[trts[j], "P"]
      iCNMA[[ests[j]]][i] <- iCNMA.i$TE.random[trts[j], "P"]
      ##
      ## Calculate MSEs for individual parameters
      ##
      aCNMA[[mses[j]]][i] <- (aCNMA[[ests[j]]][i] - lnORs[trts[j]])^2
      iCNMA[[mses[j]]][i] <- (iCNMA[[ests[j]]][i] - lnORs[trts[j]])^2
    }
    ##
    ## MSEs
    ##
    aCNMA$MSE[i] <- sum((aCNMA[i, ests] - lnORs)^2) / length(trts)
    iCNMA$MSE[i] <- sum((iCNMA[i, ests] - lnORs)^2) / length(trts)
  }
  
  
  ## Use same abbreviated names for true log odds ratios
  ##
  names(lnORs) <-
    paste0("lnOR.", gsub("+", "", names(lnORs), fixed = TRUE))
  
  
  aCNMA$network <- iCNMA$network <- "connected"
  ##
  aCNMA$model <- "aCNMA"
  iCNMA$model <- "iCNMA"
  
  
  list(aCNMA = aCNMA, iCNMA = iCNMA, lnORs = lnORs)
}


read_cp_full <- function(scenario, equal, tau2, p.baseline,
                         min.grpsize, max.grpsize, M, seq = NULL) {
  
  if (missing(seq))
    seq <- seq(1, M)
  
  path <- simpath(scenario, equal, tau2,
                  p.baseline, min.grpsize, max.grpsize)
  load(paste0(path, "settings.rda"))
  ##
  lnORs <- log(ORs)
  trts <- names(lnORs)
  
  ## Variable names for CPs and CI levels
  ##
  cps <- paste0("CP.", gsub("+", "", trts, fixed = TRUE))
  lows <- paste0("low.", gsub("+", "", trts, fixed = TRUE))
  upps <- paste0("upp.", gsub("+", "", trts, fixed = TRUE))
  
  ## Sceletions of data sets for results
  ## (variable CP contains the overall CP)
  ##
  NMA <- aCNMA <- iCNMA <-
    data.frame(network = "", model = "",
               scenario, equal, tau2, p.baseline,
               min.grpsize, max.grpsize,
               run = seq,
               CP = NA)
  
  ## Add variables for CPs and log odds ratio estimates
  ##
  for (j in seq_along(trts)) {
    NMA[[cps[j]]] <- aCNMA[[cps[j]]] <- iCNMA[[cps[j]]] <- NA
  }
  for (j in seq_along(trts)) {
    NMA[[upps[j]]] <- NMA[[lows[j]]] <- NA
    aCNMA[[upps[j]]] <- aCNMA[[lows[j]]] <- NA
    iCNMA[[upps[j]]] <- iCNMA[[lows[j]]] <- NA
  }
  
  ## Calculate overall CP
  ##
  for (i in seq) {
    path <- simpath(scenario, equal, tau2,
                    p.baseline, min.grpsize, max.grpsize,
                    i, M)
    load(paste0(path, "full-modsel.rda"))
    ##
    if (ms.full$selected == "aCNMA")
      iCNMA.i <- ms.full$aCNMA
    else
      iCNMA.i <- ms.full$iCNMA.selected
    ##
    for (j in seq_along(trts)) {
      ##
      ## Lower CI limits from NMA, aCNMA and selected iCNMA
      ##
      NMA[[lows[j]]][i] <- NMA.full$lower.random[trts[j], "P"]
      aCNMA[[lows[j]]][i] <- aCNMA.full$lower.random[trts[j], "P"]
      iCNMA[[lows[j]]][i] <- iCNMA.i$lower.random[trts[j], "P"]
      ##
      ## Upper CI limits from NMA, aCNMA and selected iCNMA
      ##
      NMA[[upps[j]]][i] <- NMA.full$upper.random[trts[j], "P"]
      aCNMA[[upps[j]]][i] <- aCNMA.full$upper.random[trts[j], "P"]
      iCNMA[[upps[j]]][i] <- iCNMA.i$upper.random[trts[j], "P"]
      ##
      ## Calculate CPs for individual parameters
      ##
      NMA[[cps[j]]][i] <-
        inCI(lnORs[trts[j]], NMA[[lows[j]]][i], NMA[[upps[j]]][i])
      aCNMA[[cps[j]]][i] <-
        inCI(lnORs[trts[j]], aCNMA[[lows[j]]][i], aCNMA[[upps[j]]][i])
      iCNMA[[cps[j]]][i] <-
        inCI(lnORs[trts[j]], iCNMA[[lows[j]]][i], iCNMA[[upps[j]]][i])
    }
    ##
    ## CPs
    ##
    NMA$CP[i] <- sum(NMA[i, cps]) / length(trts)
    aCNMA$CP[i] <- sum(aCNMA[i, cps]) / length(trts)
    iCNMA$CP[i] <- sum(iCNMA[i, cps]) / length(trts)
  }
  
  
  ## Use same abbreviated names for true log odds ratios
  ##
  names(lnORs) <-
    paste0("lnOR.", gsub("+", "", names(lnORs), fixed = TRUE))
  
  
  NMA$network <- aCNMA$network <- iCNMA$network <- "connected"
  ##
  NMA$model <- "NMA"
  aCNMA$model <- "aCNMA"
  iCNMA$model <- "iCNMA"
  
  
  list(NMA = NMA, aCNMA = aCNMA, iCNMA = iCNMA, lnORs = lnORs)
}


read_cp_disc <- function(scenario, equal, tau2, p.baseline,
                         min.grpsize, max.grpsize, M, seq = NULL) {
  
  if (missing(seq))
    seq <- seq(1, M)
  
  path <- simpath(scenario, equal, tau2,
                  p.baseline, min.grpsize, max.grpsize)
  load(paste0(path, "settings.rda"))
  ##
  lnORs <- log(ORs)
  trts <- names(lnORs)
  
  ## Variable names for CPs and CI levels
  ##
  cps <- paste0("CP.", gsub("+", "", trts, fixed = TRUE))
  lows <- paste0("low.", gsub("+", "", trts, fixed = TRUE))
  upps <- paste0("upp.", gsub("+", "", trts, fixed = TRUE))
  
  ## Sceletions of data sets for results
  ## (variable CP contains the overall CP)
  ##
  NMA <- aCNMA <- iCNMA <-
    data.frame(network = "", model = "",
               scenario, equal, tau2, p.baseline,
               min.grpsize, max.grpsize,
               run = seq,
               CP = NA)
  
  ## Add variables for CPs and log odds ratio estimates
  ##
  for (j in seq_along(trts)) {
    NMA[[cps[j]]] <- aCNMA[[cps[j]]] <- iCNMA[[cps[j]]] <- NA
  }
  for (j in seq_along(trts)) {
    NMA[[upps[j]]] <- NMA[[lows[j]]] <- NA
    aCNMA[[upps[j]]] <- aCNMA[[lows[j]]] <- NA
    iCNMA[[upps[j]]] <- iCNMA[[lows[j]]] <- NA
  }
  
  ## Calculate overall CP
  ##
  for (i in seq) {
    path <- simpath(scenario, equal, tau2,
                    p.baseline, min.grpsize, max.grpsize,
                    i, M)
    load(paste0(path, "disc-modsel.rda"))
    ##
    if (ms.disc$selected == "aCNMA")
      iCNMA.i <- ms.disc$aCNMA
    else
      iCNMA.i <- ms.disc$iCNMA.selected
    ##
    for (j in seq_along(trts)) {
      ##
      ## Lower CI limits from NMA, aCNMA and selected iCNMA
      ##
      aCNMA[[lows[j]]][i] <- aCNMA.disc$lower.random[trts[j], "P"]
      iCNMA[[lows[j]]][i] <- iCNMA.i$lower.random[trts[j], "P"]
      ##
      ## Upper CI limits from NMA, aCNMA and selected iCNMA
      ##
      aCNMA[[upps[j]]][i] <- aCNMA.disc$upper.random[trts[j], "P"]
      iCNMA[[upps[j]]][i] <- iCNMA.i$upper.random[trts[j], "P"]
      ##
      ## Calculate CPs for individual parameters
      ##
      aCNMA[[cps[j]]][i] <-
        inCI(lnORs[trts[j]], aCNMA[[lows[j]]][i], aCNMA[[upps[j]]][i])
      iCNMA[[cps[j]]][i] <-
        inCI(lnORs[trts[j]], iCNMA[[lows[j]]][i], iCNMA[[upps[j]]][i])
    }
    ##
    ## CPs
    ##
    aCNMA$CP[i] <- sum(aCNMA[i, cps]) / length(trts)
    iCNMA$CP[i] <- sum(iCNMA[i, cps]) / length(trts)
  }
  
  
  ## Use same abbreviated names for true log odds ratios
  ##
  names(lnORs) <-
    paste0("lnOR.", gsub("+", "", names(lnORs), fixed = TRUE))
  
  
  aCNMA$network <- iCNMA$network <- "disconnected"
  ##
  aCNMA$model <- "aCNMA"
  iCNMA$model <- "iCNMA"
  
  
  list(aCNMA = aCNMA, iCNMA = iCNMA, lnORs = lnORs)
}


by_full <- function(x, varname, model = unique(x$model)) {
  
  model <- meta:::setchar(model, c("NMA", "aCNMA", "iCNMA"))
  ##
  for (i in seq_along(model)) {
    dat.i <- x[x$model == model[i], ]
    ##
    res.i <- data.frame(network = "connected", model = model[i],
                        scenario = rep(c("A", "B1", "B2", "C1", "C2"), 3),
                        het = factor(rep(c("no", "low", "moderate"), rep(5, 3)),
                                     levels = c("no", "low", "moderate")))
    ##
    res.i[[varname]] <-
      as.vector(by(dat.i[[varname]], list(dat.i$scenario, dat.i$tau2),
                   mean, na.rm = TRUE))
    res.i$sumNA <-
      as.vector(by(dat.i[[varname]], list(dat.i$scenario, dat.i$tau2),
                   sumNA))
    if (i == 1)
      res <- res.i
    else
      res <- rbind(res, res.i)
  }
  ##
  res
}


by_disc <- function(x, varname, model = unique(x$model)) {
  
  model <- meta:::setchar(model, c("aCNMA", "iCNMA"))
  ##
  for (i in seq_along(model)) {
    dat.i <- x[x$model == model[i], ]
    ##
    res.i <- data.frame(network = "disconnected", model = model[i],
                        scenario = rep(c("A", "B1", "B2", "C1", "C2"), 3),
                        het = factor(rep(c("no", "low", "moderate"), rep(5, 3)),
                                     levels = c("no", "low", "moderate")))
    ##
    res.i[[varname]] <-
      as.vector(by(dat.i[[varname]], list(dat.i$scenario, dat.i$tau2),
                   mean, na.rm = TRUE))
    res.i$sumNA <-
      as.vector(by(dat.i[[varname]], list(dat.i$scenario, dat.i$tau2),
                   sumNA))
    if (i == 1)
      res <- res.i
    else
      res <- rbind(res, res.i)
  }
  ##
  res
}


sumNA <- function(x)
  sum(is.na(x))


inCI <- function(x, lower, upper)
  1L * (lower < x & x < upper)


tolong <- function(x, varname,
                   trts = c("A", "B", "C", "D", "AB", "AC", "CD")) {
  
  vars <- paste(varname, trts, sep = ".")
  ##
  M <- nrow(x)
  n <- length(vars)
  ##
  res <- data.frame(network = "", model = "", scenario = "",
                    equal = NA, tau2 = NA, p.baseline = NA,
                    min.grpsize = NA, max.grpsize = NA, run = NA,
                    treat = rep("", n * M))
  res[[varname]] <- NA
  ##
  j <- 1
  for (i in seq_len(n)) {
    idx <- j:(j + M - 1)
    ##
    res[idx, names(res)[1:9]] <- x[, names(res)[1:9]]
    ##
    res$treat[idx] <- trts[i]
    res[[varname]][idx] <- x[[vars[i]]]
    ##
    j <- j + M
  }

  trts.long <-
    ifelse(nchar(trts) != 2, trts,
           paste(substring(trts, 1, 1), substring(trts, 2, 2), sep = "+"))
  ##
  res$treat <- factor(res$treat, levels = trts,
                      labels = trts.long)
  res
}


xyplot_full <- function(x)
    xyplot(MSE ~ I(as.numeric(het)) | I(tolower(scenario)),
           data = x, groups = model, type = "b",
           layout = c(5, 1),
           scales = list(
             x = list(at = 1:3, labels = c("no", "low", "moderate")),
             y = list(at = seq(0, 0.1, by = 0.02))),
           ylim = c(-0.005, 0.07),
           xlab = "", ylab = "Average Mean Squared Error",
           col = c("blue", "green", "red"), lty = 3:1, lwd = 2)
