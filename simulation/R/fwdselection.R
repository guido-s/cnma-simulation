fwdselection <- function(x, sel.pval = 0.157,
                         sep.interactions = "*",
                         reorder = FALSE,
                         keepall = FALSE, verbose = FALSE) {
  
  
  meta:::chkclass(x, c("netcomb", "discomb"))
  meta:::chklevel(sel.pval, ci = FALSE)
  meta:::chklogical(keepall)
  meta:::chklogical(verbose)
  meta:::chkchar(sep.interactions, length = 1)
  ##
  subnets <-
    netconnection(x$treat1, x$treat2, x$studlab)$n.subnets
  ##
  if (verbose) {
    cat("Number of subnetworks:", subnets, "\n")
    cat("Number of studies:", x$k,"\n")
    cat("Number of pairwise comparisons:", x$m,"\n\n")
  }
  
  
  ## Determine and estimate all interaction CNMA models with one, two
  ## or three 2-way interactions
  ##
  iCNMAs <- icnma_2way(x, sep.interactions = sep.interactions)
  ##
  ## CNMA model selection
  ##
  ms <- icnma_fwdsel(iCNMAs, x, sel.pval = sel.pval, reorder = reorder)
  
  
  res <- list(aCNMA = x, subnets = subnets, k = x$k, m = x$m,
              ##
              iCNMA.selected = ms$iCNMA.selected[[1]],
              selected = ms$selected,
              Q = ms$Q, df.Q = ms$df.Q, pval.Q = ms$pval.Q,
              ##
              selected1 = ms$selected1,
              selected2 = ms$selected2,
              selected3 = ms$selected3,
              ##
              hetstats1 = iCNMAs$hetstats1,
              hetstats2 = iCNMAs$hetstats2,
              hetstats3 = iCNMAs$hetstats3)
  ##
  if (keepall) {
    res$iCNMA1 <- iCNMAs$iCNMA1
    res$iCNMA2 <- iCNMAs$iCNMA2
    res$iCNMA3 <- iCNMAs$iCNMA3
  }
  ##
  res
}


icnma_fwdsel <- function(x, additive, sel.pval, reorder) {
  
  
  ms1 <- x[["hetstats1"]]
  ms2 <- x[["hetstats2"]]
  ms3 <- x[["hetstats3"]]
  ##
  Q.diff <- df.Q.diff <- pval.Q.diff <- rep(NA, 3)
  
  
  ## Q and df for additive CNMA model
  ##
  Q.add <- additive$Q.additive
  df.Q.add <- additive$df.Q.additive 
  
  
  ##
  ## Models with one 2-way interaction
  ##
  if (all(is.na(ms1$Q))) {
    selected <- "aCNMA"
    iCNMA.selected <- additive
  }
  else {
    ms1 <- ms1[!is.na(ms1$Q), ]
    ##
    sel1 <- ms1$Q == min(ms1$Q)
    ##
    ## Randomly select interaction if p-values are identical
    ##
    ms1 <- ms1[sel1, ]
    if (sum(sel1) > 1)
      ms1 <- ms1[sample(seq_len(sum(sel1)), 1), ]
    ##
    ## Additive vs one 2-way interaction
    ##
    if (Q.add > ms1$Q) {
      Q.diff[1] <- Q.add - ms1$Q
      df.Q.diff[1] <- df.Q.add - ms1$df
      pval.Q.diff[1] <- pchisq(Q.diff[1], df.Q.diff[1], lower = FALSE)
    }
    ##
    ## Step 1 of model selection (one 2-way interaction vs additive)
    ##
    if (is.na(pval.Q.diff[1]) ||
        df.Q.diff[1] == 0 ||
        pval.Q.diff[1] >= sel.pval) {
      selected <- "aCNMA"
      iCNMA.selected <- additive
    }
    else {
      ##
      ## Only consider models with first selected 2-way interaction
      ##
      ms2 <- ms2[ms2$interaction1 == ms1$interaction |
                 ms2$interaction2 == ms1$interaction, ]
      ##
      if (all(is.na(ms2$Q))) {
        selected <- "iCNMA1"
        iCNMA.selected <- x$iCNMA1[ms1$id]
      }
      else {
        ms2 <- ms2[!is.na(ms2$Q), ]
        ##
        sel2 <- ms2$Q == min(ms2$Q)
        ##
        ## Randomly select interaction if Qs are identical
        ##
        ms2 <- ms2[sel2, ]
        if (sum(sel2) > 1)
          ms2 <- ms2[sample(seq_len(sum(sel2)), 1), ]
        ##
        ## Step 2 of model selection (two vs one 2-way interactions)
        ##
        Q.diff[2] <- ms1$Q - ms2$Q
        df.Q.diff[2] <- ms1$df - ms2$df
        pval.Q.diff[2] <- pchisq(Q.diff[2], df.Q.diff[2], lower = FALSE)
        ##
        if (is.na(pval.Q.diff[2]) ||
            df.Q.diff[2] == 0 ||
            pval.Q.diff[2] >= sel.pval) {
          selected <- "iCNMA1"
          iCNMA.selected <- x$iCNMA1[ms1$id]
        }
        else {
          if (all(is.na(ms3$Q))) {
            selected <- "iCNMA2"
            iCNMA.selected <- x$iCNMA2[ms2$id]
          }
          ##
          ## Only consider models with first two selected 2-way
          ## interactions
          ##
          ms3 <- ms3[(ms3$interaction1 == ms2$interaction1 |
                      ms3$interaction2 == ms2$interaction1 |
                      ms3$interaction3 == ms2$interaction1) &
                     (ms3$interaction1 == ms2$interaction2 |
                      ms3$interaction2 == ms2$interaction2 |
                      ms3$interaction3 == ms2$interaction2), ]
          ##
          ms3 <- ms3[!is.na(ms3$Q), ]
          ##
          sel3 <- ms3$Q == min(ms3$Q)
          ##
          ## Randomly select interaction if Qs are identical
          ##
          ms3 <- ms3[sel3, ]
          if (sum(sel3) > 1)
            ms3 <- ms3[sample(seq_len(sum(sel3)), 1), ]
          ##
          ## Step 3 of model selection (three vs two 2-way interactions)
          ##
          Q.diff[3] <- ms2$Q - ms3$Q
          df.Q.diff[3] <- ms2$df - ms3$df
          pval.Q.diff[3] <- pchisq(Q.diff[3], df.Q.diff[3], lower = FALSE)
          ##
          if (is.na(pval.Q.diff[3]) ||
              df.Q.diff[3] == 0 ||
              pval.Q.diff[3] >= sel.pval) {
            selected <- "iCNMA2"
            iCNMA.selected <- x$iCNMA2[ms2$id]
          }
          else {
            selected <- "iCNMA3"
            iCNMA.selected <- x$iCNMA3[ms3$id]
          }
        }
      }
    }
  }
  
  
  if (reorder) {
    ## Change order of interaction terms
    ##
    if (is.data.frame(ms2) && ms2$interaction1 != ms1$interaction) {
      ms2$interaction2 <- ms2$interaction1
      ms2$interaction1 <- ms1$interaction
    }
    ##
    if (is.data.frame(ms3)) {
      keep <- ""
      if (ms3$interaction1 != ms2$interaction1) {
        if (ms3$interaction1 != ms2$interaction2)
          keep <- ms3$interaction1
        ms3$interaction1 <- ms2$interaction1
      }
      if (ms3$interaction2 != ms2$interaction2) {
        if (ms3$interaction2 != ms2$interaction1)
          keep <- ms3$interaction2
        ms3$interaction2 <- ms2$interaction2
      }
      if (keep != "")
        ms3$interaction3 <- keep
    }
  }
  
  
  res <- list(iCNMA.selected = iCNMA.selected,
              selected = selected,
              Q = Q.diff, df.Q = df.Q.diff, pval.Q = pval.Q.diff,
              selected1 = ms1, selected2 = ms2, selected3 = ms3)
  ##
  res
}


icnma_2way <- function(TE, seTE,
                       treat1, treat2, studlab,
                       data = NULL, subset = NULL,
                       sm,
                       fixed = FALSE,
                       random = TRUE,
                       ##
                       reference.group = "",
                       ##
                       sep.comps = "+",
                       sep.interactions = "*",
                       ##
                       n1 = NULL,
                       n2 = NULL,
                       event1 = NULL,
                       event2 = NULL,
                       ##
                       keepdata = gs("keepdata"),
                       warn = TRUE,
                       verbose = FALSE,
                       debug = FALSE
                       ) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chkchar(sep.comps, nchar = 1)
  ##
  meta:::chklogical(fixed)
  meta:::chklogical(random)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch TE, treat1, treat2, seTE, studlab from data:
  ##
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  if (inherits(TE, "pairwise")) {
    is.pairwise <- TRUE
    ##
    sm <- attr(TE, "sm")
    ##
    seTE <- TE$seTE
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    ##
    if (!is.null(TE$n1))
      n1 <- TE$n1
    if (!is.null(TE$n2))
      n2 <- TE$n2
    if (!is.null(TE$event1))
      event1 <- TE$event1
    if (!is.null(TE$event2))
      event2 <- TE$event2
    ##
    pairdata <- TE
    data <- TE
    ##
    TE <- TE$TE
  }
  else if (inherits(TE, c("netmeta", "netcomb", "discomb"))) {
    is.pairwise <- TRUE
    ##
    sm <- TE$sm
    ##
    if (missing(reference.group))
      reference.group <- TE$reference.group
    ##
    if (missing(fixed))
      fixed <- TE$fixed
    ##
    if (missing(random))
      random <- TE$random
    ##
    if (missing(random))
      random <- TE$random
    ##
    if (missing(sep.comps))
      sep.comps <- TE$sep.comps
    ##
    seTE <- TE$seTE
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    ##
    if (!is.null(TE$n1))
      n1 <- TE$n1
    if (!is.null(TE$n2))
      n2 <- TE$n2
    if (!is.null(TE$event1))
      event1 <- TE$event1
    if (!is.null(TE$event2))
      event2 <- TE$event2
    ##
    pairdata <- TE$data
    data <- TE$data
    ##
    TE <- TE$TE
  }
  else {
    is.pairwise <- FALSE
    if (missing(sm))
      if (!is.null(data) && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
    ##
    seTE <- eval(mf[[match("seTE", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
    ##
    treat1 <- eval(mf[[match("treat1", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    treat2 <- eval(mf[[match("treat2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    studlab <- eval(mf[[match("studlab", names(mf))]],
                    data, enclos = sys.frame(sys.parent()))
    ##
    n1 <- eval(mf[[match("n1", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    n2 <- eval(mf[[match("n2", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
    ##
    event1 <- eval(mf[[match("event1", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
    ##
    event2 <- eval(mf[[match("event2", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
  }
  ##
  meta:::chknumeric(TE)
  meta:::chknumeric(seTE)
  ##
  if (!any(!is.na(TE) & !is.na(seTE)))
    stop("Missing data for estimates (argument 'TE') and ",
         "standard errors (argument 'seTE') in all studies.\n  ",
         "No network meta-analysis possible.",
         call. = FALSE)
  ##
  k.Comp <- length(TE)
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  ##
  if (length(studlab) == 0) {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.",
              call. = FALSE)
    studlab <- seq(along = TE)
  }
  studlab <- as.character(studlab)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  ##
  if (!is.null(event1) & !is.null(event2))
    available.events <- TRUE
  else
    available.events <- FALSE
  ##
  if (!is.null(n1) & !is.null(n2))
    available.n <- TRUE
  else
    available.n <- FALSE
  
  
  ##
  ##
  ## (3) Store complete dataset in list object data
  ##     (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
    if (nulldata & !is.pairwise)
      data <- data.frame(.studlab = studlab, stringsAsFactors = FALSE)
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.studlab <- studlab
    }
    else
      data$.studlab <- studlab
    ##
    data$.order <- seq_along(studlab)
    ##
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    ##
    data$.TE <- TE
    data$.seTE <- seTE
    ##
    data$.event1 <- event1
    data$.n1 <- n1
    data$.event2 <- event2
    data$.n2 <- n2
    ##
    ## Check for correct treatment order within comparison
    ##
    wo <- data$.treat1 > data$.treat2
    ##
    if (any(wo)) {
      data$.TE[wo] <- -data$.TE[wo]
      ttreat1 <- data$.treat1
      data$.treat1[wo] <- data$.treat2[wo]
      data$.treat2[wo] <- ttreat1[wo]
      ##
      if (meta:::isCol(data, ".n1") & meta:::isCol(data, ".n2")) {
        tn1 <- data$.n1
        data$.n1[wo] <- data$.n2[wo]
        data$.n2[wo] <- tn1[wo]
      }
      ##
      if (meta:::isCol(data, ".event1") & meta:::isCol(data, ".event2")) {
        tevent1 <- data$.event1
        data$.event1[wo] <- data$.event2[wo]
        data$.event2[wo] <- tevent1[wo]
      }
    }
    ##
    if (!missing.subset) {
      if (length(subset) == nrow(data))
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
  }
  
  
  ##
  ##
  ## (4) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    if ((is.logical(subset) & (sum(subset) > k.Comp)) ||
        (length(subset) > k.Comp))
      stop("Length of subset is larger than number of studies.",
           call. = FALSE)
    ##
    TE <- TE[subset]
    seTE <- seTE[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
    ##
    if (!is.null(n1))
      n1 <- n1[subset]
    if (!is.null(n2))
      n2 <- n2[subset]
    if (!is.null(event1))
      event1 <- event1[subset]
    if (!is.null(event2))
      event2 <- event2[subset]
  }
  ##
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (reference.group != "")
    reference.group <- netmeta:::setref(reference.group, labels)
  ##
  ##
  meta:::chkchar(sep.interactions, length = 1)
  compmatch <- netmeta:::compmatch
  ##
  if (compmatch(labels, sep.interactions)) {
    if (!missing(sep.interactions))
      warning("Interaction separator '", sep.interactions,
              "' used in at least one treatment label. ",
              "Try to use predefined separators: ",
              "':', '-', '_', '/', '+', '.', '|', '*'.",
              call. = FALSE)
    ##
    if (!compmatch(labels, ":"))
      sep.interactions <- ":"
    else if (!compmatch(labels, "-"))
      sep.interactions <- "-"
    else if (!compmatch(labels, "_"))
      sep.interactions <- "_"
    else if (!compmatch(labels, "/"))
      sep.interactions <- "/"
    else if (!compmatch(labels, "+"))
      sep.interactions <- "+"
    else if (!compmatch(labels, "."))
      sep.interactions <- "-"
    else if (!compmatch(labels, "|"))
      sep.interactions <- "|"
    else if (!compmatch(labels, "*"))
      sep.interactions <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', ",
           "'.', '|', '*') are used in at least one treatment label.",
           "\n   Please specify a different character that should be used ",
           "as separator (argument 'sep.interactions').",
           call. = FALSE)
  }



  ##
  ##
  ## (5) Check network connectivity and conduct CNMA
  ##
  ##
  n.subsets <- netconnection(treat1, treat2, studlab)$n.subnets
  connected <- n.subsets == 1
  ##
  if (!connected) {
    ##
    ## Additive component network meta-analysis model (CNMA) for
    ## disconnected network
    ##
    net.additive <- discomb(TE, seTE, treat1, treat2, studlab,
                            sm = sm,
                            fixed = fixed,
                            random = random)
  }
  else {
    ##
    ## (C)NMA for connected network
    ##

    ## Standard NMA model
    ##
    NMA <- netmeta(TE, seTE, treat1, treat2, studlab,
                    reference.group = reference.group,
                    sm = sm,
                    fixed = fixed,
                    random = random)

    ## Additive component network meta-analysis model
    ##
    net.additive <- netcomb(NMA)
  }
  
  
  ##
  ##
  ## (6) CNMA models with one 2-way interaction
  ##
  ##
  selcomp <- function(x, n = 1, sep = "+") {
    res <- lapply(strsplit(x, sep, fixed = TRUE),
                  gsub, pattern = " ", replacement = "")
    res <- sapply(res, "[", n)
    res
  }
  ##
  comb2.1 <- combinations(net.additive, n = 2)
  if (verbose)
    print(comb2.1)
  ##  
  C.matrix <- net.additive$C.matrix
  ##
  if (length(comb2.1) == 0)
    hetstats1 <- NULL
  else {
    ##
    ## Compute Q values for all observed 2-way interactions
    ##
    hetstats1 <-
      data.frame(interaction = rep_len("", length(comb2.1)),
                 Q = NA, df = NA, pval = NA,
                 id = seq_len(length(comb2.1)))
    ##
    ## Conduct CNMAs
    ##
    iCNMA1 <- C1 <- list()
    j <- 0
    ##
    for (i in comb2.1) {
      if (debug)
        print(i)
      ##
      j <- j + 1
      comp1.i <- selcomp(i, 1)
      comp2.i <- selcomp(i, 2)
      comp12 <- paste(comp1.i, comp2.i, sep = sep.interactions)
      ##
      C1.i <- cbind(C.matrix,
                     C.matrix[, comp1.i] * C.matrix[, comp2.i])
      colnames(C1.i) <- c(colnames(C.matrix), comp12)
      ##
      if (debug)
        print(colnames(C1.i))
      ##
      if (!connected)
        net.i <- discomb(TE, seTE, treat1, treat2, studlab,
                         C.matrix = C1.i, sm = sm)
      else
        net.i <- netcomb(NMA, C.matrix = C1.i)
      ##
      hetstats1$interaction[j] <- comp12
      ##
      if ((net.i$df.Q.additive < net.additive$df.Q.additive)) {
        hetstats1[j, 2:4] <-
          c(net.i$Q.additive, net.i$df.Q.additive, net.i$pval.Q.additive)
        iCNMA1[[j]] <- net.i
        C1[[j]] <- C1.i
        ##
        if (verbose)
          print(subset(hetstats1, hetstats1$interaction != ""))
      }
    }
  }
  ##
  ## Find 2-way interactions
  ##
  int2.1 <- hetstats1$interaction[which(!is.na(hetstats1$Q))]
  
  
  ##
  ##
  ## (7) CNMA models with two 2-way interactions
  ##
  ##
  selint <- function(x, n = 1, sep = "*") {
    res <- lapply(strsplit(x, sep, fixed = TRUE),
                  gsub, pattern = " ", replacement = "")
    res <- sapply(res, "[", n)
    res
  }
  ##  
  if (length(int2.1) < 2) {
    hetstats2 <- NA
    iCNMA2 <- NA
  }
  else{
    comb2.2 <- combn(int2.1, 2) # combinations of two 2-way interactions
    ##
    C2 <- iCNMA2 <- list()
    ##
    hetstats2 <-
      data.frame(interaction1 = rep("", ncol(comb2.2)),
                 interaction2 = rep("", ncol(comb2.2)),
                 Q = NA, df = NA, pval = NA,
                 id = seq_len(ncol(comb2.2)))
    ##
    for (i in seq_len(ncol(comb2.2))) {
      j <- comb2.2[2, i]
      comp1.j <- selint(j, 1)
      comp2.j <- selint(j, 2)
      comp12 <- paste(comp1.j, comp2.j, sep = sep.interactions)
      ##
      k <- comb2.2[1, i]
      index.int <- which(hetstats1$interaction == k)
      ##
      C2[[i]] <-
        cbind(C1[[index.int]],
              C.matrix[, comp1.j] * C.matrix[, comp2.j])
      colnames(C2[[i]]) <-
        c(colnames(C1[[index.int]]), comp12)
      ##
      if (debug)
        print(c(colnames(C1[[index.int]]), comp12))
      ##
      ## Conduct CNMA
      ##
      if (!connected)
        iCNMA2[[i]] <-
          discomb(TE, seTE, treat1, treat2, studlab,
                  C.matrix = C2[[i]], sm = sm)

      else
        iCNMA2[[i]] <-
          netcomb(NMA,
                  C.matrix = C2[[i]],
                  fixed = fixed, random = random)
      ##
      hetstats2$interaction1[i] <- k
      hetstats2$interaction2[i] <- j
      ##
      if ((iCNMA2[[i]]$df.Q.additive < net.additive$df.Q.additive)) {
        hetstats2[i, 3:5] <-
          c(iCNMA2[[i]]$Q.additive,
            iCNMA2[[i]]$df.Q.additive,
            iCNMA2[[i]]$pval.Q.additive)
      }
      ##
      if (verbose)
        print(subset(hetstats2, hetstats2$interaction1 != ""))
    }
  }
  
  
  ##
  ##
  ## (8) CNMA models with three 2-way interactions
  ##
  ##
  if (length(int2.1) < 3) {
    hetstats3 <- NA
    iCNMA3 <- NA
  }
  else {
    comb2.3 <- combn(int2.1, 3) # combinations of three 2-way interactions
    ##
    C3 <- iCNMA3 <- list()
    ##
    hetstats3 <-
      data.frame(interaction1 = rep("", ncol(comb2.3)),
                 interaction2 = rep("", ncol(comb2.3)),
                 interaction3 = rep("", ncol(comb2.3)),
                 Q = NA, df = NA, pval = NA,
                 id = seq_len(ncol(comb2.3)))
    ##
    for (i in seq_len(ncol(comb2.3))) {
      j <- comb2.3[2, i]
      comp1.j <- selint(j, 1)
      comp2.j <- selint(j, 2)
      comp12j <- paste(comp1.j, comp2.j, sep = sep.interactions)
      ##
      m <- comb2.3[3, i]
      comp1.m <- selint(m, 1)
      comp2.m <- selint(m, 2)
      comp12m <- paste(comp1.m, comp2.m, sep = sep.interactions)
      ##
      k <- comb2.3[1, i]
      index.int <- which(hetstats1$interaction == k)
      ##
      C3[[i]] <-
        cbind(C1[[index.int]],
              C.matrix[, comp1.j] * C.matrix[, comp2.j],
              C.matrix[, comp1.m] * C.matrix[, comp2.m])
      colnames(C3[[i]]) <-
        c(colnames(C1[[index.int]]), comp12j, comp12m)
      ##
      if (debug)
        print(c(colnames(C1[[index.int]]), comp12j, comp12m))
      ##
      ## Conduct CNMA
      ##
      if (!connected)
        iCNMA3[[i]] <-
          discomb(TE, seTE, treat1, treat2, studlab,
                  C.matrix = C3[[i]], sm = sm)
      else
        iCNMA3[[i]] <-
          netcomb(NMA,
                  C.matrix = C3[[i]],
                  fixed = fixed, random = random)
      ##
      hetstats3$interaction1[i] <- k
      hetstats3$interaction2[i] <- j
      hetstats3$interaction3[i] <- m
      ##
      if ((iCNMA3[[i]]$df.Q.additive < net.additive$df.Q.additive)) {
        hetstats3[i, 4:6] <-
          c(iCNMA3[[i]]$Q.additive,
            iCNMA3[[i]]$df.Q.additive,
            iCNMA3[[i]]$pval.Q.additive)
      }
      ##
      if (verbose)
        print(subset(hetstats3, hetstats3$interaction1 != ""))
    }
  }
  
  
  res <- list(hetstats1 = hetstats1, iCNMA1 = iCNMA1,
              hetstats2 = hetstats2, iCNMA2 = iCNMA2,
              hetstats3 = hetstats3, iCNMA3 = iCNMA3)
  ##
  res
}
