##
##
## Generate additional disconnected networks
## (starting from minimal set)
##
##

disconnect_additional <- function(TE, seTE, treat1, treat2, studlab,
                                  data = NULL, subset = NULL,
                                  minimal.set,
                                  verbose = FALSE) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chklogical(verbose)
  ##
  if (missing(minimal.set))
      stop("Argument 'minimal.set' is mandatory.",
           call. = FALSE)
  
  
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
    seTE <- TE$seTE
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    ##
    pairdata <- TE
    data <- TE
    ##
    TE <- TE$TE
  }
  else {
    is.pairwise <- FALSE
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
  ##
  ## (3) Use subset for analysis
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
  }
  
  
  ##
  ##
  ## (4) Additional checks
  ##
  ##
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').",
         call. = FALSE)
  
  
  ##
  ##
  ## (5) Generate addional disconnected data sets
  ##
  ##
  trts <- unique(c(treat1, treat2))
  n <- length(trts)
  ##
  diff <- setdiff(trts, minimal.set)
  
  data <- list()
  fulldata <- list()
  potential.set <- list()
  set <- list()
  ##
  s <- 0
  ##
  for (i in 1:(length(diff) - 2)) {
    for (j in 1:ncol(combn(diff, i))) {
      ##
      potential.set[[j]] <- c(minimal.set, combn(diff, i)[, j])
      disc <- disconnect(TE, seTE, treat1, treat2, studlab,
                         main.trts = potential.set[[j]])
      class(disc) <- c("pairwise", class(disc))
      ##
      disc.ma <- subset(disc, subnet %in% c("main", "auxiliary"))
      ##
      ## Only keep networks which (i) include all treatments and (ii)
      ## are disconnected
      ##
      nc.ij <- netconnection(disc.ma)
      ##
      if (nc.ij$n == n && nc.ij$n.subnets > 1) {
        s <- s + 1
        ##
        if (verbose) {
          if (s > 1)
            cat("\n")
          cat(paste("**********", s, "**********\n"))
          ## cat(paste0("*** i = ", i, ", j = ", j, " ***\n\n"))
          cat("Potential set:\n")
          print(potential.set[[j]])
          cat("\n")
          print(nc.ij, details.disconnected = TRUE)
        }
        ##
        attr(disc.ma, "k") <- nc.ij$k
        attr(disc.ma, "m") <- nc.ij$m
        attr(disc.ma, "n") <- nc.ij$n
        attr(disc.ma, "n.subnets") <- nc.ij$n.subnets
        ##
        disc.ma$id <- s + 1
        disc$id <- s + 1
        ##
        data[[s]] <- disc.ma           # dataset of disconnected network
        set[[s]] <- potential.set[[j]] # treatments in main subnetwork
        fulldata[[s]] <- disc          # full dataset of disconnected network
      }
    }
  }
  
  
  res <- list(data = data, set = set, fulldata = fulldata)
  res
}
