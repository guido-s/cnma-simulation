disconnect <- function(TE, seTE, treat1, treat2, studlab,
                       data = NULL, subset = NULL,
                       main.trts,
                       warn = TRUE) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  meta:::chklogical(warn)
  ##
  if (missing(main.trts))
      stop("Argument 'main.trts' is mandatory.",
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
  main.trts <-
    meta:::setchar(main.trts, sort(unique(c(treat1, treat2))),
                   paste("must contain treatment names provided in",
                         "arguments 'treat1' and 'treat2'"))
  
  ## Check NAs and zero standard errors
  ##
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  ##
  if (any(excl)) {
    ##
    studlab <- studlab[!(excl)]
    treat1  <- treat1[!(excl)]
    treat2  <- treat2[!(excl)]
    TE      <- TE[!(excl)]
    seTE    <- seTE[!(excl)]
  }
  ##
  ##
  ## (5) Define main and auxiliary subnetworks
  ##
  ##
  sel1 <- treat1 %in% main.trts
  sel2 <- treat2 %in% main.trts
  ##
  res <- data.frame(studlab, TE, seTE, treat1, treat2)
  ##
  res$comp <-
    ifelse(sel1 & sel2, "main",
           ifelse(!sel1 & !sel2, "auxiliary", "ignore"))
  ##
  ## Ignore studies with comparisons in more than one group
  ##
  sel.main <-
    res$studlab %in% unique(res$studlab[res$comp == "main"])
  sel.ignore <-
    res$studlab %in% unique(res$studlab[res$comp == "ignore"])
  sel.auxiliary <-
    res$studlab %in% unique(res$studlab[res$comp == "auxiliary"])
  ##
  res$subnet <- res$comp
  ##
  res$subnet[sel.main & (sel.ignore | sel.auxiliary)] <- "ignore"
  res$subnet[sel.auxiliary & (sel.main | sel.ignore)] <- "ignore"
  

  sel.main <- res$subnet == "main"
  attr(res, "main.trts") <-
    unique(sort(c(res$treat1[sel.main], res$treat2[sel.main])))
  ##
  sel.aux <- res$subnet == "auxiliary"
  attr(res, "aux.trts") <-
    unique(sort(c(res$treat1[sel.aux], res$treat2[sel.aux])))
  ##
  class(res) <- c("pairwise", class(res))
  ##
  nc <- netconnection(res)
  attr(res, "k") <- nc$k
  attr(res, "m") <- nc$m
  attr(res, "n") <- nc$n
  attr(res, "n.subnets") <- nc$n.subnets
  ##
  res
}
