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
    if (!meta:::is.wholenumber(n))
      stop("Argument 'n' must be a whole number.")
    sel <- unlist(lapply(strsplit(combs, x$sep.comps,
                                  fixed = TRUE), length)) == n
    combs <- combs[sel]
  }
  
  combs
}
