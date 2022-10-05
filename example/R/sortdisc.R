sortdisc <- function(init, add, verbose = FALSE) {
  
  max.m <- function(x, var = "subnet") max(table(x[[var]]))
  ##
  getattr <- function(x, attr = "n") attributes(x)[[attr]]
  ##
  sumnets <- function(x) {
    nc <- netconnection(x)
    data.frame(m = nc$m, k = nc$k, n = nc$n, n.subnets = nc$n.subnets)
  }
  
  
  init.ma <- subset(init, subnet %in% c("main", "auxiliary"))
  nc.init <- netconnection(init.ma)
  attr(init.ma, "m") <- nc.init$m
  attr(init.ma, "k") <- nc.init$k
  attr(init.ma, "n.subnets") <- nc.init$n.subnets
  
  disc.m <-
    c(getattr(init.ma, "m"),
      unlist(lapply(add$data, getattr, "m")))
  disc.k <-
    c(getattr(init.ma, "k"),
      unlist(lapply(add$data, getattr, "k")))
  disc.max.m.subnet <-
    c(max.m(init.ma),
      unlist(lapply(add$data, max.m)))
  ##
  disc.d <- data.frame(disc.m, disc.k, disc.max.m.subnet)
  if (verbose)
    print(disc.d)
  ##
  o <- with(disc.d, rev(order(disc.m, disc.k, disc.max.m.subnet)))
  if (verbose)
    print(disc.d[o, ])
  ##
  list.discdata.unsorted <- c(list(init.ma), add$data)
  list.discdata <- vector(mode = "list", length = length(o))
  j <- 0
  ##
  for (i in o) {
    j <- j + 1
    list.discdata[[j]] <- list.discdata.unsorted[[i]]
  }
  ##
  list.set.unsorted <- c(list(minset), add$set)
  list.set <- vector(mode = "list", length = length(o))
  j <- 0
  ##
  for (i in o) {
    j <- j + 1
    list.set[[j]] <- list.set.unsorted[[i]]
  }
  ##
  list.fulldata.unsorted <- c(list(init), add$fulldata)
  list.fulldata <- vector(mode = "list", length = length(o))
  j <- 0
  ##
  for (i in o) {
    j <- j + 1
    list.fulldata[[j]] <- list.fulldata.unsorted[[i]]
  }
  
  
  ##
  ## Print new order of disconnected networks
  ## (initial disconnected network with minimal set is 8th disconnected
  ##  network)
  ##
  if (verbose)
    print(lapply(list.discdata, sumnets))
  
  
  ## Data sets for disconnected networks
  ##
  res <- list(data = list.discdata,
              set = list.set, fulldata = list.fulldata)
  ##
  names(res$data) <- names(res$set) <-
    names(res$fulldata) <- seq_len(length(res$data))
  ##
  res
}
