separate <- function(x) {
  dat <- subset(x, subnet == "main")
  nc <- netconnection(treat1, treat2, studlab, data = dat)
  ##
  if (nc$n.subnets == 1)
    return(netmeta(TE, seTE, treat1, treat2, studlab, 
                   data = dat, sm = "RR", common = FALSE, 
                   reference.group = "plac"))
  else {
    dat <- subset(dat, treat1 != "apre")
    return(netmeta(TE, seTE, treat1, treat2, studlab, 
                   data = dat, sm = "RR", common = FALSE, 
                   reference.group = "plac"))
  }
}
