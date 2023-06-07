getvar <- function(x, var1 = "selected", var2 = NULL) {
  if (is.null(var2))
    return(x[[var1]])
  else
    return(x[[var1]][[var2]])
}
