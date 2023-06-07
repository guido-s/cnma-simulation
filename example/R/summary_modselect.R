summary_modselect <- function(x) {
  res <-
    with(x,
         data.frame(
           model = c("aCNMA", "iCNMA1", "iCNMA2", "iCNMA3", "iCNMA4"),
           ##
           Q = round(c(x$aCNMA$Q.additive, x$selected1$Q, x$selected2$Q,
                       x$selected3$Q, x$selected4$Q), 2),
           df.Q = c(x$aCNMA$df.Q.additive, x$selected1$df, x$selected2$df,
                    x$selected3$df, x$selected4$df),
           ##
           Q.diff = meta:::formatN(c(NA, round(Q, 2)), text.NA = "."),
           df.Q.diff = meta:::formatN(c(NA, df.Q), text.NA = "."),
           pval.Q.diff = meta:::formatPT(c(NA, pval.Q), lab.NA = "."),
           ##
           AIC = c(".", pval.Q < 0.157)))
  ##
  selmod <-
    as.numeric(substr(x$selected, nchar(x$selected), nchar(x$selected)))
  selmod <- if (is.na(selmod)) 1 else selmod + 1
  ##
  res$selected <- ""
  res$selected[selmod] <- "***"
  ##
  res
}
