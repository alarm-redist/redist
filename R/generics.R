####################################
# S3 generics

#' Extract the redistricting matrix from a \code{redist} object
#' @method as.matrix redist
#' @param x redist object
#' @param \dots additional arguments
#' @export
as.matrix.redist = function(x, ...) {
  x$plans
}

#' @method print redist
#' @importFrom utils str
#' @export
print.redist = function(x, ...) {
  cat(x$nsims, " sampled plans with ",
      ifelse(x$algorithm == 'mcmc', max(x$plans[,1]) + 1, max(x$plans[,1])),
      " districts from a ",
      length(x$adj), "-unit map, drawn\n using ",
      c(mcmc="Markov chain Monte Carlo",
        smc="Sequential Monte Carlo")[x$algorithm], sep="")
  if (x$pct_dist_parity < 1)
    cat(" and a ", 100*x$pct_dist_parity, "% population constraint.\n", sep="")
  else
    cat(".\n")
  cat(str(x$plans))
}
