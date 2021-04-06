#' Simulation Methods for Legislative Redistricting
#'
#' Enables researchers to sample redistricting plans from a pre-specified target
#' distribution using Sequential Monte Carlo and Markov Chain Monte Carlo
#' algorithms.  The package allows for the implementation of various constraints
#' in the redistricting process such as geographic compactness and population
#' parity requirements. Tools for analysis such as computation of various
#' summary statistics and plotting functionality are also included. The package
#' implements methods described in Fifield, Higgins, Imai and Tarr (2020)
#' <doi:10.1080/10618600.2020.1739532>, Fifield, Imai, Kawahara, and Kenny (2020)
#' <doi: 10.1080/2330443X.2020.1791773>, and McCartan and Imai (2020)
#' <arXiv:2008.06131>.
#'
#'
#' @name redist-package
#' @useDynLib redist, .registration = TRUE
#' @aliases redist-package redist
#' @docType package
#'
#' @references
#' Barbu, Adrian and Song-Chun Zhu. (2005) "Generalizing Swendsen-Wang to
#' Sampling Arbitrary Posterior Probabilities." IEEE Transactions on
#' Pattern Analysis and Machine Intelligence.
#'
#' Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2020) "Automated Redistricting Simulation Using Markov
#' Chain Monte Carlo."  Available at
#' \url{https://imai.fas.harvard.edu/research/files/redist.pdf}.
#'
#' Swendsen, Robert and Jian-Sheng Wang. (1987) "Nonuniversal Critical
#' Dynamics in Monte Carlo Simulations." Physical Review Letters.
#'
#' @keywords package
#' @importFrom Rcpp evalCpp
#' @importFrom coda mcmc mcmc.list traceplot autocorr.plot densplot cumuplot gelman.plot
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom grDevices dev.off pdf
#' @importFrom stats median runif quantile qnorm IQR
#' @importFrom dplyr n dplyr_row_slice dplyr_col_modify dplyr_reconstruct .data
#' @importFrom rlang :=
#' @importFrom stringr str_c str_glue
NULL

# for dplyr
utils::globalVariables(c("where", "."))
