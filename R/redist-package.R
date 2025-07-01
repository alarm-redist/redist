#' @keywords internal
#' @aliases redist-package redist
"_PACKAGE"

#' @useDynLib redist, .registration = TRUE
#'
#' @import redistmetrics
#' @importFrom Rcpp evalCpp
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom grDevices dev.off pdf
#' @importFrom stats sd var na.omit median runif quantile qnorm IQR optim splinefun qt
#' @importFrom utils str head tail packageVersion
#' @importFrom dplyr n dplyr_row_slice dplyr_col_modify dplyr_reconstruct .data
#' @importFrom cli cli_text cli_abort cli_warn cli_inform
#' @importFrom rlang :=
#' @importFrom stringr str_c str_glue
NULL

# for dplyr
utils::globalVariables(c("where", "."))
