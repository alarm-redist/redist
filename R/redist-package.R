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


# package constants
# Algorithm Types
SMC_ALG_TYPE <- "smc"
MCMC_ALG_TYPE <- "mergesplit"
MS_SMC_ALG_TYPE <- "smc_ms"
# Sampling Spaces
GRAPH_PLAN_SPACE_SAMPLING <- "graph_plan"
FOREST_SPACE_SAMPLING <- "spanning_forest"
LINKING_EDGE_SPACE_SAMPLING <- "linking_edge"
VALID_SAMPLING_SPACES <- c(GRAPH_PLAN_SPACE_SAMPLING, FOREST_SPACE_SAMPLING, LINKING_EDGE_SPACE_SAMPLING)
## Forward Kernel Splitting types
NAIVE_K_SPLITTING <- "top_k" # pick top k unif at random
UNIF_VALID_EDGE_SPLITTING <- "unif_valid" # pick valid edge unif
EXP_BIGGER_ABS_DEV_SPLITTING <- "exp_abs_dev" # exp(- alpha * abs max dev)
CONSTRAINT_SPLITTING <- "constraint"

VALID_FORWARD_KERNEL_TYPES <- c(NAIVE_K_SPLITTING, UNIF_VALID_EDGE_SPLITTING, EXP_BIGGER_ABS_DEV_SPLITTING, CONSTRAINT_SPLITTING)
VALID_FOREST_SPLITTING_METHODS <- c(UNIF_VALID_EDGE_SPLITTING, EXP_BIGGER_ABS_DEV_SPLITTING, CONSTRAINT_SPLITTING)

