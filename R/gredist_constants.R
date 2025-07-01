#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 01/2025
# Purpose: Important constants for package
####################################################


## Forward Kernel Splitting types
NAIVE_K_SPLITTING <- "top_k" # pick top k unif at random
UNIF_VALID_EDGE_SPLITTING <- "unif_valid" # pick valid edge unif
EXP_BIGGER_ABS_DEV_SPLITTING <- "exp_abs_dev"

VALID_FORWARD_KERNEL_TYPES <- c(NAIVE_K_SPLITTING, UNIF_VALID_EDGE_SPLITTING, EXP_BIGGER_ABS_DEV_SPLITTING)


VALID_FOREST_SPLITTING_METHODS <- c(UNIF_VALID_EDGE_SPLITTING, EXP_BIGGER_ABS_DEV_SPLITTING)


# Sampling Space
GRAPH_PLAN_SPACE_SAMPLING <- "graph_plan"
FOREST_SPACE_SAMPLING <- "spanning_forest"
LINKING_EDGE_SPACE_SAMPLING <- "linking_edge"

VALID_SAMPLING_SPACES <- c(GRAPH_PLAN_SPACE_SAMPLING, FOREST_SPACE_SAMPLING, LINKING_EDGE_SPACE_SAMPLING)


# Algorithm Types
SMC_ALG_TYPE <- "smc"
MCMC_ALG_TYPE <- "mergesplit"
MS_SMC_ALG_TYPE <- "smc_ms"
