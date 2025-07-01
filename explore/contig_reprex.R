devtools::load_all(".")
library(geomander)
library(ggplot2)
library(ggredist)

# don't share this data please
map <- readRDS('explore/ex.rds')
map$geometry = rmapshaper::ms_simplify(map$geometry, 0.05, keep_shapes=TRUE)
init_plan = vctrs::vec_group_id(map$test)
V = nrow(map)

# init plan is contiguous
ccm(map$adj, map$test)
#>
#>    1
#> 9401

# run 30 steps MS
set.seed(2)
ms <- redist_mergesplit(map, warmup=0, nsims=200, counties = county, silent=T, silly_adj_fix=T)

# last plan is not contiguous
ccm(map$adj, last_plan(ms))
#> [1] 2

broken = which(check_contiguity(map$adj, last_plan(ms))$component == 2)
# plot(map, check_contiguity(map$adj, last_plan(ms))$component) + ggplot2::guides(fill="none")
plot(map, (as.matrix(ms)[, 22] == 11) + 2*(as.matrix(ms)[, 22] == 77)) +
    geom_district(aes(group=county), data=map, fill=NA, color="white", inherit.aes=FALSE) +
    guides(fill="none")
plot(map, (last_plan(ms) == 11) + 2*(last_plan(ms) == 77)) +
    geom_district(aes(group=county), data=map, fill=NA, color="white", inherit.aes=FALSE) +
    guides(fill="none")
plot(map, (init_plan == 11) + 2*(init_plan == 77)) +
    geom_district(aes(group=county), data=map, fill=NA, color="white", inherit.aes=FALSE) +
    guides(fill="none")
# plot(map, seq_len(V) %in% unique(unlist(map$adj[broken]) + 1) + 2*(seq_len(V) %in% broken)) + ggplot2::guides(fill="none")
plot(map, county)

# get components per column
mat <- get_plans_matrix(ms)
dim(mat)
# all(mat[, ncol(mat) - 1] == mat[, ncol(mat)])
# all(mat[, 1] == vctrs::vec_group_id(map$test))

vapply(seq_len(ncol(mat)), function(i) {
    ccm(map$adj, mat[, i])
}, numeric(1))

