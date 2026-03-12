localize_tree <- function(tree, keep) {
    old_to_new <- integer(length(tree))
    old_to_new[keep] <- seq_along(keep)

    lapply(keep, function(old_u) {
        kids_old <- tree[[old_u]] + 1L
        kids_new <- old_to_new[kids_old]
        as.integer(kids_new[kids_new > 0] - 1L)
    })
}

bruteforce_tree_partitions <- function(tree, pop, l, lower, upper) {
    n <- length(tree)
    if (l == 1L) {
        total_pop <- sum(pop)
        return(as.integer(lower <= total_pop && total_pop <= upper))
    }

    edges <- do.call(rbind, lapply(seq_len(n), function(i) {
        kids <- tree[[i]] + 1L
        if (length(kids) == 0L) return(NULL)
        cbind(i, kids)
    }))

    undirected_components <- function(cut_rows) {
        adj <- vector("list", n)
        keep_rows <- setdiff(seq_len(nrow(edges)), cut_rows)
        for (row in keep_rows) {
            u <- edges[row, 1]
            v <- edges[row, 2]
            adj[[u]] <- c(adj[[u]], v)
            adj[[v]] <- c(adj[[v]], u)
        }

        seen <- rep(FALSE, n)
        comps <- list()
        for (start in seq_len(n)) {
            if (seen[start]) next
            stack <- start
            seen[start] <- TRUE
            comp <- integer()
            while (length(stack)) {
                u <- stack[[length(stack)]]
                stack <- stack[-length(stack)]
                comp <- c(comp, u)
                for (v in adj[[u]]) {
                    if (!seen[v]) {
                        seen[v] <- TRUE
                        stack <- c(stack, v)
                    }
                }
            }
            comps[[length(comps) + 1L]] <- comp
        }
        comps
    }

    combos <- combn(nrow(edges), l - 1L, simplify = FALSE)
    sum(vapply(combos, function(cut_rows) {
        comps <- undirected_components(cut_rows)
        if (length(comps) != l) return(FALSE)
        pops <- vapply(comps, function(idx) sum(pop[idx]), numeric(1))
        all(pops >= lower & pops <= upper)
    }, logical(1)))
}

test_that("single-tree partition counter matches brute force on the 4x4 MMSS region", {
    keep <- which(grid$init %in% 1:3)
    region_pop <- grid$pop[keep]
    target <- sum(region_pop) / 3
    delta <- get_pop_tol(grid)
    lower <- (1 - delta) * target
    upper <- (1 + delta) * target
    ignore <- !(seq_len(nrow(grid)) %in% keep)

    set.seed(1)
    for (i in 1:5) {
        tree <- sample_ust(
            get_adj(grid),
            grid[[attr(grid, "pop_col")]],
            0,
            sum(region_pop),
            rep(1L, nrow(grid)),
            ignore
        )
        local_tree <- localize_tree(tree, keep)

        brute <- bruteforce_tree_partitions(local_tree, region_pop, 3L, lower, upper)
        enum_count <- redist:::count_single_tree_partitions_impl(
            local_tree, region_pop, 3L, lower, upper, "enumeration"
        )$count
        dp_count <- redist:::count_single_tree_partitions_impl(
            local_tree, region_pop, 3L, lower, upper, "tree_dp"
        )$count

        expect_equal(enum_count, brute)
        expect_equal(dp_count, brute)
    }
})

test_that("diag_single_tree_partitions runs on the shared 4x4 grid fixture", {
    set.seed(1)
    out <- diag_single_tree_partitions(grid, l = 3, n_plans = 4, n_trees_per_region = 8)

    expect_named(out, c("hit_rate", "mean_m", "max_m", "dist_m", "cv2_Z"))
    expect_true(is.numeric(out$hit_rate) && is.finite(out$hit_rate))
    expect_true(is.numeric(out$mean_m) && is.finite(out$mean_m))
    expect_true(is.numeric(out$max_m) && is.finite(out$max_m))
    expect_true(is.integer(out$dist_m))
    expect_true(length(out$dist_m) >= 1)
    expect_true(all(out$dist_m >= 0L))
    expect_true(is.numeric(out$cv2_Z))
})

test_that("diag_single_tree_partitions runs on Iowa", {
    set.seed(1)
    out <- diag_single_tree_partitions(ia, l = 3, n_plans = 2, n_trees_per_region = 4)

    expect_named(out, c("hit_rate", "mean_m", "max_m", "dist_m", "cv2_Z"))
    expect_true(is.numeric(out$hit_rate) && is.finite(out$hit_rate))
    expect_true(is.numeric(out$mean_m) && is.finite(out$mean_m))
    expect_true(is.numeric(out$max_m) && is.finite(out$max_m))
    expect_true(is.integer(out$dist_m))
    expect_true(length(out$dist_m) >= 1)
})
