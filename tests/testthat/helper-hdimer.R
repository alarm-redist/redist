matching_oracle <- function(n_vertices, edge_u, edge_v, weights) {
    recurse <- function(vertices, available) {
        if (length(vertices) == 0L) {
            return(list(list(edges = integer(), weight = 1)))
        }
        first <- vertices[[1L]]
        incident <- available[edge_u[available] == first | edge_v[available] == first]
        unlist(
            lapply(incident, function(edge) {
                other <- edge_v[[edge]]
                if (edge_u[[edge]] != first) {
                    other <- edge_u[[edge]]
                }
                if (!(other %in% vertices)) {
                    return(list())
                }
                remaining_vertices <- setdiff(vertices, c(first, other))
                remaining_edges <- available[
                    !(edge_u[available] %in% c(first, other)) &
                        !(edge_v[available] %in% c(first, other))
                ]
                lapply(recurse(remaining_vertices, remaining_edges), function(rest) {
                    list(
                        edges = c(edge, rest$edges),
                        weight = weights[[edge]] * rest$weight
                    )
                })
            }),
            recursive = FALSE
        )
    }

    recurse(seq_len(n_vertices), seq_along(weights))
}

plan_partition_key <- function(plan) {
    same <- outer(plan, plan, '==')
    paste(as.integer(same[upper.tri(same)]), collapse = '')
}

vertices_connected_by_edges <- function(vertices, edge_u, edge_v) {
    if (length(vertices) <= 1L) {
        return(TRUE)
    }
    reached <- vertices[[1L]]
    frontier <- reached
    while (length(frontier) > 0L) {
        incident <- edge_u %in% frontier | edge_v %in% frontier
        neighbors <- unique(c(edge_u[incident], edge_v[incident]))
        neighbors <- intersect(neighbors, vertices)
        neighbors <- setdiff(neighbors, reached)
        reached <- c(reached, neighbors)
        frontier <- neighbors
    }
    length(reached) == length(vertices)
}

hierarchical_tree_count <- function(adj, vertices, hierarchy) {
    if (length(vertices) <= 1L) {
        return(1)
    }
    edge_u <- rep(seq_along(adj) - 1L, lengths(adj))
    edge_v <- unlist(adj, use.names = FALSE)
    keep <- edge_u < edge_v & edge_u %in% vertices & edge_v %in% vertices
    edge_u <- edge_u[keep]
    edge_v <- edge_v[keep]
    if (length(edge_u) < length(vertices) - 1L) {
        return(0)
    }
    candidates <- combn(
        seq_along(edge_u),
        length(vertices) - 1L,
        simplify = FALSE
    )
    sum(vapply(
        candidates,
        \(edges) {
            tree_u <- edge_u[edges]
            tree_v <- edge_v[edges]
            if (!vertices_connected_by_edges(vertices, tree_u, tree_v)) {
                return(FALSE)
            }
            all(vapply(
                seq_len(ncol(hierarchy)),
                \(level) {
                    units <- unique(hierarchy[vertices + 1L, level])
                    all(vapply(
                        units,
                        \(unit) {
                            unit_vertices <- vertices[hierarchy[vertices + 1L, level] == unit]
                            vertices_connected_by_edges(unit_vertices, tree_u, tree_v)
                        },
                        logical(1)
                    ))
                },
                logical(1)
            ))
        },
        logical(1)
    ))
}
