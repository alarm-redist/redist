#library(igraph)

min_max <- function(e0, e1) {
    matrix(c(pmin(e0, e1), pmax(e0, e1)), ncol = 2)
}

get_farthest_two_vertices <- function(G) {
    ecc <- igraph::eccentricity(G)
    max_dist <- igraph::diameter(G, ecc)
    s <- names(which(ecc == max_dist))[1]
    sp <- igraph::shortest_paths(G, from = s, mode = "out")$vpath
    e <- names(which(sapply(sp, length) == max_dist))[1]
    c(s, e)
}

choose_next <- function(G, vertex_order_dict, order_v_list) {
    scores <- rep(0, length(order_v_list))
    for (i in seq_along(order_v_list)) {
        vs <- igraph::neighbors(G, order_v_list[i])
        count <- sum(vertex_order_dict[vs] >= 1)
        scores[i] <- count
    }
    max_score <- max(scores)
    n_scores <- rep(9999999, length(order_v_list))
    for (i in seq_along(order_v_list)) {
        if (scores[i] == max_score) {
            n_scores[i] <- 0
            vs <- igraph::neighbors(G, order_v_list[i])
            n_scores[i] <- sum(vertex_order_dict[vs] >= 1)
        }
    }
    min_n_score <- min(n_scores)
    order_v_list[which(n_scores == min_n_score)]
}

order_by_NDS <- function(G, vertex_order_dict, order_vertex_set, result_edge_list) {
    order_v_list <- as.character(order_vertex_set)
    v_count <- 1 + sum(vertex_order_dict >= 1)
    while (length(order_v_list) > 0) {
        v <- choose_next(G, vertex_order_dict, order_v_list)
        order_v_list <- setdiff(order_v_list, v)
        vertex_order_dict[v] <- v_count
        v_count <- v_count + 1
        vs <- igraph::neighbors(G, v)
        ws <- sort(vs, decreasing = FALSE, index.return = TRUE)$ix
        for (w in ws) {
            for (i in seq_len(igraph::gsize(G, v, w))) {
                if (vertex_order_dict[w] >= 1) {
                    result_edge_list <- c(result_edge_list, c(v, w))
                }
            }
        }
    }
}

split_graph <- function(G, s, t, left_vertex_set, right_vertex_set, cut_set,
                        vertex_order_dict, result_edge_list) {
    H <- G
    for (v in left_vertex_set) {
        if (!(v %in% c(s, t))) {
            H <- igraph::delete_vertices(H, v)
        }
    }
    for (v in right_vertex_set) {
        if (!(v %in% c(s, t))) {
            H <- igraph::delete_vertices(H, v)
        }
    }
    if (length(igraph::V(H)) <= 10 || igraph::edge(s, t) %in% igraph::E(H)) {
        order_vertex_set <- igraph::V(G) - left_vertex_set - right_vertex_set + cut_set
        return(order_by_NDS(G, vertex_order_dict, order_vertex_set, result_edge_list))
        # return() CK: check this
    }
    cut <- igraph::min_cut(H, s, t)$cut
    Hc <- igraph::delete_vertices(H, cut)
    cc_list <- igraph::clusters(Hc)$membership
    ccs <- NULL
    for (cc in seq_along(cc_list)) {
        if (s %in% cc_list[[cc]]) {
            ccs <- cc_list[[cc]]
            cc_list <- cc_list[-cc]
            break
        }
    }
    new_right_vertex_set <- union(right_vertex_set, ccs, cut)
    # split_graph(G, s, t, left_vertex_set, new_right_vertex_set, cut,
    #             vertex_order_dict, result_edge_list)
    # CK: check this
    new_left_vertex_set <- union(left_vertex_set, ccs, cut)
    split_graph(G, s, t, new_left_vertex_set, right_vertex_set, cut_set,
                vertex_order_dict, result_edge_list)
}

remove_deg12 <- function(G) {
    deg1_edges <- list()
    deg2_dict <- list()
    deg2_cycle <- list()

    found <- TRUE
    while (found) {
        found <- FALSE
        for (n in igraph::V(G)) {
            if (igraph::degree(G, n) == 1) {
                ns <- igraph::neighbors(G, n)
                G <- igraph::delete_vertices(G, n)
                deg1_edges <- c(deg1_edges, list(c(igraph::as_ids(ns[[1]]), n)))
                found <- TRUE
                break
            }
        }
    }

    if (length(igraph::E(G)) == 0) {
        e <- deg1_edges[length(deg1_edges)]
        G <- igraph::add_edges(G, igraph::as_ids(e)[1:2])
        deg1_edges <- deg1_edges[-length(deg1_edges)]
        return(list(G = G, deg1_edges = deg1_edges, deg2_dict = deg2_dict,
                    deg2_cycle = deg2_cycle, is_cycle = FALSE, is_tree = TRUE))
    }

    found <- TRUE
    cycle <- FALSE
    while (found) {
        found <- FALSE
        for (n in igraph::V(G)) {
            if (igraph::degree(G, n) == 2) {
                walk <- list(n)
                ns <- igraph::as_ids(igraph::neighbors(G, n))
                c <- ns[[1]]
                walk <- c(c, walk)
                while (igraph::degree(G, c) == 2) {
                    if (n == c) {
                        return(list(G = G, deg1_edges = deg1_edges,
                                    deg2_dict = deg2_dict,
                                    deg2_cycle = deg2_cycle,
                                    is_cycle = TRUE, is_tree = FALSE))
                    }
                    nc <- igraph::as_ids(igraph::neighbors(G, c))
                    if (nc[[1]] != walk[[2]]) {
                        c <- nc[[1]]
                    } else {
                        c <- nc[[2]]
                    }
                    walk <- c(c, walk)
                }
                c <- ns[[2]]
                walk <- c(walk, c)
                while (igraph::degree(G, c) == 2) {
                    nc <- igraph::neighbors(G, c)
                    if (nc[[1]] != walk[length(walk) - 1]) {
                        c <- nc[[1]]
                    } else {
                        c <- nc[[2]]
                    }
                    walk <- c(walk, c)
                }

                G <- igraph::delete_vertices(G, walk[-c(1, length(walk))])
                if (walk[[1]] == walk[length(walk)]) {  # cycle
                    key <- as.character(walk[[1]])
                    deg2_cycle[[key]] <- c(deg2_cycle[[key]], walk)
                } else {
                    has_e <- igraph::edge(G, walk[[1]], walk[length(walk)]) %in% igraph::E(G)
                    G <- igraph::add_edges(G, c(walk[[1]], walk[length(walk)]))
                    key <- paste(sort(c(walk[[1]], walk[length(walk)])), collapse = ",")
                    deg2_dict[[key]] <- c(deg2_dict[[key]], list(c(walk, has_e)))
                }
                found <- TRUE
                break
            }
        }
    }

    list(G = G, deg1_edges = deg1_edges, deg2_dict = deg2_dict, deg2_cycle,
         is_cycle = FALSE, is_tree = FALSE)
}

get_cycle <- function(edge_list) {
    result_edge_list <- list()
    v <- edge_list[1, 1]
    prev_v <- -1
    while (length(result_edge_list) < nrow(edge_list)) {
        vs <- edge_list[edge_list[, 1] == v | edge_list[, 2] == v, ]
        next_v <- setdiff(union(vs[1, ], vs[2, ]), c(prev_v, v))[1]
        result_edge_list <- c(result_edge_list, list(min_max(v, next_v)))
        prev_v <- v
        v <- next_v
    }
    do.call(rbind, result_edge_list)
}

check_connected_order <- function(edge_list) {
    touched_vertices <- character(0)
    for (i in seq_len(nrow(edge_list))) {
        if (i > 1 && !(edge_list[i, 1] %in% touched_vertices) && !(edge_list[i, 2] %in% touched_vertices)) {
            # cat(paste(i, edge_list[i, 1], edge_list[i, 2]), file = stderr())
            return(FALSE)
        }
        touched_vertices <- union(touched_vertices, edge_list[i, ])
    }
    TRUE
}

recover_deg12 <- function(deg1_edges, deg2_dict, deg2_cycle, result_edge_list) {
    found <- TRUE
    while (found) {
        found <- FALSE
        for (i in seq_len(nrow(result_edge_list))) {
            e <- min_max(result_edge_list[i, 1], result_edge_list[i, 2])
            if (toString(e) %in% names(deg2_dict)) {
                found <- TRUE
                x <- deg2_dict[[toString(e)]][[1]]
                if (toString(e) %in% result_edge_list) {
                    result_edge_list <- result_edge_list[!(toString(e) %in% result_edge_list), , drop = FALSE]
                } else {
                    result_edge_list <- result_edge_list[!(toString(c(e[2], e[1])) %in% result_edge_list), , drop = FALSE]
                }
                touched_vertices <- unique(as.character(unlist(result_edge_list[seq_len(i), ])))
                if (toString(x[[1]][1]) %in% touched_vertices) {
                    for (j in rev(seq_len(length(x[[1]]) - 1))) {
                        result_edge_list <- rbind(result_edge_list[seq_len(i), ], min_max(x[[1]][j], x[[1]][j + 1]), result_edge_list[(i + 1):nrow(result_edge_list), ])
                    }
                } else {
                    for (j in seq_len(length(x[[1]]) - 1)) {
                        result_edge_list <- rbind(result_edge_list[seq_len(i), ], min_max(x[[1]][j], x[[1]][j + 1]), result_edge_list[(i + 1):nrow(result_edge_list), ])
                    }
                }
                deg2_dict[[toString(e)]] <- NULL
                break
            }
        }
        if (!found && length(deg2_cycle) > 0) {
            for (v in names(deg2_cycle)) {
                for (i in seq_len(length(result_edge_list))) {
                    if (v %in% result_edge_list[i, ]) {
                        found <- TRUE
                        cy <- deg2_cycle[[v]][[1]]
                        for (j in seq_len(length(cy) - 1)) {
                            c0 <- cy[j]
                            c1 <- cy[j + 1]
                            result_edge_list <- rbind(result_edge_list[seq_len(i + 1), ], min_max(c0, c1), result_edge_list[(i + 2):nrow(result_edge_list), ])
                        }
                        deg2_cycle[[v]] <- NULL
                        break
                    }
                }
                if (found) {
                    break
                }
            }
        }
    }

    found <- TRUE
    while (length(deg1_edges) > 0) {
        found <- FALSE
        for (e in deg1_edges) {
            for (i in seq_len(length(result_edge_list))) {
                if (e[1] %in% result_edge_list[i, ]) {
                    result_edge_list <- rbind(result_edge_list[seq_len(i + 1), ], e, result_edge_list[(i + 2):nrow(result_edge_list), ])
                    found <- TRUE
                    deg1_edges <- deg1_edges[!identical(e, deg1_edges)]
                    break
                }
            }
            if (found) {
                break
            }
        }
        if (!found) {
            cli::cli_abort("not found")
        }
    }
    result_edge_list
}




get_order_by_cut <- function(edge_list) {
    G <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    res <- remove_deg12(G)
    if (res$is_cycle) {
        result_edge_list <- get_cycle(edge_list)
        recover_deg12(res$deg1_edges, res$deg2_dict, res$deg2_cycle, result_edge_list)
    } else if (res$is_tree) {
        result_edge_list <- edge_list
        recover_deg12(res$deg1_edges, res$deg2_dict, res$deg2_cycle, result_edge_list)
    } else {
        s_t <- get_farthest_two_vertices(G)
        result_edge_list <- character(0)
        vertex_order_dict <- rep(-1, igraph::vcount(G))
        vertex_order_dict[s_t[1]] <- 1
        left <- as.character(c(s_t[1]))
        right <- as.character(c(s_t[2]))
        cut_set <- right
        split_graph(G, s_t[1], s_t[2], left, right, cut_set, vertex_order_dict, result_edge_list)
        recover_deg12(res$deg1_edges, res$deg2_dict, res$deg2_cycle, result_edge_list)
    }
    result_edge_list
}

get_order_by_cut_with_check <- function(edge_list) {
    result_edge_list <- get_order_by_cut(edge_list)
    G1 <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    G2 <- igraph::graph_from_edgelist(result_edge_list, directed = FALSE)
    if (!igraph::isomorphic(G1, G2)) {
        cli::cli_abort("Graphs are not isomorphic!")
    }
    if (!check_connected_order(result_edge_list)) {
        cli::cli_abort("Not connected order!")
    }
    result_edge_list
}

# Example usage:
# edge_list <- matrix(c(1, 2, 2, 3, 3, 1, 3, 4), ncol = 2, byrow = TRUE)
# result_edge_list <- get_order_by_cut_with_check(edge_list)
# print(result_edge_list)
