# Helper function to ensure edges are ordered (min, max)
edge_minmax <- function(e0, e1) {
    c(min(e0, e1), max(e0, e1))
}

# Find the two vertices that are farthest apart in the graph
get_farthest_two_vertices <- function(G) {
    # eccentricity: max distance from a vertex to any other
    ecc <- igraph::eccentricity(G)

    # Find the diameter of the graph
    max_dist <- igraph::diameter(G)

    # find a vertex with max eccentricity (diameter)
    s_idx <- which(ecc == max_dist)[1]
    s <- as.numeric(igraph::V(G)$name[s_idx])

    # find the vertex farthest from s
    sp <- igraph::distances(G, v = s_idx, to = igraph::V(G))
    e_idx <- which(sp[1, ] == max_dist)[1]
    e <- as.numeric(igraph::V(G)$name[e_idx])

    list(s = s, e = e)
}



# Choose the next vertex in order based on scoring criteria
choose_next <- function(G, vertex_order, order_v_list) {

    neigh_ids <- function(v) {
        as.numeric(unique(igraph::neighbors(G, as.character(v))$name))
    }

    vertex_order_vec <- unlist(vertex_order)

    scores <- vapply(
        order_v_list,
        function(v) {
            vs <- neigh_ids(v)
            sum(vertex_order_vec[as.character(vs)] >= 1)
        },
        numeric(1)
    )

    max_score <- max(scores)
    candidates <- order_v_list[scores == max_score]

    n_scores <- vapply(
        candidates,
        function(v) {
            vs <- neigh_ids(v)
            ords <- vertex_order_vec[as.character(vs)]
            sum(ords[ords >= 1])
        },
        numeric(1)
    )

    candidates[which.min(n_scores)]
}

# Order vertices by NDS criterion
order_by_NDS <- function(G, vertex_order, order_vertex_set, result_edge_list) {
    # Convert once, keep names identical
    vertex_order_vec <- unlist(vertex_order)

    order_v_list <- as.numeric(order_vertex_set)

    # Count already-ordered vertices
    v_count <- sum(vertex_order_vec >= 1) + 1

    while (length(order_v_list) > 0) {

        v <- choose_next(G, vertex_order_vec, order_v_list)

        order_v_list <- order_v_list[order_v_list != v]
        vertex_order_vec[as.character(v)] <- v_count
        v_count <- v_count + 1

        # order neighbors by vertex order, keeping only processed ones
        vs <- as.numeric(unique(igraph::neighbors(G, as.character(v))$name))

        # Order neighbors by their vertex order
        vs <- vs[order(vertex_order_vec[as.character(vs)])]

        # Keep only neighbors already ordered
        vs_done <- vs[vertex_order_vec[as.character(vs)] >= 1]

        if (length(vs_done) > 0) {
            edge_counts <- vapply(
                vs_done,
                function(w) count_multiple(G, as.character(v), as.character(w)),
                integer(1)
            )

            new_edges <- Map(
                function(w, k) replicate(k, c(v, w), simplify = FALSE),
                vs_done,
                edge_counts
            )

            result_edge_list <- c(
                result_edge_list,
                unlist(new_edges, recursive = FALSE)
            )
        }
    }

    list(
        vertex_order = vertex_order_vec,
        result_edge_list  = result_edge_list
    )
}

# count multiple edges between two vertices
count_multiple <- function(G, v1, v2) {
    # Get all edges and their endpoints
    el <- igraph::ends(G, igraph::E(G))
    # Count edges between v1 and v2 (in either direction for undirected)
    sum((el[, 1] == v1 & el[, 2] == v2) | (el[, 1] == v2 & el[, 2] == v1))
}

# contract node v into u, moving v's edges to u
contracted_nodes <- function(G, u, v) {

    u_chr <- as.character(u)
    v_chr <- as.character(v)

    # Neighbors of v (character names)
    v_neighbors <- igraph::neighbors(G, v_chr)$name

    # Exclude self-loops (u)
    v_neighbors <- v_neighbors[v_neighbors != u_chr]

    if (length(v_neighbors) > 0) {

        # count multiplicities of edges connected to v
        edge_counts <- vapply(
            v_neighbors,
            function(w) count_multiple(G, v_chr, w),
            integer(1)
        )

        # build new edge vector to attach to u
        new_edges <- unlist(
            Map(
                function(w, k) rep(c(u_chr, w), k),
                v_neighbors,
                edge_counts
            )
        )

        if (length(new_edges) > 0) {
            G <- igraph::add_edges(G, new_edges)
        }
    }

    # Remove v
    igraph::delete_vertices(G, v_chr)
}

# recursively split graph with minimum node cuts
split_graph <- function(G, s, t,
                        left_vertex_set,
                        right_vertex_set,
                        cut_set,
                        vertex_order,
                        result_edge_list) {

    # base case: order vertices if graph is small
    base_case <- function(G,
                          left_vertex_set,
                          right_vertex_set,
                          cut_set,
                          vertex_order,
                          result_edge_list) {

        order_vertex_set <- igraph::V(G)$name
        order_vertex_set <- setdiff(order_vertex_set, as.character(left_vertex_set))
        order_vertex_set <- setdiff(order_vertex_set, as.character(right_vertex_set))
        order_vertex_set <- union(order_vertex_set, as.character(cut_set))
        order_vertex_set <- as.numeric(order_vertex_set)

        res <- order_by_NDS(G, vertex_order, order_vertex_set, result_edge_list)

        list(
            vertex_order     = res$vertex_order,
            result_edge_list = res$result_edge_list
        )
    }

    # operate on a contracted copy of the graph
    H <- G

    # contract all left-side vertices into s
    left_to_contract <- setdiff(left_vertex_set, c(s, t))
    left_to_contract <- left_to_contract[
        as.character(left_to_contract) %in% igraph::V(H)$name
    ]

    for (v in left_to_contract) {
        H <- contracted_nodes(H, s, v)
    }

    # contract all right-side vertices into t
    right_to_contract <- setdiff(right_vertex_set, c(s, t))
    right_to_contract <- right_to_contract[
        as.character(right_to_contract) %in% igraph::V(H)$name
    ]

    for (v in right_to_contract) {
        H <- contracted_nodes(H, t, v)
    }

    # if graph is small or s-t are connected, run base case
    if (igraph::vcount(H) <= 10 ||
        igraph::are_adjacent(H, as.character(s), as.character(t))) {

        return(base_case(
            G,
            left_vertex_set,
            right_vertex_set,
            cut_set,
            vertex_order,
            result_edge_list
        ))
    }

    # find the minimum s-t node cut
    cut <- minimum_st_node_cut(H, s, t)
    cut <- cut[as.character(cut) %in% igraph::V(H)$name]

    if (length(cut) == 0) {
        return(base_case(
            G,
            left_vertex_set,
            right_vertex_set,
            cut_set,
            vertex_order,
            result_edge_list
        ))
    }

    # remove cut vertices and find connected components
    Hc <- igraph::delete_vertices(H, as.character(cut))
    cc_list <- igraph::decompose(Hc)

    cc_vertices <- lapply(cc_list, function(cc) as.numeric(igraph::V(cc)$name))

    # identify the component containing the start vertex s
    s_idx <- which(vapply(
        cc_vertices,
        function(vs) s %in% vs,
        logical(1)
    ))

    ccs <- cc_vertices[[s_idx]]
    other_ccs <- cc_vertices[-s_idx]

    # first recursive call, expanding the right side
    new_right_vertex_set <- union(right_vertex_set, cut)
    if (length(other_ccs) > 0) {
        new_right_vertex_set <- union(
            new_right_vertex_set,
            unlist(other_ccs)
        )
    }

    res <- split_graph(
        G, s, t,
        left_vertex_set,
        new_right_vertex_set,
        cut,
        vertex_order,
        result_edge_list
    )

    # second recursive call, expanding the left side
    new_left_vertex_set <- union(left_vertex_set, ccs)
    new_left_vertex_set <- union(new_left_vertex_set, cut)

    split_graph(
        G, s, t,
        new_left_vertex_set,
        right_vertex_set,
        cut_set,
        res$vertex_order,
        res$result_edge_list
    )
}


# Find minimum s-t node cut (wrapper for igraph's st_min_cuts)
minimum_st_node_cut <- function(G, s, t) {

    # igraph's st_min_cuts requires a directed graph
    if (igraph::is_directed(G)) {
        G_dir <- G
    } else {
        G_dir <- igraph::as_directed(G, mode = "mutual")
    }

    vnames <- igraph::V(G_dir)$name
    s_chr <- as.character(s)
    t_chr <- as.character(t)

    # Match vertex names to ids
    s_id <- match(s_chr, vnames)
    t_id <- match(t_chr, vnames)

    if (is.na(s_id) || is.na(t_id)) {
        return(numeric(0))
    }

    cuts <- igraph::st_min_cuts(G_dir, s_id, t_id)

    if (length(cuts$cuts) == 0) {
        return(numeric(0))
    }

    # find the smallest cut by number of vertices
    cut_sizes <- lengths(cuts$cuts)
    min_cut <- cuts$cuts[[which.min(cut_sizes)]]

    as.numeric(vnames[as.integer(min_cut)])
}


# simplify graph by removing degree 1 & 2 vertices, tracking for later recovery
remove_deg12 <- function(G) {

    deg1_edges <- list()
    deg2_dict  <- list()
    deg2_cycle <- list()

    # ---- Helper: current degrees as named integer vector
    get_deg <- function(G) {
        igraph::degree(G, v = igraph::V(G), mode = "all")
    }

    # iteratively remove all degree-1 vertices
    repeat {
        deg <- get_deg(G)
        deg1 <- names(deg)[deg == 1]

        if (length(deg1) == 0) break

        n <- deg1[1]
        ns <- igraph::neighbors(G, n)$name

        G <- igraph::delete_vertices(G, n)
        deg1_edges[[length(deg1_edges) + 1]] <-
            c(as.numeric(ns[1]), as.numeric(n))
    }

    # if graph is a tree, all edges will be removed
    if (igraph::ecount(G) == 0) {

        e <- deg1_edges[[length(deg1_edges)]]
        deg1_edges <- deg1_edges[-length(deg1_edges)]

        existing <- igraph::V(G)$name
        to_add <- setdiff(as.character(e), existing)

        if (length(to_add) > 0) {
            G <- igraph::add_vertices(G, length(to_add), name = to_add)
        }

        G <- igraph::add_edges(G, as.character(e))

        return(list(
            G = G,
            deg1_edges = deg1_edges,
            deg2_dict  = deg2_dict,
            deg2_cycle = deg2_cycle,
            is_cycle   = FALSE,
            is_tree    = TRUE
        ))
    }

    # iteratively remove chains of degree-2 vertices
    repeat {
        deg <- get_deg(G)
        deg2 <- names(deg)[deg == 2]

        if (length(deg2) == 0) break

        n <- deg2[1]
        walk <- n
        ns <- igraph::neighbors(G, n)$name

        # trace the chain in one direction
        c <- ns[1]
        walk <- c(c, walk)

        while (igraph::degree(G, c) == 2) {

            if (c == n) {
                return(list(
                    G = G,
                    deg1_edges = deg1_edges,
                    deg2_dict  = deg2_dict,
                    deg2_cycle = deg2_cycle,
                    is_cycle   = TRUE,
                    is_tree    = FALSE
                ))
            }

            nc <- igraph::neighbors(G, c)$name
            c <- if (nc[1] != walk[2]) nc[1] else nc[2]
            walk <- c(c, walk)
        }

        # trace the chain in the other direction
        c <- ns[2]
        walk <- c(walk, c)

        while (igraph::degree(G, c) == 2) {
            nc <- igraph::neighbors(G, c)$name
            c <- if (nc[1] != walk[length(walk) - 1]) nc[1] else nc[2]
            walk <- c(walk, c)
        }

        # remove the internal vertices of the chain
        internal <- walk[2:(length(walk) - 1)]
        G <- igraph::delete_vertices(G, internal)

        walk_num <- as.numeric(walk)

        # store chain info to recover it later
        if (walk[1] == walk[length(walk)]) {

            key <- walk[1]
            deg2_cycle[[key]] <- c(deg2_cycle[[key]], list(walk_num))

        } else {

            has_e <- igraph::are_adjacent(G, walk[1], walk[length(walk)])
            G <- igraph::add_edges(G, c(walk[1], walk[length(walk)]))

            key <- paste(
                edge_minmax(as.numeric(walk[1]), as.numeric(walk[length(walk)])),
                collapse = "-"
            )

            deg2_dict[[key]] <- c(
                deg2_dict[[key]],
                list(list(walk = walk_num, has_e = has_e))
            )
        }
    }

    list(
        G = G,
        deg1_edges = deg1_edges,
        deg2_dict  = deg2_dict,
        deg2_cycle = deg2_cycle,
        is_cycle   = FALSE,
        is_tree    = FALSE
    )
}


# Get cycle edges in order
get_cycle <- function(edge_list) {
    # Convert to a numeric 2‑column matrix for fast vector ops
    edges <- do.call(rbind, lapply(edge_list, as.numeric))

    # Start at the first vertex of the first edge
    v <- edges[1, 1]
    prev_v <- NA

    # Preallocate result
    n <- nrow(edges)
    result <- vector("list", n)

    # walk the cycle and record edges in order
    for (i in seq_len(n)) {
        # Find all edges touching v
        idx <- which(edges[,1] == v | edges[,2] == v)

        # Extract the two vertices from those edges
        candidates <- unique(c(edges[idx, 1], edges[idx, 2]))

        # Next vertex is whichever is not prev_v or v
        next_v <- setdiff(candidates, c(prev_v, v))[1]

        # Store edge in min–max order
        result[[i]] <- sort(c(v, next_v))

        prev_v <- v
        v <- next_v
    }

    result
}


# re-insert the degree 1 & 2 vertices that were removed
recover_deg12 <- function(deg1_edges, deg2_dict, deg2_cycle, result_edge_list) {

    # Helper: normalize edge
    norm_edge <- function(e) sort(as.numeric(e))

    # Helper: compute key
    edge_key <- function(e) paste(norm_edge(e), collapse = "-")

    # Precompute keys for result edges
    get_keys <- function(edges) vapply(edges, edge_key, character(1))

    # Insert edges at a given position
    insert_edges <- function(lst, pos, new_edges) {
        if (pos <= 1) {
            c(new_edges, lst)
        } else if (pos >= length(lst)) {
            c(lst, new_edges)
        } else {
            c(lst[1:(pos - 1)], new_edges, lst[pos:length(lst)])
        }
    }

    # first, re-insert degree-2 vertex chains
    repeat {
        keys <- get_keys(result_edge_list)
        found <- FALSE

        # Look for any edge whose key appears in deg2_dict
        for (i in seq_along(result_edge_list)) {
            e <- result_edge_list[[i]]
            key <- edge_key(e)

            if (key %in% names(deg2_dict)) {
                found <- TRUE
                walks <- deg2_dict[[key]]
                deg2_dict[[key]] <- NULL

                for (walk_info in walks) {
                    # Remove the edge (search for it each time)
                    keys <- get_keys(result_edge_list)
                    j <- match(key, keys, nomatch = 0)
                    if (j == 0) next

                    result_edge_list <- result_edge_list[-j]

                    # Determine touched vertices using the ORIGINAL position i (before removal)
                    # but slicing result_edge_list AFTER the removal
                    touched <- if (i > 1)
                        unique(unlist(result_edge_list[1:(i - 1)]))
                    else
                        numeric(0)

                    walk <- walk_info$walk

                    # Build edges to insert
                    if (walk[1] %in% touched) {
                        seqs <- (length(walk) - 1):1
                    } else {
                        seqs <- 1:(length(walk) - 1)
                    }

                    edges_to_insert <- lapply(seqs, function(k) norm_edge(c(walk[k], walk[k + 1])))

                    # Reverse because Python inserts one at a time (which reverses), but we insert as a block
                    edges_to_insert <- rev(edges_to_insert)

                    # Insert at position i (the original position)
                    result_edge_list <- insert_edges(result_edge_list, i, edges_to_insert)
                }
                break # for i
            }
        }

        if (!found && length(deg2_cycle) > 0) {
            # handle special case for degree-2 cycles
            keys <- get_keys(result_edge_list)

            for (v in names(deg2_cycle)) {
                v_num <- as.numeric(v)

                # Find first edge containing v
                idx <- vapply(result_edge_list, function(e) v_num %in% e, logical(1))
                i <- match(TRUE, idx, nomatch = 0)

                if (i > 0) {
                    found <- TRUE

                    # Build cycle edges
                    edges_to_insert <- unlist(
                        lapply(deg2_cycle[[v]], function(cy) {
                            lapply(seq_len(length(cy) - 1),
                                   function(j) norm_edge(c(cy[j], cy[j + 1])))
                        }),
                        recursive = FALSE
                    )

                    result_edge_list <- insert_edges(result_edge_list, i + 1, edges_to_insert)
                    deg2_cycle[[v]] <- NULL
                    break
                }
            }
        }

        if (!found) break
    }

    # next, re-insert degree-1 vertices
    while (length(deg1_edges) > 0) {
        keys <- get_keys(result_edge_list)
        placed <- FALSE

        for (e_idx in seq_along(deg1_edges)) {
            e <- deg1_edges[[e_idx]]

            # Find first edge touching e[1]
            idx <- vapply(result_edge_list, function(x) e[1] %in% x, logical(1))
            i <- match(TRUE, idx, nomatch = 0)

            if (i > 0) {
                result_edge_list <- insert_edges(result_edge_list, i + 1, list(e))
                deg1_edges <- deg1_edges[-e_idx]
                placed <- TRUE
                break
            }
        }

        if (!placed) {
            cli::cli_abort("Unable to find connection point for degree 1 edge")
        }
    }

    result_edge_list
}


# Main function to order edges by cut
get_order_by_cut <- function(edge_list) {

    # Convert to igraph
    edges_chr <- do.call(rbind, lapply(edge_list, as.character))
    G <- igraph::graph_from_edgelist(edges_chr, directed = FALSE)
    G <- igraph::simplify(G, remove.multiple = FALSE, remove.loops = FALSE)

    # Remove degree‑1 and degree‑2 vertices
    r <- remove_deg12(G)
    G <- r$G

    # Helper: convert igraph edges to numeric list
    as_numeric_edges <- function(g) {
        ed <- igraph::as_edgelist(g)
        lapply(seq_len(nrow(ed)), function(i) as.numeric(ed[i, ]))
    }

    # Case 1: the simplified graph is a cycle
    if (r$is_cycle) {
        base_edges <- as_numeric_edges(G)
        result <- get_cycle(base_edges)
        result <- recover_deg12(r$deg1_edges, r$deg2_dict, r$deg2_cycle, result)
        return(lapply(result, as.numeric))
    }

    # Case 2: the simplified graph is a tree
    if (r$is_tree) {
        base_edges <- as_numeric_edges(G)
        result <- recover_deg12(r$deg1_edges, r$deg2_dict, r$deg2_cycle, base_edges)
        return(lapply(result, as.numeric))
    }

    # Case 3: general graph, use recursive splitting
    # Find farthest pair
    ft <- get_farthest_two_vertices(G)
    s <- ft$s
    t <- ft$e

    # Initialize vertex order
    verts <- as.numeric(igraph::V(G)$name)
    vertex_order <- setNames(rep(-1, length(verts)), as.character(verts))
    vertex_order[[as.character(s)]] <- 1

    # Initialize left/right/cut sets
    left <- s
    right <- t
    cut_set <- right

    # Split graph
    split_res <- split_graph(
        G, s, t,
        left, right, cut_set,
        vertex_order,
        result_edge_list = list()
    )

    result <- recover_deg12(
        r$deg1_edges,
        r$deg2_dict,
        r$deg2_cycle,
        split_res$result_edge_list
    )

    lapply(result, as.numeric)
}


# Check the output and validate
get_order_by_cut_with_check <- function(edge_list) {

    result_edge_list <- get_order_by_cut(edge_list)

    # Convert list-of-edges → 2‑column numeric matrix
    edges1 <- do.call(rbind, lapply(edge_list, as.numeric))
    edges2 <- do.call(rbind, lapply(result_edge_list, as.numeric))

    # Build graphs
    G1 <- igraph::graph_from_edgelist(edges1, directed = FALSE)
    G2 <- igraph::graph_from_edgelist(edges2, directed = FALSE)

    # Keep multiple edges, keep loops
    G1 <- igraph::simplify(G1, remove.multiple = FALSE, remove.loops = FALSE)
    G2 <- igraph::simplify(G2, remove.multiple = FALSE, remove.loops = FALSE)

    # Structural check
    if (!igraph::isomorphic(G1, G2)) {
        cli::cli_abort("Output graph is not isomorphic to input graph.")
    }

    # Connectivity check
    if (!igraph::is_connected(G2)) {
        cli::cli_abort("Output edge order is not connected.")
    }

    result_edge_list
}



# convert adjacency list to edge matrix
adj_to_edges <- function(adj) {
    lens <- lengths(adj)
    edges <- cbind(rep(seq_along(adj), times = lens), unlist(adj) + 1L)
    edges[edges[, 1] < edges[, 2], , drop = FALSE]
}

ndscut <- function(edges) {
    # Convert input to a list of edge pairs if it's not already
    edge_list <- list()

    if (is.matrix(edges) || is.data.frame(edges)) {
        # Check for empty input first
        if (nrow(edges) == 0) {
            return(NULL)
        }
        # If input is a matrix or data frame, convert to a list of pairs
        for (i in seq_len(nrow(edges))) {
            edge_list[[i]] <- c(edges[i, 1], edges[i, 2])
        }
    } else if (is.list(edges)) {
        # Convert adjacency list to edge matrix, then to edge list
        edge_mat <- adj_to_edges(edges)
        if (is.null(edge_mat) || nrow(edge_mat) == 0) {
            return(NULL)
        }
        for (i in seq_len(nrow(edge_mat))) {
            edge_list[[i]] <- c(edge_mat[i, 1], edge_mat[i, 2])
        }
    }

    if (length(edge_list) == 0) {
        return(NULL)
    }

    new_edge_list <- get_order_by_cut(edge_list)

    # Create output
    result <- list()
    for (i in seq_along(new_edge_list)) {
        e <- new_edge_list[[i]]
        # Skip NULL or invalid edges
        if (!is.null(e) && length(e) >= 2) {
            c <- edge_minmax(e[1], e[2])
            result[[length(result) + 1]] <- c(c[1], c[2])
        }
    }

    # final check to ensure output is a single connected component
    edges_df <- do.call(rbind, result)
    G <- igraph::graph_from_edgelist(edges_df, directed = FALSE)
    if (!igraph::is_connected(G)) {
        cli::cli_warn('Output edge order is not connected')
    }

    result
}
