# Helper function to ensure edges are ordered (min, max)
minmax <- function(e0, e1) {
  c(min(e0, e1), max(e0, e1))
}

# Find the two vertices that are farthest apart in the graph
get_farthest_two_vertices <- function(G) {
  # Calculate eccentricity (maximum distance from each vertex to any other)
  ecc <- igraph::eccentricity(G)

  # Find the diameter of the graph
  max_dist <- igraph::diameter(G)

  # Find a vertex with maximum eccentricity
  s <- names(ecc)[which(ecc == max_dist)][1]

  # Find the vertex farthest from s
  sp <- igraph::distances(G, v = s, to = igraph::V(G))
  e <- names(which(sp == max_dist))[1]

  list(s = as.numeric(s), e = as.numeric(e))
}

# Choose the next vertex in order based on scoring criteria
choose_next <- function(G, vertex_order_dict, order_v_list) {
  scores <- rep(0, length(order_v_list))

  for (i in 1:length(order_v_list)) {
    # Get neighbors of the vertex
    vs <- as.numeric(igraph::neighbors(G, order_v_list[i]))

    # Count neighbors that are already ordered
    for (v in vs) {
      if (vertex_order_dict[[as.character(v)]] >= 1) {
        scores[i] <- scores[i] + 1
      }
    }
  }

  max_score <- max(scores)
  n_scores <- rep(9999999, length(order_v_list))

  for (i in 1:length(order_v_list)) {
    if (scores[i] == max_score) {
      n_scores[i] <- 0
      vs <- as.numeric(igraph::neighbors(G, order_v_list[i]))

      for (v in vs) {
        if (vertex_order_dict[[as.character(v)]] >= 1) {
          n_scores[i] <- n_scores[i] + vertex_order_dict[[as.character(v)]]
        }
      }
    }
  }

  min_n_score <- min(n_scores)
  order_v_list[which(n_scores == min_n_score)[1]]
}

# Order vertices by NDS criterion
order_by_NDS <- function(G, vertex_order_dict, order_vertex_set, result_edge_list) {
  order_v_list <- as.numeric(order_vertex_set)

  v_count <- 1
  for (v in names(vertex_order_dict)) {
    if (vertex_order_dict[[v]] >= 1) {
      v_count <- v_count + 1
    }
  }

  while (length(order_v_list) > 0) {
    v <- choose_next(G, vertex_order_dict, order_v_list)
    order_v_list <- order_v_list[order_v_list != v]
    vertex_order_dict[[as.character(v)]] <- v_count
    v_count <- v_count + 1

    vs <- as.numeric(igraph::neighbors(G, v))

    # Sort neighbors by their order
    vs_order <- sapply(as.character(vs), function(w) vertex_order_dict[[w]])
    vs <- vs[order(vs_order)]

    for (w in vs) {
      # Handle multiple edges between vertices
      edge_count <- count_multiple(G, v, w)

      for (i in 1:edge_count) {
        if (vertex_order_dict[[as.character(w)]] >= 1) {
          result_edge_list[[length(result_edge_list) + 1]] <- c(v, w)
        }
      }
    }
  }

  list(vertex_order_dict = vertex_order_dict, result_edge_list = result_edge_list)
}

# Helper function to count multiple edges between vertices
count_multiple <- function(G, v1, v2) {
  edges <- igraph::get.edge.ids(G, c(rbind(v1, v2)), directed = FALSE)
  length(edges)
}

# Split graph recursively using minimum node cuts
split_graph <- function(G, s, t, left_vertex_set, right_vertex_set, cut_set, vertex_order_dict, result_edge_list) {
  H <- G

  # Contract vertices to s and t
  for (v in left_vertex_set) {
    if (v != s && v != t) {
      H <- igraph::contract(H, c(s, v), 'first')
    }
  }

  for (v in right_vertex_set) {
    if (v != s && v != t) {
      H <- igraph::contract(H, c(t, v), 'first')
    }
  }

  # Base case: small graph or direct edge
  if (igraph::vcount(H) <= 10 || igraph::are.adjacent(H, s, t)) {
    order_vertex_set <- as.character(igraph::V(G))
    order_vertex_set <- setdiff(order_vertex_set, as.character(left_vertex_set))
    order_vertex_set <- setdiff(order_vertex_set, as.character(right_vertex_set))
    order_vertex_set <- union(order_vertex_set, as.character(cut_set))
    order_vertex_set <- as.numeric(order_vertex_set)

    result <- order_by_NDS(G, vertex_order_dict, order_vertex_set, result_edge_list)
    vertex_order_dict <- result$vertex_order_dict
    result_edge_list <- result$result_edge_list

    return(list(vertex_order_dict = vertex_order_dict, result_edge_list = result_edge_list))
  }

  # Find minimum s-t node cut
  cut <- minimum_st_node_cut(H, s, t)

  # Remove cut vertices
  Hc <- H - cut

  # Get connected components
  cc_list <- igraph::decompose(Hc)
  ccs <- NULL

  for (i in 1:length(cc_list)) {
    cc_vertices <- as.numeric(igraph::V(cc_list[[i]]))
    if (s %in% cc_vertices) {
      ccs <- cc_vertices
      cc_list <- cc_list[-i]
      break
    }
  }

  new_right_vertex_set <- right_vertex_set
  new_right_vertex_set <- union(new_right_vertex_set, cut)

  for (cc in cc_list) {
    new_right_vertex_set <- union(new_right_vertex_set, as.numeric(igraph::V(cc)))
  }

  result <- split_graph(G, s, t, left_vertex_set, new_right_vertex_set, cut, vertex_order_dict, result_edge_list)
  vertex_order_dict <- result$vertex_order_dict
  result_edge_list <- result$result_edge_list

  new_left_vertex_set <- left_vertex_set
  new_left_vertex_set <- union(new_left_vertex_set, ccs)
  new_left_vertex_set <- union(new_left_vertex_set, cut)

  result <- split_graph(G, s, t, new_left_vertex_set, right_vertex_set, cut_set, vertex_order_dict, result_edge_list)

  result
}

# Find minimum s-t node cut (wrapper for igraph's st_min_cuts)
minimum_st_node_cut <- function(G, s, t) {
  # Convert to vertex ids if needed
  s_id <- ifelse(is.character(s), igraph::V(G)[s], s)
  t_id <- ifelse(is.character(t), igraph::V(G)[t], t)

  # Find all s-t cuts and get the minimum one
  cuts <- igraph::st_min_cuts(G, s_id, t_id)
  if (length(cuts$cut) == 0) {
    return(numeric(0))
  }

  min_cut <- cuts$cut[[which.min(sapply(cuts$cut, length))]]
  as.numeric(igraph::V(G)[min_cut])
}

# Remove degree 1 and 2 vertices and track them
remove_deg12 <- function(G) {
  deg1_edges <- list()
  deg2_dict <- list()
  deg2_cycle <- list()

  # Remove degree 1 vertices
  found <- TRUE
  while (found) {
    found <- FALSE
    for (n in as.numeric(igraph::V(G))) {
      if (igraph::degree(G, n) == 1) {
        ns <- as.numeric(igraph::neighbors(G, n))
        G <- G - n
        deg1_edges[[length(deg1_edges) + 1]] <- c(ns[1], n)
        found <- TRUE
        break
      }
    }
  }

  # Check if it's a tree (no edges left)
  if (igraph::ecount(G) == 0) {
    e <- deg1_edges[[length(deg1_edges)]]
    G <- G + igraph::edge(e)
    deg1_edges <- deg1_edges[-length(deg1_edges)]
    return(list(
      G = G, deg1_edges = deg1_edges, deg2_dict = deg2_dict,
      deg2_cycle = deg2_cycle, is_cycle = FALSE, is_tree = TRUE
    ))
  }

  # Handle degree 2 vertices
  found <- TRUE
  cycle <- FALSE

  while (found) {
    found <- FALSE
    for (n in as.numeric(igraph::V(G))) {
      if (igraph::degree(G, n) == 2) {
        walk <- c(n)
        ns <- as.numeric(igraph::neighbors(G, n))

        c <- ns[1]
        walk <- c(c, walk)

        while (igraph::degree(G, c) == 2) {
          if (n == c) {
            return(list(
              G = G, deg1_edges = deg1_edges, deg2_dict = deg2_dict,
              deg2_cycle = deg2_cycle, is_cycle = TRUE, is_tree = FALSE
            ))
          }

          nc <- as.numeric(igraph::neighbors(G, c))
          if (nc[1] != walk[2]) {
            c <- nc[1]
          } else {
            c <- nc[2]
          }
          walk <- c(c, walk)
        }

        c <- ns[2]
        walk <- c(walk, c)

        while (igraph::degree(G, c) == 2) {
          nc <- as.numeric(igraph::neighbors(G, c))
          if (nc[1] != walk[length(walk) - 1]) {
            c <- nc[1]
          } else {
            c <- nc[2]
          }
          walk <- c(walk, c)
        }

        # Remove internal vertices of walk
        G <- G - as.numeric(walk[2:(length(walk) - 1)])

        if (walk[1] == walk[length(walk)]) {
          # It's a cycle
          key <- as.character(walk[1])
          if (!(key %in% names(deg2_cycle))) {
            deg2_cycle[[key]] <- list()
          }
          deg2_cycle[[key]][[length(deg2_cycle[[key]]) + 1]] <- walk
        } else {
          # Add edge between endpoints
          has_e <- igraph::are.adjacent(G, walk[1], walk[length(walk)])
          G <- G + igraph::edge(c(walk[1], walk[length(walk)]))

          key <- paste(minmax(walk[1], walk[length(walk)]), collapse = '-')
          if (!(key %in% names(deg2_dict))) {
            deg2_dict[[key]] <- list()
          }
          deg2_dict[[key]][[length(deg2_dict[[key]]) + 1]] <- list(walk = walk, has_e = has_e)
        }

        found <- TRUE
        break
      }
    }
  }

  list(
    G = G, deg1_edges = deg1_edges, deg2_dict = deg2_dict,
    deg2_cycle = deg2_cycle, is_cycle = FALSE, is_tree = FALSE
  )
}

# Get cycle edges in order
get_cycle <- function(edge_list) {
  result_edge_list <- list()
  v <- edge_list[[1]][1]
  prev_v <- -1

  while (length(result_edge_list) < length(edge_list)) {
    vs <- Filter(function(x) x[1] == v || x[2] == v, edge_list)
    next_v <- setdiff(union(vs[[1]], vs[[2]]), c(prev_v, v))[1]
    result_edge_list[[length(result_edge_list) + 1]] <- minmax(v, next_v)
    prev_v <- v
    v <- next_v
  }

  result_edge_list
}

# Check if edge list forms a connected order
check_connected_order <- function(edge_list) {
  touched_vertices <- c()

  for (i in 1:length(edge_list)) {
    if (i > 1 && !(edge_list[[i]][1] %in% touched_vertices) &&
      !(edge_list[[i]][2] %in% touched_vertices)) {
      return(FALSE)
    }
    touched_vertices <- union(touched_vertices, edge_list[[i]])
  }

  TRUE
}

# Recover degree 1 and 2 vertices
recover_deg12 <- function(deg1_edges, deg2_dict, deg2_cycle, result_edge_list) {
  # Recover degree 2 vertices
  found <- TRUE

  while (found) {
    found <- FALSE

    for (i in 1:length(result_edge_list)) {
      e <- minmax(result_edge_list[[i]][1], result_edge_list[[i]][2])
      key <- paste(e, collapse = '-')

      if (key %in% names(deg2_dict)) {
        found <- TRUE

        for (x in deg2_dict[[key]]) {
          # Remove the edge
          result_edge_list <- result_edge_list[-i]

          # Get vertices touched so far
          touched_vertices <- unique(unlist(result_edge_list[1:i - 1]))

          walk <- x$walk

          if (walk[1] %in% touched_vertices) {
            # Insert walk edges in reverse order
            for (j in (length(walk) - 1):1) {
              result_edge_list <- c(
                result_edge_list[1:i - 1],
                list(c(walk[j], walk[j + 1])),
                result_edge_list[i:length(result_edge_list)]
              )
            }
          } else {
            # Insert walk edges in order
            for (j in 1:(length(walk) - 1)) {
              result_edge_list <- c(
                result_edge_list[1:i - 1],
                list(c(walk[j], walk[j + 1])),
                result_edge_list[i:length(result_edge_list)]
              )
            }
          }
        }

        # Remove processed entry
        deg2_dict[[key]] <- NULL
        break
      }
    }

    # Handle cycles if no degree 2 paths were found
    if (!found && length(deg2_cycle) > 0) {
      for (v in names(deg2_cycle)) {
        v_num <- as.numeric(v)

        for (i in 1:length(result_edge_list)) {
          if (v_num %in% result_edge_list[[i]]) {
            found <- TRUE

            for (cy in deg2_cycle[[v]]) {
              for (j in 1:(length(cy) - 1)) {
                c0 <- cy[j]
                c1 <- cy[j + 1]
                result_edge_list <- c(
                  result_edge_list[1:i],
                  list(minmax(c0, c1)),
                  result_edge_list[(i + 1):length(result_edge_list)]
                )
              }
            }

            # Remove processed entry
            deg2_cycle[[v]] <- NULL
            break
          }
        }

        if (found) break
      }
    }
  }

  # Recover degree 1 vertices
  while (length(deg1_edges) > 0) {
    found <- FALSE

    for (e_idx in 1:length(deg1_edges)) {
      e <- deg1_edges[[e_idx]]

      for (i in 1:length(result_edge_list)) {
        if (e[1] %in% result_edge_list[[i]]) {
          result_edge_list <- c(
            result_edge_list[1:i],
            list(e),
            result_edge_list[(i + 1):length(result_edge_list)]
          )
          found <- TRUE
          deg1_edges <- deg1_edges[-e_idx]
          break
        }
      }

      if (found) break
    }

    if (!found) {
      stop('Error: unable to find connection point for degree 1 edge')
    }
  }

  result_edge_list
}

# Main function to order edges by cut
get_order_by_cut <- function(edge_list) {
  # Create igraph from edge list
  edges_df <- do.call(rbind, edge_list)
  G <- igraph::graph_from_edgelist(edges_df, directed = FALSE)

  # Handle multiple edges
  G <- igraph::simplify(G, remove.multiple = FALSE, remove.loops = FALSE)

  # Remove degree 1 and 2 vertices
  result <- remove_deg12(G)
  G <- result$G
  deg1_edges <- result$deg1_edges
  deg2_dict <- result$deg2_dict
  deg2_cycle <- result$deg2_cycle
  is_cycle <- result$is_cycle
  is_tree <- result$is_tree

  if (is_cycle) {
    # Handle cycle case
    edges_vec <- igraph::as_edgelist(G)
    edges_list <- lapply(1:nrow(edges_vec), function(i) edges_vec[i, ])
    result_edge_list <- get_cycle(edges_list)
    result_edge_list <- recover_deg12(deg1_edges, deg2_dict, deg2_cycle, result_edge_list)
  } else if (is_tree) {
    # Handle tree case
    edges_vec <- igraph::as_edgelist(G)
    result_edge_list <- lapply(1:nrow(edges_vec), function(i) edges_vec[i, ])
    result_edge_list <- recover_deg12(deg1_edges, deg2_dict, deg2_cycle, result_edge_list)
  } else {
    # Handle general graph case
    G <- result$G
    vertices <- get_farthest_two_vertices(G)
    s <- vertices$s
    t <- vertices$e

    result_edge_list <- list()
    vertex_order_dict <- list()

    for (v in as.numeric(igraph::V(G))) {
      vertex_order_dict[[as.character(v)]] <- -1
    }

    vertex_order_dict[[as.character(s)]] <- 1
    left <- c(s)
    right <- c(t)
    cut_set <- right

    result <- split_graph(G, s, t, left, right, cut_set, vertex_order_dict, result_edge_list)
    vertex_order_dict <- result$vertex_order_dict
    result_edge_list <- result$result_edge_list

    result_edge_list <- recover_deg12(deg1_edges, deg2_dict, deg2_cycle, result_edge_list)
  }

  result_edge_list
}

# Check the output and validate
get_order_by_cut_with_check <- function(edge_list) {
  result_edge_list <- get_order_by_cut(edge_list)

  # Create graphs for comparison
  edges_df1 <- do.call(rbind, edge_list)
  edges_df2 <- do.call(rbind, result_edge_list)

  G1 <- igraph::graph_from_edgelist(edges_df1, directed = FALSE)
  G2 <- igraph::graph_from_edgelist(edges_df2, directed = FALSE)

  # Simplify for isomorphism check (keep multiple edges)
  G1 <- igraph::simplify(G1, remove.multiple = FALSE, remove.loops = FALSE)
  G2 <- igraph::simplify(G2, remove.multiple = FALSE, remove.loops = FALSE)

  if (!igraph::isomorphic(G1, G2)) {
    stop('Error: output graph is not isomorphic to input graph')
  }

  if (!check_connected_order(result_edge_list)) {
    stop('Error: output edge order is not connected')
  }

  result_edge_list
}

ndscut <- function(edges) {
  # Convert input to a list of edge pairs if it's not already
  edge_list <- list()

  if (is.matrix(edges) || is.data.frame(edges)) {
    # If input is a matrix or data frame, convert to a list of pairs
    for (i in 1:nrow(edges)) {
      edge_list[[i]] <- c(edges[i, 1], edges[i, 2])
    }
  } else if (is.list(edges)) {
    # If already a list, use it directly
    edge_list <- edges
  } else {
    stop('Input must be a matrix, data frame, or list of edge pairs')
  }

  if (length(edge_list) == 0) {
    message('The input graph is empty.')
    return(NULL)
  }

  new_edge_list <- get_order_by_cut(edge_list)

  # Create output
  result <- list()
  for (i in 1:length(new_edge_list)) {
    e <- new_edge_list[[i]]
    c <- minmax(e[1], e[2])
    result[[i]] <- c(c[1], c[2])
  }

  # Final connectivity check
  if (!check_connected_order(new_edge_list)) {
    warning('Output edge order is not connected')
  }

  result
}
