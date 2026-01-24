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
  s_idx <- which(ecc == max_dist)[1]
  s <- as.numeric(igraph::V(G)$name[s_idx])

  # Find the vertex farthest from s
  sp <- igraph::distances(G, v = s_idx, to = igraph::V(G))
  e_idx <- which(sp[1, ] == max_dist)[1]
  e <- as.numeric(igraph::V(G)$name[e_idx])

  list(s = s, e = e)
}

# Choose the next vertex in order based on scoring criteria
choose_next <- function(G, vertex_order_dict, order_v_list) {
  scores <- rep(0, length(order_v_list))

  for (i in 1:length(order_v_list)) {
    # Get neighbors of the vertex (unique only, as multigraphs return duplicates)
    vs <- unique(igraph::neighbors(G, as.character(order_v_list[i]))$name)
    vs <- as.numeric(vs)

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
      vs <- unique(igraph::neighbors(G, as.character(order_v_list[i]))$name)
      vs <- as.numeric(vs)

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

    # Get neighbors using vertex name (unique only, as multigraphs return duplicates)
    vs <- unique(igraph::neighbors(G, as.character(v))$name)
    vs <- as.numeric(vs)

    # Sort neighbors by their order
    vs_order <- sapply(as.character(vs), function(w) vertex_order_dict[[w]])
    vs <- vs[order(vs_order)]

    for (w in vs) {
      # Handle multiple edges between vertices (use vertex names)
      edge_count <- count_multiple(G, as.character(v), as.character(w))

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
  # Get all edges between v1 and v2
  length(igraph::E(G)[v1 %--% v2])
}

# Contract node v into node u (NetworkX-style contraction)
# Keeps u, removes v, adds edges from u to all neighbors of v
contracted_nodes <- function(G, u, v) {
  # Get neighbors of v (as character names)
  v_neighbors <- igraph::neighbors(G, as.character(v))$name

  # Add edges from u to each neighbor of v (except u itself)
  for (neighbor in v_neighbors) {
    if (neighbor != as.character(u)) {
      # Count existing edges between v and neighbor
      n_edges <- count_multiple(G, as.character(v), neighbor)
      # Add same number of edges between u and neighbor
      for (i in seq_len(n_edges)) {
        G <- igraph::add_edges(G, c(as.character(u), neighbor))
      }
    }
  }

  # Remove vertex v
  G <- igraph::delete_vertices(G, as.character(v))

  G
}

# Split graph recursively using minimum node cuts
split_graph <- function(G, s, t, left_vertex_set, right_vertex_set, cut_set, vertex_order_dict, result_edge_list) {
  H <- G

  # Contract vertices to s and t (NetworkX-style contraction)
  # Contract left vertices into s
  for (v in left_vertex_set) {
    if (v != s && v != t && as.character(v) %in% igraph::V(H)$name) {
      H <- contracted_nodes(H, s, v)
    }
  }

  # Contract right vertices into t
  for (v in right_vertex_set) {
    if (v != s && v != t && as.character(v) %in% igraph::V(H)$name) {
      H <- contracted_nodes(H, t, v)
    }
  }

  # Base case: small graph or direct edge
  if (igraph::vcount(H) <= 10 || igraph::are_adjacent(H, as.character(s), as.character(t))) {
    order_vertex_set <- igraph::V(G)$name
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

  # Filter cut to only vertices that exist in H
  cut <- cut[as.character(cut) %in% igraph::V(H)$name]

  # If no valid cut vertices, fall back to base case
  if (length(cut) == 0) {
    order_vertex_set <- igraph::V(G)$name
    order_vertex_set <- setdiff(order_vertex_set, as.character(left_vertex_set))
    order_vertex_set <- setdiff(order_vertex_set, as.character(right_vertex_set))
    order_vertex_set <- union(order_vertex_set, as.character(cut_set))
    order_vertex_set <- as.numeric(order_vertex_set)

    result <- order_by_NDS(G, vertex_order_dict, order_vertex_set, result_edge_list)
    return(list(vertex_order_dict = result$vertex_order_dict, result_edge_list = result$result_edge_list))
  }

  # Remove cut vertices (convert to character names)
  Hc <- igraph::delete_vertices(H, as.character(cut))

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
  # Convert undirected graph to directed (both directions) for st_min_cuts
  # st_min_cuts only works with directed graphs
  if (!igraph::is_directed(G)) {
    G_dir <- igraph::as_directed(G, mode = 'mutual')
  } else {
    G_dir <- G
  }

  # Convert to vertex ids in the directed graph
  s_id <- which(igraph::V(G_dir)$name == as.character(s))
  t_id <- which(igraph::V(G_dir)$name == as.character(t))

  if (length(s_id) == 0 || length(t_id) == 0) {
    return(numeric(0))
  }

  # Find all s-t cuts and get the minimum one
  cuts <- igraph::st_min_cuts(G_dir, s_id, t_id)
  if (length(cuts$cuts) == 0) {
    return(numeric(0))
  }

  min_cut <- cuts$cuts[[which.min(sapply(cuts$cuts, length))]]
  # Convert vertex sequence back to numeric vertex names
  as.numeric(igraph::V(G_dir)$name[as.numeric(min_cut)])
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
    for (n in igraph::V(G)$name) {
      if (igraph::degree(G, n) == 1) {
        ns <- igraph::neighbors(G, n)$name
        G <- igraph::delete_vertices(G, n)
        deg1_edges[[length(deg1_edges) + 1]] <- c(as.numeric(ns[1]), as.numeric(n))
        found <- TRUE
        break
      }
    }
  }

  # Check if it's a tree (no edges left)
  if (igraph::ecount(G) == 0) {
    e <- deg1_edges[[length(deg1_edges)]]
    # Add vertices back if they don't exist
    existing_vertices <- igraph::V(G)$name
    vertices_to_add <- setdiff(as.character(e), existing_vertices)
    if (length(vertices_to_add) > 0) {
      G <- igraph::add_vertices(G, length(vertices_to_add), name = vertices_to_add)
    }
    G <- igraph::add_edges(G, as.character(e))
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
    for (n in igraph::V(G)$name) {
      if (igraph::degree(G, n) == 2) {
        walk <- c(n)
        ns <- igraph::neighbors(G, n)$name

        c <- ns[1]
        walk <- c(c, walk)

        while (igraph::degree(G, c) == 2) {
          if (n == c) {
            return(list(
              G = G, deg1_edges = deg1_edges, deg2_dict = deg2_dict,
              deg2_cycle = deg2_cycle, is_cycle = TRUE, is_tree = FALSE
            ))
          }

          nc <- igraph::neighbors(G, c)$name
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
          nc <- igraph::neighbors(G, c)$name
          if (nc[1] != walk[length(walk) - 1]) {
            c <- nc[1]
          } else {
            c <- nc[2]
          }
          walk <- c(walk, c)
        }

        # Remove internal vertices of walk
        G <- igraph::delete_vertices(G, walk[2:(length(walk) - 1)])

        # Convert walk back to numeric for storage
        walk_numeric <- as.numeric(walk)

        if (walk[1] == walk[length(walk)]) {
          # It's a cycle
          key <- as.character(walk[1])
          if (!(key %in% names(deg2_cycle))) {
            deg2_cycle[[key]] <- list()
          }
          deg2_cycle[[key]][[length(deg2_cycle[[key]]) + 1]] <- walk_numeric
        } else {
          # Add edge between endpoints
          has_e <- igraph::are_adjacent(G, walk[1], walk[length(walk)])
          G <- igraph::add_edges(G, c(walk[1], walk[length(walk)]))

          key <- paste(minmax(as.numeric(walk[1]), as.numeric(walk[length(walk)])), collapse = '-')
          if (!(key %in% names(deg2_dict))) {
            deg2_dict[[key]] <- list()
          }
          deg2_dict[[key]][[length(deg2_dict[[key]]) + 1]] <- list(walk = walk_numeric, has_e = has_e)
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
  v <- as.numeric(edge_list[[1]][1])
  prev_v <- -1

  while (length(result_edge_list) < length(edge_list)) {
    vs <- Filter(function(x) as.numeric(x[1]) == v || as.numeric(x[2]) == v, edge_list)
    next_v <- setdiff(union(as.numeric(vs[[1]]), as.numeric(vs[[2]])), c(prev_v, v))[1]
    result_edge_list[[length(result_edge_list) + 1]] <- as.numeric(minmax(v, next_v))
    prev_v <- v
    v <- next_v
  }

  result_edge_list
}

# Check if edge list forms a connected order
check_connected_order <- function(edge_list) {
  touched_vertices <- c()

  for (i in 1:length(edge_list)) {
    # Skip NULL or invalid edges
    if (is.null(edge_list[[i]]) || length(edge_list[[i]]) < 2) {
      return(FALSE)
    }

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

    for (i in seq_along(result_edge_list)) {
      e <- minmax(result_edge_list[[i]][1], result_edge_list[[i]][2])
      key <- paste(e, collapse = '-')

      if (key %in% names(deg2_dict)) {
        found <- TRUE

        # Process each walk for this edge
        # Each walk corresponds to one occurrence of this edge in a multigraph
        insert_pos <- i
        for (x in deg2_dict[[key]]) {
          # Find and remove this edge occurrence
          edge_found <- FALSE
          for (j in insert_pos:length(result_edge_list)) {
            e_j <- minmax(result_edge_list[[j]][1], result_edge_list[[j]][2])
            if (all(e_j == e)) {
              # Remove this occurrence
              result_edge_list <- result_edge_list[-j]
              insert_pos <- j
              edge_found <- TRUE
              break
            }
          }

          if (!edge_found) {
            # Edge not found, skip this walk
            next
          }

          # Get vertices touched so far
          if (insert_pos > 1) {
            touched_vertices <- unique(unlist(result_edge_list[1:(insert_pos - 1)]))
          } else {
            touched_vertices <- c()
          }

          walk <- x$walk

          # Build list of edges to insert
          edges_to_insert <- list()
          if (walk[1] %in% touched_vertices) {
            # Insert walk edges in reverse order
            for (j in (length(walk) - 1):1) {
              edges_to_insert[[length(edges_to_insert) + 1]] <- c(walk[j], walk[j + 1])
            }
          } else {
            # Insert walk edges in order
            for (j in 1:(length(walk) - 1)) {
              edges_to_insert[[length(edges_to_insert) + 1]] <- c(walk[j], walk[j + 1])
            }
          }

          # Insert all edges at position insert_pos-1
          if (insert_pos == 1) {
            result_edge_list <- c(edges_to_insert, result_edge_list)
            insert_pos <- insert_pos + length(edges_to_insert)
          } else {
            result_edge_list <- c(
              result_edge_list[1:(insert_pos - 1)],
              edges_to_insert,
              result_edge_list[insert_pos:length(result_edge_list)]
            )
            insert_pos <- insert_pos + length(edges_to_insert)
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

        for (i in seq_along(result_edge_list)) {
          if (v_num %in% result_edge_list[[i]]) {
            found <- TRUE

            # Build list of cycle edges to insert
            edges_to_insert <- list()
            for (cy in deg2_cycle[[v]]) {
              for (j in 1:(length(cy) - 1)) {
                c0 <- cy[j]
                c1 <- cy[j + 1]
                edges_to_insert[[length(edges_to_insert) + 1]] <- minmax(c0, c1)
              }
            }

            # Insert all edges at once after position i
            if (i == length(result_edge_list)) {
              result_edge_list <- c(result_edge_list, edges_to_insert)
            } else {
              result_edge_list <- c(
                result_edge_list[1:i],
                edges_to_insert,
                result_edge_list[(i + 1):length(result_edge_list)]
              )
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

    for (e_idx in seq_along(deg1_edges)) {
      e <- deg1_edges[[e_idx]]

      for (i in seq_along(result_edge_list)) {
        if (e[1] %in% result_edge_list[[i]]) {
          # Insert edge after position i
          if (i == length(result_edge_list)) {
            result_edge_list <- c(result_edge_list, list(e))
          } else {
            result_edge_list <- c(
              result_edge_list[1:i],
              list(e),
              result_edge_list[(i + 1):length(result_edge_list)]
            )
          }
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
  # Convert to character to preserve vertex names across deletions
  edges_df <- do.call(rbind, lapply(edge_list, as.character))
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
    edges_list <- lapply(1:nrow(edges_vec), function(i) as.numeric(edges_vec[i, ]))
    result_edge_list <- get_cycle(edges_list)
    result_edge_list <- recover_deg12(deg1_edges, deg2_dict, deg2_cycle, result_edge_list)
  } else if (is_tree) {
    # Handle tree case
    edges_vec <- igraph::as_edgelist(G)
    result_edge_list <- lapply(1:nrow(edges_vec), function(i) as.numeric(edges_vec[i, ]))
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

  # Ensure all edges are numeric
  result_edge_list <- lapply(result_edge_list, as.numeric)

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

# Helper: detect if input is an adjacency list (list of neighbor vectors)
# vs a list of edge pairs (each element has exactly 2 elements)
is_adjacency_list <- function(x) {
  if (!is.list(x) || length(x) == 0) {
    return(FALSE)
  }
  # Adjacency lists typically have varying lengths or lengths != 2
  # Edge lists have all elements of length 2
  lengths <- vapply(x, length, integer(1))
  # If all lengths are 2, it's likely an edge list

  # If lengths vary or any length != 2, it's an adjacency list
  !all(lengths == 2)
}

# Helper: convert adjacency list to edge matrix
# adj is 0-indexed (as from redist.adjacency)
adj_to_edges <- function(adj) {
  # Convert to 1-indexed
  adj <- lapply(adj, function(x) x + 1L)
  # Build edge list (only include edges where j > i to avoid duplicates)
  edges <- do.call(rbind, lapply(seq_along(adj), function(i) {
    neighbors <- adj[[i]]
    neighbors <- neighbors[neighbors > i]
    if (length(neighbors) > 0) cbind(i, neighbors) else NULL
  }))
  edges
}

ndscut <- function(edges) {
  # Convert input to a list of edge pairs if it's not already
  edge_list <- list()

  if (is.matrix(edges) || is.data.frame(edges)) {
    # Check for empty input first
    if (nrow(edges) == 0) {
      message('The input graph is empty.')
      return(NULL)
    }
    # If input is a matrix or data frame, convert to a list of pairs
    for (i in seq_len(nrow(edges))) {
      edge_list[[i]] <- c(edges[i, 1], edges[i, 2])
    }
  } else if (is.list(edges) && is_adjacency_list(edges)) {
    # Convert adjacency list to edge matrix, then to edge list
    edge_mat <- adj_to_edges(edges)
    if (is.null(edge_mat) || nrow(edge_mat) == 0) {
      message('The input graph is empty.')
      return(NULL)
    }
    for (i in seq_len(nrow(edge_mat))) {
      edge_list[[i]] <- c(edge_mat[i, 1], edge_mat[i, 2])
    }
  } else if (is.list(edges)) {
    # If already a list of edge pairs, use it directly
    edge_list <- edges
  } else {
    stop('Input must be a matrix, data frame, adjacency list, or list of edge pairs')
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
    # Skip NULL or invalid edges
    if (!is.null(e) && length(e) >= 2) {
      c <- minmax(e[1], e[2])
      result[[length(result) + 1]] <- c(c[1], c[2])
    }
  }

  # Final connectivity check
  if (!check_connected_order(new_edge_list)) {
    warning('Output edge order is not connected')
  }

  result
}
