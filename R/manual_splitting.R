#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2022/12/29
# Purpose: R wrapper for manual map splitting code
####################################################



#' Build New Tree from Vertex
#'
#' Takes a tree stored in a 0-indexed, directed graph form and returns the subtree
#' created by only traversing the tree from `root`.
#'
#' @param tree_adj_list The adjacency list of the tree. It should be zero indexed and
#' directed.
#' @param root The 0-indexed vertex to start building the new tree from.
#'
#' @return The subtree created by traversing the tree only starting from `root`
#'
#'
get_subtree_from_root <- function(tree_adj_list, root) {
    # Initialize variables
    visited <- rep(FALSE, length(tree_adj_list))
    new_tree_adj_list <- rep(list(integer(0)), length(tree_adj_list)) # New adjacency list to be returned

    stack <- list(root)  # Initialize stack with the root node

    # Perform iterative DFS
    while (length(stack) > 0) {
        # Pop the top element from the stack
        node <- stack[[length(stack)]]
        stack <- stack[-length(stack)]

        if (!visited[node + 1]) {
            # Mark node as visited
            visited[node + 1] <- TRUE

            # Add neighbors to the stack and update the new adjacency list
            for (neighbor in rev(tree_adj_list[[node + 1]])) {  # rev() ensures DFS goes deeper first
                if (!visited[neighbor + 1]) {
                    stack <- append(stack, neighbor)  # Push the neighbor to the stack
                    new_tree_adj_list[[node + 1]] <- c(new_tree_adj_list[[node + 1]], neighbor)  # Store the edge
                }
            }
        }
    }

    # Return the new adjacency list containing only the visited edges
    return(new_tree_adj_list)

}




#' Draw a spanning tree on a region of a plan
#'
#' Draws a spanning tree uniformly at random on a region of a plan.
#'
#' @inheritParams redist_smc
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided,
#' the algorithm will only generate spanning trees where removing an edge splits
#' at most one county.
#' @param vertex_region_ids A vector or matrix of plan id labels. It should be
#' zero indexed. If its a matrix then only the first column is used.
#' @param region_sizes A vector or matrix of region sizes. If its a matrix then
#' only the first column is used.
#' @param region_id_to_draw_tree_on The id of the region to draw a tree on
#' @param verbose Whether or not to print the plan before a tree is drawn on it
#'
#' @return A list with the following
#' \itemize{
#'   \item{uncut_tree}{ - The spanning tree drawn stored as a 0-indexed directed
#'   adjacency list.}
#'   \item{root}{ - The 0-indexed root of the tree.}
#'   \item{num_attempts}{ - The number of attempts it took to draw the tree.}
#' }
#'
#' @export
draw_tree_on_region <- function(
        map, counties, vertex_region_ids, region_sizes,
        region_id_to_draw_tree_on, verbose
){
    ndists <- attr(state_map, "ndists")



    # get the map parameters
    map_params <- get_map_parameters(map, counties)

    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds


    if (!is.matrix(vertex_region_ids)) {
        # Convert to a matrix
        vertex_region_ids <- matrix(vertex_region_ids, nrow = length(vertex_region_ids))
    }
    num_regions <- unique(vertex_region_ids[,1]) |> length()
    validate_initial_region_id_mat(vertex_region_ids, V, 1, num_regions)

    if(region_id_to_draw_tree_on >= num_regions){
        cli_abort("Input {.arg region_id_to_draw_tree_on} must be between zero and the number of regions ({num_regions}).")
    }


    if (!is.matrix(region_sizes)) {
        # Convert to a matrix
        region_sizes <- matrix(region_sizes, nrow = length(region_sizes))
    }
    #
    if(nrow(region_sizes) != ndists){
        region_sizes <- rbind(
            region_sizes,
            matrix(0L, nrow = ndists - nrow(region_sizes), ncol = ncol(region_sizes))
            )
    }
    validate_initial_region_sizes_mat(region_sizes, ndists, 1, num_regions)
    num_districts <- sum(region_sizes[,1] == 1)


    return_list <- draw_a_tree_on_a_region(
        adj_list=adj_list, counties=counties, pop=pop,
        ndists=ndists, num_regions=num_regions, num_districts=num_districts,
        region_id_to_draw_tree_on=region_id_to_draw_tree_on,
        lower=pop_bounds[1], upper=pop_bounds[3],
        region_ids=vertex_region_ids, region_dvals=region_sizes,
        verbose
    )

    return(return_list)

}



manually_split_a_multidistrict <- function(){

}
