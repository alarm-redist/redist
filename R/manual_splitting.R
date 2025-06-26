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
#' @param partial_plan_labels A vector or matrix of plan id labels. It should be
#' zero indexed. If its a matrix then only the first column is used.
#' @param region_sizes A vector or matrix of region sizes. If its a matrix then
#' only the first column is used.
#' @param region_id_to_draw_tree_on The id of the region to draw a tree on
#' @param verbose Whether or not to print the plan before a tree is drawn on it
#'
#' @returns A list with the following
#' * `uncut_tree` - he spanning tree drawn stored as a 0-indexed directed
#'   adjacency list.
#' * `root` - The 0-indexed root of the tree.
#' * `num_attempts` -  The number of attempts it took to draw the tree.
#' * `pop_below` - The population below each vertex in `uncut_tree` ie
#'    the population induced by removing the edge terminating in that vertex
#' * `uncut_tree_vertex_parents`: The parents of each zero-indexed vertex in
#'    the tree. A value of -1 means that vertex is the root and -2 means that
#'    vertex is not in the tree.
#'
#' @export
draw_tree_on_region <- function(
        map, counties, partial_plan_labels, region_sizes,
        region_id_to_draw_tree_on, verbose
){
    ndists <- attr(map, "ndists")



    # get the map parameters
    map_params <- get_map_parameters(map, counties)

    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds


    if (!is.matrix(partial_plan_labels)) {
        # Convert to a matrix
        partial_plan_labels <- matrix(partial_plan_labels, nrow = length(partial_plan_labels))
    }
    num_regions <- unique(partial_plan_labels[,1]) |> length()
    validate_initial_region_id_mat(partial_plan_labels, V, 1, num_regions)

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
        region_ids=partial_plan_labels, region_sizes=region_sizes,
        verbose
    )

    return(return_list)

}




#' Manually Split A Multidistrict within a Plan
#'
#' Splits a multidistrict in a plan according to the procedure in ???.
#'
#' @inheritParams redist_smc
#' @param counties A vector containing county (or other administrative or
#' geographic unit) labels for each unit, which may be integers ranging from 1
#' to the number of counties, or a factor or character vector.  If provided,
#' the algorithm will only generate spanning trees where removing an edge splits
#' at most one county.
#' @param partial_plan_labels A vector or matrix of plan id labels. It should be
#' zero indexed. If its a matrix then only the first column is used.
#' @param region_sizes A vector or matrix of region sizes. If its a matrix then
#' only the first column is used.
#' @param region_id_to_split The id of the region to split
#' @param verbose Whether or not to print the number of iterations until a
#' successful split is made and print the split plan
#'
#' @returns A list with the following elements
#' * `num_attempts` - The number of attempts it took to draw the tree.
#' * `uncut_tree` - The spanning tree drawn stored as a 0-indexed directed
#'   adjacency list.
#' * `region_id_that_was_split` - The label of the region that was split.
#' * `region_id_that_was_split` - The label of the region that was split.
#' * `region_sizes` - 1-indexed region sizes.
#' * `partial_plan_labels` - 0-indexed region labels of the split partial
#'   plan.
#' * `region_pops` - 1-indexed region populations.
#' * `num_regions` - The number of regions in the split plan.
#' * `num_districts` - The number of districts in the split plan.
#' * `uncut_tree` - The uncut spanning tree drawn stored as a 0-indexed directed
#'    adjacency list.
#' * `uncut_tree_root` - The root of the uncut spanning tree drawn.
#' * `cut_tree` - The cut spanning tree drawn stored as a 0-indexed directed
#'   adjacency list.
#' * `pop_below` - The population below each vertex in `uncut_tree` ie
#'    the population induced by removing the edge terminating in that vertex
#' * `uncut_tree_vertex_parents`: The parents of each zero-indexed vertex in
#'    the tree. A value of -1 means that vertex is the root and -2 means that
#'    vertex is not in the tree.
#' * `new_region1_id` - The label of the first of the two new split
#'   regions.
#' * `new_region1_tree_root` - The root of the cut tree associated with
#'   the first of the two new split regions.
#' * `new_region1_size` - The size of the first of the two new split
#'   regions.
#' * `new_region1_pop` - The population of the first of the two new split
#'   regions.
#' * `new_region2_id` - The label of the second of the two new split
#'   regions.
#' * `new_region2_tree_root` - The root of the cut tree associated with
#'   the second of the two new split regions.
#' * `new_region2_size` - The size of the second of the two new split
#'   regions.
#' * `new_region2_pop` - The population of the second of the two new split
#'   regions.
#' * `cut_edge` - The two zero-indexed vertices making up the edge that
#'   was cut from `uncut_tree`.
#'
#' @export
manually_split_a_multidistrict <- function(
        map, counties, partial_plan_labels, region_sizes,
        region_id_to_split, split_dval_min=1, split_dval_max=NULL,
        split_district_only, k_param=1, verbose=F){

    ndists <- attr(map, "ndists")

    # get the map parameters
    map_params <- get_map_parameters(map, counties)

    # check the plan inputs
    if (!is.matrix(partial_plan_labels)) {
        # Convert to a matrix
        partial_plan_labels <- matrix(partial_plan_labels, nrow = length(partial_plan_labels))
    }
    num_regions <- unique(partial_plan_labels[,1]) |> length()
    validate_initial_region_id_mat(partial_plan_labels, map_params$V, 1, num_regions)

    if(region_id_to_split >= num_regions){
        cli_abort("Input {.arg region_id_to_split} must be between zero and the number of regions ({num_regions}).")
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

    if(is.null(split_dval_max)){
        if(split_district_only){
            split_dval_max <- 1
        }else{
            split_dval_max <- region_sizes[region_id_to_split+1,1]
        }
    }


    split_algout <- perform_a_valid_multidistrict_split(
        adj_list=map_params$adj_list,
        counties=map_params$counties,
        pop=map_params$pop,
        ndists=ndists,
        num_regions=num_regions,
        num_districts=num_districts,
        region_id_to_split=region_id_to_split,
        target=map_params$pop_bounds[2],
        lower=map_params$pop_bounds[1],
        upper=map_params$pop_bounds[3],
        region_ids=partial_plan_labels,
        region_sizes=region_sizes,
        split_dval_min=split_dval_min,
        split_dval_max=split_dval_max,
        split_district_only=split_district_only,
        verbose=verbose,
        k_param=k_param
    )


    # find the edge that was cut between the two trees
    cut_edge <- list()
    # for each vertex check if they have any neighbors not shared with the
    # cut tree
    for (v in seq_len(length(split_algout$uncut_tree))) {
        # find vertices in the precut tree vs the cut tree
        vertex_nbor_diff <- setdiff(split_algout$uncut_tree[[v]], split_algout$cut_tree[[v]])
        if(length(vertex_nbor_diff) > 1){
            cli_abort("More than 1 edge was cut from the tree!")
        }else if(length(vertex_nbor_diff) == 1){
            cut_edge <- c(v-1, vertex_nbor_diff[1])
        }
    }

    split_algout$cut_edge <- cut_edge

    return(split_algout)
}







#' Splits a map into `ndist` districts
#'
#' Splits a map all the way into `ndist` districts and saves the splitting
#' information from each step including the pre-cut and cut spanning tree drawn
#' at each step that was successfully split. Very useful for debugging and
#' visual inspection of internal results.
#'
#'
#' see \code{\link{manually_split_a_multidistrict}}.
#'
#' @inheritParams redist_smc
#'
#' @returns A list of length `ndists - 1` where each element of the list is
#' itself a list where each elements are the return values of
#' \code{\link{manually_split_a_multidistrict}} along with two new additional
#' elements of the list:
#'
#' * `cut_edges` - A list of 2-element vectors (u,v) where u,v are the zero-indexed
#' vertices of the edges removed from the drawn spanning trees at each step to
#' create the two new regions
#' * `spanning_forest` - The spanning forest of the cut spanning trees drawn up
#' to that point.
#' @export
split_entire_map_into_plan <- function(
        map, counties,
        min_region_cut_sizes = NULL, max_region_cut_sizes = NULL,
        split_district_only=F, k_param=1, verbose=F
){
    # get N
    ndists <- attr(map, "ndists")

    if(is.null(min_region_cut_sizes)){
        min_region_cut_sizes <- rep(1, ndists-1)
    }
    if(is.null(max_region_cut_sizes)){
        if(split_district_only){
            max_region_cut_sizes <- rep(1, ndists-1)
        }else{
            max_region_cut_sizes <- rev(seq_len(ndists-1))
        }
    }

    # validate the cut sizes
    validate_cut_sizes(ndists, ndists-1, min_region_cut_sizes, max_region_cut_sizes)

    partial_plan_labels <- matrix(0L, ncol = 1, nrow = nrow(map))
    region_sizes <- matrix(0L, ncol = 1, nrow = ndists)
    region_sizes[1,1] <- ndists

    region_id_to_split <- 0

    num_splits <- ndists-1

    split_output_list <- vector("list", length = num_splits)


    for (split_num in seq_len(num_splits) ) {
        # split a multidistrict
        split_output <- manually_split_a_multidistrict(
            map=map,
            counties=counties,
            partial_plan_labels=partial_plan_labels,
            region_sizes=region_sizes,
            region_id_to_split=region_id_to_split,
            split_dval_min=min_region_cut_sizes[split_num],
            split_dval_max=max_region_cut_sizes[split_num],
            split_district_only=split_district_only,
            k_param=k_param,
            verbose=verbose
            )

        # a list of all the split trees for each step
        a_split_trees_list <-vector("list", length = split_num+1)
        # a list with the previously split trees and the one drawn
        # before splitting
        a_presplit_trees_list <-vector("list", length = split_num)

        # for the first split we just use the two trees
        if(split_num == 1){
            # for first one forest is the just two cut trees
            split1_tree1 <- get_subtree_from_root(
                split_output$cut_tree,
                split_output$new_region1_tree_root
            )
            split1_tree2 <- get_subtree_from_root(
                split_output$cut_tree,
                split_output$new_region2_tree_root
            )


            a_split_trees_list[[(split_output$new_region1_id+1)]] <- split1_tree1
            a_split_trees_list[[(split_output$new_region2_id+1)]] <- split1_tree2

            # presplit tree is entire thing for this
            a_presplit_trees_list[[1]] <- split_output$uncut_tree
        }else{ # else we use previous one

            # forest will be trees from previous step but replace the split one
            # with the two new trees
            old_split_region_id <- split_output$region_id_that_was_split + 1
            # add one bc R is 1 indexed
            new_region1_id <- split_output$new_region1_id+1
            new_region2_id <- split_output$new_region2_id+1

            # check new region 1 has same id as old split id
            assertthat::assert_that(
                old_split_region_id == new_region1_id,
                msg = "Splitting procedure was changed. We expected the new region 1
            id to be the same as the id of the old split region."
            )
            # check new region 2 is just the end of the list
            assertthat::assert_that(
                split_num+1 == new_region2_id,
                msg = "Splitting procedure was changed. We expected the new region 2
            id to be the number of regions minus 1, so the split number"
            )

            # copy over all the trees from regions that didn't change
            unchanged_old_regions <- unique(
                split_output_list[[split_num-1]]$partial_plan_labels + 1
            ) |> setdiff(old_split_region_id)


            for (i in unchanged_old_regions) {
                # make the first split_num trees the same as previous
                a_split_trees_list[[i]] <- split_output_list[[split_num-1]]$post_split_trees_list[[i]]
                a_presplit_trees_list[[i]] <- split_output_list[[split_num-1]]$post_split_trees_list[[i]]
            }
            a_presplit_trees_list[[old_split_region_id]] <- split_output$uncut_tree

            # replace two new ones
            split_tree1 <- get_subtree_from_root(
                split_output$cut_tree,
                split_output$new_region1_tree_root
            )
            split_tree2 <- get_subtree_from_root(
                split_output$cut_tree,
                split_output$new_region2_tree_root
            )
            a_split_trees_list[[new_region1_id]] <- split_tree1
            a_split_trees_list[[new_region2_id]] <- split_tree2

        }

        # store the list of trees
        split_output$post_split_trees_list <- a_split_trees_list
        split_output$presplit_trees_list <- a_presplit_trees_list


        # combine into forest
        split_forest <- lapply(
            1:length(split1_tree1),
            function(i) lapply(a_split_trees_list, function(l) l[[i]]) |> unlist()
        )

        split_output$post_split_forest <- split_forest

        presplit_forest <- lapply(
            1:length(split1_tree1),
            function(i) lapply(a_presplit_trees_list, function(l) l[[i]]) |> unlist()
        )
        split_output$presplit_forest <- presplit_forest


        # store output
        split_output_list[[split_num]] <- split_output

        # add cut edges as all the cut_edge up to this point
        cut_edges <- lapply(
            seq_len(split_num),
            function(i) split_output_list[[i]]$cut_edge)

        split_output_list[[split_num]]$cut_edges <- cut_edges

        # now update the plan info
        partial_plan_labels <- split_output$partial_plan_labels
        region_sizes <- split_output$region_sizes


        num_regions <- split_output$num_regions
        num_districts <- split_output$num_districts

        # pick a new region to split
        # get indices
        if(num_regions - num_districts == 1){
            region_id_to_split <- which(region_sizes > 1) - 1
        }else if(num_regions - num_districts > 1){
            region_id_to_split <- sample(which(region_sizes > 1), 1) - 1
        }else if(verbose){
            print("Done!")
        }

    }


    return(split_output_list)


}


