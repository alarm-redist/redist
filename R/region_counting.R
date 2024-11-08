# Creates a hash map that counts the number of times a specific region appears
# in a matrix of plans or partial plans
#
# Only counts regions in district_nums. If null counts them all
get_region_count_map <- function(plan_mat, num_regions, region_ids_to_count = NULL){
    # check region ids only go from 1:num_regions
    assertthat::assert_that(
        assertthat::are_equal(
            1:num_regions,
            sort(unique(plan_mat[,1]))
        )
    )

    region_ids_to_count <- unique(region_ids_to_count)

    # check not counting an invalid region id
    assertthat::assert_that(
        all(region_ids_to_count %in% 1:num_regions)
    )

    region_to_count_vec <- rep(FALSE, num_regions)

    if(is.null(region_ids_to_count)){
        checking_region_ids <- FALSE
    }else{
        checking_region_ids <- TRUE
        region_to_count_vec[region_ids_to_count] <- TRUE
    }

    # create hash map
    region_count_map <- hash::hash()

    V <- nrow(plan_mat)

    # now iterate over each plan
    for (col_num in 1:ncol(plan_mat)) {

        # list mapping district number to precincts with district
        region_list <- vector("list", length = num_regions)

        # now check which district each vertex is in
        for (vertex_num in 1:V) {

            region_id <- plan_mat[vertex_num, col_num]
            # now add this vertex to the set of vertices associated with that
            # district id

            # skip if its not one to count
            if(checking_region_ids && !region_to_count_vec[region_id]){
                next
            }


            # Check if
            if(is.null(region_list[[region_id]])){
                region_list[[region_id]] <- c(vertex_num)
            }else{ # else if already there just append
                region_list[[region_id]] <- append(
                    region_list[[region_id]], vertex_num
                )
            }
        }

        region_list <- region_list[!sapply(region_list, is.null)]

        # now add count of districts to hash map
        for(region_vertices in region_list){
            # convert set of precincts to string
            str_region_vertices <- toString(region_vertices)
            # check if district already counted
            if(hash::has.key(str_region_vertices, region_count_map)){
                region_count_map[[str_region_vertices]] <- region_count_map[[str_region_vertices]] + 1
            }else{
                region_count_map[[str_region_vertices]] <- 1
            }
        }
    }

    return(region_count_map)
}


get_region_count_df <- function(plan_mat, N){
    region_count_map <- get_region_count_map(plan_mat, N)

    region_count_df <- data.frame(
        district=names(hash::values(region_count_map, USE.NAMES=TRUE)),
        d_count=hash::values(region_count_map, USE.NAMES=FALSE),
        row.names=NULL
    )  %>%
        arrange(desc(d_count))

    return(region_count_df)
}

#' @export
get_district_count_df <- function(plans, chain_nums, district_ids_to_count = NULL){

    # if not specific districts specified then just do all
    district_count_map <- plans %>%
        filter(chain %in% chain_nums) %>%
        as.matrix() %>%
        get_region_count_map(attr(plans, "ndist"), district_ids_to_count)



    district_count_df <- data.frame(
        district=names(hash::values(district_count_map, USE.NAMES=TRUE)),
        d_count=hash::values(district_count_map, USE.NAMES=FALSE),
        row.names=NULL
    )  %>%
        arrange(desc(d_count))

    return(district_count_df)
}







