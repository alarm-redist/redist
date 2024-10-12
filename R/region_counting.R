# Creates a hash map that counts the number of times a specific region appears
# in a matrix of plans or partial plans
get_region_count_map <- function(plan_mat, N){
    # create hash map
    region_count_map <- hash::hash()

    V <- nrow(plan_mat)

    # now iterate over each column
    for (col_num in 1:ncol(plan_mat)) {

        # list mapping district number to precincts with district
        district_list <- vector("list", length = N)

        # now check which district each vertex is in
        for (vertex_num in 1:V) {

            district_id <- plan_mat[vertex_num, col_num]
            # now add this vertex to the set of vertices associated with that
            # district id

            # Check if
            if(is.null(district_list[[district_id]])){
                district_list[[district_id]] <- c(vertex_num)
            }else{ # else if already there just append
                district_list[[district_id]] <- append(
                    district_list[[district_id]], vertex_num
                )
            }
        }

        # now add count of districts to hash map
        for(district_vertices in district_list){
            # convert set of precincts to string
            str_district_vertices <- toString(district_vertices)
            # check if district already counted
            if(hash::has.key(str_district_vertices, region_count_map)){
                region_count_map[[str_district_vertices]] <- region_count_map[[str_district_vertices]] + 1
            }else{
                region_count_map[[str_district_vertices]] <- 1
            }
        }
    }

    return(region_count_map)
}


get_region_count_df <- function(plan_mat, N){
    district_count_map <- get_region_count_map(plan_mat, N)

    district_count_df <- data.frame(
        district=names(hash::values(district_count_map, USE.NAMES=TRUE)),
        d_count=hash::values(district_count_map, USE.NAMES=FALSE),
        row.names=NULL
    )  %>%
        arrange(desc(d_count))

    return(district_count_df)
}

#' @export
get_district_count_df <- function(plans, chain_num){

    district_count_map <- plans %>%
            filter(chain == chain_num) %>%
            as.matrix() %>%
            get_region_count_map(attr(plans, "ndist"))

    district_count_df <- data.frame(
        district=names(hash::values(district_count_map, USE.NAMES=TRUE)),
        d_count=hash::values(district_count_map, USE.NAMES=FALSE),
        row.names=NULL
    )  %>%
        arrange(desc(d_count))

    return(district_count_df)
}





