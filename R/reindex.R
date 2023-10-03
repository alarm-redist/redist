#' Reorders district numbers
#'
#' Ensures that for each column in the plans object, the first
#' district listed is 1, the second is 2, up to n districts. Assumes that all
#' columns have the same number of districts as the first.
#' @param plans A numeric vector (if only one map) or
#' matrix with one row for each precinct and one column for each map.
#'
#' @return integer matrix
#' @export
#'
#' @examples
#' cds <- matrix(c(rep(c(4L, 5L, 2L, 1L, 3L), 5),
#'     rep(c(5L, 4L, 3L, 2L, 1L), 2), rep(c(4L, 5L, 2L, 1L, 3L), 3)), nrow = 25)
#' redist.reorder(cds)
#'
redist.reorder <- function(plans) {
    # Check inputs
    if (!is.matrix(plans)) {
        if (is.numeric(plans)) {
            plans <- as.matrix(plans, ncol = 1)
        } else {
            cli_abort("{.arg plans} must be a matrix or integer vector.")
        }
    }

    apply(plans, 2, vctrs::vec_group_id)
}


#' Sink Plans to 1:ndists
#'
#' Takes a plan and renumbers it to be from 1:ndists
#'
#' @param plan vector of assignments, required.
#'
#' @return A vector with an ID that corresponds from 1:ndists, and attribute `n` indicating the number of districts.
#'
#' @concept prepare
#' @export
#' @examples
#' data(fl25_enum)
#' plan <- fl25_enum$plans[, 5118]
#' # Subset based on something:
#' plan <- plan[plan != 2]
#' plan <- redist.sink.plan(plan)
#' # Now plan can be used with redist.flip()
#' plan
#'
redist.sink.plan <- function(plan) {
    .Deprecated("vctrs::vec_group_id")
    vctrs::vec_group_id(plan)
}

