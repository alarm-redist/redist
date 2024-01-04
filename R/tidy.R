# tidy accessor/helper functions ----


#' Access the Current `redist_plans()` Object
#'
#' Useful inside piped expressions and `dplyr` functions.
#'
#' @returns A `redist_plans` object, or `NULL` if not called from inside a
#'   `dplyr` function.
#'
#' @examples
#' pl()
#'
#' @export
pl <- function() {
    cur_plans()
}

#' Helper function to get current plans/map objects
#' Traverses call stack to find plans object passed to dplyr verbs
#' @noRd
get_cur_df <- function(dplyr_funcs) {
    calls <- sys.calls()
    frames <- sys.frames()
    for (i in rev(seq_along(calls))) {
        call <- calls[[i]]
        frame <- frames[[i]]
        if (is.null(rlang::call_name(call))) next
        if (any(vapply(dplyr_funcs,
            function(x) identical(x, rlang::frame_fn(frame)),
            logical(1)))) {
            return(rlang::env_get(frame, ".data"))
        }
    }
    NULL
}

#' Helper function to get current map object
#' @noRd
cur_map <- function(verbs = c("mutate", "summarize", "merge_by",
                        "filter", "arrange", "transmute")) {
    get_cur_df(list(mutate = mutate.redist_map,
        transmute = transmute.redist_map,
        summarize = summarise.redist_map,
        merge_by = merge_by,
        filter = filter.redist_map,
        arrange = arrange.redist_map)[verbs])
}

#' Helper function to get current plans object
#' @noRd
cur_plans <- function(verbs = c("mutate", "summarize", "filter",
                          "arrange", "transmute")) {
    get_cur_df(list(mutate = mutate.redist_plans,
        transmute = transmute.redist_plans,
        summarize = summarise.redist_plans,
        filter = filter.redist_plans,
        arrange = arrange.redist_plans)[verbs])
}
#

