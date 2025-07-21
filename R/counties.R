#' Relabel Discontinuous Counties
#'
#' @param adj adjacency list
#' @param counties character vector of county names
#' @param simplify boolean - TRUE returns a numeric vector of ids, while FALSE
#' appends a number when there are multiple connected components.
#'
#' @return character vector of county names
#'
#' @importFrom dplyr group_by mutate ungroup rowwise
#'
#' @concept prepare
#' @export
#' @examples
#' set.seed(2)
#' data(fl25)
#' data(fl25_adj)
#' counties <- sample(c(rep("a", 20), rep("b", 5)))
#' redist.county.relabel(fl25_adj, counties)
#'
redist.county.relabel <- function(adj, counties,  simplify = TRUE) {
    if (length(adj) != length(counties)) {
        cli::cli_abort("{.arg adj} and {.arg counties} must have the same length.")
    }

    if ("numeric" %in% class(counties) || "integer" %in% class(counties)) {
        counties <- as.character(counties)
    }

    sorted <- sort(unique(counties))
    groups <- match(counties, sorted)

    component <- data.frame(counties = counties, comp = contiguity(adj, groups)) %>%
        group_by(counties) %>%
        mutate(comps = max(.data$comp)) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(countiescomp = ifelse(.data$comps > 1, paste0(counties, "-", .data$comp), counties)) %>%
        ungroup()
    if (simplify) {
        redist.county.id(component$countiescomp)
    } else {
        component$countiescomp
    }
}


#' Create County IDs
#'
#' @param counties vector of counties, required.
#'
#' @return A vector with an ID that corresponds from 1:n counties
#'
#' @concept prepare
#' @export
#' @examples
#' set.seed(2)
#' counties <- sample(c(rep("a", 20), rep("b", 5)))
#' redist.county.id(counties)
#'
redist.county.id <- function(counties) {
    if (class(counties) %in% c("character", "numeric", "integer")) {
        uc <- unique(sort(counties))
        county_id <- rep(0, length(counties))
        for (i in 1:length(counties)) {
            county_id[i] <- which(uc == counties[i])
        }
    } else {
        cli::cli_abort("{.arg counties} must be a character, numeric, or integer vector.")
    }

    county_id
}
