#' Convert perim_df to borderlength_mat
#'
#' @param n number of rows in the original object
#' @param perim_df a perim_df from \code{redist.prep.polsbypopper}
#'
#' @noRd
#' @return matrix for
perim_df_2_borderlength_mat <- function(n, perim_df) {
  blm <- matrix(0, n, n)

  for (i in 1:nrow(perim_df)) {
    if (perim_df$origin[i] == -1) {
      blm[ perim_df$touching[i], perim_df$touching[i]] <- blm[perim_df$touching[i], perim_df$touching[i]] + perim_df$edge[i]
    } else {
      blm[perim_df$origin[i], perim_df$touching[i]] <-blm[perim_df$origin[i], perim_df$touching[i]] + perim_df$edge[i]

      blm[perim_df$touching[i], perim_df$origin[i]] <- blm[perim_df$touching[i], perim_df$origin[i]] + perim_df$edge[i]
    }
  }

  return(blm)
}
