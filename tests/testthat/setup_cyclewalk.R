# Shared setup for cyclewalk tests

# Iowa map (iowa data loaded in setup.R)
ia <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)

# 4x4 grid with unit population for distribution tests
grid_sf <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 4, ymax = 4)) |>
    sf::st_as_sfc() |>
    sf::st_make_grid(n = c(4, 4)) |>
    sf::st_sf(geometry = _)
grid_sf <- grid_sf[c(13, 14, 9, 10, 15, 16, 11, 12, 5, 6, 1, 2, 7, 8, 3, 4), ]
grid_sf$adj <- redist.adjacency(grid_sf)
grid_sf$pop <- rep(1L, 16)
grid_sf$init <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L)
grid <- redist_map(grid_sf, existing_plan = init, pop_bounds = c(3.99, 4, 4.01),
                   adj = grid_sf$adj) |> suppressWarnings()

# Tolerance checker for distribution tests
# Allows 10% relative tolerance for large probabilities (>1%)
# Allows 40% relative tolerance for small probabilities
is_close <- function(obs, exp) {
    ratio <- obs / exp
    if (exp > 0.01) {
        ratio >= 0.9 && ratio <= 1.1
    } else {
        ratio >= 0.6 && ratio <= 1.4
    }
}

# Stricter tolerance for longer chains
is_close_tight <- function(obs, exp) {
    ratio <- obs / exp
    ratio >= 0.95 && ratio <= 1.05
}
