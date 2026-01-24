skip_if_not_installed('igraph')

result_to_matrix <- function(result) do.call(rbind, result)

check_isomorphic <- function(edges1, edges2) {
  G1 <- igraph::graph_from_edgelist(edges1, directed = FALSE)
  G2 <- igraph::graph_from_edgelist(edges2, directed = FALSE)
  igraph::isomorphic(G1, G2)
}

is_connected_order <- function(result) {
  touched <- result[[1]]
  for (i in seq_along(result)[-1]) {
    if (!any(result[[i]] %in% touched)) {
      return(FALSE)
    }
    touched <- union(touched, result[[i]])
  }
  TRUE
}

test_that('ndscut preserves graph structure', {
  cases <- list(
    path = matrix(c(1, 2, 2, 3, 3, 4), ncol = 2, byrow = TRUE),
    triangle = matrix(c(1, 2, 2, 3, 3, 1), ncol = 2, byrow = TRUE),
    star = matrix(c(1, 2, 1, 3, 1, 4, 1, 5), ncol = 2, byrow = TRUE),
    grid = matrix(c(1, 2, 2, 3, 1, 4, 2, 5, 3, 6, 4, 5, 5, 6), ncol = 2, byrow = TRUE)
  )

  for (name in names(cases)) {
    edges <- cases[[name]]
    result <- ndscut(edges)
    result_mat <- result_to_matrix(result)

    expect_equal(length(result), nrow(edges), info = name)
    expect_true(check_isomorphic(edges, result_mat), info = name)
    expect_true(is_connected_order(result), info = name)
  }
})

test_that('ndscut handles empty input', {
  expect_message(result <- ndscut(matrix(ncol = 2, nrow = 0)), 'empty')
  expect_null(result)
})

test_that('ndscut accepts adjacency list (0-indexed)', {
  # Simple path as adjacency list
  adj_path <- list(c(1L), c(0L, 2L), c(1L, 3L), c(2L))
  result <- ndscut(adj_path)

  expect_equal(length(result), 3)
  expect_true(is_connected_order(result))

  # Triangle
  adj_tri <- list(c(1L, 2L), c(0L, 2L), c(0L, 1L))
  expect_equal(length(ndscut(adj_tri)), 3)
})

test_that('ndscut produces equivalent frontier size for FL25', {
  # Python output (ordered.dat) for reference
  python_ordered <- matrix(c(
    6, 23, 3, 23, 3, 6, 21, 23, 3, 21, 1, 6, 1, 3, 3, 4, 4, 21, 1, 4,
    2, 21, 1, 2, 2, 4, 6, 9, 1, 9, 18, 21, 2, 18, 1, 14, 9, 14, 1, 15,
    14, 15, 1, 13, 13, 15, 2, 16, 16, 18, 1, 12, 12, 13, 1, 11, 11, 12, 1, 10,
    10, 11, 2, 17, 16, 17, 1, 8, 8, 10, 1, 7, 7, 8, 2, 19, 17, 19, 2, 20,
    19, 20, 1, 5, 5, 7, 2, 22, 20, 22, 5, 22, 5, 24, 22, 24, 5, 25, 22, 25, 24, 25
  ), ncol = 2, byrow = TRUE)

  # Get R result from adjacency list
  result <- ndscut(adj)
  result_mat <- result_to_matrix(result)

  # Basic checks
  expect_equal(nrow(result_mat), 51)
  expect_true(is_connected_order(result))
  expect_true(check_isomorphic(python_ordered, result_mat))

  # Write both to temp files and compare frontier sizes
  dir <- tempdir()
  py_path <- file.path(dir, 'py_ordered')
  r_path <- file.path(dir, 'r_ordered')

  utils::write.table(data.frame(python_ordered), paste0(py_path, '.dat'),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  utils::write.table(data.frame(result_mat), paste0(r_path, '.dat'),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  py_frontier <- redist.calc.frontier.size(py_path)
  r_frontier <- redist.calc.frontier.size(r_path)

  # R should produce equivalent or better frontier size
  expect_equal(r_frontier$max, py_frontier$max)
  expect_equal(r_frontier$average, py_frontier$average, tolerance = 0.1)
})
