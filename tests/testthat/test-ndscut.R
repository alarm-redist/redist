skip_if_not_installed('igraph')

result_to_matrix <- function(result) {
    do.call(rbind, result)
}

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
  result <- ndscut(matrix(ncol = 2, nrow = 0))
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

test_that('ndscut produces equivalent frontier size for Iowa', {
  # old python implementation
  python_ordered <- matrix(
    c(23L, 49L, 23L, 53L, 49L, 53L, 16L, 23L, 16L, 53L, 31L, 49L,
      31L, 53L, 23L, 82L, 16L, 82L, 16L, 70L, 70L, 82L, 53L, 57L, 16L,
      57L, 28L, 53L, 28L, 31L, 28L, 57L, 16L, 52L, 52L, 70L, 52L, 57L,
      22L, 31L, 22L, 28L, 10L, 57L, 10L, 28L, 58L, 70L, 52L, 58L, 6L,
      57L, 6L, 10L, 22L, 33L, 10L, 33L, 52L, 92L, 58L, 92L, 48L, 52L,
      6L, 48L, 48L, 92L, 7L, 10L, 6L, 7L, 3L, 22L, 3L, 96L, 33L, 96L,
      44L, 58L, 44L, 92L, 6L, 86L, 7L, 86L, 29L, 58L, 29L, 44L, 54L,
      92L, 48L, 54L, 51L, 92L, 44L, 51L, 51L, 54L, 48L, 79L, 79L, 86L,
      54L, 79L, 44L, 56L, 29L, 56L, 44L, 89L, 51L, 89L, 56L, 89L, 54L,
      90L, 51L, 90L, 54L, 62L, 62L, 79L, 62L, 90L, 50L, 79L, 50L, 62L,
      26L, 89L, 26L, 90L, 68L, 90L, 62L, 68L, 62L, 63L, 50L, 63L, 63L,
      68L, 4L, 26L, 4L, 68L, 59L, 68L, 59L, 63L, 4L, 93L, 59L, 93L,
      27L, 93L, 9L, 33L, 7L, 9L, 19L, 33L, 19L, 96L, 9L, 19L, 7L, 38L,
      38L, 86L, 64L, 86L, 50L, 64L, 38L, 64L, 7L, 12L, 9L, 12L, 12L,
      38L, 45L, 96L, 19L, 45L, 50L, 77L, 63L, 77L, 63L, 91L, 59L, 91L,
      77L, 91L, 50L, 85L, 64L, 85L, 77L, 85L, 20L, 59L, 20L, 27L, 20L,
      91L, 38L, 42L, 42L, 64L, 42L, 85L, 19L, 34L, 12L, 34L, 12L, 35L,
      35L, 42L, 61L, 91L, 20L, 61L, 45L, 66L, 34L, 66L, 20L, 88L, 61L,
      88L, 27L, 80L, 80L, 88L, 80L, 87L, 17L, 34L, 17L, 35L, 17L, 66L,
      8L, 77L, 8L, 85L, 40L, 85L, 40L, 42L, 8L, 40L, 25L, 77L, 25L,
      61L, 8L, 25L, 1L, 61L, 1L, 88L, 2L, 88L, 2L, 87L, 1L, 2L, 66L,
      98L, 17L, 98L, 35L, 99L, 40L, 99L, 73L, 87L, 8L, 94L, 40L, 94L,
      94L, 99L, 8L, 37L, 25L, 37L, 37L, 94L, 25L, 39L, 1L, 39L, 37L,
      39L, 17L, 41L, 41L, 99L, 1L, 15L, 2L, 15L, 2L, 69L, 69L, 73L,
      15L, 69L, 46L, 99L, 46L, 94L, 95L, 98L, 41L, 95L, 13L, 94L, 13L,
      37L, 14L, 37L, 14L, 39L, 13L, 14L, 5L, 39L, 5L, 15L, 5L, 14L,
      36L, 73L, 36L, 65L, 65L, 69L, 15L, 78L, 69L, 78L, 65L, 78L, 15L,
      83L, 5L, 83L, 78L, 83L, 76L, 94L, 46L, 76L, 13L, 76L, 41L, 55L,
      46L, 55L, 55L, 95L, 13L, 81L, 14L, 81L, 14L, 24L, 24L, 83L, 24L,
      81L, 43L, 78L, 43L, 83L, 24L, 43L, 74L, 76L, 55L, 74L, 11L, 76L,
      11L, 81L, 47L, 81L, 24L, 47L, 32L, 55L, 32L, 74L, 24L, 67L, 43L,
      67L, 21L, 74L, 11L, 21L, 11L, 18L, 18L, 47L, 47L, 97L, 67L, 97L,
      18L, 97L, 30L, 32L, 21L, 30L, 21L, 71L, 18L, 71L, 18L, 75L, 75L,
      97L, 30L, 72L, 71L, 72L, 71L, 84L, 75L, 84L, 60L, 72L, 60L, 84L
    ), ncol = 2, byrow = TRUE)

  # Get R result from adjacency list
  adj_ia <- redist.adjacency(iowa)
  result <- ndscut(adj_ia)
  result_mat <- result_to_matrix(result)

  # Basic checks
  expect_equal(nrow(result_mat), 222)
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
