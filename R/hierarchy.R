normalize_hierarchy_levels <- function(adj, levels) {
  if (is.null(levels)) {
    return(NULL)
  }

  normalized <- lapply(levels, function(level) {
    if (anyNA(level)) {
      cli::cli_abort(
        "{.arg counties} hierarchy levels must not contain missing values."
      )
    }

    group <- vctrs::vec_group_id(level)
    component <- contiguity(adj, group)
    n_components <- tapply(component, group, max)
    split_unit <- n_components[group] > 1L
    component_label <- component
    component_label[!split_unit] <- 1L
    vctrs::vec_group_id(
      data.frame(level = level, component = component_label)
    )
  })

  if (length(normalized) > 1L) {
    for (level_index in seq_len(length(normalized) - 1L)) {
      nested <- tapply(
        normalized[[level_index]],
        normalized[[level_index + 1L]],
        function(x) length(unique(x)) == 1L
      )
      if (!all(unlist(nested))) {
        cli::cli_abort(
          "{.arg counties} hierarchy levels must be nested from coarser to finer."
        )
      }
    }
  }

  do.call(cbind, normalized)
}

project_plan_hierarchy <- function(adj, plan, levels) {
  if (is.null(levels)) {
    return(NULL)
  }

  if (is.null(dim(plan))) {
    if (length(plan) != length(adj)) {
      cli::cli_abort('{.arg plan} must have one value per unit in {.arg adj}.')
    }
    plan <- data.frame(plan = plan)
  } else {
    if (nrow(plan) != length(adj)) {
      cli::cli_abort('{.arg plan} must have one row per unit in {.arg adj}.')
    }
    plan <- as.data.frame(plan)
  }
  names(plan) <- paste0('plan_', seq_len(ncol(plan)))
  if (anyNA(plan)) {
    cli::cli_abort('{.arg plan} must not contain missing values.')
  }

  if (!is.list(levels)) {
    if (is.null(dim(levels))) {
      levels <- list(levels)
    } else {
      levels <- lapply(seq_len(ncol(levels)), function(index) levels[, index])
    }
  }
  if (length(levels) == 0L) {
    return(NULL)
  }
  if (any(lengths(levels) != length(adj))) {
    cli::cli_abort('Each hierarchy level must have one value per unit in {.arg adj}.')
  }
  if (any(vapply(levels, anyNA, logical(1)))) {
    cli::cli_abort('Hierarchy levels must not contain missing values.')
  }

  projected <- vector('list', length(levels))
  for (level_index in seq_along(levels)) {
    key <- plan
    key$unit <- levels[[level_index]]
    group <- vctrs::vec_group_id(key)
    key$component <- contiguity(adj, group)
    if (level_index > 1L) {
      key$parent <- projected[[level_index - 1L]]
    }
    projected[[level_index]] <- vctrs::vec_group_id(key)
  }

  do.call(cbind, projected)
}
