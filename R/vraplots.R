#' Majority Minority Plots
#'
#' @param grouppercent output from redist.group.percent
#' @param type string in 'hist', 'toptwo', or 'box'
#' @param title ggplot title
#'
#' @importFrom ggplot2 geom_histogram theme_bw geom_point geom_vline geom_hline annotate geom_boxplot
#' @importFrom dplyr arrange desc first nth summarise
#'
#' @return ggplot
#' @concept plot
#' @export
redist.plot.majmin <- function(grouppercent, type = 'hist', title = ''){

  if(type == 'hist'){
    mm <- colSums(grouppercent > 0.5)

    p <- tibble(mm = mm) %>%
      ggplot(aes(x = mm)) +
      geom_histogram() +
      theme_bw() +
      labs(x = 'Minority Majority Districts', y = 'Count', title = title)

    return(p)

  } else if(type == 'toptwo'){
    tibble(blk_pct = c(grouppercent), district = rep(1:nrow(grouppercent), ncol(grouppercent)),
                 nloop = rep(1:ncol(grouppercent), each = nrow(grouppercent) )) %>%
      group_by(nloop) %>%
      arrange(desc(blk_pct), .by_group = TRUE) %>%
      mutate(group_id = row_number()) %>%
      filter(group_id <= 2) %>%
      summarise(first_blk = first(tot_blk_pct), second_blk = nth(tot_blk_pct, n = 2)) %>%
      ggplot(aes(x = first_blk, y = second_blk)) +
      geom_point() +
      geom_vline(xintercept = .5, color = 'blue', linetype = 'dotted') +
      geom_hline(yintercept = .5, color = 'blue', linetype = 'dotted') +
      labs(x = 'Larger Minority Percent', y = 'Smaller Minority Percent') +
      theme_bw() +
      annotate("rect", xmin = 0.5, ymin = 0.5, xmax = 1, ymax = 1, alpha = 0.2)


  } else if(type == 'box'){
    tibble(blk_pct = c(grouppercent), district = rep(1:nrow(grouppercent), ncol(grouppercent)),
           nloop = rep(1:ncol(grouppercent), each = nrow(grouppercent))) %>%
      group_by(district) %>%
      mutate(blk_pct = sort(blk_pct)) %>%
      ungroup() %>%
      ggplot(aes(x = district, y = blk_pct, group = district)) %>%
      geom_boxplot() +
      theme_bw() +
      labs(x = 'Districts, Sorted by Minority Percent', y = 'Minority Percent')

  }


  stop('No available type specified.')
}


globalVariables(c('nloop', 'blk_pct', 'group_id', 'tot_blk_pct', 'first_blk', 'second_blk', 'district'))
