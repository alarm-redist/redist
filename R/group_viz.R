#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2020/11/12
# Purpose: Redistricting plan group population visualization
####################################################

#' Visualize Plans and Their Treatment of Groups
#'
#' Launches a \code{shiny} app which allows a user to explore a set of plans.
#'
#' The app plots plans as points on a two-dimensional plane, based on
#' multidimensional scaling applied to a variation of information distance
#' matrix (see \code{\link{redist.distances}}).  Users can click any point on the
#' plane to see the corresponding plan.  Points are colored by their share of
#' the group population (where the group can be minority voters, Democratic
#' voters, or some other group).  The app also clusters plans according to their
#' similarity, and plots a representative plan for each cluster, allowing users
#' to quickly get a sense of the high-level structure of the space of possible
#' redistricting plans.
#'
#' @param prep The output of the \code{\link{redist.prep_group_viz}} function,
#'   which handles all the computations and preparatory work.
#' @param browser if \code{TRUE} (the default), launch the Shiny app in the
#'   default browser. Otherwise, just start the Shiny server.
#'
#' @examples \dontrun{
#' data("fl25")
#' data(algdat.p10)
#' # Sample
#' sampled_plans = redist.smc(algdat.p10$adjlist, fl25$TotPop,
#'                            nsims=200, ndists=3, popcons=0.1)
#'
#' # Prep
#' viz_prep = redist.prep_group_viz(sampled_plans, fl25, fl25$VAP, fl25$BlackVAP)
#'
#' # Run the app!
#' redist.group_viz(viz_prep)
#' }
#'
#' @concept visualize
#' @export
redist.group_viz = function(prep, browser=TRUE) {
    if (!requireNamespace("shiny", quietly=T)
            | !requireNamespace("rmarkdown", quietly=T)
            | !requireNamespace("cluster", quietly=T)
            | !requireNamespace("tidyr", quietly=T)
            | !requireNamespace("ggiraph", quietly=T)
            | !requireNamespace("leaflet", quietly=T)
            | !requireNamespace("rmapshaper", quietly=T)
            | !requireNamespace("gridExtra", quietly=T)) {
        stop('The "shiny", "rmarkdown", "cluster", "tidyr", "ggiraph", "leaflet", and "gridExtra"',
        "packages are required for this function to work.", call.=F)
    }
    if (!methods::is(prep, "redist_prep_viz"))
        stop("`prep` must be of class `redist_prep_viz`. ",
             "You must call `redist.prep_group_viz` first.", call.=F)

    app = list(ui = viz_ui(), server = do.call(viz_server, prep))

    shiny::runApp(app, launch.browser=browser, port=1965)
}

#' Prepare to Visualize Plans
#'
#' This helper function performs all necessary computations and preparation to
#' run the Shiny app in \code{\link{redist.group_viz}}
#'
#' @param result a \code{redist} object produced by \code{\link{redist.smc}}
#'   or \code{\link{redist.mcmc}}, or a matrix of district assignments.
#' @param geometry  an \code{sf} data frame containing shapefile information for
#'   the current map.
#' @param full_pop a numeric vector with the population of the group for every precinct.
#' @param group_pop a numeric vector with the population of every precinct.
#' @param current a numeric vector containing the current plan, which will be plotted separately.
#' @param comp_meas the compactness measure desired.  See \code{\link{redist.compactness}} for details.
#' @param renumber whether to renumber districts in ascending order by their group share, so
#'   district 1 has the smallest fraction of the group, district 2 has the second-smallest,
#'   and so on.
#' @param dist_by the way to compute the variation of information distance, by \code{"people"},
#'   the recommended default, by \code{"group"} or \code{"minority"}, or by
#'   \code{"land"}, which treats each precinct equally.
#' @param clusters the number of clusters to group the plans into. If left unspecified,
#'   this will be determined automatically using the silhouette method.
#' @param simp_geom if \code{TRUE}, the geometry will be slightly simplified to
#'   allow for faster plotting.
#'
#' @return A list of class \code{redist_prep_viz} containing the preprocessed
#'   inputs and calculations.
#'
#' @examples \dontrun{
#' data("fl25")
#' data(algdat.p10)
#' # Sample
#' sampled_basic = redist.smc(algdat.p10$adjlist, fl25$TotPop,
#'                            nsims=200, ndists=3, popcons=0.1)
#'
#' # Prep
#' viz_prep = redist.prep_group_viz(sampled_plans, fl25, fl25$VAP, fl25$BlackVAP)
#'
#' # Run the app!
#' redist.group_viz(viz_prep)
#' }
#'
#' @concept visualize
#' @export
redist.prep_group_viz = function(result, geometry, full_pop, group_pop, current=NULL,
                                 comp_meas="EdgesRemoved", renumber=TRUE, dist_by="people",
                                 clusters=NULL, simp_geom=TRUE) {
    # Type checks
    if (!methods::is(result, "redist"))
        stop("`result` must be of class `redist`.", call.=F)
    if (min(result$cdvec[,1]) == 0) # 1-index
        result$cdvec = 1 + result$cdvec
    if (!any(class(full_pop) %in% c('numeric', 'integer')))
        stop('Please provide "full_pop" as a numeric vector.')
    if (!any(class(group_pop) %in% c('numeric', 'integer')))
        stop('Please provide "group_pop" as a numeric vector.')

    N = max(result$cdvec[,1])
    full_pop = as.numeric(full_pop)
    group_pop = as.numeric(group_pop)

    # set up geometry
    cat("Reprojecting geometry\n")
    if (is.na(sf::st_crs(geometry)))
        sf::st_crs(geometry) = 4326
    geometry = sf::st_transform(geometry, '+proj=longlat +datum=WGS84')
    if (simp_geom) {
        cat("Simplifying geometry\n")
        geometry = rmapshaper::ms_simplify(geometry, keep_shapes=T)
    }
    geometry$.id = as.character(1:nrow(geometry))

    centroids = geometry %>%
        sf::st_centroid() %>%
        suppressWarnings() %>%
        sf::st_coordinates() %>%
        dplyr::as_tibble()

    no_current = is.null(current)
    if (no_current) current = rep(N, nrow(result$cdvec))
    m_plot = cbind(current, result$cdvec)

    # basic vitals
    tot_edges = sum(sapply(result$aList, length))/2
    maxdev = c(redist.parity(result$cdvec[,1:2,drop=F], full_pop)[1], result$maxdev)
    comp = redist.compactness(geometry, m_plot, measure=comp_meas,
                              population=full_pop, adjacency=result$aList) %>%
        dplyr::group_by(.data$nloop) %>%
        dplyr::summarize_all(mean) %>%
        dplyr::pull()


    # different units to use
    dist_by = dplyr::case_when(
        dist_by == "group" ~ group_pop,
        dist_by == "minority" ~ group_pop,
        dist_by == "land" ~ scale(st_area(geometry), center=F),
        T ~ full_pop
    )
    cat("Calculating pairwise distances\n")
    dists = redist.distances(m_plot, "info", pop=dist_by)$VI

    together = which(dists <= 0.001, arr.ind=T) %>%
        dplyr::as_tibble() %>%
        dplyr::filter(row > col, !col %in% row)
    size_t = table(together$col)
    sizes = rep(1, nrow(dists))
    sizes[as.integer(names(size_t))] = size_t + 1

    cat("Performing multidimensional scaling\n")
    mds = stats::cmdscale(dists)

    cat("Computing group percentages\n")
    pct_min = redist.group.percent(m_plot, group_pop, full_pop)
    if (renumber)
        pct_min = apply(pct_min, 2, sort, na.last=T)
    if (no_current)
        pct_min[,1] = rowMeans(pct_min[,-1])
    rownames(pct_min) = 1:nrow(pct_min)

    cat("Clustering plans\n")
    cl_m = cbind(mds, t(pct_min))
    if (is.null(clusters)) {
        try_k = 12
        silh = rep(0, try_k)
        for (k in 3:try_k) {
            #cl = pam(cl_m, k, diss=F, pamonce=T)
            cl = cluster::pam(dists, k, diss=T)
            silh[k] = cl$silinfo$avg.width
        }
        clusters = which.max(silh)
    }
    cat("  (Using", clusters, "clusters)\n")
    cl = cluster::pam(dists, clusters, diss=T)

    cat("Finding typical plans for each cluster\n")
    N = nrow(pct_min)
    typical = 1:clusters %>%
        lapply(function(i) m_plot[,cl$medoids[i]])
    #typical = 1:k %>%
    #    lapply(function(i) prec_cooccur(m_plot, which(cl$clustering == i))) %>%
    #    lapply(function(cm) stats::kmeans(cm, N, iter.max=20, nstart=3)$cluster)

    out = list(m=m_plot, dev=maxdev, comp=comp, geom=geometry,
               mds=mds, sizes=sizes, hide=c(1, together$row),
               pct_min=pct_min, pop=full_pop, grp_pop=group_pop, centroids=centroids,
               cluster=cl$clustering, centers=mds[cl$medoids,], typical=typical, k=clusters)
    class(out) = c("redist_prep_viz", "list")
    out
}



viz_ui = function() {
    shiny::fluidPage(
        # add in methods from https://github.com/rstudio/leaflet/pull/598
        shiny::tags$head(
            shiny::tags$script(shiny::HTML(viz_script)),
            shiny::tags$style(type="text/css", shiny::HTML(
            ".container-fluid {
                max-width: 1200px;
            }
            h2 {
                font-weight: bold;
            }
            .container-fluid > h2:first-child {
                font-size: 2.8em;
                border-bottom: 4px solid #777;
                padding-bottom: 4px;
            }
            #mds svg {
                cursor: crosshair;
            }
            #mindistr .shiny-options-group {
                column-count: 2;
            }
            .leaflet-pane.leaflet-overlay-pane {
                mix-blend-mode: multiply;
            }
            .girafe_container_std g > text {
                pointer-events: none;
            }
            "))
        ),

        shiny::titlePanel("Redistricting Ensemble Explorer"),

        shiny::h2("Plan Vitals"),
        shiny::p("The first set of plots show how the set of sampled redistricting
                 plans measures up, according to basic metrics, and how the
                 current plan, if any, compares."),
        shiny::fluidRow(
            shiny::column(6, shiny::verticalLayout(
                shiny::h3("Population deviation"),
                shiny::plotOutput("popdev", width="100%", height=300)
            )),
            shiny::column(6, shiny::verticalLayout(
                shiny::h3("Compactness"),
                shiny::plotOutput("comp", width="100%", height=300)
            ))
        ),
        shiny::br(), shiny::br(),

        shiny::h2("Types of Plans"),
        shiny::p(
            "The scatterplot below lays out plans according to the shared information
                distance between them.",
            shiny::tags$ul(
                shiny::tags$li("The current plan, if there is one, is marked by a red square."),
                shiny::tags$li(shiny::strong("Hover over a plan to learn more,
                    or click a plan to see it on the map.")),
                shiny::tags$li("Click again to show the fraction of group voters in each precinct."),
                shiny::tags$li("Districts are shaded in ascending order by group share
                    (dark blue is the smallest group fraction, bright green is the most)")
            )
        ),

        shiny::fluidRow(
            shiny::column(6, shiny::verticalLayout(
                shiny::h3("All plans"),
                ggiraph::girafeOutput("mds", width="100%", height="100%"),
                shiny::uiOutput("mindistr")
            )),
            shiny::column(6, shiny::verticalLayout(
                shiny::h3(shiny::textOutput("map_title")),
                leaflet::leafletOutput("map", width="100%", height=640))
            )
        ),

        shiny::h3(shiny::HTML("&lsquo;Typical&rsquo; plans for each cluster")),
        shiny::p("These maps show the plan which is at the center of each cluster,
                 and can be considered representative."),
        shiny::plotOutput("typical", width="100%", height="640px"),
        shiny::br(), shiny::br(),


        shiny::h2("Seats and Votes"),
        shiny::fluidRow(
            shiny::column(6, shiny::verticalLayout(
                shiny::h3("Group-controlled seats"),
                shiny::plotOutput("numseats", width="100%", height=400),
                shiny::sliderInput("cutoff", "Group control threshold", width="100%",
                                   min=0.4, max=0.65, value=0.5, step=0.005, round=T)
            )),
            shiny::column(6, shiny::verticalLayout(
                shiny::h3("Fraction of group population by district"),
                shiny::plotOutput("distrpct", width="100%", height=480)
            ))
        ),
        shiny::br(), shiny::br(), shiny::br(), shiny::br()#,

        #downloadButton("report", "Generate report")
    )
}

viz_server = function(m, dev, comp, geom, mds, sizes, hide, pct_min, pop,
                      grp_pop, centroids, cluster, centers, typical, k) {

    N = max(m)
    colnames(centers) = c("x", "y")
    centers = dplyr::as_tibble(centers) %>%
      dplyr::mutate(id = 1:n())
    show_current = diff(range(m[,1])) > 0

    function(input, output, session) {
        make_hist = function(x, show_current) {
            nbins = floor(ncol(m)/20)
            p = ggplot2::ggplot(NULL, ggplot2::aes(x=x[-1])) +
              ggplot2::geom_histogram(bins=nbins) +
              ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0, 0.05))) +
              ggplot2::theme_bw()

            if (show_current) {
                ymax = max(graphics::hist(x, plot=F)$counts)
                p = p +
                  ggplot2::geom_vline(xintercept=x[1], color="red", size=1) +
                  ggplot2::annotate("text", y=ymax, x=x[1] + 0.01*diff(range(x)),
                                    label="Current plan", hjust="left")
            }

            p
        }

        # POPULATION DEVIATION PLOT
        plot_popdev =  make_hist(dev, show_current) +
            ggplot2::labs(x="Maximum population deviation", y="Number of plans") +
            ggplot2::scale_x_continuous(labels=scales::percent)
        output$popdev = shiny::renderPlot(plot_popdev)

        # COMPACTNESS PLOT
        plot_comp = make_hist(comp, show_current) +
            ggplot2::labs(x="Average compactness measure", y="Number of plans")
        output$comp = shiny::renderPlot(plot_comp)

        # SEATS PLOT
        output$numseats = shiny::renderPlot({
            num_seats = apply(pct_min, 2, function(x) sum(x > input$cutoff))
            breaks = seq(min(num_seats)-1/4, max(num_seats)+1/4, 1/2)
            p = ggplot2::ggplot(NULL, ggplot2::aes(x=num_seats[-1])) +
              ggplot2::geom_histogram(breaks=breaks) +
              ggplot2::labs(x="Seats controlled by group", y="Number of plans") +
              ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0, 0.05))) +
              ggplot2::theme_bw()

            if (show_current) {
                p = p +
                  ggplot2::geom_vline(xintercept=num_seats[1], color="red", size=1) +
                  ggplot2::annotate("text", y=max(table(num_seats)),
                                    x=num_seats[1]+0.01*diff(range(num_seats)),
                                    hjust="left", label="Current plan")
            }

            p
        })

        # DISTRICT SHARE BOXPLOTS
        output$distrpct = shiny::renderPlot({
            pct_d = dplyr::as_tibble(pct_min, .name_repair=function(x) as.character(1:ncol(pct_min))) %>%
                dplyr::mutate(district=factor(as.character(1:N), as.character(1:N))) %>%
                tidyr::pivot_longer(-tidyselect::starts_with("district"), names_to="plan", values_to="pct")
            p = ggplot2::ggplot(pct_d[pct_d$plan != 1,],
                                ggplot2::aes(.data$district, .data$pct)) +
                ggplot2::geom_boxplot(fill="#bbbbbb", outlier.size=0.2) +
                ggplot2::scale_y_continuous(labels=scales::percent) +
                ggplot2::labs(x="District", y="Group percentage") +
                ggplot2::theme_bw()

            if (show_current) {
                p = p + ggplot2::geom_point(data=pct_d[pct_d$plan == 1,],
                                            color="red", size=3.5, shape=18) +
                  ggplot2::annotate("text", min(pct_min), x=N, label="Current plan in red")
            }

            p
        })

        output$mindistr = shiny::renderUI({
            shiny::radioButtons("mindistr", label="Color points by", width="100%",
                                choiceValues=c(1:N, -1), selected=N,
                                choiceNames = c(str_c("Group pct. in district ", 1:N), "Cluster"))
        })

        # MDS PLOT
        output$mds = ggiraph::renderGirafe({
            mindistr = if (is.null(input$mindistr)) N else as.integer(input$mindistr)
            plot_col = if (mindistr > 0) str_c("min_", mindistr) else "cluster"
            labels = apply(rbind(1:ncol(m), cluster, pct_min), 2, function(x) {
                pct_labels = str_glue("District {1:N}: {scales::percent(x[-1:-2], 0.1)}")
                map_label = if (x[1] == 1) "Current plan\n" else str_c("Plan #", x[1]-1, "\n")
                cl_label = str_c("Cluster #", x[2], "\n")
                str_c(map_label, cl_label, paste(pct_labels, collapse="\n"))
            })

            plot_d = cbind(mds, cluster, t(pct_min)) %>%
                dplyr::as_tibble(.name_repair=function(x) {
                    c("x", "y",  "cluster", str_c("min_", x[-1:-3]))
                }) %>%
                dplyr::select(x, y, col=!!plot_col) %>%
                dplyr::mutate(id = 1:n(), label = .data$labels, size = .data$sizes)

            if (mindistr > 0)
                legend = str_glue("Group pct. in\ndistrict {mindistr}")
            else
                legend = "Cluster"

            p = ggplot2::ggplot(plot_d[-hide,],
                                ggplot2::aes(x, y, color=col, size=.data$size,
                                             data_id=id, tooltip=.data$label)) +
                ggiraph::geom_point_interactive(alpha=0.85) +
                ggplot2::scale_color_viridis_c(labels=if (mindistr > 0) scales::percent else scales::comma) +
                {if (show_current)
                    ggiraph::geom_point_interactive(data=plot_d[1,], color="#dd0000",
                                                    shape=15, size=2.5)} +
                {if (length(hide) > 1)
                    ggplot2::scale_size_continuous(range=c(1.2, 5.5))
                 else
                    ggplot2::scale_size_continuous(range=c(1.5, 1.5))} +
                ggplot2::labs(color=legend, x="MDS 1", y="MDS 2") +
                ggplot2::theme_bw(base_size=11) +
                ggplot2::guides(color=ggplot2::guide_colourbar(barwidth=20), size=F) +
                ggplot2::theme(legend.position="bottom")

            if (mindistr < 0)
                p = p + ggplot2::geom_text(ggplot2::aes(x, y, label=id), data=centers, inherit.aes=F, size=5)

            tooltip_css = "
            background: #000000c0;
            padding: 7px;
            border-radius: 2px;
            color: white;
            backdrop-filter: blur(4px);
            font-family: Arial, sans-serif"
            select_css = "stroke: black; stroke-width:2px; z-index: 100;"

            ggiraph::girafe(ggobj=p) %>%
                ggiraph::girafe_options(ggiraph::opts_hover(""),
                                        ggiraph::opts_zoom(max=3),
                                        ggiraph::opts_tooltip(tooltip_css),
                                        ggiraph::opts_toolbar(saveaspng=F),
                                        ggiraph::opts_selection(select_css, type="single"))
        })

        selected_map = shiny::reactive({
            id = as.numeric(input$mds_selected)
            if (is.null(id) || length(id) == 0) -1 else id
        })

        # MAP
        output$map = leaflet::renderLeaflet({
            leaflet::leaflet(data=geom) %>%
                leaflet::addProviderTiles(leaflet::providers$Stamen.Toner) %>%
                leaflet::addPolygons(layerId=~.id, weight=1, fillOpacity=0.85, opacity=0.6)
        })

        output$map_title = shiny::renderText({
            id = selected_map()
            if (id > 1) str_c("Plan #", id-1)
            else if (id==1) "Current plan"
            else "Fraction of group population"
        })

        pal = scales::col_numeric("viridis", range(0.5 + 0:N))
        mean_pct = sum(grp_pop)/sum(pop)
        pct = scales::rescale((grp_pop+mean_pct)/(pop+1))
        col_min = scales::gradient_n_pal(c("#00000000", "#440154ff"))(pct)
        label_opts = leaflet::labelOptions(noHide=T, textOnly=T,
                                           style=list("font-size"="13pt",
                                                      "font-weight"="bold",
                                                      "color"="#ffffffaa"))
        # MAP UPDATE
        shiny::observe({
            id = selected_map()
            if (id >= 0) {
                distr_pops = tapply(pop, m[,id], sum)
                centers_x = tapply(centroids$X*pop, m[,id], sum) / distr_pops
                centers_y = tapply(centroids$Y*pop, m[,id], sum) / distr_pops
                pcts = tapply(grp_pop, m[,id], sum) / distr_pops
                cols = pal(rank(pcts))[m[,id]]

                leaflet::leafletProxy("map", data=geom) %>%
                    viz_setShapeStyle(layerId=~.id, fillColor=cols, color=cols) %>%
                    leaflet::clearGroup("distr_numbers") %>%
                    leaflet::addLabelOnlyMarkers(centers_x, centers_y, label=1:N,
                                                 group="distr_numbers", labelOptions=label_opts)
            } else {
                cols = col_min

                leaflet::leafletProxy("map", data=geom) %>%
                    viz_setShapeStyle(layerId=~.id, fillColor=cols, color=cols) %>%
                    leaflet::clearGroup("distr_numbers")
            }
        })

        # CLUSTER MEDOIDS
        output$typical = shiny::renderPlot({
            grobs = lapply(1:length(typical), function(i) {
                plan = typical[[i]]
                pcts = tapply(grp_pop, plan, sum) / tapply(pop, plan, sum)
                cols = pal(rank(pcts))[plan]

                dplyr::mutate(geom, d = plan) %>%
                ggplot2::ggplot(ggplot2::aes(fill=cols)) +
                    ggplot2::geom_sf(size=0, color="transparent") +
                    ggplot2::scale_fill_identity() +
                    ggplot2::guides(fill=F) +
                    ggplot2::labs(title=str_c("Cluster #", i)) +
                    ggplot2::theme_void()
            })
            cols = round(sqrt(k/0.7))
            gridExtra::marrangeGrob(grobs, ncol=cols, nrow=ceiling(k/cols), top=NULL)
        })


        ### REPORT
        output$report = shiny::downloadHandler(
            filename="redistricting_report.pdf",
            content = function(file) {
                # Copy the report file to a temporary directory before processing it, in
                # case we don't have write permissions to the current working directory
                tempReport = file.path(tempdir(), "viz_report.Rmd")
                file.copy("R/viz_report.Rmd", tempReport, overwrite=T)

                params = list(plot_popdev=plot_popdev)
                rmarkdown::render(tempReport, output_file=file, params=params,
                                  envir=globalenv())
                                  #envir=new.env(parent=globalenv()))
            })
    }
}

# To update map colors quickly
viz_setShapeStyle = function(map, data=leaflet::getMapData(map), layerId, stroke=NULL,
                         color=NULL, weight=NULL, opacity=NULL, fill=NULL,
                         fillColor=NULL, fillOpacity=NULL, dashArray=NULL,
                         smoothFactor=NULL, noClip=NULL, options=NULL) {
    options = c(list(layerI=layerId), options,
                 leaflet::filterNULL(list(stroke=stroke, color=color, weight=weight,
                                 opacity=opacity, fill=fill, fillColor=fillColor,
                                 fillOpacity=fillOpacity, dashArray=dashArray,
                                 smoothFactor=smoothFactor, noClip=noClip)))
    # evaluate all options
    options = leaflet::evalFormula(options, data=data)
    # make them the same length (by building a data.frame)
    options = do.call(data.frame, c(options, list(stringsAsFactors=F)))

    layerId = options[[1]]
    style = options[-1]

    leaflet::invokeMethod(map, data, "setStyle", "shape", layerId, style)
}

# To update map colors quickly
viz_script = "
window.LeafletWidget.methods.setStyle = function(category, layerId, style) {
  var map = this;
  if (!layerId){
    return;
  } else if (!(typeof(layerId) === 'object' && layerId.length)) { // in case a single layerid is given
    layerId = [layerId];
  }

  //convert columnstore to row store
  style = HTMLWidgets.dataframeToD3(style);
  //console.log(style);

  layerId.forEach(function(d,i) {
    var layer = map.layerManager.getLayer(category, d);
    if (layer){ // or should this raise an error?
      layer.setStyle(style[i]);
    }
  });
};
"

utils::globalVariables(c("x", "y")) # for R CHECK
