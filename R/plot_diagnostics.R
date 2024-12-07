#' Diagnostic plotting functionality for MCMC redistricting.
#'
#' \code{gredist.diagplot} generates several common MCMC diagnostic plots.
#'
#' @usage gredist.diagplot(sumstat,
#' plot = c("trace", "autocorr", "densplot", "mean", "gelmanrubin"),
#' logit = FALSE, savename = NULL)
#'
#' @param sumstat A vector, list, \code{mcmc} or \code{mcmc.list} object
#' containing a summary statistic of choice.
#' @param plot The type of diagnostic plot to generate: one of "trace",
#' "autocorr", "densplot", "mean", "gelmanrubin". If \code{plot = "gelmanrubin"},
#' the input \code{sumstat} must be of class \code{mcmc.list} or \code{list}.
#' @param logit Flag for whether to apply the logistic transformation for the
#' summary statistic. The default is \code{FALSE}.
#' @param savename Filename to save the plot. Default is \code{NULL}.
#'
#' @details This function allows users to generate several standard diagnostic
#' plots from the MCMC literature, as implemented by Plummer et. al (2006).
#' Diagnostic plots implemented include trace plots, autocorrelation plots,
#' density plots, running means, and Gelman-Rubin convergence diagnostics
#' (Gelman & Rubin 1992).
#'
#' @return Returns a plot of file type \code{.pdf}.
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte
#' Carlo." Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/gredist.pdf}.
#'
#' Gelman, Andrew and Donald Rubin. (1992) "Inference from iterative simulations
#' using multiple sequences (with discussion)." Statistical Science.
#'
#' Plummer, Martin, Nicky Best, Kate Cowles and Karen Vines. (2006) "CODA:
#' Convergence Diagnosis and Output Analysis for MCMC." R News.
#'
#' @examples
#' \donttest{
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#'
#' ## Get an initial partition
#' init_plan <- fl25_enum$plans[, 5118]
#' fl25$init_plan <- init_plan
#'
#' ## 25 precinct, three districts - no pop constraint ##
#' fl_map <- redist_map(fl25, existing_plan = 'init_plan', adj = fl25_adj)
#' alg_253 <- redist_flip(fl_map, nsims = 10000)
#'
#' ## Get Republican Dissimilarity Index from simulations
#' rep_dmi_253 <- redistmetrics::seg_dissim(alg_253, fl25, mccain, pop) |>
#'     redistmetrics::by_plan(ndists = 3)
#'
#' ## Generate diagnostic plots
#' gredist.diagplot(rep_dmi_253, plot = "trace")
#' gredist.diagplot(rep_dmi_253, plot = "autocorr")
#' gredist.diagplot(rep_dmi_253, plot = "densplot")
#' gredist.diagplot(rep_dmi_253, plot = "mean")
#'
#' ## Gelman Rubin needs two chains, so we run a second
#' alg_253_2 <- redist_flip(fl_map, nsims = 10000)
#'
#' rep_dmi_253_2 <- redistmetrics::seg_dissim(alg_253_2, fl25, mccain, pop) |>
#'     redistmetrics::by_plan(ndists = 3)
#'
#' ## Make a list out of the objects:
#' rep_dmi_253_list <- list(rep_dmi_253, rep_dmi_253_2)
#'
#' ## Generate Gelman Rubin diagnostic plot
#' gredist.diagplot(sumstat = rep_dmi_253_list, plot = "gelmanrubin")
#'
#' }
#' @concept plot
#' @export
gredist.diagplot <- function(sumstat, plot = c("trace", "autocorr", "densplot",
                                "mean", "gelmanrubin"),
                            logit = FALSE, savename = NULL) {
    if (!requireNamespace("coda", quietly = TRUE))
        cli_abort(c("{.fn gredist.diagplot} requires the {.pkg coda} package.",
            ">" = 'Install it with {.code install.packages("coda")}'))

    ##############
    ## Warnings ##
    ##############
    if (!inherits(sumstat, c("integer", "numeric", "list", "mcmc", "mcmc.list"))) {
        cli_abort("{.arg sumstat} should be either a numeric vector, list, or {.cls mcmc} object.")
    }
    if (!(plot %in% c("trace", "autocorr", "densplot",
        "mean", "gelmanrubin"))) {
        cli_abort("Sorry. We don't currently support the {.value {plot}} diagnostic.")
    }
    if (plot == "gelmanrubin" & !inherits(sumstat, c("list", "mcmc.list"))) {
        cli_abort("If generating a Gelman-Rubin plot, please provide an object of class list or mcmc.list")
    }

    ########################
    ## Create mcmc object ##
    ########################
    if (is.numeric(sumstat)) {
        segout <- coda::mcmc(sumstat)
    } else if (is.list(sumstat)) {
        for (i in 1:length(sumstat)) {
            sumstat[[i]] <- coda::mcmc(sumstat[[i]])
        }
        segout <- coda::mcmc.list(sumstat)
    } else if (inherits(sumstat, c("mcmc", "mcmc.list"))) {
        segout <- sumstat
    }

    ## Logit transform
    if (logit) {
        if (inherits(segout, "mcmc")) {
            segout <- log(segout/(1 - segout))
        } else if (inherits(segout, "mcmc.list")) {
            for (i in 1:length(segout)) {
                segout[[i]] <- log(segout[[i]]/(1 - segout[[i]]))
            }
        }
    }

    ##################
    ## Create plots ##
    ##################
    if (plot == "trace") {
        if (!is.null(savename)) {
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        coda::traceplot(segout)
        if (!is.null(savename)) {
            dev.off()
        }
    }
    if (plot == "autocorr") {
        if (!is.null(savename)) {
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        coda::autocorr.plot(segout, lag.max = 50)
        if (!is.null(savename)) {
            dev.off()
        }
    }
    if (plot == "densplot") {
        if (!is.null(savename)) {
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        coda::densplot(segout)
        if (!is.null(savename)) {
            dev.off()
        }
    }
    if (plot == "mean") {
        if (!is.null(savename)) {
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        coda::cumuplot(segout, probs = .5, type = "l", lty = 1)
        if (!is.null(savename)) {
            dev.off()
        }
    }
    if (plot == "gelmanrubin" && inherits(segout, "mcmc.list")) {
        if (!is.null(savename)) {
            pdf(file = paste(savename, ".pdf", sep = ""))
        }
        coda::gelman.plot(segout, transform = FALSE)
        if (!is.null(savename)) {
            dev.off()
        }
    }

}
