\name{redist.diagplot}
\alias{redist.diagplot}
\title{
Diagnostic plotting functionality for MCMC redistricting.
}
\description{
\code{redist.diagplot} generates several common MCMC diagnostic plots. 
}
\usage{
redist.diagplot(sumstat,
                plot = c("trace", "autocorr", "densplot",
                    "mean", "gelmanrubin"),
                logit = FALSE, 
                savename = NULL)
}
\arguments{
\item{sumstat}{A vector, list, \code{mcmc} or \code{mcmc.list} object containing a summary statistic of choice.}
\item{plot}{The type of diagnostic plot to generate: one of "trace", "autocorr", "densplot", "mean",
"gelmanrubin". If \code{plot = "gelmanrubin"}, the input \code{sumstat} must be of class \code{mcmc.list} or \code{list}.}
\item{logit}{Flag for whether to apply the logistic transformation for the summary statistic. The default is \code{FALSE}.}
\item{savename}{Filename to save the plot. Default is \code{NULL}.}
}

\details{
This function allows users to generate several standard diagnostic plots from the MCMC literature, as implemented by Plummer et. al (2006). Diagnostic plots implemented include trace plots, autocorrelation plots, density plots, running means, and Gelman-Rubin convergence diagnostics (Gelman & Rubin 1992).
}
\value{
Returns a plot of file type \code{.pdf}.
}
\references{
Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr. (2015) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo." Working 
  Paper. Available at \url{http://imai.princeton.edu/research/files/redist.pdf}.

Gelman, Andrew and Donald Rubin. (1992) "Inference from iterative simulations using multiple sequences
(with discussion)." Statistical Science.

Plummer, Martin, Nicky Best, Kate Cowles and Karen Vines. (2006) "CODA: Convergence Diagnosis and Output Analysis for MCMC." R News.
}

\examples{
\dontrun{
data(algdat.pfull)

## Get an initial partition
set.seed(1)
initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]

## 25 precinct, three districts - no pop constraint ##
alg_253 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                       popvec = algdat.pfull$precinct.data$pop,
                       initcds = initcds,
                       nsims = 10000)

## Get Republican Dissimilarity Index from simulations
rep_dmi_253 <- redist.segcalc(alg_253,
                              algdat.pfull$precinct.data$repvote,
                              algdat.pfull$precinct.data$pop)
                              
## Generate diagnostic plots
redist.diagplot(rep_dmi_253, plot = "trace")
redist.diagplot(rep_dmi_253, plot = "autocorr")
redist.diagplot(rep_dmi_253, plot = "densplot")
redist.diagplot(rep_dmi_253, plot = "mean")

}
}

\keyword{diagplot}
