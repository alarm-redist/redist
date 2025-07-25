% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redist_ms_parallel.R
\name{redist_mergesplit_parallel}
\alias{redist_mergesplit_parallel}
\title{Parallel Merge-Split/Recombination MCMC Redistricting Sampler}
\usage{
redist_mergesplit_parallel(
  map,
  nsims,
  chains = 1,
  warmup = if (is.null(init_plan)) 10 else max(100, nsims\%/\%5),
  thin = 1L,
  init_plan = NULL,
  counties = NULL,
  compactness = 1,
  constraints = list(),
  constraint_fn = function(m) rep(0, ncol(m)),
  adapt_k_thresh = 0.99,
  k = NULL,
  ncores = NULL,
  cl_type = "PSOCK",
  return_all = TRUE,
  init_name = NULL,
  silly_adj_fix = FALSE,
  verbose = FALSE,
  silent = FALSE
)
}
\arguments{
\item{map}{A \code{\link{redist_map}} object.}

\item{nsims}{The number of samples to draw, including warmup.}

\item{chains}{the number of parallel chains to run. Each chain will have
\code{nsims} draws. If \code{init_plan} is sampled, each chain will be initialized
with its own sampled plan.}

\item{warmup}{The number of warmup samples to discard. Recommended to be at
least the first 20\% of samples, and in any case no less than around 100
samples, unless initializing from a random plan.}

\item{thin}{Save every \code{thin}-th sample. Defaults to no thinning (1).}

\item{init_plan}{The initial state of the map, provided as a single vector
to be shared across all chains, or a matrix with \code{chains} columns.
If not provided, will default to the reference map of the map object, or if
none exists, will sample a random initial state using redist_smc. You can
also request a random initial state for each chain by setting
init_plan="sample".}

\item{counties}{A vector containing county (or other administrative or
geographic unit) labels for each unit, which may be integers ranging from 1
to the number of counties, or a factor or character vector.  If provided,
the algorithm will generate maps tend to follow county lines. There is no
strength parameter associated with this constraint. To adjust the number of
county splits further, or to constrain a second type of administrative
split, consider using \code{add_constr_splits()}, \code{add_constr_multisplits()},
and \code{add_constr_total_splits()}.}

\item{compactness}{Controls the compactness of the generated districts, with
higher values preferring more compact districts. Must be nonnegative. See the
'Details' section for more information, and computational considerations.}

\item{constraints}{A list containing information on constraints to implement.
See the 'Details' section for more information.}

\item{constraint_fn}{A function which takes in a matrix where each column is
a redistricting plan and outputs a vector of log-weights, which will be
added the the final weights.}

\item{adapt_k_thresh}{The threshold value used in the heuristic to select a
value \code{k_i} for each splitting iteration. Set to 0.9999 or 1 if
the algorithm does not appear to be sampling from the target distribution.
Must be between 0 and 1.}

\item{k}{The number of edges to consider cutting after drawing a spanning
tree. Should be selected automatically in nearly all cases.}

\item{ncores}{the number of parallel processes to run. Defaults to the
maximum available.}

\item{cl_type}{the cluster type (see \code{\link[=makeCluster]{makeCluster()}}). Safest is \code{"PSOCK"},
but \code{"FORK"} may be appropriate in some settings.}

\item{return_all}{if \code{TRUE} return all sampled plans; otherwise, just return
the final plan from each chain.}

\item{init_name}{a name for the initial plan, or \code{FALSE} to not include
the initial plan in the output.  Defaults to the column name of the
existing plan, or "\code{<init>}" if the initial plan is sampled.}

\item{silly_adj_fix}{Heuristic for fixing weird inputs.}

\item{verbose}{Whether to print out intermediate information while sampling.
Recommended.}

\item{silent}{Whether to suppress all diagnostic information.}
}
\value{
A \code{\link{redist_plans}} object with all of the simulated plans, and an
additional \code{chain} column indicating the chain the plan was drawn from.
}
\description{
\code{redist_mergesplit_parallel()} runs \code{\link[=redist_mergesplit]{redist_mergesplit()}} on several
chains in parallel.
}
\details{
This function draws samples from a specific target measure, controlled by the
\code{map}, \code{compactness}, and \code{constraints} parameters.

Key to ensuring good performance is monitoring the acceptance rate, which
is reported at the sample level in the output.
Users should also check diagnostics of the sample by running
\code{summary.redist_plans()}.

Higher values of \code{compactness} sample more compact districts;
setting this parameter to 1 is computationally efficient and generates nicely
compact districts.
}
\examples{
\dontrun{
data(fl25)
fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
sampled <- redist_mergesplit_parallel(fl_map, nsims = 100, chains = 100)
}

}
\references{
Carter, D., Herschlag, G., Hunter, Z., and Mattingly, J. (2019). A
merge-split proposal for reversible Monte Carlo Markov chain sampling of
redistricting plans. arXiv preprint arXiv:1911.01503.

McCartan, C., & Imai, K. (2023). Sequential Monte Carlo for Sampling
Balanced and Compact Redistricting Plans. \emph{Annals of Applied Statistics} 17(4).
Available at \doi{10.1214/23-AOAS1763}.

DeFord, D., Duchin, M., and Solomon, J. (2019). Recombination: A family of
Markov chains for redistricting. arXiv preprint arXiv:1911.05725.
}
\concept{simulate}
