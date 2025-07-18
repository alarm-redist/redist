% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redist-package.R
\docType{package}
\name{redist-package}
\alias{redist-package}
\alias{redist}
\title{redist: Simulation Methods for Legislative Redistricting}
\description{
\if{html}{\figure{logo.png}{options: style='float: right' alt='logo' width='120'}}

Enables researchers to sample redistricting plans from a pre-specified target distribution using Sequential Monte Carlo and Markov Chain Monte Carlo algorithms. The package allows for the implementation of various constraints in the redistricting process such as geographic compactness and population parity requirements. Tools for analysis such as computation of various summary statistics and plotting functionality are also included. The package implements the SMC algorithm of McCartan and Imai (2023) \doi{10.1214/23-AOAS1763}, the enumeration algorithm of Fifield, Imai, Kawahara, and Kenny (2020) \doi{10.1080/2330443X.2020.1791773}, the Flip MCMC algorithm of Fifield, Higgins, Imai and Tarr (2020) \doi{10.1080/10618600.2020.1739532}, the Merge-split/Recombination algorithms of Carter et al. (2019) \doi{10.48550/arXiv.1911.01503} and DeFord et al. (2021) \doi{10.1162/99608f92.eb30390f}, and the Short-burst optimization algorithm of Cannon et al. (2020) \doi{10.48550/arXiv.2011.02288}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://alarm-redist.org/redist/}
  \item Report bugs at \url{https://github.com/alarm-redist/redist/issues}
}

}
\author{
\strong{Maintainer}: Christopher T. Kenny \email{ctkenny@proton.me} (\href{https://orcid.org/0000-0002-9386-6860}{ORCID})

Authors:
\itemize{
  \item Cory McCartan \email{mccartan@psu.edu} (\href{https://orcid.org/0000-0002-6251-669X}{ORCID})
  \item Ben Fifield \email{benfifield@gmail.com} (\href{https://orcid.org/0000-0002-2247-0201}{ORCID})
  \item Kosuke Imai \email{imai@harvard.edu} (\href{https://orcid.org/0000-0002-2748-1022}{ORCID})
}

Other contributors:
\itemize{
  \item Jun Kawahara \email{jkawahara@i.kyoto-u.ac.jp} [contributor]
  \item Alexander Tarr \email{atarr@princeton.edu} [contributor]
  \item Michael Higgins \email{mikehiggins@k-state.edu} [contributor]
}

}
\keyword{internal}
