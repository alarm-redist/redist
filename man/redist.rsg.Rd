\name{redist.rsg}
\alias{redist.rsg}
\title{ Redistricting via Random Seed and Grow Algorithm }
\description{ 
  \code{redist.rsg} generates redistricting plans using a random seed a grow algorithm.  This
  is the non-compact districting algorithm described in Chen and Rodden (2013).  The algorithm
  can provide start values for the other redistricting routines in this package.
}

\usage{
  redist.rsg(adj.list, population, ndists, thresh, verbose = TRUE, maxiter=5000)
}

\arguments{
  \item{adj.list}{  list of length N, where N is the number of precincts.  Each list element is an integer
    vector indicating which precincts that precinct is adjacent to.  It is assumed that precinct numbers
    start at 0.}
  \item{population}{ numeric vector of list N, where N is the number of precincts. Each element lists the
    population total of the corresponding precinct, and is used to enforce population constraints.
  }
  \item{ndists}{ integer, the number of districts we want to partition the precincts into.
  }
  \item{thresh}{ numeric, indicating how close district population targets have to be to the target population
    before algorithm converges.  thresh=0.05 for example means that all districts must be between 0.95 and 1.05
    times the size of target.pop in population size.
  }
  \item{verbose}{ boolean, indicating whether the time to run the algorithm is printed.}
  \item{maxiter}{ integer, indicating maximum number of iterations to attempt before convergence to population
    constraint fails.  If it fails once, it will use a different set of start values and try again.  If it
    fails again, redist.rsg() returns an object of all NAs, indicating that use of more iterations may be advised.}
}

\value{
  list, containing three objects containing the completed redistricting plan.
    \itemize{
    \item{\code{district_membership}}{  A vector of length N, indicating the district membership of each precinct.}
    \item{\code{district_list}}{  A list of length Ndistrict.  Each list contains a vector of the precincts in the respective district.}
    \item{\code{district_pop}}{  A vector of length Ndistrict, containing the population totals of the respective districts.}
    }                 

}


\author{ 
  Benjamin Fifield, Department of Politics, Princeton University 
  \email{bfifield@princeton.edu}, \url{http://www.benfifield.com}

  Michael Higgins, Department of Politics, Princeton University
  \email{mjh5@princeton.edu}, \url{http://www.princeton.edu/~mjh5}
  
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@princeton.edu}, \url{http://imai.princeton.edu}

  James Lo, \email{jameslo@princeton.edu}
  
  Alexander Tarr, Department of Electrical Engineering, Princeton University
  \email{atarr@princeton.edu}

}

\references{

Jowei Chen and Jonathan Rodden (2013) ``Unintentional Gerrymandering: Political Geography and Electoral Bias in Legislatures.''
Quarterly Journal of Political Science. 8(3): 239-269.

}


\examples{

### Real data example from test set

data("algdat.pfull")
res <- redist.rsg(algdat.pfull$adjlist, algdat.pfull$precinct.data$pop, 3, 0.05)


\dontrun{

### Example that generates test data from a square map with equal population districts
### Number of precincts is Nrows*Ncols
### getTest() outputs an adjacency list out of specified rows and columns

genTest <- function(Nrows,Ncols){

NN <- Nrows * Ncols

geog <- matrix(NA,nrow=Nrows+2, ncol=Ncols+2)
geog[2:(Nrows+1), 2:(Ncols+1)] <- 0:(NN-1)

adj.list <- vector("list", NN)

for(i in 2:(Nrows+1)){
	for(j in 2:(Ncols+1)){
		adj.list[[ geog[i,j] + 1 ]] <- c(geog[i-1,j],geog[i+1,j],geog[i,j-1],geog[i,j+1])
	}
}
adj.list <- lapply(adj.list, na.omit)
adj.list <- lapply(adj.list, as.numeric)
return(adj.list)
}

### Generate a 100x100 precinct map and redistrict it into 10 districts
adj.list <- genTest(100,100)
population <- rep(300,length(adj.list))
tmp <- redist.rsg(adj.list, population, 10, 0.05)

}

}

\keyword{ redistricting }
