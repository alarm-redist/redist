#!/bin/bash

STARTDIR=$(pwd)
cd redist
Rscript -e "library('Rcpp') ; compileAttributes(verbose = TRUE)"
cd ..
R CMD INSTALL redist


