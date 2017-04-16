#!/bin/bash

# Build documentation, compile C++ attributes
# R -e 'sink("src/redist_init.c");tools::package_native_routine_registration_skeleton(".");sink()'
R -e 'devtools::document()'
R -e 'library(Rcpp);compileAttributes(verbose = TRUE)'

# Clean up src folder before build
cd src/
rm -rf *.o
rm -rf *.so
rm -rf *.rds

cd ../..

# Build and run CRAN checks
R CMD BUILD redist --resave-data 
R CMD CHECK redist_*.tar.gz --as-cran
