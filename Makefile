default : install

attrs :
	cd redist ; Rscript -e "library(Rcpp) ; compileAttributes(verbose=TRUE)"

build : attrs clean
	R CMD build redist

install : attrs
	R CMD INSTALL redist

fullinstall : build
	R CMD INSTALL redist*.tar.gz

remove :
	R CMD REMOVE redist

clean :
	rm -f redist*.tar.gz
	rm -f ./redist/src/*.o
	rm -f ./redist/src/*.so
	rm -f ./redist/src/*.rds
