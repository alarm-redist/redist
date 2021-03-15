#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector talisman(NumericMatrix dvs, double nd, double alpha = 1, double beta = 1) {
  NumericVector comp(dvs.ncol());
  double Tp, Te, Br, curr;
  for(int c = 0; c < dvs.ncol(); c++){
    Br = 0.0;
    Tp = 0.0;
    for(int r = 0; r < dvs.nrow(); r++){
      curr = dvs(r,c);
      // Handle Tp
      Tp += fabs(0.5 - curr);
      // Handle Te
      if(curr < 0.5){
        Br += 1.0;
      }
    }
    Te = fabs(Br/nd - 0.5);
    Tp = Tp/nd;
    
    comp(c) = Tp * (1.0 + alpha * Te); 
  }
  comp = comp * beta;
  return comp;
}


/*** R
#dvs <- matrix(c(0.4,0.5,0.51,0.53,0.7,0.2,0.3,0.4,0.8,0.9), nrow = 5)
#nd <- 5
#talisman(dvs, nd)
# should return 0.0884 and 0.2860
*/


