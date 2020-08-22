#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix agg_p2d(IntegerMatrix dm, IntegerVector vote, int nd) {
  IntegerMatrix mat = IntegerMatrix(nd, dm.ncol());
  for(int j = 0; j < dm.ncol(); j++){
    for(int i = 0; i < dm.nrow(); i++){
      mat(dm(i,j)-1,j) += vote[i];
    }
  }
  return mat;
}


// [[Rcpp::export]]
IntegerVector dseats(IntegerMatrix dm, IntegerMatrix dcounts, IntegerMatrix rcounts, int nd){
  IntegerVector dseats(dm.ncol());
  
  for(int c = 0; c < dcounts.ncol(); c++){
    for(int r = 0; r < dcounts.nrow(); r++){
      if(dcounts(r,c) >= rcounts(r,c)){
        dseats(c) += 1;
      }
    }
  }
  
  return dseats;
}

// [[Rcpp::export]]
IntegerVector dseatsDVS(NumericMatrix dvs){
  IntegerVector dseats = IntegerVector(dvs.ncol());
  for(int c = 0; c < dvs.ncol(); c++){
    for(int r = 0; r < dvs.nrow(); r++){
      if(dvs(r,c) > .5){
        dseats[c] += 1;
      }
    }
  }
  return dseats;
}

// [[Rcpp::export]]
NumericMatrix DVS(IntegerMatrix dcounts, IntegerMatrix rcounts){
  NumericMatrix mat = NumericMatrix(dcounts.nrow(), dcounts.ncol());
  
  for(int c = 0; c < mat.ncol(); c++){
    for(int r = 0; r < mat.nrow(); r++){
      mat(r,c) = (double)dcounts(r,c)/(dcounts(r,c)+rcounts(r,c));
    }
  }
  return mat;
}

// [[Rcpp::export]]
NumericVector effgapEP(NumericMatrix dvs, IntegerVector dseat_vec, int nd){
  NumericVector V = colMeans(dvs);
  NumericVector S(dseat_vec.size());
  for(int s = 0; s < dseat_vec.size(); s++){
    S[s] = dseat_vec[s]/(double) nd;
  }
  NumericVector eg(dseat_vec.size());
  for(int s = 0; s < dseat_vec.size(); s++){
    eg(s) = -2*V[s] + S[s] + .5;
  }
  return eg;
}

// [[Rcpp::export]]
NumericVector effgap(IntegerMatrix dcounts, IntegerMatrix rcounts, int totvote){
  NumericVector eg(dcounts.ncol());
  
  IntegerMatrix dwaste(dcounts.nrow(), dcounts.ncol());
  IntegerMatrix rwaste(rcounts.nrow(), rcounts.ncol());
  int minwin;
  for(int c = 0; c < dcounts.ncol(); c++){
    for(int r = 0; r < dcounts.nrow(); r++){
      minwin = floor((dcounts(r,c) + rcounts(r,c))/2.0)+1;
      if(dcounts(r,c) > rcounts(r,c)){
        dwaste(r,c) += (dcounts(r,c) - minwin);
        rwaste(r,c) += rcounts(r,c);
      } else{
        dwaste(r,c) += dcounts(r,c);
        rwaste(r,c) += (rcounts(r,c) - minwin);
      }
    }
  }
  
  IntegerVector netwaste(dcounts.ncol());
  netwaste = colSums(dwaste) - colSums(rwaste);
  
  for(int i = 0; i < netwaste.size(); i++){
    eg[i] = netwaste[i]/(double)totvote;
  }
  
  return eg;
}


// [[Rcpp::export]]
NumericVector taugap(double tau, NumericMatrix dvs, IntegerVector dseat_vec, int nd){
  NumericMatrix ai_mat = NumericMatrix(dvs.nrow(), dvs.ncol());
  IntegerMatrix ei_mat = IntegerMatrix(dvs.nrow(), dvs.ncol());
  NumericMatrix expr = NumericMatrix(ai_mat.nrow(), ai_mat.ncol());
  for(int c = 0; c < ai_mat.ncol(); c++){
    for(int r = 0; r < ai_mat.nrow(); r++){
      ai_mat(r,c) = 2*dvs(r,c) - 1;
      if(ai_mat(r,c) >= 0){
        ei_mat(r,c) = 1;
      } else{
        ei_mat(r,c) = -1;
      }
      expr(r,c) = ei_mat(r,c)*pow(ai_mat(r,c)*ei_mat(r,c), tau+1);
    }
  }
  
  NumericVector temp = colSums(expr)/2;
  NumericVector dseat_share = NumericVector(dseat_vec.size());
  for(int i =0; i < dseat_share.size(); i++){
    dseat_share(i) = dseat_vec(i)/(double) nd;
  }

  return -2*(temp +.5 - dseat_share);
}

// [[Rcpp::export]]
NumericVector meanmedian(NumericMatrix dvs){
  NumericVector mm = NumericVector(dvs.ncol());
  NumericVector med = NumericVector(dvs.ncol());
  NumericVector col = dvs(_,1);
  for(int c = 0; c < dvs.ncol(); c++){
    col = dvs(_,c);
   med(c) = median(col); 
  }
  mm = colMeans(dvs) - med;
  return mm;
}
// [[Rcpp::export]]
NumericVector bias(NumericMatrix dvs, int nd){
  NumericVector sw = .5 - colMeans(dvs);
  NumericMatrix dvs_sw =  clone(dvs);
  for(int c = 0; c < dvs_sw.ncol(); c++){
    for(int r = 0; r < dvs_sw.nrow(); r++){
      dvs_sw(r,c) += sw(c);
    }
  }
  
  IntegerVector newseats = dseatsDVS(dvs_sw);
  NumericVector seatshare = (NumericVector) newseats/nd;
  NumericVector bias = seatshare - 0.5;
  
  return bias;
}

// [[Rcpp::export]]
NumericVector declination(NumericMatrix dvs, IntegerVector dseat_vec, int nd){
  NumericVector Dwin = NumericVector(dvs.ncol()); 
  NumericVector Rwin = NumericVector(dvs.ncol());
  
  for(int c = 0; c < dvs.ncol(); c++){
    for(int r = 0; r < dvs.nrow(); r++){
      if(dvs(r,c) >= .5){
        Dwin(c) += dvs(r,c);
      } else{
        Rwin(c) += dvs(r,c);
      }
    }
  }
  for(int i =0; i < Dwin.size(); i++){
    Dwin(i) = Dwin(i)/dseat_vec(i);
    Rwin(i) = Rwin(i)/(nd-dseat_vec(i));
  }
  
  NumericVector dseatshare = (NumericVector)dseat_vec/(double)nd;
  
  return ((Dwin-.5)/dseatshare)-((0.5-Rwin)/(1-dseatshare));
}

// [[Rcpp::export]]
NumericVector lopsidedwins(NumericMatrix dvs, IntegerVector dseat_vec, int nd){
  NumericVector Dwin = NumericVector(dvs.ncol()); 
  NumericVector Rwin = NumericVector(dvs.ncol());
  
  for(int c = 0; c < dvs.ncol(); c++){
    for(int r = 0; r < dvs.nrow(); r++){
      if(dvs(r,c) >= .5){
        Dwin(c) += dvs(r,c);
      } else{
        Rwin(c) += dvs(r,c);
      }
    }
  }
  for(int i =0; i < Dwin.size(); i++){
    Dwin(i) = Dwin(i)/dseat_vec(i);
    Rwin(i) = Rwin(i)/(nd-dseat_vec(i));
  }
  
  return Dwin + Rwin - 1.0;
}


// [[Rcpp::export]]
NumericVector responsiveness(NumericMatrix dvs, double v, int nd, double bandwidth = .01){
  NumericVector right = (v + (bandwidth/2.0)) - colMeans(dvs);
  NumericVector left = (v - (bandwidth/2.0)) - colMeans(dvs);
  NumericMatrix dvs_right = clone(dvs);
  NumericMatrix dvs_left = clone(dvs);
  
  for(int c = 0; c < dvs.ncol(); c++){
    for(int r = 0; r < dvs.nrow(); r++){
      dvs_right(r,c) += right(c);
      dvs_left(r,c) += left(c);
    }
  }
  
  NumericVector seat_right = (NumericVector)dseatsDVS(dvs_right)/(double)nd;
  NumericVector seat_left = (NumericVector)dseatsDVS(dvs_left)/(double)nd;

  return (seat_right - seat_left)/bandwidth;
}

// [[Rcpp::export]]
NumericVector biasatv(NumericMatrix dvs, double v, int nd){
  NumericVector dshift = (v) - colMeans(dvs);
  NumericVector rshift = (1-v) - colMeans(dvs);
  NumericMatrix dvs_dshift = clone(dvs);
  NumericMatrix dvs_rshift = clone(dvs);
  
  for(int c = 0; c < dvs.ncol(); c++){
    for(int r = 0; r < dvs.nrow(); r++){
      dvs_dshift(r,c) += dshift(c);
      dvs_rshift(r,c) += rshift(c);
    }
  }
  
  NumericVector seat_dshift = (NumericVector)dseatsDVS(dvs_dshift)/(double)nd;
  NumericVector seat_rshift = 1.0 - (NumericVector)dseatsDVS(dvs_rshift)/(double)nd;
  
  return (seat_dshift - seat_rshift)/2;
}

