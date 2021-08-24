#include "smc_base.h"

#ifndef DISTANCE_H
#define DISTANCE_H

double distance(double x1, double x2,
                double y1, double y2);

#endif


#ifndef DISTANCE_MATRIX_H
#define DISTANCE_MATRIX_H

NumericMatrix distance_matrix(NumericVector x,
                              NumericVector y);

#endif

#ifndef CLOSEST_ADJ_H
#define CLOSEST_ADJ_H

int closest_adj(IntegerVector adj,
                int i_dist,
                NumericVector x,
                NumericVector y);

#endif

#ifndef DIST_DIST_DIFF_H
#define DIST_DIST_DIFF_H

double dist_dist_diff(int p,
                      int i_dist,
                      int j_dist,
                      NumericVector x_center,
                      NumericVector y_center,
                      NumericVector x,
                      NumericVector y);

#endif
