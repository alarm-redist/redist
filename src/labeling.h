#include <map>
#include "smc_base.h"
#include "tree_op.h"

#define LABELING_H
#ifndef LABELING_H

double log_labelings_exact(const Graph &g);

double log_labelings_IS(const Graph &g, int n_eff=1000);


#endif
