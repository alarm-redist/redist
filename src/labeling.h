#include <map>
#include "smc_base.h"
#include "tree_op.h"

#ifndef LABELING_H
#define LABELING_H

double log_labelings_exact(const Graph &g);

double log_labelings_IS(const Graph &g, int n=1000);

#endif
