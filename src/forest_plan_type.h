#include "base_plan_type.h"


class TreePlan : public Plan {


public:
    // We now need to keep track of trees as undirected graphs
    Graph forest_graph;
    std::vector<int> tree_roots;

    // implementation of the pure virtual function
    void fun() { cout << "fun() called"; }
};