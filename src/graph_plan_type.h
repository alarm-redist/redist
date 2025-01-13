
#include "base_plan_type.h"


class GraphPlan : public Plan {
    // private member variable
    

public:
    using Plan::Plan;
    // implementation of the pure virtual function
    

    void update_vertex_info_from_cut(
            Tree &ust, EdgeCut cut_edge, 
            const int split_region1_id, const int split_region2_id,
            bool split_district_only
    );
};