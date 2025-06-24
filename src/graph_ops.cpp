/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/3
* Purpose: Graph functions
********************************************************/

#include "graph_ops.h"



RegionMultigraphCount build_region_multigraph(
    Graph const &g, 
    PlanVector const &region_ids,
    int const num_regions
){
    RegionMultigraphCount region_multigraph(num_regions);
    int const V = g.size();

    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        int v_region_num = region_ids[v];

        // now iterate over its neighbors
        for (int v_nbor : g[v]) {
            // find which region neighbor corresponds to
            int v_nbor_region_num = region_ids[v_nbor];

            // to avoid double counting only count when v less u
            if(v_region_num < v_nbor_region_num){
                // we increase the count of edges
                region_multigraph[v_region_num][v_nbor_region_num]++;
                region_multigraph[v_nbor_region_num][v_region_num]++;
            }
        }
    }

    return region_multigraph;
}





double get_log_merged_region_multigraph_spanning_tree(
    RegionMultigraphCount const &region_multigraph,
    std::vector<int> &merge_index_reshuffle,
    int region1_id, int region2_id
){
    // TODO make the merged reindex the biggest one so we can skip it later
    int num_regions = region_multigraph.size();
    int smaller_id = std::min(region1_id, region2_id);
    int bigger_id = std::max(region1_id, region2_id);
    int merged_reindex = num_regions-2;
    for (int current_reindex = 0, i = 0; i < num_regions; i++){
        if(i == region1_id || i == region2_id){
            merge_index_reshuffle[i] = merged_reindex;
        }else{
            merge_index_reshuffle[i] = current_reindex;
            ++current_reindex;
        }
    }
    
    // Rprintf("Reindex is: ");
    // for (size_t i = 0; i < num_regions; i++)
    // {
    //     Rprintf("%d maps to %d\n", i, merge_index_reshuffle[i]);
    // }
    


    arma::mat laplacian_mat(
        num_regions-2, 
        num_regions-2, 
        arma::fill::zeros);
    // iterate over the multigraph but ignore last column and row 
    for (size_t region_id = 0; region_id < num_regions; region_id++)
    {
        int reindexed_id = merge_index_reshuffle[region_id];
        // Rprintf("%d remapped to %d!\n", region_id, reindexed_id);
        // skip 
        if(reindexed_id == num_regions-2) continue;
        int vertex_degree = 0;
        // iterate over neighbors
        for (auto const &it : region_multigraph[region_id])
        {
            int reindexed_nbor = merge_index_reshuffle[it.first];
            // ignore loops within the merged region 
            if(reindexed_id == reindexed_nbor) continue;
            // Rprintf("%d remapped to %d!\n", it.first, reindexed_nbor);
            
            // Rprintf("(%d,%d) Edges %d \n", region_id, it.first, it.second);

            // add number of edges to degree
            vertex_degree += it.second;
            // ignore if num_regions-1 since we're skipping last one
            if(reindexed_nbor == num_regions-2) continue;
            // need minus equals here to add for both merged regions
            laplacian_mat(reindexed_id, reindexed_nbor) -= it.second;
            // laplacian_mat(reindexed_nbor, reindexed_id) -= it.second;
        }
        // need plus equals here to add for both merged regions
        laplacian_mat(reindexed_id, reindexed_id) += vertex_degree;
    }

    
    // laplacian_mat.print();
    
    double lst, sign;
    arma::log_det(lst, sign, laplacian_mat);

    return lst;

}

arma::mat build_region_laplacian(
    RegionMultigraphCount const &region_multigraph
){
    int num_regions = region_multigraph.size();
    arma::mat laplacian_mat(num_regions, num_regions, arma::fill::zeros);
    // iterate over the multigraph
    for (size_t region_id = 0; region_id < num_regions; region_id++)
    {
        int vertex_degree = 0;
        // iterate over neighbors
        for (auto const &it : region_multigraph[region_id])
        {
            // add number of edges to degree
            vertex_degree += it.second;
            laplacian_mat(region_id, it.first) = -it.second;
            laplacian_mat(it.first, region_id) = -it.second;
        }
        laplacian_mat(region_id, region_id) = vertex_degree;
    }

    return(laplacian_mat);
}

// Can call from R
RegionMultigraphCount get_region_multigraph(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
){
    std::unordered_set<int> uniqueElements;
    for (int element : region_ids) {
        uniqueElements.insert(element);
    }

    AllPlansVector underlying_id_vec(
        region_ids.begin(),
        region_ids.end()
    );

    PlanVector region_id_vec(underlying_id_vec, 0, underlying_id_vec.size()-1);


    int num_regions = uniqueElements.size();
    return(build_region_multigraph(
        list_to_graph(adj_list), 
        region_id_vec,
        num_regions
    ));
}

arma::mat get_region_laplacian(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
){
    return(
        build_region_laplacian(get_region_multigraph(
            adj_list,
            region_ids
        ))
    );
}


// avoid submat because it copies 
double compute_log_region_multigraph_spanning_tree(
    RegionMultigraphCount const &region_multigraph
){
    int num_regions = region_multigraph.size();
    arma::mat laplacian_mat(num_regions, num_regions, arma::fill::zeros);
    // iterate over the multigraph
    for (size_t region_id = 0; region_id < num_regions; region_id++)
    {
        int vertex_degree = 0;
        // iterate over neighbors
        for (auto const &it : region_multigraph[region_id])
        {
            // add number of edges to degree
            vertex_degree += it.second;
            laplacian_mat(region_id, it.first) = -it.second;
            laplacian_mat(it.first, region_id) = -it.second;
        }
        laplacian_mat(region_id, region_id) = vertex_degree;
    }

    // laplacian_mat.submat(0, 0, num_regions-2, num_regions-2).print();

    return arma::log_det_sympd(laplacian_mat.submat(0, 0, num_regions-2, num_regions-2));
}


double get_log_number_linking_edges(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
){
    return(compute_log_region_multigraph_spanning_tree(
        get_region_multigraph(
            adj_list,
            region_ids
        )
    ));
}




double get_merged_log_number_linking_edges(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids,
    int const region1_id, int const region2_id
){
    RegionMultigraphCount region_multigraph = get_region_multigraph(
        adj_list,
        region_ids
    );

    std::vector<int> merge_index_reshuffle(region_multigraph.size());
    return get_log_merged_region_multigraph_spanning_tree(
        region_multigraph, merge_index_reshuffle,
        region1_id, region2_id
    );
}

