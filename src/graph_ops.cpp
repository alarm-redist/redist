/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/3
* Purpose: Graph functions
********************************************************/

#include "graph_ops.h"

RegionMultigraph build_region_multigraph(
    Graph const &g, 
    arma::subview_col<arma::uword> const &region_ids,
    int const num_regions
){
    RegionMultigraph region_multigraph(num_regions);
    int const V = g.size();

    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        int v_region_num = region_ids(v);

        // now iterate over its neighbors
        for (int v_nbor : g[v]) {
            // find which region neighbor corresponds to
            int v_nbor_region_num = region_ids(v_nbor);

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

arma::imat build_region_laplacian(
    RegionMultigraph const &region_multigraph
){
    int num_regions = region_multigraph.size();
    arma::imat laplacian_mat(num_regions, num_regions, arma::fill::zeros);
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
RegionMultigraph get_region_multigraph(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
){
    std::unordered_set<int> uniqueElements;
    for (int element : region_ids) {
        uniqueElements.insert(element);
    }

    int num_regions = uniqueElements.size();
    return(build_region_multigraph(
        list_to_graph(adj_list), 
        region_ids.head(region_ids.size()),
        num_regions
    ));
}

arma::imat get_region_laplacian(
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



double compute_log_region_multigraph_spanning_tree(
    RegionMultigraph const &region_multigraph
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

    double lst, sign;
    arma::log_det(lst, sign, 
        laplacian_mat.submat(0, 0, num_regions-2, num_regions-2) );

    return lst;
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