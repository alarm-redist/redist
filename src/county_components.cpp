/********************************************************
* Author: Philip O'Sullivan
* Institution: Harvard University
* Date Created: 2025/4
* Purpose: Class for Managing county component aspect of weight calculations
********************************************************/

#include "county_components.h"

CountyComponents::CountyComponents(
    MapParams const &map_params, int const num_regions
) :
map_params(map_params), 
num_regions(num_regions),
max_possible_num_componets(num_regions-1+map_params.num_counties),
county_district_lookup_table(map_params.num_counties > 1 ? num_regions * (map_params.num_counties) : 0, MAX_SUPPORTED_COUNTYREGION_VALUE),
vertices_visited(map_params.num_counties > 1 ? map_params.V : 0),
component_pairs_visited(map_params.num_counties > 1 ? max_possible_num_componets*max_possible_num_componets : 0),
county_component_graph(map_params.num_counties > 1 ? max_possible_num_componets : 0),
county_component_forest(map_params.num_counties > 1 ? max_possible_num_componets : 0),
region_component_counts(map_params.num_counties > 1 ? num_regions : 0){
    // check not too many regions or counties
    if(num_regions > MAX_SUPPORTED_NUM_DISTRICTS){
        REprintf("The maximum number of supported districts is %d!\n", MAX_SUPPORTED_NUM_DISTRICTS);
        throw Rcpp::exception("Number of districts is too large!\n");
    }
    if(map_params.num_counties > MAX_SUPPORTED_NUM_DISTRICTS){
        REprintf("The maximum number of supported counties is %d!\n", MAX_SUPPORTED_NUM_DISTRICTS);
        throw Rcpp::exception("Number of counties is too large!\n");
    }
}



// This builds the county component graph and tree 
void CountyComponents::build_component_graph_and_tree(
    PlanVector const &region_ids
){
    // first reset visited vectors
    std::fill(vertices_visited.begin(), vertices_visited.end(), false);
    std::fill(component_pairs_visited.begin(), component_pairs_visited.end(), false);
    // reset graph and tree
    for (size_t i = 0; i < max_possible_num_componets; i++)
    {
        county_component_graph[i].clear();
        county_component_forest[i].clear();
    }
    // reset component lookup table
    std::fill(
        county_district_lookup_table.begin(),
        county_district_lookup_table.end(),
        MAX_SUPPORTED_COUNTYREGION_VALUE
    );

    // max possible number of components is num_counties + num_regions-1
    int max_num_components = map_params.num_counties + num_regions - 1;
    // 
    int num_splits = 0;
    num_components = 0; // counts the total number of connected components in all county intersect region
    // We iterate over all vertices. For each vertex we explore all
    // connected neighbors with the same district and county 


    // counter for new vertices in the graph
    // each vertex represents a component 
    CountyRegion new_component_graph_index = 0;

    // visited every vertices in the graph
    for (int u = 0; u < map_params.V; u++)
    {
        // skip if we've already visited 
        if(vertices_visited[u]) continue;
        
        int current_county = map_params.counties(u)-1;
        int current_region = region_ids[u];
        
        int current_component_lookup_index = (current_county * num_regions) + current_region;
        CountyRegion current_component_vertex_id;
        // check if we've already declared this component 
        // if default value then we haven't encountered this component yet 
        if(county_district_lookup_table[current_component_lookup_index] == MAX_SUPPORTED_COUNTYREGION_VALUE){
            // increase number of connected components
            ++num_components;
            // we add this as a vertex 
            county_district_lookup_table[current_component_lookup_index] = new_component_graph_index;
            current_component_vertex_id = new_component_graph_index;
            // increment the counter for next one 
            new_component_graph_index++;            
        }else{ // else we've already seen this before 
            current_component_vertex_id = county_district_lookup_table[current_component_lookup_index];
        }
        // declare the vertex which is component ID, region, county
        std::tuple<CountyRegion, RegionID, RegionID> current_component_vertex(
            current_component_vertex_id, current_region, current_county
        );

        
        // Now we traverse all neighbors in the same component 
        std::queue<int> vertex_queue;
        vertex_queue.push(u);

        while(!vertex_queue.empty()){
            // get from queue
            int v = vertex_queue.front(); vertex_queue.pop();
            int v_region = region_ids[v];
            int v_county = map_params.counties(v);

            if(current_county != v_county || current_region != v_region){
                REprintf("BIG TIME ERROR IN COUNT SPLITS!!\n");
                throw Rcpp::exception("WOAHHHHH\n");
            }
            // mark this as visited 
            vertices_visited[v] = true;
            // go through children 
            for(auto const child_vertex: map_params.g[v]){
                // ignore if we already visited 
                if(vertices_visited[child_vertex]) continue;

                // check if this vertex is in the same component as our current one
                int child_region = region_ids[child_vertex];
                int child_county = map_params.counties(child_vertex)-1;
                // if same then same component 
                if(child_region == current_region && child_county == current_county){
                    // mark as visited to avoid being added later  
                    vertices_visited[child_vertex] = true;
                    vertex_queue.push(child_vertex);
                }else{ // else different components 
                    // check if we've already defined this component 
                    int child_component_lookup_index = (child_county * num_regions) + child_region;
                    CountyRegion child_component_vertex_id;
                    if(county_district_lookup_table[child_component_lookup_index] == MAX_SUPPORTED_COUNTYREGION_VALUE){
                        // this means we haven't so we need to give this vertex a component id
                        county_district_lookup_table[child_component_lookup_index] = new_component_graph_index;
                        child_component_vertex_id = new_component_graph_index;
                        // increment the counter for next one 
                        ++new_component_graph_index;
                        // also increase number of connected components
                        ++num_components;
                    }else{ // else we've already defined this component 
                        child_component_vertex_id = county_district_lookup_table[child_component_lookup_index];
                    }
                    // check if we've already added this edge to the graph if so then continue 
                    auto component_pair_index1 = (current_component_vertex_id * max_num_components) + child_component_vertex_id;
                    if(component_pairs_visited[component_pair_index1]) continue;
                    // else add this edge to the graph and mark as added 
                    auto component_pair_index2 = (child_component_vertex_id * max_num_components) + current_component_vertex_id;
                    std::tuple<CountyRegion, RegionID, RegionID> child_component_vertex(
                        child_component_vertex_id, child_region, child_county
                    );
                    // add current to child and vice verse
                    county_component_graph[child_component_vertex_id].push_back(current_component_vertex);
                    county_component_graph[current_component_vertex_id].push_back(child_component_vertex);
                    // now mark this component pair as added 
                    component_pairs_visited[component_pair_index1] = component_pairs_visited[component_pair_index2] = true;
                }
            }
        }
    }

    // REprintf("%d components and %d regions = %d \n", num_connected_components, map_params.num_counties,
    //     num_connected_components - map_params.num_counties);

    return;
    //return num_connected_components - map_params.num_counties;
}