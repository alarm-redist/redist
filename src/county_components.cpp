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
counties_on(map_params.num_counties > 1),
num_regions(num_regions),
max_possible_num_componets(num_regions-1+map_params.num_counties),
county_district_lookup_table(counties_on ? num_regions * (map_params.num_counties) : 0, MAX_SUPPORTED_COUNTYREGION_VALUE),
vertices_visited(counties_on ? map_params.V : 0),
components_visited(counties_on ? max_possible_num_componets : 0),
merged_components_visited(counties_on ? max_possible_num_componets : 0),
component_pairs_visited(counties_on ? max_possible_num_componets*max_possible_num_componets : 0),
counties_component_adj(counties_on ? (map_params.num_counties * (map_params.num_counties-1)) / 2 : 0),
county_component_graph(counties_on ? max_possible_num_componets : 0),
region_vertices(counties_on ? num_regions : 0),
region_component_counts(counties_on ? num_regions : 0),
component_graph_vertices(counties_on ? max_possible_num_componets : 0){
    // check not too many regions or counties
    if(num_regions > MAX_SUPPORTED_NUM_DISTRICTS){
        REprintf("The maximum number of supported districts is %d!\n", MAX_SUPPORTED_NUM_DISTRICTS);
        throw Rcpp::exception("Number of districts is too large!\n");
    }
    if(map_params.num_counties > MAX_SUPPORTED_NUM_COUNTIES){
        REprintf("The maximum number of supported counties is %d!\n", MAX_SUPPORTED_NUM_DISTRICTS);
        throw Rcpp::exception("Number of counties is too large!\n");
    }
}



// This builds the county component graph and tree 
void CountyComponents::build_component_graph_and_tree(
    PlanVector const &region_ids
){
    // if no counties then just return 
    if(!counties_on) return;
    // first reset visited vectors
    std::fill(vertices_visited.begin(), vertices_visited.end(), false);
    std::fill(component_pairs_visited.begin(), component_pairs_visited.end(), false);
    // reset graph and tree
    for (size_t i = 0; i < max_possible_num_componets; i++)
    {
        county_component_graph[i].clear();
    }
    // reset component lookup table
    std::fill(
        county_district_lookup_table.begin(),
        county_district_lookup_table.end(),
        MAX_SUPPORTED_COUNTYREGION_VALUE
    );

    // REprintf("Component pair visited size %u\n",
    //     component_pairs_visited.size());

    num_components = 0; // counts the total number of connected components in all county intersect region
    // We iterate over all vertices. For each vertex we explore all
    // connected neighbors with the same district and county 


    // counter for new vertices in the graph
    // each vertex represents a component 
    CountyRegion new_component_graph_index = 0;

    // visited every vertices in the graph to build component graph 
    for (int u = 0; u < map_params.V; u++)
    {
        // skip if we've already visited 
        if(vertices_visited[u]) continue;
        
        CountyID current_county = map_params.counties(u)-1;
        auto current_region = region_ids[u];
        // Rprintf("Component (%u, %u)\n", current_region, current_county);
        // for lookup we pretend table is num_counties x num_regions
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
        CountyComponentVertex current_component_vertex(
            current_component_vertex_id, current_region, current_county
        );
        component_graph_vertices[current_component_vertex_id] = current_component_vertex;
        // also make this the vertex for the region
        region_vertices[current_region] = current_component_vertex;

        
        // Now we traverse all neighbors in the same component 
        std::queue<int> vertex_queue;
        vertex_queue.push(u);

        while(!vertex_queue.empty()){
            // get from queue
            // we don't actually need to save this. Can delete in future just for error checking now
            int v = vertex_queue.front(); vertex_queue.pop();
            int v_region = region_ids[v];
            int v_county = map_params.counties(v)-1;

            if(current_county != v_county || current_region != v_region){
                REprintf("BIG TIME ERROR IN COUNT SPLITS!!\n");
                throw Rcpp::exception("WOAHHHHH\n");
            }
            // if(v_county == 31){
            //     REprintf("v=%d, region=%d, county=%u\n", v, v_region, v_county);
            // }
            // mark this as visited 
            vertices_visited[v] = true;
            // go through children 
            for(auto const child_vertex: map_params.g[v]){
                // ignore if we already visited 
                if(vertices_visited[child_vertex]) continue;
            
                auto child_region = region_ids[child_vertex];
                CountyID child_county = map_params.counties(child_vertex)-1;

                // REprintf("Child %d - (%u, %u) ---- ", child_vertex, child_region, child_county);
                // check if this vertex is in the same component as our current one
                if(child_region == current_region && child_county == current_county){
                    // if same then same component  mark as visited to avoid being added later  
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
                    // REprintf(" --- %d, VID=(%u, %u) ", child_component_lookup_index,
                    //     child_component_vertex_id, current_component_vertex_id);
                    // check if we've already added this edge to the graph if so then continue 
                    auto component_pair_index1 = (current_component_vertex_id * max_possible_num_componets) + child_component_vertex_id;
                    if(component_pairs_visited[component_pair_index1]){
                        // REprintf("Skipped!\n");
                        continue;
                    } 
                    // else add this edge to the graph and mark as added 
                    auto component_pair_index2 = (child_component_vertex_id * max_possible_num_componets) + current_component_vertex_id;
                    // REprintf(" %d, %d ---- ", component_pair_index1, component_pair_index2);
                    CountyComponentVertex child_component_vertex(
                        child_component_vertex_id, child_region, child_county
                    );
                    // REprintf("-- Size = %d --", county_component_graph.size());
                    // add current to child and vice verse
                    county_component_graph[child_component_vertex_id].push_back(current_component_vertex);
                    county_component_graph[current_component_vertex_id].push_back(child_component_vertex);
                    // now mark this component pair as added 
                    component_pairs_visited[component_pair_index1] = component_pairs_visited[component_pair_index2] = true;
                }
                // REprintf("Done!\n");
            }
        }
        // REprintf("Done with component!\n");
    }

    // if(num_regions == 2){
    // Rprintf("There are %d components and size %u\n",
    //     num_components, county_component_graph.size());
    
    // for (size_t i = 0; i < county_component_graph.size(); i++)
    // {
    //     auto component_vtx = component_graph_vertices[i];
    //     Rprintf("(%d, %d): ",
    //         std::get<1>(component_vtx),
    //         std::get<2>(component_vtx));
    //     for(auto const &a_vtx: county_component_graph[i]){
    //         Rprintf("(%d, %d), ",
    //         std::get<1>(a_vtx),
    //         std::get<2>(a_vtx));
    //         // better be connected
    //         auto component_pair_index1 = (std::get<0>(component_vtx) * max_possible_num_componets) + std::get<0>(a_vtx);
    //         if(!component_pairs_visited[component_pair_index1]){
    //             REprintf(" [Index %u but false!] ", 
    //                 component_pair_index1);
    //         }
    //     }
    //     Rprintf("\n");
    // }
    // REprintf("Now printing the plan!\nc(");
    // // for (size_t i = 0; i < size.region_ids; i++)
    // // {
    // //     /* code */
    // // }
    
    // for(auto const &v: region_ids){
    //     REprintf("%u, ", v);
    //     // Rcpp::Rcerr << v << ", ";
    // }
    // REprintf("\n");
    // }
    



    return;
}




std::pair<bool, int> CountyComponents::count_county_splits(
    Plan const &plan
){
    // if 1 county then its just num_regions - 1
    if(map_params.num_counties == 1) return std::make_pair(false, plan.num_regions - 1);
    // reset internal lookup table 
    std::fill(
        county_district_lookup_table.begin(),
        county_district_lookup_table.end(),
        MAX_SUPPORTED_COUNTYREGION_VALUE
    );
    // reset visit vector 
    std::fill(
        vertices_visited.begin(), 
        vertices_visited.end(), 
        false
    );
    

    

    bool more_than_one_district_intersect_county = false;
    // 
    int num_splits = 0;
    int num_connected_components = 0; // counts the total number of connected components in all county intersect region
    // We iterate over all vertices. For each vertex we explore all
    // connected neighbors with the same district and county 


    for (int u = 0; u < map_params.V; u++)
    {
        // skip if we've already visited 
        if(vertices_visited[u]) continue;
        // else we've encountered a new connected component so increase the count
        ++num_connected_components;
        // Now we traverse all neighbors with the same county and region 
        int current_county = map_params.counties(u)-1;
        int current_region = plan.region_ids[u];
        // check if we've already visited this district x county component 
        int component_lookup_index = (current_county * num_regions) + current_region;
        if(county_district_lookup_table[component_lookup_index] == MAX_SUPPORTED_COUNTYREGION_VALUE){
            // means we haven't so we'll set this to dummy value of 0
            county_district_lookup_table[component_lookup_index] = 0;
        }else{
            // means we've seen this before so more than one 
            more_than_one_district_intersect_county = true;
        }

        std::queue<int> vertex_queue;
        vertex_queue.push(u);

        while(!vertex_queue.empty()){
            // get from queue
            int v = vertex_queue.front(); vertex_queue.pop();
            int v_region = plan.region_ids[v];
            int v_county = map_params.counties(v)-1;

            if(current_county != v_county || current_region != v_region){
                REprintf("BIG TIME ERROR IN COUNT SPLITS!!\n");
                throw Rcpp::exception("WOAHHHHH\n");
            }
            // mark this as visited 
            vertices_visited[v] = true;
            // go through children 
            for(auto const child_vertex: map_params.g[v]){
                // add the children if same region and county and not visited before 
                if(map_params.counties(child_vertex)-1 == current_county && 
                    plan.region_ids[child_vertex] == current_region &&
                    !vertices_visited[child_vertex]){
                    // mark as visited to avoid being added later  
                    vertices_visited[child_vertex] = true;
                    vertex_queue.push(child_vertex);
                }                
            }
        }

    }

    int num_county_splits = num_connected_components - map_params.num_counties;
    return std::make_pair(more_than_one_district_intersect_county, num_county_splits);
}



bool CountyComponents::check_merging_regions_is_ok(
    int const region1_id, int const region2_id
){
    // reset visited vector 
    std::fill(
        components_visited.begin(), 
        components_visited.end(), 
        false
    );
    std::fill(
        merged_components_visited.begin(), 
        merged_components_visited.end(), 
        false
    );
    int removed_splits = 0;
    // we start at some component vertex in region 1
    std::queue<CountyComponentVertex> vertex_queue;
    vertex_queue.push(region_vertices[region1_id]);

    while(!vertex_queue.empty()){
        // get the vertex component id
        CountyComponentVertex actual_component_vertex = vertex_queue.front(); 
        vertex_queue.pop();

        CountyRegion actual_component_vertex_id = std::get<0>(actual_component_vertex);
        RegionID actual_region = std::get<1>(actual_component_vertex);
        RegionID other_region = actual_region == region1_id ? region2_id : region1_id;
        CountyID current_county = std::get<2>(actual_component_vertex);
        // REprintf("Visiting (%u, %u)\n", actual_region, current_county);
        // get the lookup of this component 
        CountyRegion actual_current_component_lookup_index = (current_county * num_regions) + actual_region;
        // get the lookup of the component of the other region (if it exists)
        CountyRegion actual_other_component_lookup_index = (current_county * num_regions) + other_region;
        CountyRegion other_component_vertex_id = county_district_lookup_table[actual_other_component_lookup_index];
        // see if other component exists meaning lookup index isn't dummy value
        if(other_component_vertex_id != MAX_SUPPORTED_COUNTYREGION_VALUE){
            // if it exists then we need to mark the merged component id as visited 
            // if region1 the id is just itself but if region2 then its the id for region1
            // if component is in region 2 and region 1 has a component in this county 
            // we pretend its just region 1
            CountyRegion merged_component_vertex_id = actual_region == region1_id ?
                actual_component_vertex_id: other_component_vertex_id;

            // check if already visited 
            if(merged_components_visited[merged_component_vertex_id]){
                // REprintf("Already visited (%u, %u)\n", actual_region, current_county);
                // if already visited then we need to see if there is an adjacent vertex
                // in the same county but the other region. If not then this is an illegal merge
                auto other_region_component_vertex_id = county_district_lookup_table[actual_other_component_lookup_index];
                auto component_pair_index = (other_region_component_vertex_id * max_possible_num_componets) + actual_component_vertex_id;
                if(component_pairs_visited[component_pair_index]){
                    // this means they're adjacent and in the same county so increase removed
                    // splits by 1
                    ++removed_splits;
                }else{
                    // this means they're not adjacent despite being in the same county so 
                    // not allowed meaning we immideately return false!
                    // REprintf("\n ALERT: (%u, %u) Disconnected in County %u!\n",
                    //     region1_id, region2_id, current_county);
                    return false;
                }
            }else{
                // mark merged as visited and move on 
                merged_components_visited[merged_component_vertex_id] = true;
            }
        }
        // mark this actual component as visited 
        components_visited[actual_component_vertex_id] = true;

        // Now we add any unvisited adjacent components in either region 1 or 2
        for(auto const &child_vertex: county_component_graph[actual_component_vertex_id]){
            auto child_vertex_id = std::get<0>(child_vertex);
            auto child_region = std::get<1>(child_vertex);
            // only add if not visited and in region 1 or 2
            if(!components_visited[child_vertex_id] &&
                (child_region == region1_id || child_region == region2_id)
            ){
                // mark as visited to avoid being added later  
                components_visited[child_vertex_id] = true;
                // REprintf("Adding (%u, %u)\n", child_region, std::get<2>(child_vertex));
                vertex_queue.push(child_vertex);
            }
        }
    }
    // now the number of splits for the plan with the two regions merged is
    // num_components - removed_splits - map_params.num_counties
    // Now return whether or not thats less than or equal to number of regions minus 2
    // if(num_components - removed_splits - map_params.num_counties > num_regions - 2){
    //     REprintf("Comparing %u vs %d!\n",
    //         num_components - removed_splits - map_params.num_counties, num_regions - 2);
    // }
    
    return num_components - removed_splits - map_params.num_counties <= num_regions - 2;
}



// only returns true if 
// - NO edge between (county_A, region1_id), (county_B, region1_id)
// - NO edge between (county_A, region2_id), (county_B, region2_id)
bool CountyComponents::count_county_boundary(
    RegionID region1_id, CountyID county_A,
    RegionID region2_id, CountyID county_B
) const{
    auto county_A_region1_lookup = mat_index_from_pair(county_A, region1_id, num_regions);
    // check if this component even exists 
    auto county_B_region1_lookup = mat_index_from_pair(county_B, region1_id, num_regions);
    // check if county B intersect region 1 exists and they share an edge
    if(county_district_lookup_table[county_B_region1_lookup] != MAX_SUPPORTED_COUNTYREGION_VALUE &&
        component_pairs_visited[mat_index_from_pair(
            county_district_lookup_table[county_B_region1_lookup], 
            county_district_lookup_table[county_A_region1_lookup], 
            max_possible_num_componets)]
        ){  
            // if yes then this edge can't be counted
            return false;
    }


    auto county_A_region2_lookup = mat_index_from_pair(county_A, region2_id, num_regions);
    // check if this component even exists 
    auto county_B_region2_lookup = mat_index_from_pair(county_B, region2_id, num_regions);
    // check if county B intersect region 2 exists and they share an edge
    if(county_district_lookup_table[county_A_region2_lookup] != MAX_SUPPORTED_COUNTYREGION_VALUE &&
        component_pairs_visited[mat_index_from_pair(
            county_district_lookup_table[county_B_region2_lookup], 
            county_district_lookup_table[county_A_region2_lookup], 
            max_possible_num_componets)]
        ){  
            // if yes then this edge can't be counted
            return false;
    }

    return true;
}


bool CountyComponents::check_is_county_component_multigraph_valid(
    Graph &county_graph
){
    // reset graph 
    for (size_t i = 0; i < map_params.num_counties; i++)
    {
        county_graph[i].clear();
    }
    
    // reset component adj
    std::fill(
        counties_component_adj.begin(), 
        counties_component_adj.end(), 
        false
    ); 
    // reset visited vector 
    std::fill(
        components_visited.begin(), 
        components_visited.end(), 
        false
    );
    // now we walk through the graph 
    std::queue<CountyComponentVertex> vertex_queue;
    // start at region 0
    vertex_queue.push(region_vertices[0]);
    


    while(!vertex_queue.empty()){
        // 
        CountyComponentVertex current_vertex = vertex_queue.front(); 
        vertex_queue.pop();
        // get the county of this component 
        CountyID current_county = std::get<2>(current_vertex);
        RegionID current_region = std::get<1>(current_vertex);
        auto current_vertex_id = std::get<0>(current_vertex);
        // mark as visited
        components_visited[current_vertex_id] = true;
        // now we check its children 
        for(auto const &child_vertex: county_component_graph[current_vertex_id]){
            auto child_vertex_id = std::get<0>(child_vertex);
            auto child_region = std::get<1>(child_vertex);
            auto child_county = std::get<2>(child_vertex);
            // check if the child shares the same region but not the same county 
            // and current component index less than child to avoid double counting 
            if(current_region == child_region && 
                current_county != child_county &&
                current_vertex_id < child_vertex_id
            ){
                    // REprintf("Regions (%u, %u), County (%u, %u)\n",
                        // current_region, child_region, current_county, child_county);
                    // REprintf("(Region %u, County %u) and (Region %u, County %u)\n",
                    //     current_region, current_county, child_region, child_county
                    // );
                    // REprintf("(County %u) and (County %u)\n",
                    //     current_county, child_county
                    // );
                    // check if this county pair has already been visited 
                    auto pair_index= index_from_ordered_pair(
                        std::min(child_county, current_county), 
                        std::max(child_county, current_county), 
                        map_params.num_counties);
                    // if true it means these counties already have an index between them 
                    // and so this is illegal
                    if(counties_component_adj[pair_index]){
                        // REprintf("Caught!\n");
                        return false;
                    }else{
                        // else mark as true 
                        counties_component_adj[pair_index] = true;
                        // add this edge to the county graph 
                        county_graph[child_county].push_back(current_county);
                        county_graph[current_county].push_back(child_county);
                    } 
                }
            // add if not visited 
            if(!components_visited[child_vertex_id]){
                // mark as visited to avoid being added later  
                components_visited[child_vertex_id] = true;
                // REprintf("Adding (%u, %u)\n", child_region, std::get<2>(child_vertex));
                vertex_queue.push(child_vertex);
            }
        }

    }
    // now reuse this for visited county graph 
    std::fill(
        counties_component_adj.begin(), 
        counties_component_adj.end(), 
        false
    );
    // print graph 
    // for (size_t i = 0; i < map_params.num_counties; i++)
    // {
    //     REprintf("%u: ", i);
    //     for(const auto &el: county_graph[i]){
    //         REprintf("%u, ", el);
    //     }
    //     REprintf("\n");
    // }
    
    // return true;

    // now we check the graph for cycles 
    for (int u = 0; u < map_params.num_counties; u++)
    {
        // skip if we've already visited somehwere else 
        if(counties_component_adj[u]) continue;

        // else we look through this connected component
        std::queue<std::pair<int,int>> vertex_queue;
        vertex_queue.push({u, -1});

        while(!vertex_queue.empty()){
            
            auto current_entry = vertex_queue.front(); 
            vertex_queue.pop();
            int current_vertex = current_entry.first;
            int parent = current_entry.second;
            // Rprintf("v=%d, parent=%d\n",current_vertex, parent);
            // if we already visited this then that means there's a cycle
            if(counties_component_adj[current_vertex]){
                // Rprintf("Cycle!\n");
                return false;
            } 
            // else mark as visited 
            counties_component_adj[current_vertex] = true;
            // now we add all the children
            for(auto const &child_vertex: county_graph[current_vertex]){
                // ignore if the parent
                if(child_vertex == parent){
                    continue;
                }else{
                    vertex_queue.push({child_vertex, current_vertex});
                }
            }
        }
    }
    return true;

}


bool CountyComponents::check_valid_hiearchical_plan(Plan const &plan,
    Graph &county_graph){
    if(!counties_on) return true;
    // first we check if the number of splits and connected components is ok
    std::pair<bool, int> component_results = count_county_splits(plan);
    if(component_results.first || component_results.second > num_regions-1) return false;

    // if thats ok then build the component tree 
    build_component_graph_and_tree(plan.region_ids);

    // now check the restricted county graph of that 
    return check_is_county_component_multigraph_valid(county_graph);
}