
// Given a map and a plan, find precincts in dist_1 that border dist_2

// [[Rcpp::export]]
LogicalVector walnuts_find_boundary_prec(List map, IntegerVector plan, int dist_1, int dist_2, int n_rows) {
    LogicalVector boundary(n_rows, false);

    for (int i = 0; i < n_rows; i++) {
        IntegerVector adj = map[i];
        int adj_length = adj.size();
        if (plan[i] == dist_1) {
            for (int j = 0; j < adj_length; j++) {
                if (plan[adj[j]] == dist_2) {
                    boundary[i] = true;
                    break;
                }
            }
        }
    }

    return boundary;
}

// Given a map and a plan, find precincts in dist_1 that border dist_2 that are in the gpp

// [[Rcpp::export]]
LogicalVector walnuts_find_boundary_blk(List map, IntegerVector plan, int dist_1, int dist_2, int n_rows, CharacterVector geoids, std::string gpp) {
    LogicalVector boundary(n_rows, false);

    for (int i = 0; i < n_rows; i++) {
        IntegerVector adj = map[i];
        int adj_length = adj.size();
        if (as<std::string>(geoids[i]) == gpp && plan[i] == dist_1) {
            for (int j = 0; j < adj_length; j++) {
                if (plan[adj[j]] == dist_2) {
                    boundary[i] = true;
                    break;
                }
            }
        }
    }

    return boundary;
}
