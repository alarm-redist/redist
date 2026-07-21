// Generic helper functions 
#include "utils.h"


Rcpp::List cli_config(bool clear, const char *fmt) {
    return Rcpp::List::create(
        Rcpp::_["clear"] = clear, 
        Rcpp::_["show_after"] = 0.25, 
        Rcpp::_["format"] = fmt
    );
}


void tree_size_check(
    MapParams const &map_params, 
    int const proposed_tree_size, 
    std::vector<int> const &tree_sizes,
    std::string_view const where
){
    if (proposed_tree_size <= 0 ||
        proposed_tree_size > static_cast<int>(tree_sizes.size())
    ) {
        std::ostringstream oss;
        oss << "Invalid tree size.\n" << where << "\n"; 
        oss << "proposed_tree_size=" << proposed_tree_size << "\n";
        oss << "tree sizes vector.size()=" << tree_sizes.size() << "\n";
        throw std::runtime_error(oss.str());
    }
}