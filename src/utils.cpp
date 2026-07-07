// Generic helper functions 
#include "utils.h"





Rcpp::List cli_config(bool clear, const char *fmt) {
    return Rcpp::List::create(
        Rcpp::_["clear"] = clear, 
        Rcpp::_["show_after"] = 0.25, 
        Rcpp::_["format"] = fmt
    );
}
