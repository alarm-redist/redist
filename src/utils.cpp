// Generic helper functions 
#include "utils.h"


Rcpp::List cli_config(bool clear, const char *fmt) {
    return Rcpp::List::create(
        Rcpp::_["clear"] = clear, 
        Rcpp::_["show_after"] = 0.25, 
        Rcpp::_["format"] = fmt
    );
}



RcppThread::ThreadPool get_thread_pool(int const num_threads) {
    if (num_threads == 1) {
        return RcppThread::ThreadPool(0);
    } else if (num_threads > 1) {
        return RcppThread::ThreadPool(num_threads);
    } else {
        return RcppThread::ThreadPool(std::thread::hardware_concurrency());
    }
}


int get_num_threads(const RcppThread::ThreadPool &pool){
    auto const pool_threads = pool.getNumThreads();
    if (pool_threads <= 1) {
        return 1;
    } else {
        return pool_threads;
    }
}