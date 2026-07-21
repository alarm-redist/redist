#ifndef THREADING_HELPERS_H
#define THREADING_HELPERS_H



#include <RcppThread.h>


/*
 * Creates a Rcpp Threadpool
 *
 *
 */
RcppThread::ThreadPool get_thread_pool(int const num_threads);

// Returns the number of threads being used. 
// When single threading it returns as one thread
int get_num_threads(const RcppThread::ThreadPool &pool);





#endif
