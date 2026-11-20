#ifndef THREADING_HELPERS_H
#define THREADING_HELPERS_H



#include <RcppThread.h>


/*
 * Threading convention used throughout the package:
 *
 *   num_threads == 1  run serially, with no worker threads
 *   num_threads  > 1  use exactly that many threads
 *   num_threads <= 0  use every available core
 *
 * Creates a Rcpp Threadpool following that convention.
 */
RcppThread::ThreadPool get_thread_pool(int const num_threads);

/*
 * The same convention, for the free `RcppThread::parallelFor` rather than a
 * pool. Note that RcppThread itself reads 0 as "no worker threads", which is
 * the opposite of what 0 means in our own arguments, so this translation is
 * needed wherever a user supplied thread count reaches parallelFor.
 */
int resolve_num_threads(int const num_threads);

// Returns the number of threads being used. 
// When single threading it returns as one thread
int get_num_threads(const RcppThread::ThreadPool &pool);





#endif
