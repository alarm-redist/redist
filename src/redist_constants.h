
#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <chrono>

// Important Constants
constexpr bool DEBUG_VERBOSE = false; // Compile-time constant
constexpr double SELECTION_ALPHA = 0.6321206;

// Performance configurations. These can be turned on locally for more granular tracking
// and to help with debugging by performing potentially onerous extra checks
// Useful for anyone forking the code or doing in depth performance analysis
namespace perf_config {
    inline constexpr bool track_granular_times = false; 
    // Turns on granular tracking of each part of SMC algorithm
    // No solid estimates of performance hit but could cause anywhere from 
    // .1 to 5% slowdown depending on how expensive a call to CPU clock is versus
    // the SMC functions    
    inline constexpr bool supposedly_safe_input_checks = true; 
    // Turns on input checking for functions that should theoretically be unneccesary. 
    // Helpful to turn on when checking for bugs as it will help catch weird things sooner
    inline constexpr bool check_threadpool_integrity = true;
    // Will perform extra checking on RcppThread::Pool object calls like ensuring 
    // thread ids are generated properly. 
    inline constexpr bool bounds_checking = true; 
    // Performs bounds checking that should be unnccesary, useful when debugging
    // to spot bounds errors as soon as they appear. 
    inline constexpr bool object_integrity_checking = true; 
    // Performs integrity checks on key objects in the code 
    // Can be very expensive but good for spotting strange errors
    // Currently checks the following types of objects
    //  - Spanning trees returned from Wilson's algorithm
    //  - Integrity of plan objects 
}


using Clock = std::chrono::steady_clock;



// Debugging constants
constexpr bool TREE_SPLITTING_DEBUG_VERBOSE =
    false;                                    //  Turns on verbose debugging for tree stuff
constexpr bool WEIGHTS_DEBUG_VERBOSE = false; //  Turns on verbose debugging for weight stuff

#endif