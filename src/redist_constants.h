
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
    inline constexpr bool track_granular_times = true; 
    // Turns on granular tracking of each part of SMC algorithm
    // No solid estimates of performance hit but could cause anywhere from 
    // .1 to 5% slowdown depending on how expensive a call to CPU clock is versus
    // the SMC functions 
    inline constexpr bool supposedly_safe_input_checks = false; 
    // Performs checks that shouldn't be neccessary if the input is as expected
    // Helpful to turn on when checking for bugs as it will help catch problems 
    // where something is corrupted somewhere without segfaulting and it gets passed on
    // to later functions. Helps catch weird things sooner
    inline constexpr bool check_threadpool_integrity = false;
    // Will perform extra checking on RcppThread::Pool object calls like ensuring 
    // thread ids are generated properly. 
    inline constexpr bool bounds_checking = false; 
    // Performs bounds checking that should be unnccesary, useful when debugging
    // to spot bounds errors as soon as they appear. 
    inline constexpr bool object_integrity_checking = false; 
    // Performs integrity checks on key objects in the code 
    // Can be very expensive but good for spotting strange errors
    // Currently checks the following types of objects
    //  - Spanning trees returned from Wilson's algorithm
    //  - Integrity of plan objects
    inline constexpr bool redundancy_checks = false; 
    // Performs checks in a function that should be redundant if the function is 
    // doing what its supposed to. e.g., imagine a function should walk through every
    // vertex in a county and we know the size of the county. Then this would turn on
    // a check that all county vertices were indeed visitied
    inline constexpr bool special_timing = false;
    // For performance profiling of specific code. This uses a global SpecialTimes
    // struct and can be dropped into any function. It assumes you are single threading.
    // Simple struct for tracking granular time
    struct SpecialTimes {
        double sample_sub_ust_time = 0.0;
        int sub_ust_num_calls = 0;
        double other_stuff_time = 0.0;
    };

    inline SpecialTimes SPECIAL_TIMES;

}


using Clock = std::chrono::steady_clock;


// Debugging constants
constexpr bool TREE_SPLITTING_DEBUG_VERBOSE =
    false;                                    //  Turns on verbose debugging for tree stuff
constexpr bool WEIGHTS_DEBUG_VERBOSE = false; //  Turns on verbose debugging for weight stuff



#endif