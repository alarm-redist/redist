
#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

// Important Constants

constexpr bool DEBUG_VERBOSE = false; // Compile-time constant
constexpr double SELECTION_ALPHA = 0.6321206;

// Performance configurations. These can be turned on locally for more granular tracking
// Useful for anyone forking the code or doing in depth performance analysis
namespace perf_config {
    inline constexpr bool track_granular_times = false; 
    // Turns on granular tracking of each part of SMC algorithm
    // No solid estimates of performance hit but could cause anywhere from 
    // .1 to 5% slowdown depending on how expensive a call to CPU clock is versus
    // the SMC functions    
}


using Clock = std::chrono::steady_clock;



// Debugging constants
constexpr bool TREE_SPLITTING_DEBUG_VERBOSE =
    false;                                    //  Turns on verbose debugging for tree stuff
constexpr bool WEIGHTS_DEBUG_VERBOSE = false; //  Turns on verbose debugging for weight stuff

#endif