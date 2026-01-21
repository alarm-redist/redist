# C++ Redist Code Overview 
# Author: Philip O'Sullivan

## Description
This document provides an overview of the c++ written or edited by Philip O'Sullivan from 2024-2025. As of writing it is probably about half of the c++ code.


# Important Naming convenctions 

    - `ndists` from the R code still means the number of districts
    - `total_seats` is equivalent to `nseats` in the R code. IE `total_seats` is the total number of seats in a map.


# Classes
Most of the code I wrote tries to use object oriented programming to manage things.
    - `MapParams` Class - This class is largely just a wrapper for all the important map related variables that frequently get passed around to functions. By storing it all in this class it helps limit inputs needed to functions as well as store useful information like the number of counties and help avoid needing to recompute it. 
    - `TreeSplitter` Class - This is the abstract base class used for choosing how to split an edge in a spanning tree. Specifically, given a vector of potential `EdgeCuts` it attempts to select one according to some rule. This abstract class makes it easy to add new forward kernels. 
    - `Plan` Class - This is the abstract base class used for representing a plan. At its core it is a lightweight wrapper for a vector mapping vertex ids to the region they are assigned to. For each of the different sampling spaces there are concrete derived classes. In general all the splitting and calculation code is designed to work on an abstract `Plan` object and the sampling space specifics 

    The derived classes are 
        - `GraphPlan` 
        - `ForestPlan`
        - `LinkingEdgePlan`

    This class has methods that require the following classes as inputs 
        - `MapParams`
        - `USTSampler`
        - `PlanMultigraph`
        - `TreeSplitter`