# GENERAL INFORMATION

The hope of this document is manage the list of tasks to do and also keep a record of what was done. 

# ---- IMPORTANT TO KEEP IN MIND -----
Note all of the diagnostic information is accurately updated to account for a final resampling. When the weights are really low variance it doesn't really matter but just take that with a grain of salt for everything. 


# ---- Possible Bugs and Other Issues -----

**Random seed bug** 
Since `state_xo` and `state_sr` are global variables when you run without multiprocessing this seems to cause issues (validation sometimes fails.) This needs to be addressed as when running on the cluster multiprocessing is very inefficient (since memory is more expensive than cpu cores). 


# ----- ACTIVE TASKS -----
 
**Use Sparse Matrix LDL for log tau**

Use this: https://github.com/samuel-watson/SparseChol for computing log tau terms since graph laplacians are symmetric and 
sparse for large number of districts 

**Do Resampling in c++**
For `run_gsmc_plans` make the resampling happen in c++ instead of R. This will save on memory overhead since you don't need to do reindexing in R and can use the dummy matrix for space. 

**Make immutable class members constant**
For all classes where a variable doesn't change after initialization (think things like `ndists`, `V`, etc.) make them constant. This might help a little bit with thread safety but the main reason is to avoid modifying things that shouldn't be changed in code. 

**Create new SpanningTree class**
Create a new `SpanningTree` class that is a wrapper for the directed spanning tree itself along with the vectors that always
have to be created for it. Namely the `visited, ignore, parent, pop_below` vectors. Right now we allocate `parent, pop_below` vectors everytime and I imagine that is making things expensive. 

**Make get_all_valid_edges_in_directed_tree work directly on tree**
Right now the function to get all edges in a directed tree requires a population above and parent vector but I think this should be possible to be run on the tree itself starting at the root without needing to compute the parent or pop above vectors. 

**Fix merge split bug**
Right now when running smc+merge split the maps somehow become disconnected for Oregon with the county constraint turned on. Investigate this further, as I'm not sure what is causing the bug. 

**Make tree_dev Support Generalized Splits**
Right now `tree_dev` only works for 1-district splits but it can easily be generalized to support arbitrary region splits.

**Create Diagnostic Levels**
Make it so the code supports multiple different diagnostic levels. I am thinking 
    - Level 0: Bare Minimum information needed so
        - Log weight standard deviations for each step
        - Acceptance rates for each step 
        - Number of unique parents at each step
        - The k values used
        - Effective sample size 
        - Step type (smc or merge split)
    - Level 1: (Default) 
        - Log weights for every step
        - Draw tries matrix for every step
        - Parent index Matrix

    - Level 3: (Very memory intensive)
        - The plan information (vertex region ids and dvals) for every step


**Pass Integer results back in Rcpp Integer data types**
To make things even more memory efficient change it so all integer data is passed back as Rcpp Integer vectors or matrices. This ensures it is returned to R as integer data which takes up half the space of numeric. 


**Make Parent Tries Atomic**
For the parent tries stuff need to change it to a vector of atomic integers since current implementation is not thread safe.
Code is 
```
std::vector<std::atomic<uint>> new_parent_unsuccessful_tries_mat(M);

// Initialize the vector elements
for (size_t i = 0; i < M; ++i) {
    new_parent_unsuccessful_tries_mat.at(i).store(static_cast<uint>(0)); // Explicitly store values
}
```

**Do Resampling in c++**
Make it so the resampling step can be done in c++, not R.

**Fix merge split on district only splits**
Right now I think there's a bug where if you do merge split after the final smc step (so N districts) the district split only version still uses the remainder meaning it only merges regions adjacent to the remainder. Since this is unlabeled in that case it should just pick an arbitrary pair

**Change diagnostics.R to be better**
Make the display for the new algorithm types better, also figure out a better way to calculate Rhats. For the ones that are by district it should compute for all districts and display the maximum. It also needs to have some way to store these results maybe .

**Change redist_plans storage scheme**
Right now I lump everything into diagnostics but there is probably a smarter way to do this. There should be another things called like internal diagnostics or something that only contains things saved during diagnostic mode.

**Consolidate Redist gsmc R files and c++ code**
Right now the redist_gsmc and redist_gsmc with merge split are in seperate c++ and R files but this seems unneccesary. It should be possible to combine them and just add a flag for whether or not merge split steps need to be run. 

**Shrink Incremental Weight Mat Size and other related variables**
For the MH step the weights don't actually change so we don't need to store the output for every step, just every smc step. That means we also don't need to track the effective sample size.

**Remove unneccesary variables when not in diagnostics**
In order to save memory make it so that some things like the log incremental weights matrix are not created or are not big when its not diagnostic mode. This will be really important for saving memory. 

Also need to decide if things like the original ancestors matrix should be created when its not diagnostic mode. Probably not but we will see.



**Create Function for Diagnostics**
Since both the smc and smc with merge split share some diagnostics its probably better to write a function to create the vectors shared between them to avoid duplicate code.


**Create Diagnostics Object for cpp code**
To avoid needing to pass so many different parameters in the function header consider creating a new diagnostic class in cpp which serves as a wrapper for all the diagnostic things that are collected in the `split_maps` function. 



**Create More Internal Visualization/Examination Stuff**
Already created a new file, `splitting_inspection.cpp` to help with this but general idea is create code that allows someone in R to pass in a graph in adjacency list form (and potentially other stuff) and see things like the cut tree output by the splitting procedure or an internal region level graph. This will be helpful for making figures for the paper and also for deep level debugging/diagnostic things.  

**Deal with final sample vs diagnostics**
Need to more cleanly seperate diagnostic information from stuff in the final sample. This is because the resampling step if done changes things slightly and just re-indexing based on the final re-sample might destroy diagnostic information. 




**For Previous tries track previous index not the current new particle**
- `draw_tries_vec` is the one. Make one that tracks the previous sampled index

**Ask why log compactness takes district matrix despite only using column**
The compute log tau term here: https://github.com/alarm-redist/redistmetrics/blob/5f7b36d8a7f9c7bc3c9098a7b6c6aa561d9074c7/inst/include/kirchhoff_inline.h#L16 takes a matrix of district values but only actually uses a single column. That can probably be changed.

# ----- GENERAL THOUGHTS -----

- For trying to predict if something can never be split
    - Sort population in order and see if gap is too large. This doesn't check every case but removes some
    - Partition split problem 

# ----- COMPLETED TASKS -----

**Add reordering Function - Done December 2024**
Task: Add a function that reorders the regions in a `Plan` object to be in order of most recently split. 

Comments after completion: Right now this function does a redundant copy. In the future think about adding a flag for whether or not the dummy plan needs to be copied first. 

**Pass results back in arma data types - DONE 11/8/2024**
Task: To make things more memory efficient change it so all results are passed back as arma vectors, matrices, or cubes for things that are multiple matrices. 

Comments after completion: Realized that all arma matrix types (even unsigned integers) are passed back to R with numeric data. So to make things as memory efficient as possible need to pass things back as Rcpp::IntegerMatrix or Rcpp::IntegerVector whenever possible as integer data in R takes up about half the space of numeric. Not a big deal for most applications but for the FAS cluster this will be a big help as you get heavily penalized for requesting more memory. 

**Create MH Ratio Calculator - DONE 11/6/2024**
Task: Create a function (along with helpers as needed) that computes the MH ratio for valid new proposed plans.

Comments after completion: Initially got the probability of merging pairs flipped in the ratio but fixed it.

**Create Full Pass through - DONE 11/2/2024**
Task: In the `smc_and_mcmc.cpp` file make a function analagous to `optimal_gsmc_plans` that runs it for the whole thing. For now just stick with doing MCMC moves with a set frequency with the amount set to 1/acceptance rate. Need to think about what kind of diagnostics to keep. Right now I think it should just be the number of attempts made and how many were successful along with storing the plans and weights at the end 

Comments after completion: None


**Make get_edge_to_cut only need vertex region id vector* - DONE 11/2/2024**
Task: Change `get_edge_to_cut` to only take in a vector of vertex region ids and a max dval to try instead of taking in a plan. This will make merge split easier by allowing you to just pass in a vertex region id vector with two regions merged without having to bother editing the actual plan.

Comments after completion: Successfully changed that for merge split stuff.

**Seperate New Region ID Creation from Update Edges - DONE 11/1/2024**
Task: To make code more resuable for merge split stuff make it so that creating IDs/resizing attributes in the plan is done outside the `update_plan_from_cut` function and instead that function only handles updating the plan. 

Comments after completion: Changed `update_plan_from_cut` so it just does updating. Resizing plan attributes and coming up with region ids now happen outside the function so it should be fully usable for merge split stuff later on.

**Remove Plan String Label and Replace with Split Order Number - DONE 11/1/2024**
Task: Remove the string label attribute from plan objects and instead just have a vector which tracks the relative order regions were added. This will be a vector mapping region id to a number where the order of the number relative to the others indicates the order it was added. Will also need to add a region_order_max attribute so each time a new region is added it is set to that plus one. At the end you can use this to recover the order districts were created but it can also be used to inspect intermediate results as well. 

Comments after completion: Successfully redesigned cpp and R code so it is now possible to track the relative order regions were added and the R function now automatically makes it so that region ID is numbered oldest to most recent added.

**Track Remainder Better - DONE 11/12024**
Task: Right now the remainder region in a plan is blindly set to region2_id but this will be a problem for the merge split because the remainder might get changed. 

Comments after completion: Added a new option to the `Plan` constructor for specifying whether the plan is for district only splits or not. If its district only splits remainder is initialized to 0 and if not set to -1. That should mean if any future code tries to use remainder region when it shouldn't that will throw an error.

**Rewrite gsmc and basic_smc into one function - DONE 11/1/2024**
Task: Use the code from the new `smc_and_mcmc.cpp` file to rewrite the gsmc and basic_smc stuff. It should be possible to make one version of most of the functions and just have a single input control whether or not it just splits districts. This should cut down on the code and eliminate the big redundancy problem. 

Comments after completion: The functionality previously in `gsmc.cpp` and `basic_smc.cpp` has now been moved to `optimal_gsmc.cpp` and some of the functions in there (namely the splitting and weight ones) have been moved to their own files. There is now a general purpose function that handles this with a new flag to control whether or not its doing generalized splits or one-district splits only. I also deleted the `redist_gsmc.R` and `redist_basic_smc.R` and consolidated things to the `redist_optimal_gmsc.R` which also now has a flag for generalized vs one district splits.

**Create Separate Updating Function - DONE 11/1/2024**
Create a function which takes a cut tree and information on the new regions and updates them. Previously this was in the `cut_regions` function but it has now been seperated out since sometimes for the MCMC step we want to know what a cut region would look like but don't neccesarily always want to actually change it. 


**Compress dval storage in `plan` object - DONE FORGOT TO RECORD DATE**
-   We don't actually need to store an `V` length vector of the dvals associated with each vertex. Instead all you need to do is keep a vector mapping region id values to the dvalue 
-   Change the attribute in `Plan` linking region id to d value from a `map` to just an array since the regions are already 0 indexed we don't need to bother with a map
-   In the final output then the dval vectors don't need to be length `V` but instead just the number of regions 


**Add population tempering - DONE 9/11/2024 **
- Its the `log_temper` stuff right now. 
    - It can be added to every term except subtract from the final one 
    - Just create a seperate weights computation function and have it done in there
    - Remember the value must be calculated for each of the three regions
