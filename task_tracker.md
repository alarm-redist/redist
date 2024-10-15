# GENERAL INFORMATION

The hope of this document is manage the list of tasks to do and also keep a record of what was done. 

# ---- IMPORTANT TO KEEP IN MIND -----
Note all of the diagnostic information is accurately updated to account for a resampling. When the weights are really low variance it doesn't really matter but just take that with a grain of salt for everything. 


# ----- ACTIVE TASKS -----

**Create More Internal Visualization/Examination Stuff**
Already created a new file, `splitting_inspection.cpp` to help with this but general idea is create code that allows someone in R to pass in a graph in adjacency list form (and potentially other stuff) and see things like the cut tree output by the splitting procedure or an internal region level graph. This will be helpful for making figures for the paper and also for deep level debugging/diagnostic things.  

**Deal with final sample vs diagnostics**
Need to more cleanly seperate diagnostic information from stuff in the final sample. This is because the resampling step if done changes things slightly and just re-indexing based on the final re-sample might destroy diagnostic information. 

**Try storing output as arma matrix instead of c++**
Instead of making vectors of vectors or vector of vectors of vectors try just making `arma::Matrix` instead. This should help with
optionally setting the storage size 

**Compress dval storage in `plan` object:**
-   We don't actually need to store an `V` length vector of the dvals associated with each vertex. Instead all you need to do is keep a vector mapping region id values to the dvalue 
-   Change the attribute in `Plan` linking region id to d value from a `map` to just an array since the regions are already 0 indexed we don't need to bother with a map
-   In the final output then the dval vectors don't need to be length `V` but instead just the number of regions 


**For Previous tries track previous index not the current new particle**
- `draw_tries_vec` is the one. Make one that tracks the previous sampled index



**Add optimal weights to current smc**
- `get_wgts`

# ----- GENERAL THOUGHTS -----

- For trying to predict if something can never be split
    - Sort population in order and see if gap is too large. This doesn't check every case but removes some
    - Partition split problem 

# ----- COMPLETED TASKS -----

**Add population tempering - DONE 9/11/2024 **
- Its the `log_temper` stuff right now. 
    - It can be added to every term except subtract from the final one 
    - Just create a seperate weights computation function and have it done in there
    - Remember the value must be calculated for each of the three regions
