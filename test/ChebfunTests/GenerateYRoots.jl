using Pkg
include("ChebfunVars.jl")
include("../../../Julia-Rootfinding/src/CombinedSolver.jl")

using DelimitedFiles

# fs and gs are the input functions
# as are the lower bounds of the interval, bs are the upper bounds.
# file_names are the names for the output CSV files.
for (f, g, a, b, file_name) in zip(fs, gs, as, bs, file_names)
    roots = solve([f, g], a, b) # returns a vector of vectors of floats
    
    # Save the results to a CSV file without column names
    writedlm("test/ChebfunTests/YRoots_results/" * file_name * ".csv", roots, ',')
end



