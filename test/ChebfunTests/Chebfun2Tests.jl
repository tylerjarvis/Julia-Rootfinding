include("ChebfunVars.jl")

# Get the list of all files in the directory
result_files_paths = readdir("ChebfunTests/Chebfun_results/")

# Include each file
for path in result_files_paths
    include("Chebfun_results/" * path)
end

print("Hello")

