using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")

# Get the list of all files in the directory
cheb_result_files_paths = readdir("ChebfunTests/Chebfun_results/")
YRoots_result_files_paths = readdir("ChebfunTests/YRoots_results/")

# Include each file
for path in cheb_result_files_paths
    include("Chebfun_results/" * path)
end

for path in YRoots_result_files_paths
    include("YRoots_results/" * path)
end

using CSV
using DataFrames

# Function to load data from a CSV file
function load_data(filename)
    df = CSV.read(filename, DataFrame)
    return df
end

# Call the function to load the CSV data
df = load_data("data.csv")

# Display the DataFrame
println(df)
