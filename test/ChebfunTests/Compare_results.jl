# Get the list of all files in the directory
cheb_result_files_paths = readdir("test/ChebfunTests/Chebfun_results/")
YRoots_result_files_paths = readdir("test/ChebfunTests/YRoots_results/")

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
    return CSV.read(filename, DataFrame)
end

# Function to compare two DataFrames within a given tolerance
function compare_dataframes(file1, file2; tolerance=1e-9)
    # Load the data from the CSV files
    df1 = load_data(file1)
    df2 = load_data(file2)

    # Check if DataFrames have the same structure
    if size(df1) != size(df2)
        println("DataFrames have different shapes.")
        return false
    end

    # Convert DataFrames to matrices
    matrix1 = Matrix(df1)
    matrix2 = Matrix(df2)

    # Compare element-wise with tolerance
    for i in 1:size(matrix1, 1)
        for j in 1:size(matrix1, 2)
            if abs(matrix1[i, j] - matrix2[i, j]) > tolerance
                println("Difference at position ($i, $j) is greater than tolerance: ", abs(matrix1[i, j] - matrix2[i, j]))
                return false
            end
        end
    end

    return true
end

all_match = true
for (chebFun_file, YRoots_file) in zip(sort!(cheb_result_files_paths), sort!(YRoots_result_files_paths))
    # println("comparing files: " * chebFun_file * " \t and \t " * YRoots_file)
    if compare_dataframes("test/ChebfunTests/Chebfun_results/" * chebFun_file, "test/ChebfunTests/YRoots_results/" * YRoots_file; tolerance=1e-9)
        continue
    else
        println("This discrepancy is in files: " * chebFun_file * " \t and \t " * YRoots_file * "\n")
        global all_match = false
    end
end

if all_match
    println("All results match!")
end