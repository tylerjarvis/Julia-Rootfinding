using Pkg
using Test

# Activate the project environment
Pkg.activate(".") # MAKE SURE YOU ARE IN THE /Julia-Rootfinding folder when you run this script in the terminal
Pkg.instantiate()

print("Testing...")
include("ChebyshevApproximatorTest.jl")  # Adjust the path as needed

# Uncomment the lines below to run specific test sets
# test_transformPoints()
# test_getFinalDegree()
# test_startedConverging()
# test_checkConstantInDimension()
# test_hasConverged()
# test_createMeshgrid()
# test_getApproxError()
# test_intervalApproximateND()
# test_getChebyshevDegrees()
test_all()
