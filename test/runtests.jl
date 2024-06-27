using Pkg
using Test

# Activate the project environment
Pkg.activate(".") # MAKE SURE YOU ARE IN THE /Julia-Rootfinding folder when you run this script in the terminal
Pkg.instantiate()

print("Testing...")
include("ChebyshevApproximatorTest.jl")  
include("ChebyshevSubdivisionSolverTest.jl") # Adjust the path as needed

# Uncomment the lines below to run specific test sets
# 
# ============================================ All Tests ============================================
# test_all()
#
# ============================================ ChebyshevApproximator Tests ============================================
test_all_ChebyshevApproximator()
# 
# test_transformPoints()
# test_getFinalDegree()
# test_startedConverging()
# test_checkConstantInDimension()
# test_hasConverged()
# test_createMeshgrid()
# test_getApproxError()
# test_intervalApproximateND()
# test_getChebyshevDegrees()
# test_chebApproximate()
# 
# ============================================ ChebyshevSubdivisionSolver Tests ============================================
# test_all_ChebyshevSubdivisionSolver()
# 
# test_getLinearTerms()
# test_linearCheck1()


# SETUP TO TEST ALL
function test_all()
    test_all_ChebyshevApproximator()
    test_all_ChebyshevSubdivisionSolver()
end