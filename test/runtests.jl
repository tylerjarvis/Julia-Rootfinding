using Pkg
using Test

# Activate the project environment
Pkg.activate(".") # MAKE SURE YOU ARE IN THE /Julia-Rootfinding folder when you run this script in the terminal
Pkg.instantiate()

println("Testing...")
include("ChebyshevApproximatorTest.jl")
include("ChebyshevSubdivisionSolverTest.jl")
include("TrackedIntervalTest.jl")
include("QuadraticCheckTest.jl")
include("../../Julia-Rootfinding/src/StructsWithTheirFunctions/TrackedInterval.jl")
include("../../Julia-Rootfinding/src/StructsWithTheirFunctions/SolverOptions.jl")
function test_all()
    test_all_ChebyshevApproximator()
    test_all_ChebyshevSubdivisionSolver()
    test_all_TrackedInterval()
    test_all_QuadraticCheck()
end

# Uncomment the lines below to run specific test sets
# 
# ============================================ All Tests ============================================
# test_all()
#
# ============================================ ChebyshevApproximator Tests ============================================
# test_all_ChebyshevApproximator()
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
# test_reduceSolvedDim()
# test_transformChebInPlace1D()
# test_transformChebInPlaceND()
# test_getTransformationError()
# test_transformCheb()
# test_transformChebToInterval()
# test_getSubdivisionDims()
# test_getInverseOrder()
# test_getSubdivisionIntervals()
# test_boundingIntervalLinearSystem()
# test_zoomInOnIntervalIter()
# test_isExteriorInterval()
# test_trimMs()
test_solvePolyRecursive()
# test_solveChebyshevSubdivision()

# ============================================ TrackedInterval Tests ============================================

# test_all_TrackedInterval()

# test_copyInterval()
# test_addTransform()
# test_getIntervalForCombining()
# test_isPoint()
# test_getFinalInterval()
# test_getFinalPoint()
# test_overlapsWith()
# test_startFinalStep()

# ============================================ QuadraticCheck Tests ============================================

# test_all_QuadraticCheck()


# test_quadraticCheck2D()
# test_quadraticCheck3D()
# test_get_fixed_vars()
# test_quadraticCheckND()
# test_quadraticCheck()
