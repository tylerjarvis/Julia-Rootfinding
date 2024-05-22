using Test
using HasConverged: hasConverged

function test_has_converged()
    # Test case 1: Large tolerance, converged
    coeff = [0.1, 0.0000001, -0.3]
    coeff2 = [0.11, -0.01, -0.31]
    tol = 0.01001
    @test hasConverged(coeff, coeff2, tol) == true

    # Test case 2: Small tolerance, not converged
    coeff = [0.1, 0.2, 0.3]
    coeff2 = [0.11, 0.21, 0.31]
    tol = 0.001
    @test hasConverged(coeff, coeff2, tol) == false

    # Test case 3: Zero tolerance, not converged
    coeff = [0.1, 0.2, 0.3]
    coeff2 = [0.1000000001, 0.2, 0.3]
    tol = 0
    @test hasConverged(coeff, coeff2, tol) == false

    # Test case 4: Large number of inputs, converged
    coeff = [i * 0.1 for i in 0:9]
    coeff2 = [(i * 0.1) + 0.0001 for i in 0:9]
    tol = 0.001
    @test hasConverged(coeff, coeff2, tol) == true
end

# Run the tests
@testset "test_has_converged" begin
    test_has_converged()
end

println("All tests passed!")
