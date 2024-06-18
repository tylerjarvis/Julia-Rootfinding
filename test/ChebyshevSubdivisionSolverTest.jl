include("../../Julia-Rootfinding/src/ChebyshevSubdivisionSolver.jl")
using Test

function test_all_ChebyshevSubdivisionSolver()
    @testset "All tests in ChebyshevSubdivisionSolverTest.jl" begin
        test_getLinearTerms()
    end
end

function test_getLinearTerms()
    @testset "getLinearTerms unit tests" begin
        M_1 = reshape(collect(0:15),(2,2,2,2))
        A_1_expected = [8; 4; 2; 1]
        @test isapprox(A_1_expected,getLinearTerms(M_1))

        M_2 = reshape(collect(0:15),(2,2,4))
        A_2_expected = [4, 2, 1]
        @test isapprox(A_2_expected,getLinearTerms(M_2))

        M_3 = [43.2;12.2;-9.2]
        A_3_expected = [12.2]
        @test isapprox(A_3_expected,getLinearTerms(M_3))
        
        M_4 = reshape(collect(0:(2*15*4*6-1)),(6,4,15,2))
        A_4_expected = [360, 24, 6, 1]
        @test isapprox(A_4_expected,getLinearTerms(M_4))
    end
end