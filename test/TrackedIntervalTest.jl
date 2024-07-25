include("../../Julia-Rootfinding/src/StructsWithTheirFunctions/TrackedInterval.jl")
using Test

function tests_all_TrackedInterval()
    @testset "All tests in TrackedIntervalTest.jl" begin
        test_intervalCopy()
        test_addTransform()
    end
end

function test_intervalCopy()
    @testset "intervalCopy unit tests" begin
        trackedInterval_1 = TrackedInterval([1;2;;-1;1;;-3;2])
        copiedInterval_1 = intervalCopy(trackedInterval_1)
        # BASIC INTERVAL COPYING
        @test copiedInterval_1.interval == trackedInterval_1.interval
        @test copiedInterval_1.empty == false
        @test copiedInterval_1.nextTransformPoints == fill(0.0394555475981047,3)
        trackedInterval_1.empty = true
        trackedInterval_1.transforms = [[2;3;;4;5;;1.234234;1.3212331]]
        trackedInterval_1.nextTransformPoints[2] = 4
        copiedInterval_2 = intervalCopy(trackedInterval_1)
        # TEST OTHER ATTRIBUTES COPIED
        @test copiedInterval_2.empty == true
        @test copiedInterval_2.transforms == [[2;3;;4;5;;1.234234;1.3212331]]
        @test copiedInterval_2.nextTransformPoints[1] == 0.0394555475981047
        @test copiedInterval_2.nextTransformPoints[2] == 4

        # FINAL STEP TEST
        @test copiedInterval_2.finalStep == false
        @test copiedInterval_2.possibleDuplicateRoots == []
        @test copiedInterval_2.possibleExtraRoot == false
        @test copiedInterval_2.preFinalInterval == []
        @test copiedInterval_2.preFinalTransforms == []
        trackedInterval_1.finalStep = true
        trackedInterval_1.possibleDuplicateRoots = [3,1,2]
        trackedInterval_1.possibleExtraRoot = true
        trackedInterval_1.preFinalInterval = [1;2;;-1;1;;-3;2]
        trackedInterval_1.preFinalTransforms = [[1;2;;-1;1;;-3;2]]
        copiedInterval_3 = intervalCopy(trackedInterval_1)
        @test copiedInterval_3.finalStep == true
        @test copiedInterval_3.possibleDuplicateRoots == [3,1,2]
        @test copiedInterval_3.possibleExtraRoot == true
        @test copiedInterval_3.preFinalInterval == [1;2;;-1;1;;-3;2]
        @test copiedInterval_3.preFinalTransforms == [[1;2;;-1;1;;-3;2]]


    end
end

function test_addTransform()
    @testset "addTransform unit tests" begin
        # Empty interval test 1
        trackedInterval_1 = TrackedInterval([-1.;2;;-1;1])
        trackedInterval_1.finalStep = true
        trackedInterval_1.canThrowOutFinalStep = true
        subInterval_1 = [0;-.0001;;-1;1]
        addTransform(trackedInterval_1,subInterval_1)
        @test trackedInterval_1.empty == true
        @test trackedInterval_1.transforms == []
        @test trackedInterval_1.interval == [-1.;2;;-1;1]

        # Empty interval test 2
        trackedInterval_2 = TrackedInterval([-1;2;;-1;1])
        subInterval_2 = [0;-.0001;;-1;1]
        addTransform(trackedInterval_2,subInterval_2)
        @test trackedInterval_2.empty == true
        @test trackedInterval_2.transforms == []
        @test trackedInterval_2.interval == [-1.;2;;-1;1]

        trackedInterval_3 = TrackedInterval([-1.;2;;-1;1])
        trackedInterval_3.finalStep = true
        subInterval_3 = [0;-.0001;;-1;1]
        addTransform(trackedInterval_3,subInterval_3)
        @test trackedInterval_3.empty == false
        @test isapprox(trackedInterval_3.interval,[ 0.5;  0.5;;-1.;   1. ])
    

    end
end