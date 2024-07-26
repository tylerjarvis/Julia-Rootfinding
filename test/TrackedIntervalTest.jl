include("../../Julia-Rootfinding/src/StructsWithTheirFunctions/TrackedInterval.jl")
using Test

function test_all_TrackedInterval()
    @testset "All tests in TrackedIntervalTest.jl" begin
        test_intervalCopy()
        test_addTransform()
        test_getIntervalForCombining()
        test_isPoint()
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
    
        trackedInterval_4 = TrackedInterval([-1.;2;;-1;1])
        trackedInterval_4.finalStep = true
        subInterval_4 = [.0001;0;;-1;1]
        addTransform(trackedInterval_4,subInterval_4)
        @test trackedInterval_4.empty == false
        @test isapprox(trackedInterval_4.interval,[ 0.50015;  0.50015;;-1.;       1.     ])
        @test isapprox(trackedInterval_4.transforms[1],[0.e+00; 1.e+00;;1.e-04; 0.e+00])

        trackedInterval_5 = TrackedInterval([-1.;2;;-1;1;;-1;.9])
        trackedInterval_5.finalStep = true
        subInterval_5 = [.0001;0;;-1;1;;-1;1]
        addTransform(trackedInterval_5,subInterval_5)
        @test trackedInterval_5.empty == false
        @test isapprox(trackedInterval_5.interval,[ 0.50015;  0.50015;;-1.;       1.     ;;-1.;       0.9    ])
        @test isapprox(trackedInterval_5.transforms[1],[0.e+00; 1.e+00; 1.e+00;;1.e-04; 0.e+00; 0.e+00])

        trackedInterval_6 = TrackedInterval([-1.;2;;-1;1;;-1;.9])
        trackedInterval_6.finalStep = true
        subInterval_6 = [-18;-18.01;;-3;1;;-3;1]
        addTransform(trackedInterval_6,subInterval_6)
        @test trackedInterval_6.empty == false
        @test isapprox(trackedInterval_6.interval,[-1.;  -1. ;;-1.;   1. ;;-1.;   0.9])
        @test isapprox(trackedInterval_6.transforms[1],[ 0.;  1.;  1.;;-1.;  0.;  0.])

        trackedInterval_7 = TrackedInterval([-1.;2;;-1;1;;-1;.9])
        trackedInterval_7.finalStep = true
        subInterval_7 = [3;-3.01;;3;1;;3;1]
        addTransform(trackedInterval_7,subInterval_7)
        @test trackedInterval_7.empty == false
        @test isapprox(trackedInterval_7.interval,[2.;  2. ;;1.;  1. ;;0.9; 0.9])
        @test isapprox(trackedInterval_7.transforms[1],[0.; 0.; 0.;;1.; 1.; 1.])

    end
end

function test_getIntervalForCombining()
    @testset "getIntervalForCombining unit tests" begin
        trackedInterval = TrackedInterval([-1;1;;-1.2332;1.2134;;-5;1])
        @test isapprox(getIntervalForCombining(trackedInterval),[-1;1;;-1.2332;1.2134;;-5;1])
        @test trackedInterval.preFinalInterval == []
        trackedInterval.preFinalInterval = [-5;5;;-1;1;;500;1]
        trackedInterval.finalStep = true
        @test isapprox(getIntervalForCombining(trackedInterval),[-5;5;;-1;1;;500;1])
        trackedInterval.finalStep = false
        @test isapprox(getIntervalForCombining(trackedInterval),[-1;1;;-1.2332;1.2134;;-5;1])
    end
end

function test_isPoint()
    @testset "test_isPoint unit tests" begin
        trackedInterval_1 = TrackedInterval([1;2;;3;4])
        @test !isPoint(trackedInterval_1)

        trackedInterval_2 = TrackedInterval([1;1;;3;4])
        @test !isPoint(trackedInterval_2)

        trackedInterval_3 = TrackedInterval([1;2;;3;3])
        @test !isPoint(trackedInterval_3)

        trackedInterval_4 = TrackedInterval([1;2;;1;2])
        @test !isPoint(trackedInterval_4)

        trackedInterval_5 = TrackedInterval([1;1;;3;3])
        @test isPoint(trackedInterval_5)

        trackedInterval_6 = TrackedInterval([2;2;;3;4;;5;5;;7;7])
        @test !isPoint(trackedInterval_6)

        trackedInterval_7 = TrackedInterval([1;1;;pi;pi;;5;5;;7;7])
        @test isPoint(trackedInterval_7)

        trackedInterval_8 = TrackedInterval([1+10^-20;1;;3-10^-20;3])
        @test isPoint(trackedInterval_8)
    end
end